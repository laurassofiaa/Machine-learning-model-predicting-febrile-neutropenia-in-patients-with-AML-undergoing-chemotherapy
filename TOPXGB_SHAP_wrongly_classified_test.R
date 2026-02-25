# ----------------------------
# Libraries
# ----------------------------
library(xgboost)
library(SHAPforxgboost)
library(ggplot2)
library(data.table)
library(dplyr)
library(jsonlite)
library(scales)
library(grid)
library(cli)

# ----------------------------
# Config
# ----------------------------
POINT_SIZE <- 1.25
POINT_ALPHA <- 0.75

FOLDER_TAG <- "NOimputation"
TRAIN_TAG  <- "train2_test2"

ROOT <- "/path/to/infection_model"
DATA_PATH <- file.path(ROOT, "data", TRAIN_TAG, FOLDER_TAG)

# input files (CHANGE if your filenames differ)
MODEL_FILE    <- file.path(DATA_PATH, "simple_xgb_model.json")
FEATURES_FILE <- file.path(DATA_PATH, "simple_features.json")

X_TRAIN_FILE  <- file.path(DATA_PATH, "X_train.csv")
X_TEST_FILE   <- file.path(DATA_PATH, "X_test_drop.csv")

# must contain columns: y_original, y_pred
TEST_PRED_FILE <- file.path(DATA_PATH, "eval_df.csv")

# plots output
PLOT_BASE       <- file.path(ROOT, "plots", TRAIN_TAG, FOLDER_TAG, "wrong_classified_shap")
PLOT_ONLY_WRONG <- file.path(PLOT_BASE, "only_wrong")
PLOT_HIGHLIGHT  <- file.path(PLOT_BASE, "highlight_wrong")

dir.create(PLOT_ONLY_WRONG, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_HIGHLIGHT,  recursive = TRUE, showWarnings = FALSE)

# number of top features to plot (same “top” idea as your training script)
TOP_K <- 6L

# Downsample for plotting (set NULL to disable)
MAX_POINTS_PER_PLOT <- 10000L

# Colours for errors
FP_COLOUR <- "red"
FN_COLOUR <- "deepskyblue"

# Reproducibility for sampling (plotting only)
set.seed(123)

# ----------------------------
# Small helpers
# ----------------------------
`%||%` <- function(x, y) {
  if (length(x) == 0) return(y)
  if (is.null(x)) return(y)
  if (all(is.na(x))) return(y)
  x
}

read_feature_list <- function(path) {
  txt <- readLines(path, warn = FALSE)
  txt <- paste(txt, collapse = "")
  
  # Try JSON first
  cols <- tryCatch(jsonlite::fromJSON(txt), error = function(e) NULL)
  if (!is.null(cols)) return(cols)
  
  # Fallback: Python list style -> JSON-ish
  txt2 <- gsub("'", "\"", txt)
  jsonlite::fromJSON(txt2)
}

# Adds a vertical arrow outside the plot area indicating lower/higher risk
# relative to y=0 on SHAP axis.
add_risk_arrow_outside_from_zero <- function(
    g,
    offset_mm = 8,
    gap_npc = 0.02,
    trim_npc = 0.05,
    low_lbl  = "Lower risk",
    high_lbl = "Higher risk",
    lwd = 1.2,
    head_mm = 1
) {
  gb <- ggplot2::ggplot_build(g)
  yr <- gb$layout$panel_scales_y[[1]]$range$range
  
  # map y=0 to npc [0..1]
  y0_npc <- (0 - yr[1]) / diff(yr)
  y0_npc <- min(max(y0_npc, 0), 1)
  
  # panel borders in npc
  y_bot <- 0 + trim_npc
  y_top <- 1 - trim_npc
  
  # gap around origin
  y_up_start <- min(1, y0_npc + gap_npc)
  y_dn_start <- max(0, y0_npc - gap_npc)
  
  # end points
  y_up_end <- y_top
  y_dn_end <- y_bot
  
  x_out <- grid::unit(1, "npc") + grid::unit(offset_mm, "mm")
  gp <- grid::gpar(col = "black", fill = "black", lwd = lwd)
  
  make_seg <- function(y0, y1) {
    if (abs(y1 - y0) < 1e-6) return(NULL)
    grid::segmentsGrob(
      x0 = x_out, x1 = x_out,
      y0 = grid::unit(y0, "npc"),
      y1 = grid::unit(y1, "npc"),
      gp = gp,
      arrow = grid::arrow(ends = "last", type = "closed", length = grid::unit(head_mm, "mm"))
    )
  }
  
  seg_up <- make_seg(y_up_start, y_up_end)
  seg_dn <- make_seg(y_dn_start, y_dn_end)
  
  txt_low <- grid::textGrob(low_lbl,
                            x = x_out + grid::unit(3, "mm"),
                            y = grid::unit(y_dn_end, "npc") * 1.5,
                            rot = 90, gp = grid::gpar(cex = 0.7)
  )
  txt_high <- grid::textGrob(high_lbl,
                             x = x_out + grid::unit(3, "mm"),
                             y = grid::unit(y_up_end, "npc") * 0.9,
                             rot = 90, gp = grid::gpar(cex = 0.7)
  )
  
  g +
    { if (!is.null(seg_up)) ggplot2::annotation_custom(seg_up) else NULL } +
    { if (!is.null(seg_dn)) ggplot2::annotation_custom(seg_dn) else NULL } +
    ggplot2::annotation_custom(txt_low) +
    ggplot2::annotation_custom(txt_high) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = grid::unit(c(5.5, 34, 5.5, 5.5), "pt"))
}

# Choose which columns to use for overlaying points
get_xy_cols <- function(df) {
  if (!("value" %in% names(df))) stop("Expected column `value` in shap_long.")
  y_candidates <- c("rfvalue", "phi", "shap_value")
  y_col <- y_candidates[y_candidates %in% names(df)][1]
  if (is.na(y_col)) stop("Expected one of: rfvalue, phi, shap_value in shap_long.")
  list(x = "value", y = y_col)
}

make_manual_legend_grob <- function(
    title = "Legend",
    correct_lbl = "Correct",
    fp_lbl = "False positive",
    fn_lbl = "False negative",
    mean_lbl = "Mean",
    correct_col = "black",
    fp_col = "red",
    fn_col = "blue",
    mean_col = "darkgray"
) {
  # Layout inside the legend box (0..1 NPC)
  x_symbol <- grid::unit(0.12, "npc")
  x_text   <- grid::unit(0.22, "npc")
  
  y_title  <- grid::unit(0.88, "npc")
  y1 <- grid::unit(0.68, "npc")
  y2 <- grid::unit(0.50, "npc")
  y3 <- grid::unit(0.32, "npc")
  y4 <- grid::unit(0.14, "npc")
  
# 
#   bg <- grid::rectGrob(
#     x = grid::unit(0, "npc"), y = grid::unit(0, "npc"),
#     width = grid::unit(1, "npc"), height = grid::unit(1, "npc"),
#     just = c("left", "bottom"),
#     gp = grid::gpar(fill = "white", col = "grey70", lwd = 0.7)
#   )
  bg <- grid::nullGrob()

  # Title
  ttl <- grid::textGrob(
    title, x = x_text, y = y_title,
    just = c("left", "center"),
    gp = grid::gpar(fontface = "bold", cex = 0.9, col = "black")
  )
  
  # Symbols
  dot_correct <- grid::pointsGrob(x_symbol, y1, pch = 16, size = grid::unit(3, "mm"),
                                  gp = grid::gpar(col = correct_col))
  dot_fp <- grid::pointsGrob(x_symbol, y2, pch = 16, size = grid::unit(3, "mm"),
                             gp = grid::gpar(col = fp_col))
  dot_fn <- grid::pointsGrob(x_symbol, y3, pch = 16, size = grid::unit(3, "mm"),
                             gp = grid::gpar(col = fn_col))
  
  line_mean <- grid::segmentsGrob(
    x0 = grid::unit(0.06, "npc"), x1 = grid::unit(0.18, "npc"),
    y0 = y4, y1 = y4,
    gp = grid::gpar(col = mean_col, lwd = 2)
  )
  
  # Labels
  t1 <- grid::textGrob(correct_lbl, x = x_text, y = y1, just = c("left", "center"),
                       gp = grid::gpar(cex = 0.85, col = "black"))
  t2 <- grid::textGrob(fp_lbl, x = x_text, y = y2, just = c("left", "center"),
                       gp = grid::gpar(cex = 0.85, col = "black"))
  t3 <- grid::textGrob(fn_lbl, x = x_text, y = y3, just = c("left", "center"),
                       gp = grid::gpar(cex = 0.85, col = "black"))
  t4 <- grid::textGrob(mean_lbl, x = x_text, y = y4, just = c("left", "center"),
                       gp = grid::gpar(cex = 0.85, col = "black"))
  
  grid::grobTree(bg, ttl, dot_correct, dot_fp, dot_fn, line_mean, t1, t2, t3, t4)
}

add_manual_legend_box <- function(
    g,
    x_npc = 0.72, y_npc = 0.55,   # bottom-left position inside the panel (0..1)
    w_npc = 0.26, h_npc = 0.22,   # width/height inside the panel (0..1)
    ...
) {
  gb <- ggplot2::ggplot_build(g)
  xr <- gb$layout$panel_scales_x[[1]]$range$range
  yr <- gb$layout$panel_scales_y[[1]]$range$range
  
  xmin <- xr[1] + x_npc * diff(xr)
  xmax <- xr[1] + (x_npc + w_npc) * diff(xr)
  ymin <- yr[1] + y_npc * diff(yr)
  ymax <- yr[1] + (y_npc + h_npc) * diff(yr)
  
  leg <- make_manual_legend_grob(...)
  
  g +
    ggplot2::annotation_custom(leg, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax) +
    ggplot2::coord_cartesian(clip = "off")
}


# ----------------------------
# Load model + features
# ----------------------------
cli::cli_alert_info("Loading model and feature list...")
model1 <- xgb.load(MODEL_FILE)
cols   <- read_feature_list(FEATURES_FILE)

# ----------------------------
# Load train/test X (ensure same columns/order)
# ----------------------------
cli::cli_alert_info("Loading X_train / X_test...")
dt_x_train <- data.table::fread(X_TRAIN_FILE)[, -1, with = FALSE]
dt_x_test  <- data.table::fread(X_TEST_FILE)[, -1, with = FALSE]

dt_x_train <- dt_x_train[, ..cols]
dt_x_test  <- dt_x_test[,  ..cols]

x_train <- as.matrix(dt_x_train)
x_test  <- as.matrix(dt_x_test)

# ----------------------------
# Load test preds/meta and build FP/FN flags from eval_df
# ----------------------------
cli::cli_alert_info("Loading test predictions/meta...")
dt_test_pred <- data.table::fread(TEST_PRED_FILE)
dt_test_pred <- as.data.table(dt_test_pred)

stopifnot(nrow(dt_test_pred) == nrow(x_test))

if (!all(c("y_original", "y_pred") %in% names(dt_test_pred))) {
  stop("TEST_PRED_FILE must contain columns: y_original, y_pred")
}

dt_test_pred[, y_original := as.integer(y_original)]
dt_test_pred[, y_pred     := as.integer(y_pred)]

# Derive FP/FN directly from eval_df (no need for fp/fn columns)
dt_test_pred[, fp := as.integer(y_original == 0L & y_pred == 1L)]
dt_test_pred[, fn := as.integer(y_original == 1L & y_pred == 0L)]

dt_test_pred[, err_type := fifelse(fp == 1L, "False positive",
                                   fifelse(fn == 1L, "False negative", NA_character_))]

dt_test_pred[, is_wrong := !is.na(err_type)]
wrong_flag <- dt_test_pred$is_wrong

cli::cli_alert_info("Wrong rows: {sum(wrong_flag)} / {length(wrong_flag)}")
cli::cli_alert_info("FP: {sum(dt_test_pred$fp)}  |  FN: {sum(dt_test_pred$fn)}")

# ----------------------------
# Compute SHAP
#  - use TRAIN to get top features
#  - compute TEST SHAP for plotting
# ----------------------------
cli::cli_alert_info("Computing SHAP on train (to get top features)...")
shap_train <- shap.values(xgb_model = model1, X_train = x_train)

TOP_K <- min(TOP_K, length(shap_train$mean_shap_score))
top_feats <- names(shap_train$mean_shap_score)[1:TOP_K]
cli::cli_alert_info("Top features: {paste(top_feats, collapse = ', ')}")

cli::cli_alert_info("Computing SHAP on test (for plots)...")
shap_test <- shap.values(xgb_model = model1, X_train = x_test)

shap_long_test <- shap.prep(
  shap_contrib = shap_test$shap_score,
  X_train = x_test
) %>%
  dplyr::filter(variable %in% top_feats)

# attach wrong flag + error type by row index (ID is 1..nrow(x_test))
shap_long_test$is_wrong  <- wrong_flag[shap_long_test$ID]
shap_long_test$err_type  <- dt_test_pred$err_type[shap_long_test$ID]

# ----------------------------
# Nice axis titles (edit as needed)
# ----------------------------
feature_titles <- c(
  "temperature"                   = "Temperature (°C)",
  "b_neut"                        = "ANC (10⁹/L)",
  "b_leuk"                        = "WBC (10⁹/L)",
  "temperature_trend_from_2_d"     = "Temperature\n2-day trend (°C)",
  "b_leuk_trend_from_7_d"          = "WBC 7-day trend (10⁹/L)",
  "p_crp_trend_from_3_d"           = "CRP 3-day trend (mg/L)",
  "sykli_IND"                      = "Induction phase",
  "bneut_accumulated_max_050_cycle"= "Duration of neutropenia\n(<0.5 10⁹/L) per cycle (days)",
  "bneut_accumulated_max_050_total"= "Duration of neutropenia\n(<0.5 10⁹/L) per regimen (days)",
  "countdown"                      = "Time from treatment initiation (days)"
)

# ----------------------------
# Plot loop (both methods)
# ----------------------------
cli::cli_alert_info("Starting plot loop...")

for (k in seq_along(top_feats)) {
  feat <- top_feats[k]
  cli::cli_alert_info("Plotting feature: {feat}")
  
  sl <- shap_long_test[shap_long_test$variable == feat, , drop = FALSE]
  
  sl$pt_label <- ifelse(is.na(sl$err_type), "Correct", sl$err_type)
  sl$pt_label <- factor(sl$pt_label, levels = c("Correct", "False positive", "False negative"))

  
  # (Optional) downsample for speed/clarity
  # Keep ALL wrong points; sample from correct points if needed.
  if (!is.null(MAX_POINTS_PER_PLOT) && nrow(sl) > MAX_POINTS_PER_PLOT) {
    set.seed(1000 + k)
    sl_wrong   <- sl[sl$is_wrong, , drop = FALSE]
    sl_correct <- sl[!sl$is_wrong, , drop = FALSE]
    
    remaining <- MAX_POINTS_PER_PLOT - nrow(sl_wrong)
    if (remaining <= 0) {
      # too many wrong points: sample wrong points only
      idx <- sample(nrow(sl_wrong), MAX_POINTS_PER_PLOT)
      sl <- sl_wrong[idx, , drop = FALSE]
    } else {
      idx <- sample(nrow(sl_correct), min(remaining, nrow(sl_correct)))
      sl <- rbind(sl_wrong, sl_correct[idx, , drop = FALSE])
    }
  }
  
  xy_cols <- get_xy_cols(sl)
  xcol <- xy_cols$x
  ycol <- xy_cols$y
  
  # --------------------------
  # Method B: ALL points + highlight FP/FN
  # --------------------------
  g_all <- shap.plot.dependence(
    data_long = sl,
    x = feat
  )
  
  
  pt_layers <- which(sapply(g_all$layers, function(l) inherits(l$geom, "GeomPoint")))
  if (length(pt_layers) > 0) {
    g_all$layers[[pt_layers[1]]]$aes_params$size <- POINT_SIZE
    g_all$layers[[pt_layers[1]]]$geom_params$size <- POINT_SIZE
  }
  
  # styling
  x_title <- feature_titles[feat] %||% feat
  g_all <- g_all +
    ggplot2::theme_bw() +
    ggplot2::xlab(x_title) +
    ggplot2::ylab("SHAP value") +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(size = 12, colour = "black"),
      axis.text.y  = ggplot2::element_text(size = 12, colour = "black"),
      axis.title   = ggplot2::element_text(size = 14, colour = "black"),
      axis.line    = ggplot2::element_line(colour = "black"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30", alpha = 0.5, linetype = "dashed")
  
  # push hline behind points (same trick you used)
  nL <- length(g_all$layers)
  if (nL >= 2) g_all$layers <- c(g_all$layers[nL], g_all$layers[-nL])
  
  # smooth line
  g_all <- g_all + ggplot2::geom_smooth(
    se = FALSE,
    method = "loess",
    colour = "darkgray",
    linewidth = 0.7
  )
  
  # feature-specific x-axis styling (safer with scale_x_continuous + pseudo_log_trans)
  if (feat == "temperature") {
    g_all <- g_all + ggplot2::scale_x_continuous(
      breaks = c(35, 36, 37, 38, 39, 40, 41),
      labels = scales::label_number(accuracy = 1, trim = TRUE)
    )
  }
  
  if (feat == "b_leuk") {
    g_all <- g_all +
      ggplot2::scale_x_continuous(
        trans = scales::pseudo_log_trans(base = 10, sigma = 1e-3),
        breaks = c(0, 0.01, 0.1, 1, 10, 100),
        labels = scales::label_number(accuracy = 0.01, trim = TRUE)
      ) +
      ggplot2::xlab(paste0(x_title, " (log)"))
  }
  
  if (feat == "b_neut") {
    g_all <- g_all +
      ggplot2::scale_x_continuous(
        trans = scales::pseudo_log_trans(base = 10, sigma = 1e-2),
        breaks = c(0, 0.05, 0.5, 5, 50),
        labels = scales::label_number(accuracy = 0.01, trim = TRUE)
      ) +
      ggplot2::xlab(paste0(x_title, " (log)")) +
      ggplot2::theme(
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.text = ggplot2::element_text(size = 7)
      )
  }
  
  # overlay wrong points: FP red, FN blue (use *fill* so we don't clash with any existing colour scale)
  sl_wrong <- subset(sl, is_wrong)
  xy_cols <- get_xy_cols(sl)
  xcol <- xy_cols$x   # should be "value"
  ycol <- xy_cols$y   # should be "rfvalue" (or "phi"/"shap_value")
  
  if (nrow(sl_wrong) > 0) {
    g_all <- g_all +
      ggplot2::geom_point(
        data = sl_wrong,
        mapping = ggplot2::aes_string(x = ycol, y = xcol, fill = "err_type"),
        inherit.aes = FALSE,
        shape = 21,
        colour = "black",
        stroke = 0.1,
        alpha = POINT_ALPHA,
        size = POINT_SIZE,
        show.legend = TRUE
      ) +
      ggplot2::scale_fill_manual(
        values = c("False positive" = FP_COLOUR, "False negative" = FN_COLOUR),
        name = "",
        breaks = c("False positive", "False negative")
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))
  }
  
  legend_pos_map <- list(
    temperature = list(x = 0.75, y = 0.03, w = 0.28, h = 0.24), # OK
    b_leuk      = list(x = 0.0001, y = 0.03, w = 0.005, h = 0.24),  #OK
    b_neut      = list(x = 0.0005, y = 0.03, w = 0.03, h = 0.24), # OK
    countdown   = list(x = 0.75, y = 0.03, w = 0.28, h = 0.24), # OK
    p_crp_trend_from_3_d   = list(x = 0.75, y = 0.03, w = 0.28, h = 0.24), # OK
    sykli_IND   = list(x = 0.75, y = 0.03, w = 0.28, h = 0.24) # OK
  )
  
  legend_pos_default <- list(x = 0.72, y = 0.52, w = 0.26, h = 0.24)
  lp <- legend_pos_map[[feat]]
  if (is.null(lp)) lp <- legend_pos_default
  
  # add risk arrow (outside panel)
  g_all <- add_risk_arrow_outside_from_zero(g_all, gap_npc = 0.025, trim_npc = 0.06, offset_mm = 5)
  
  # Remove any ggplot-built legend (so only your manual one shows)
  g_all <- g_all + ggplot2::theme(legend.position = "none")
  
  # Add manual legend box (tweak x_npc/y_npc to move it)
  g_all <- add_manual_legend_box(
    g_all,
    x_npc = lp$x, y_npc = lp$y,
    w_npc = lp$w, h_npc = lp$h,
    title = "",
    correct_col = "black",
    fp_col = FP_COLOUR,
    fn_col = FN_COLOUR,
    mean_col = "darkgray"
  )
  
  
  out_all <- file.path(PLOT_HIGHLIGHT, paste0("SHAP_test_highlight_wrong_", feat, ".png"))
  ggplot2::ggsave(filename = out_all, plot = g_all, width = 6, height = 5, dpi = 600)
  
  # --------------------------
  # Method A: ONLY wrong points (FP/FN colored + legend)
  # --------------------------
  sl_wrong2 <- subset(sl, is_wrong)
  
  if (nrow(sl_wrong2) == 0) {
    cli::cli_alert_warning("No wrong points for {feat} — skipping Method A plot.")
    next
  }
  
  if (!is.null(MAX_POINTS_PER_PLOT) && nrow(sl_wrong2) > MAX_POINTS_PER_PLOT) {
    set.seed(2000 + k)
    idx2 <- sample(nrow(sl_wrong2), MAX_POINTS_PER_PLOT)
    sl_wrong2 <- sl_wrong2[idx2, , drop = FALSE]
  }
  
  g_wrong <- ggplot2::ggplot(sl_wrong2, ggplot2::aes_string(x = xcol, y = ycol)) +
    ggplot2::geom_point(
      ggplot2::aes(fill = err_type),
      shape = 21,
      colour = "black",
      stroke = 0.25,
      alpha = POINT_ALPHA,
      size = POINT_SIZE
    ) +
    ggplot2::scale_fill_manual(
      values = c("False positive" = FP_COLOUR, "False negative" = FN_COLOUR),
      name = "Errors",
      breaks = c("False positive", "False negative")
    ) +
    ggplot2::theme_bw() +
    ggplot2::xlab(x_title) +
    ggplot2::ylab("SHAP value") +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(size = 12, colour = "black"),
      axis.text.y  = ggplot2::element_text(size = 12, colour = "black"),
      axis.title   = ggplot2::element_text(size = 14, colour = "black"),
      axis.line    = ggplot2::element_line(colour = "black"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30", alpha = 0.5, linetype = "dashed") +
    ggplot2::geom_smooth(
      se = FALSE,
      method = "loess",
      colour = "darkgray",
      linewidth = 0.7
    )
  
  # feature-specific x-axis styling (same as above)
  if (feat == "temperature") {
    g_wrong <- g_wrong + ggplot2::scale_x_continuous(
      breaks = c(35, 36, 37, 38, 39, 40, 41),
      labels = scales::label_number(accuracy = 1, trim = TRUE)
    )
  }
  
  if (feat == "b_leuk") {
    g_wrong <- g_wrong +
      ggplot2::scale_x_continuous(
        trans = scales::pseudo_log_trans(base = 10, sigma = 1e-3),
        breaks = c(0, 0.01, 0.1, 1, 10, 100),
        labels = scales::label_number(accuracy = 0.01, trim = TRUE)
      ) +
      ggplot2::xlab(paste0(x_title, " (log)"))
  }
  
  if (feat == "b_neut") {
    g_wrong <- g_wrong +
      ggplot2::scale_x_continuous(
        trans = scales::pseudo_log_trans(base = 10, sigma = 1e-2),
        breaks = c(0, 0.05, 0.5, 5, 50),
        labels = scales::label_number(accuracy = 0.01, trim = TRUE)
      ) +
      ggplot2::xlab(paste0(x_title, " (log)")) +
      ggplot2::theme(
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.text = ggplot2::element_text(size = 7)
      )
  }
  
  g_wrong <- add_risk_arrow_outside_from_zero(g_wrong, gap_npc = 0.025, trim_npc = 0.06, offset_mm = 9)
  
  out_wrong <- file.path(PLOT_ONLY_WRONG, paste0("SHAP_test_only_wrong_", feat, ".png"))
  ggplot2::ggsave(filename = out_wrong, plot = g_wrong, width = 6.5, height = 5, dpi = 600)
  
  gc()
}

cli::cli_alert_success("Done. Plots saved under: {PLOT_BASE}")
# =========================
# Libraries
# =========================
library(caret)
library(randomForestSRC)
library(tidyr)
library(pROC)
library(MLmetrics)
library(xgboost)
library(doParallel)
library(PRROC)
library(cutpointr)
library(SHAPforxgboost)
library(ggplot2)
library(cli)
library(data.table)
library(dplyr)
library(jsonlite)
library(grid)
library(scales)

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && a != "") a else b

# =========================
# Your arrow function
# =========================
add_risk_arrow_outside_from_zero <- function(
    g,
    offset_mm = 8,
    gap_npc   = 0.02,
    trim_npc  = 0.05,
    low_lbl   = "Lower risk",
    high_lbl  = "Higher risk",
    lwd       = 1.2,
    head_mm   = 1
) {
  gb <- ggplot2::ggplot_build(g)
  yr <- gb$layout$panel_scales_y[[1]]$range$range
  
  y0_npc <- (0 - yr[1]) / diff(yr)
  y0_npc <- min(max(y0_npc, 0), 1)
  
  y_bot <- 0 + trim_npc
  y_top <- 1 - trim_npc
  
  y_up_start <- min(1, y0_npc + gap_npc)
  y_dn_start <- max(0, y0_npc - gap_npc)
  
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
      arrow = grid::arrow(ends = "last", type = "closed",
                          length = grid::unit(head_mm, "mm"))
    )
  }
  
  seg_up   <- make_seg(y_up_start, y_up_end)
  seg_down <- make_seg(y_dn_start, y_dn_end)
  
  txt_low <- grid::textGrob(low_lbl,  x = x_out + grid::unit(3, "mm"),
                            y = grid::unit(y_dn_end, "npc") * 1.5,
                            rot = 90, gp = grid::gpar(cex = 0.7))
  txt_high <- grid::textGrob(high_lbl, x = x_out + grid::unit(3, "mm"),
                             y = grid::unit(y_up_end, "npc") * 0.9,
                             rot = 90, gp = grid::gpar(cex = 0.7))
  
  g +
    { if (!is.null(seg_up)) ggplot2::annotation_custom(seg_up) else NULL } +
    { if (!is.null(seg_down)) ggplot2::annotation_custom(seg_down) else NULL } +
    ggplot2::annotation_custom(txt_low) +
    ggplot2::annotation_custom(txt_high) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = grid::unit(c(5.5, 34, 5.5, 5.5), "pt"))
}

# =========================
# Tags / paths
# =========================
FOLDER_TAG <- "Noimputation"
TRAIN_TAG  <- "train2_test2"

ROOT <- "/path/to/infection_model"
DATA_PATH <- file.path(ROOT, "data", TRAIN_TAG, FOLDER_TAG)
PLOT_PATH <- file.path(ROOT, "plots", TRAIN_TAG, FOLDER_TAG, "SHAP_contributions_TOP_yzero")
if (!dir.exists(PLOT_PATH)) dir.create(PLOT_PATH, recursive = TRUE)

# =========================
# Load model + X_train + features list
# =========================
model1 <- xgb.load(file.path(DATA_PATH, "simple_xgb_model.json"))

txt <- readLines(file.path(DATA_PATH, "simple_features.json"), warn = FALSE)
txt <- paste(txt, collapse = "")
txt <- gsub("'", "\"", txt)
cols <- fromJSON(txt)

dt_x <- data.table::fread(file.path(DATA_PATH, "X_train.csv"))[, -1, with = FALSE]
x_train <- dt_x[, ..cols]
x_train <- as.matrix(x_train)
# =========================
# SHAP computation
# =========================
shap_values <- shap.values(xgb_model = model1, X_train = x_train)
shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = x_train)

# =========================
# Save selected features over threshold
# =========================
THRESHOLD <- 0.01
ms  <- shap_values$mean_shap_score
sel <- sort(ms[ms > THRESHOLD], decreasing = TRUE)

#fwrite(data.table(feature = names(sel)),
#       file.path(DATA_PATH, sprintf("SHAP_selected_features_over_%g.csv", THRESHOLD)))

# =========================
# Nice titles
# =========================
name_map <- c(
  "temperature"                      = "Temperature (°C)",
  "b_neut"                           = "ANC (10^9/L)",
  "b_leuk"                           = "WBC (10^9/L)",
  "temperature_trend_from_2_d"       = "Temperature\n2-day trend (°C)",
  "b_leuk_trend_from_7_d"            = "WBC 7-day trend (10^9/L)",
  "p_crp_trend_from_3_d"             = "CRP 3-day trend (mg/L)",
  "p_crp_trend_from_7_d"             = "CRP 7-day trend (mg/L)",
  "sykli_IND"                        = "Induction phase",
  "bneut_accumulated_max_050_cycle"  = "Duration of neutropenia\n(≤0.5 10^9/L) per cycle (days)",
  "bneut_accumulated_max_050_total"  = "Duration of neutropenia\n(≤0.5 10^9/L) per regimen (days)",
  "countdown"                        = "Time since\ntreatment initiation (days)",
  "b_trom"                           ="PLT (10⁹/L)"
)

# =========================
# Collect crossings for all features
# =========================
zero_table <- data.table(feature = character(), x_at_zero = character())

# =========================
# Plot top 10 dependence plots
# =========================
top10 <- names(shap_values$mean_shap_score)[1:6]

for (feat in top10) {
  message("Plotting SHAP dependence for: ", feat)
  
  sl <- shap_long[shap_long$variable == feat, ]
  
  # Optional downsample
  if (nrow(sl) > 10000) {
    set.seed(1)
    idx <- sample(nrow(sl), 10000)
    sl <- sl[idx, , drop = FALSE]
  }
  
  # ------------------------------------------------------------
  # 1) Build base dependence plot (includes the red smooth curve)
  # ------------------------------------------------------------
  g <- shap.plot.dependence(
    data_long  = sl,
    x          = feat
  )
  
  if (feat == "b_neut") {
    x_point <- 0.15
    
    g <- g +
      ggplot2::geom_vline(xintercept = x_point, linetype = "dashed", colour = "grey") +
      ggplot2::annotate(
        "text",
        x = x_point, y = -1,
        label = "x=0.15",
        angle = 90, vjust = -0.5, size = 3, colour = "grey"
      )
  }
  # ------------------------------------------------------------
  # 2) Extract the ACTUAL red-line layer data and solve for y=0
  #    (do this BEFORE any scale_x_* transforms)
  # ------------------------------------------------------------
  gb <- ggplot2::ggplot_build(g)
  
  # candidate layers: must have x,y and many unique x
  line_layer_idx <- which(vapply(gb$data, function(d) {
    is.data.frame(d) &&
      all(c("x", "y") %in% names(d)) &&
      nrow(d) > 20 &&
      length(unique(d$x)) > 20
  }, logical(1)))
  
  x_at_zero <- numeric(0)
  
  if (length(line_layer_idx) > 0) {
    # usually the smooth curve layer has fewer rows than the scatter layer
    li <- line_layer_idx[which.min(vapply(gb$data[line_layer_idx], nrow, integer(1)))]
    dline <- gb$data[[li]][, c("x", "y")]
    dline <- dline[is.finite(dline$x) & is.finite(dline$y), ]
    dline <- dline[order(dline$x), ]
    
    x <- dline$x
    y <- dline$y
    
    # crossing intervals
    idx <- which(sign(y[-1]) * sign(y[-length(y)]) < 0)
    
    if (length(idx) > 0) {
      # linear interpolation for x when y=0
      x_at_zero <- vapply(idx, function(i) {
        x[i] + (0 - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
      }, numeric(1))
      x_at_zero <- sort(unique(x_at_zero))
    }
  }
  
  # store results
  zero_table <- rbind(
    zero_table,
    data.table(
      feature = feat,
      x_at_zero = if (length(x_at_zero)) paste(signif(x_at_zero, 6), collapse = ";") else NA_character_
    )
  )
  
  # print for visibility
  if (length(x_at_zero) > 0) {
    message(feat, ": red-line y=0 at x = ", paste(signif(x_at_zero, 4), collapse = ", "))
  } else {
    message(feat, ": no red-line y=0 crossing found")
  }
  
  # ------------------------------------------------------------
  # 3) Now add your styling + mark crossings on plot
  # ------------------------------------------------------------
  g <- g +
    ggplot2::theme_bw() +
    ggplot2::xlab(feature_titles[feat] %||% feat) +
    ggplot2::ylab("SHAP value") +
    ggplot2::theme(
      axis.text.x   = ggplot2::element_text(size = 12, colour = "black"),
      axis.text.y   = ggplot2::element_text(size = 12, colour = "black"),
      axis.title    = ggplot2::element_text(size = 14, colour = "black"),
      axis.line     = ggplot2::element_line(colour = "black"),
      panel.grid    = ggplot2::element_blank(),
      panel.border  = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30",
                        alpha = 0.5, linetype = "dashed")
  
  # mark crossings (if any)
  if (length(x_at_zero) > 0) {
    g <- g +
      ggplot2::geom_vline(xintercept = x_at_zero, linetype = "dotted") +
      ggplot2::annotate(
        "text",
        x = x_at_zero, y = 0,
        label = signif(x_at_zero, 3),
        angle = 90, vjust = -0.5, size = 3
      )
  }
  
  # Optional: your ylim tweak for ANC
  if (feat == "b_neut") {
    g <- g + ggplot2::coord_cartesian(ylim = c(-0.7, 0.2))
  }
  
  # ------------------------------------------------------------
  # 4) Your x-axis transforms/labels
  # ------------------------------------------------------------
  if (feat == "temperature") {
    g <- g + ggplot2::scale_x_continuous(
      breaks = c(35, 36, 37, 38, 39, 40, 41),
      labels = scales::label_number(accuracy = 1, trim = TRUE)
    )
  }
  
  if (feat == "b_leuk") {
    g <- g + ggplot2::scale_x_log10(
      trans  = scales::pseudo_log_trans(base = 10, sigma = 1e-3),
      breaks = c(0, 0.01, 0.1, 1, 10, 100),
      labels = scales::label_number(accuracy = 0.01, trim = TRUE)
    ) +
      ggplot2::xlab(paste0((feature_titles[feat] %||% feat), " (log)"))
  }
  
  if (feat == "b_neut") {
    g <- g + ggplot2::scale_x_log10(
      trans  = scales::pseudo_log_trans(base = 10, sigma = 1e-2),
      breaks = c(0, 0.05, 0.5, 5, 50),
      labels = scales::label_number(accuracy = 0.01, trim = TRUE)
    ) +
      ggplot2::xlab(paste0((feature_titles[feat] %||% feat), " (log)")) +
      ggplot2::theme(
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.text = ggplot2::element_text(size = 10)
      )
  }
  
  # ------------------------------------------------------------
  # 5) Your risk arrow additions
  # ------------------------------------------------------------
  if (feat == "b_neut") {
    g <- add_risk_arrow_outside_from_zero(g, gap_npc = 0.025, trim_npc = 0.06, offset_mm = 9)
  } else if (feat == "b_leuk") {
    g <- add_risk_arrow_outside_from_zero(g, gap_npc = 0.025, trim_npc = 0.06, offset_mm = 9)
  } else {
    g <- add_risk_arrow_outside_from_zero(g, gap_npc = 0.02, trim_npc = 0.06, offset_mm = 9)
  }
  
  # ------------------------------------------------------------
  # 6) Save plot
  # ------------------------------------------------------------
  out_file <- file.path(PLOT_PATH, paste0("SHAP_contribution_", feat, ".png"))
  ggplot2::ggsave(
    filename = out_file,
    plot     = g,
    width    = 6,
    height   = 5,
    dpi      = 600
  )
  
  gc()
}

# =========================
# Save crossings table
# =========================
fwrite(zero_table, file.path(PLOT_PATH, "redline_y0_crossings.csv"))
message("Saved crossings to: ", file.path(PLOT_PATH, "redline_y0_crossings.csv"))


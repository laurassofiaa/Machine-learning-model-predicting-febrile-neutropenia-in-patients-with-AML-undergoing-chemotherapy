

# Libraries
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

# FOLDER_TAG <- "interpolation"
# FOLDER_TAG <- "uncertainty"
FOLDER_TAG <- "NOimputation" 

# TRAIN_TAG <- "train3_test2"
# TRAIN_TAG <- "train3_test1"
TRAIN_TAG <- "train2_test2"
# TRAIN_TAG <- "train2_test1"
# TRAIN_TAG <- "train1_test1"

# Root once:
ROOT <- "/path/to/infection_model"

# Paths derived from the tag:
DATA_PATH <- file.path(ROOT, "data", TRAIN_TAG, FOLDER_TAG)
PLOT_PATH <- file.path(ROOT, "plots", TRAIN_TAG, FOLDER_TAG, "SHAP_contributions_TOP")
if (!dir.exists(PLOT_PATH)) dir.create(PLOT_PATH, recursive = TRUE)


# Load the model
model1 <- xgb.load(file.path(DATA_PATH, "simple_xgb_model.json"))

# load top variables
txt <- readLines(file.path(DATA_PATH,"simple_features.json"), warn = FALSE)
txt <- paste(txt, collapse = "")

# turn Python-style list into JSON-style
txt <- gsub("'", "\"", txt)

cols <- fromJSON(txt)   # now it's valid JSON text
dt_x <- data.table::fread(file.path(DATA_PATH, "X_train.csv"))[, -1, with = FALSE]
x_train <- dt_x[, ..cols]
x_train = as.matrix(x_train)

# SHAP
shap_values <- shap.values(xgb_model = model1, X_train = x_train)

# The ranked features by mean |SHAP|
shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = x_train) %>%
  dplyr::filter(variable %in% names(shap_values$mean_shap_score[1:6]))

# -----------------------------------

# shap features over 0.5 to csv file
THRESHOLD <- 0.01

# 1) Use the provided mean |SHAP| vector
ms <- shap_values$mean_shap_score                 # named numeric vector
sel <- sort(ms[ms > THRESHOLD], decreasing = TRUE)

# ------------------------------------------

# Define a dictionary for nicer titles
feature_titles <- c(
  "temperature"                       = "Temperature (°C)",
  "b_neut"                            = "ANC (10⁹/L)",
  "b_leuk"                            = "WBC (10⁹/L)",
  "temperature_trend_from_2_d"        = "Temperature \n2-day trend (°C)",
  'b_leuk_trend_from_7_d'             = "WBC 7-day trend (10⁹/L)",
  'p_crp_trend_from_3_d'              = "CRP 3-day trend (mg/L)",
  'sykli_IND'                         = "Induction phase",
  "bneut_accumulated_max_050_cycle"   = "Duration of neutropenia\n(<0.5 10⁹/L) per cycle (days)",
  "bneut_accumulated_max_050_total"   = "Duration of neutropenia\n(<0.5 10⁹/L) per regimen (days)",
  "countdown"                         = "Time from \ntreatment initiation (days)"
)



add_risk_arrow_outside_from_zero <- function(
    g,
    offset_mm   = 8,      # distance to the right of the panel
    gap_npc     = 0.02,   # empty gap around y=0 (fraction of panel height)
    trim_npc    = 0.05,   # shorten arrows at top/bottom
    low_lbl     = "Lower risk",
    high_lbl    = "Higher risk",
    lwd         = 1.2,
    head_mm     = 1
) {
  gb <- ggplot2::ggplot_build(g)
  yr <- gb$layout$panel_scales_y[[1]]$range$range

  # map y=0 to npc [0..1]; clamp to panel
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
  gp    <- grid::gpar(col = "black", fill = "black", lwd = lwd)

  make_seg <- function(y0, y1) {
    # draw regardless of direction; skip only if (near-)zero length
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
  seg_up   <- make_seg(y_up_start, y_up_end)     # 0 → up
  seg_down <- make_seg(y_dn_start, y_dn_end)     # 0 → down (now works)

  txt_low  <- grid::textGrob(low_lbl,  x = x_out-grid::unit(3, "mm"), y = grid::unit(y_dn_end, "npc")*1.5,
                             rot = 90, gp = grid::gpar(cex = 0.7))
  txt_high <- grid::textGrob(high_lbl, x = x_out-grid::unit(3, "mm"), y = grid::unit(y_up_end, "npc")*0.9,
                             rot = 90, gp = grid::gpar(cex = 0.7))
  g +
    { if (!is.null(seg_up))   ggplot2::annotation_custom(seg_up)   else NULL } +
    { if (!is.null(seg_down)) ggplot2::annotation_custom(seg_down) else NULL } +
    ggplot2::annotation_custom(txt_low) +
    ggplot2::annotation_custom(txt_high) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = grid::unit(c(5.5, 34, 5.5, 5.5), "pt"))
}




for (k1 in seq_len(6)) {
  feat <- names(shap_values$mean_shap_score)[k1]
  message("Plotting SHAP dependence for: ", feat)
  
  # 1) Filter to just that feature
  sl <- shap_long[ shap_long$variable == feat, ]
  
  # 2) If >10k rows, randomly pick 10k
  if (nrow(sl) > 10000) {
    idx    <- sample(nrow(sl), 10000)
    sl     <- sl[idx, , drop = FALSE]
  }
  
  # 3) Build the dependence plot
  g <- shap.plot.dependence(
    data_long     = sl,
    x             = feat
  )

    
  g <- g + ggplot2::theme_bw() +
    ggplot2::xlab(feature_titles[feat]) +
    ggplot2::ylab("SHAP value") +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(size = 12, colour = "black"),
      axis.text.y      = ggplot2::element_text(size = 12, colour = "black"),
      axis.title       = ggplot2::element_text(size = 14, colour = "black"),
      axis.line = element_line(colour = "black"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
      
    )
  
  g <- g + geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30", alpha=0.5, linetype = "dashed")
  n <- length(g$layers)
  g$layers <- c(g$layers[n], g$layers[-n])
  
  
  
  if (feat == 'bneut_accumulated_max_050_cycle'){
    g <- g + ggplot2::coord_cartesian(ylim = c(-0.7, 0.2))
    
  }
  
  # custom log x-axis for b-leuk
  if (feat=='temperature'){
    g <- g+ggplot2::scale_x_continuous(
      breaks = c(35, 36, 37, 38, 39, 40, 41),
      labels = scales::label_number(accuracy = 1, trim = TRUE) 
    )
  }
  if (feat == "b_leuk") {
    g <- g + ggplot2::scale_x_log10(
      trans  = scales::pseudo_log_trans(base = 10, sigma = 1e-3),
      breaks = c(0, 0.01, 0.1, 1, 10, 100),
      labels = scales::label_number(accuracy = 0.01, trim = TRUE)  # "0.01", "0.05", "0.1", "10", etc.
    )+
      ggplot2::xlab(paste0(feature_titles[feat], "(log)"))
  }
  if (feat == "b_neut"){
    g <- g +
      ggplot2::scale_x_log10(
        trans  = scales::pseudo_log_trans(base = 10, sigma = 1e-2),
        breaks = c(0, 0.05, 0.5, 5, 50),
        labels = scales::label_number(accuracy = 0.01, trim = TRUE)        )+
      ggplot2::xlab(paste0(feature_titles[feat], "(log)"))+
      theme(
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 10)
      )
  }

  if (feat == "b_neut") {
    # g <- add_risk_arrow(g, x_override=100)
    g <- add_risk_arrow_outside_from_zero(g, gap_npc = 0.025, trim_npc = 0.06, offset_mm = 9)
  }
  else if (feat == "b_leuk") {
    # g <- add_risk_arrow(g, x_override=1000)
    g <- add_risk_arrow_outside_from_zero(g, gap_npc = 0.025, trim_npc = 0.06, offset_mm = 9)
    # default
  }
  else {
    # g <- add_risk_arrow(g)   # ← unconditionally add the arrow
    g <- add_risk_arrow_outside_from_zero(g, gap_npc = 0.025, trim_npc = 0.06, offset_mm = 9)
  }

  
  # 5) Save to your specific folder
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
 






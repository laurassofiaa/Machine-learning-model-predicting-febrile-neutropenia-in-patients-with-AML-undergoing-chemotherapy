# Predict AML, MDS and MF using XGBoost
# install.packages("jsonlite")

#pak::pkg_install('jsonlite', "mounts/research/Rpackage")
#pak::pkg_install('xgboost', "mounts/research/Rpackage")

# Libraries
#source("/path/to/library.R")
library(caret)
#library(randomForestSRC)
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
PLOT_PATH <- file.path(ROOT, "plots", TRAIN_TAG, FOLDER_TAG, "SHAP_contributions")
if (!dir.exists(PLOT_PATH)) dir.create(PLOT_PATH, recursive = TRUE)


# Load the model
model1 <- xgb.load(file.path(DATA_PATH, "xgb_model.json"))

x_train <- data.table::fread(file.path(DATA_PATH, "X_train.csv"))[, -1, with = FALSE]
x_train = as.matrix(x_train)

# SHAP
shap_values <- shap.values(xgb_model = model1, X_train = x_train)

# The ranked features by mean |SHAP|
shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = x_train) %>%
  dplyr::filter(variable %in% names(shap_values$mean_shap_score[1:10]))

# -----------------------------------

# shap features over 0.5 to csv file
THRESHOLD <- 0.01

# 1) Use the provided mean |SHAP| vector
ms <- shap_values$mean_shap_score                 # named numeric vector
sel <- sort(ms[ms > THRESHOLD], decreasing = TRUE)

# a) just feature names
fwrite(data.table(feature = names(sel)),
       file.path(DATA_PATH, sprintf("SHAP_selected_features_over_%g.csv", THRESHOLD)))

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
  "bneut_accumulated_max_050_cycle"   = "Duration of neutropenia\n(≤0.5 10⁹/L) per cycle (days)",
  "bneut_accumulated_max_050_total"   = "Duration of neutropenia\n(≤0.5 10⁹/L) per regimen (days)",
  "countdown"                         = "Time since \ntreatment initiation (days)"
)
  


add_risk_arrow <- function(
    g,
    outside_frac = 0.1,        # used only if x_override is NA/NaN
    x_override   = NA_real_,     # set to a number to force the x position
    low_lbl = "Lower risk",
    high_lbl = "Higher risk"
) {
  gb <- ggplot2::ggplot_build(g)
  xr <- gb$layout$panel_scales_x[[1]]$range$range
  yr <- gb$layout$panel_scales_y[[1]]$range$range
  
  # choose x position
  use_default <- is.null(x_override) || is.na(x_override) || is.nan(x_override)
  x_arrow <- if (use_default) xr[2] + outside_frac * (xr[2] - xr[1]) else x_override
  print(x_arrow)
  y_low  <- yr[1] + 0.15 * diff(yr)
  y_high <- yr[2] - 0.15 * diff(yr)
  
  g +
    ggplot2::annotate(
      "segment",
      x = x_arrow, xend = x_arrow,
      y = y_low,   yend = y_high,
      arrow = ggplot2::arrow(ends = "both", type = "closed", length = grid::unit(2, "mm"))
    ) +
    ggplot2::annotate("text", x = x_arrow, y = y_low*1,  label = low_lbl,  vjust = 1.8, angle = 90, size = 3) +
    ggplot2::annotate("text", x = x_arrow, y = y_high*0.8, label = high_lbl, vjust = 1.8, angle = 90, size = 3) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = grid::unit(c(5.5, 30, 5.5, 5.5), "pt"))
}




add_risk_arrow_outside <- function(g,
                                   offset_mm = 6,      # distance to the right of the panel
                                   y_pad_npc = 0.03,   # top/bottom padding as fraction of panel height
                                   low_lbl = "Lower risk",
                                   high_lbl = "Higher risk") {
  
  seg <- grid::segmentsGrob(
    x0 = grid::unit(1, "npc") + grid::unit(offset_mm, "mm"),
    x1 = grid::unit(1, "npc") + grid::unit(offset_mm, "mm"),
    y0 = grid::unit(y_pad_npc, "npc"),
    y1 = grid::unit(1 - y_pad_npc, "npc"),
    gp = grid::gpar(lwd = 1.2),
    arrow = ggplot2::arrow(ends = "both", type = "open", length = grid::unit(2, "mm"))
  )
  
  txt_low  <- grid::textGrob(low_lbl,
                             x = grid::unit(1, "npc") + grid::unit(offset_mm-3, "mm"),
                             y = grid::unit(y_pad_npc, "npc") + grid::unit(15, "mm"),
                             rot = 90, gp = grid::gpar(cex = 0.8))
  txt_high <- grid::textGrob(high_lbl,
                             x = grid::unit(1, "npc") + grid::unit(offset_mm-3, "mm"),
                             y = grid::unit(1 - y_pad_npc, "npc")- grid::unit(15, "mm"),
                             rot = 90, gp = grid::gpar(cex = 0.8))
  
  g +
    ggplot2::annotation_custom(seg) +
    ggplot2::annotation_custom(txt_low) +
    ggplot2::annotation_custom(txt_high) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = grid::unit(c(5.5, 30, 5.5, 5.5), "pt"))
}




for (k1 in seq_len(10)) {
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
  ) +
    ggplot2::theme_bw() +
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
  
  g <- g + geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30", alpha=0.5)
  n <- length(g$layers)
  g$layers <- c(g$layers[n], g$layers[-n])
    
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
        #theme(legend.position = "top") +
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
  
  # 4) Arrow on EVERY plot
  if (feat == "b_neut") {
    # g <- add_risk_arrow(g, x_override=100)
    g <- add_risk_arrow_outside(g, offset_mm=6)
  }
  else if (feat == "b_leuk") {
    # g <- add_risk_arrow(g, x_override=1000)
    g <- add_risk_arrow_outside(g, offset_mm=6)
  }
  else {
    # g <- add_risk_arrow(g)   # ← unconditionally add the arrow
    g <- add_risk_arrow_outside(g, offset_mm=6)
  }
  
  # 5) Save to your specific folder 
  out_file <- file.path(PLOT_PATH, paste0("SHAP_contribution_", feat, ".png"))
  ggplot2::ggsave(
    filename = out_file,
    plot     = g,
    width    = 5,
    height   = 4,
    dpi      = 600
  )
  
  gc()
}
 






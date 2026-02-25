library(data.table)
library(xgboost)
library(SHAPforxgboost)

# ----------------------------
# Config
# ----------------------------
FOLDER_TAG <- "NOimputation"
TRAIN_TAG  <- "train2_test2"

ROOT <- "/path/to/infection_model"

DATA_PATH <- file.path(ROOT, "data",  TRAIN_TAG, FOLDER_TAG)
PLOT_PATH <- file.path(ROOT, "plots", TRAIN_TAG, FOLDER_TAG, "SHAP")
if (!dir.exists(PLOT_PATH)) dir.create(PLOT_PATH, recursive = TRUE)

# ----------------------------
# Name map for prettier labels
# (extend/edit as you like)
# ----------------------------
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
  "b_trom"                           = "PLT (10⁹/L)"
  )

# ----------------------------
# Helper: align matrix to model features
# ----------------------------
align_to_model <- function(X, xgb_mod) {
  X <- as.matrix(X)
  if (is.null(colnames(X))) stop("X must have colnames().")
  
  # sanitize + unique (prevents issues with special chars/duplicates)
  colnames(X) <- make.names(colnames(X), unique = TRUE)
  
  mod_feats_str <- xgboost::xgb.attr(xgb_mod, "feature_names")
  if (!is.null(mod_feats_str) && nzchar(mod_feats_str)) {
    mod_feats <- strsplit(mod_feats_str, ",")[[1]]
    mod_feats <- make.names(mod_feats, unique = TRUE)
    
    missing <- setdiff(mod_feats, colnames(X))
    if (length(missing) > 0) {
      stop("X is missing model features: ", paste(missing, collapse = ", "))
    }
    
    # drop extras + reorder to match model
    X <- X[, mod_feats, drop = FALSE]
  } else {
    # store names safely as an attribute (NOT xgb_mod$feature_names <- ...)
    xgboost::xgb.attr(xgb_mod, "feature_names") <- paste(colnames(X), collapse = ",")
  }
  
  X
}

# ----------------------------
# SHAP analysis function (with name_map + colored barplot)
# ----------------------------
analyze_shap_group <- function(X_mat, group_label, out_dir_base, xgb_mod, top_n = 10, name_map = NULL) {
  
  # Align X to model features
  X_mat <- align_to_model(X_mat, xgb_mod)
  
  # Compute SHAP
  shap_out <- SHAPforxgboost::shap.values(xgb_model = xgb_mod, X_train = X_mat)
  
  # Mean SHAP %
  dt_mean <- data.table(
    feature   = names(shap_out$mean_shap_score),
    mean_shap = 100 * shap_out$mean_shap_score / sum(shap_out$mean_shap_score)
  )[order(-mean_shap)]
  
  top_n <- min(top_n, nrow(dt_mean))
  top_feats <- dt_mean$feature[seq_len(top_n)]
  
  # Remove BIAS column if present
  shap_score <- shap_out$shap_score
  if (!is.null(colnames(shap_score)) && "BIAS" %in% colnames(shap_score)) {
    shap_score <- shap_score[, setdiff(colnames(shap_score), "BIAS"), drop = FALSE]
  }
  
  # Slice to top features
  shap_score_sub_mat <- as.matrix(shap_score)[, top_feats, drop = FALSE]
  X_mat_sub <- X_mat[, top_feats, drop = FALSE]
  
  shap_score_sub_df <- as.data.frame(shap_score_sub_mat)
  X_df_sub <- as.data.frame(X_mat_sub)
  
  # Prepare long format for beeswarm
  shap_long <- SHAPforxgboost::shap.prep(
    shap_contrib = shap_score_sub_df,
    X_train      = X_df_sub
  )
  
  # Output dir
  out_dir <- file.path(out_dir_base, group_label)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # ----------------------------
  # Plot 1: Beeswarm
  # ----------------------------
  png(
    filename = file.path(out_dir, "SHAP_beeswarm.png"),
    width = 1600, height = 1200, res = 150
  )
  print(SHAPforxgboost::shap.plot.summary(shap_long))
  dev.off()
  
  # ----------------------------
  # Plot 2: Mean SHAP barplot (plain)
  # ----------------------------
  vals <- dt_mean$mean_shap[seq_len(top_n)]
  labs_plain <- dt_mean$feature[seq_len(top_n)]
  
  png(
    filename = file.path(out_dir, "Mean_SHAP.png"),
    width = 1600, height = 1200, res = 150
  )
  
  old_par <- par(mar = c(10, 5, 4, 2) + 0.1)
  
  bp <- barplot(
    height = vals,
    names.arg = rep("", length(vals)),
    ylab = "SHAP Contribution (%)",
    main = "Top Feature Importance (by % SHAP)",
    ylim = c(0, max(vals) * 1.15),
    xaxt = "n"
  )
  
  text(x = bp, y = vals + max(vals) * 0.03,
       labels = sprintf("%.3f", vals), pos = 3, cex = 1.3)
  
  axis(1, at = bp, labels = FALSE)
  
  text(
    x = bp,
    y = par("usr")[3] - 0.03 * max(vals),
    labels = labs_plain,
    srt = 45, adj = 1, xpd = TRUE, cex = 0.9
  )
  
  par(old_par)
  dev.off()
  
  # ----------------------------
  # Plot 3: Mean SHAP barplot (colored by direction + pretty labels)
  #   Color rule:
  #     correlation(feature_value, shap_value) < 0  -> steelblue
  #     correlation(feature_value, shap_value) >= 0 -> firebrick
  # ----------------------------
  # pretty labels
  if (!is.null(name_map)) {
    labs_pretty <- ifelse(top_feats %in% names(name_map), name_map[top_feats], top_feats)
  } else {
    labs_pretty <- top_feats
  }
  
  # direction colors
  bar_colors <- sapply(seq_along(top_feats), function(i) {
    f <- top_feats[i]
    # correlation between feature and SHAP contribution
    cor_val <- suppressWarnings(cor(X_df_sub[[f]], shap_score_sub_df[[f]], use = "complete.obs"))
    if (is.na(cor_val)) return("gray")
    if (cor_val < 0) "steelblue" else "firebrick"
  })
  
  # group-specific title
  if (group_label == "IND") {
    title_label <- "Top Feature Importance for Inductions"
  } else if (group_label == "KONS") {
    title_label <- "Top Feature Importance for Consolidations"
  } else {
    title_label <- "Top Feature Importance"
  }
  
  png(
    filename = file.path(out_dir, "Mean_SHAP_percent_colored.png"),
    width = 12, height = 10, units = "in", res = 600
  )
  
  old_par <- par(mar = c(11.5, 5, 4, 0) + 0.1)
  
  bp <- barplot(
    height = vals,
    names.arg = rep("", length(vals)),
    ylab = "SHAP Contribution (%)",
    main = title_label,
    col = bar_colors,
    xaxt = "n",
    ylim = c(0, max(vals) * 1.15),
    font.main = 1,
    cex.main = 1.5,
    cex.lab = 1.7,
    cex.axis = 1.3
  )
  
  # percent labels on bars
  text(
    x = bp,
    y = vals + max(vals) * 0.02,
    labels = sprintf("%.1f%%", vals),
    pos = 3,
    cex = 1.7
  )
  
  axis(1, at = bp, labels = FALSE)
  
  # angled pretty labels
  text(
    x = bp,
    y = par("usr")[3] - 0.04 * max(vals),
    labels = labs_pretty,
    srt = 45, adj = 0.99, xpd = TRUE, cex = 1.15
  )
  
  par(old_par)
  dev.off()
  
  # Return objects for optional further use
  list(
    group = group_label,
    mean_table = dt_mean,
    shap_long = shap_long,
    shap_raw = shap_out
  )
}

# ----------------------------
# Load model + data
# ----------------------------
xgb_mod <- xgb.load(file.path(DATA_PATH, "xgb_model.json"))

dt_x <- data.table::fread(file.path(DATA_PATH, "X_train.csv"))[, -1, with = FALSE]

stopifnot("sykli_IND" %in% colnames(dt_x))

# Split groups
dt_KONS <- dt_x[sykli_IND == 0]
dt_IND  <- dt_x[sykli_IND == 1]

# Matrices (keep sykli_IND here; align_to_model() will drop it if model doesn't use it)
X_mat  <- as.matrix(dt_x)
X_KONS <- as.matrix(dt_KONS)
X_IND  <- as.matrix(dt_IND)

# OPTIONAL: if your model already has feature names stored, keep them.
# If not, this will set them safely once, based on X_mat after alignment:
X_mat_aligned <- align_to_model(X_mat, xgb_mod)
xgboost::xgb.attr(xgb_mod, "feature_names") <- paste(colnames(X_mat_aligned), collapse = ",")

# ----------------------------
# Run SHAP
# ----------------------------
out_dir_base <- PLOT_PATH

res_all  <- analyze_shap_group(X_mat,  "ALL",  out_dir_base, xgb_mod, top_n = 10, name_map = name_map)
res_kons <- analyze_shap_group(X_KONS, "KONS", out_dir_base, xgb_mod, top_n = 10, name_map = name_map)
res_ind  <- analyze_shap_group(X_IND,  "IND",  out_dir_base, xgb_mod, top_n = 10, name_map = name_map)

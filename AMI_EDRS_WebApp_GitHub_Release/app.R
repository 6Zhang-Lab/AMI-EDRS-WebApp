############################################################
## AMI-EDRS Shiny App v13.0 User-input State Reset Dashboard
## Early Dynamic Risk Stratification for Acute Myocardial Infarction
##
## Improvements over v1:
## v7 updates:
## 1) Public-facing labels no longer mention internal pipeline step names.
## 2) Demo buttons automatically run prediction and display results.
## 3) C1-C5 phenotype can be assigned from complete organ-domain dynamics;
##    if insufficient, a clear data-requirement message is shown.
## 4) Demo/example data are downloadable from the interface.
## 5) Results page removes model-version display and uses a deeper executive clinical style.
## 6) Reactive-context bug fixed by isolating reactive reads in button/download callbacks.
## 7) Public interface uses only one patient-data upload entry point; model-object uploads are hidden.
## 8) Result plots use larger text and improved executive color contrast.
## 9) v8: improved left quick-start panel, equal-width large workflow tabs, and differentiated action/demo buttons.
## 10) v9: raw clinical measurements are entered by the user; organ-domain scores, C1-C5 phenotype, and AMI-EDRS risk are calculated automatically.
############################################################

options(stringsAsFactors = FALSE)

## Explicit library calls help shinyapps.io detect package dependencies.
library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(tibble)
library(ggplot2)
library(scales)
library(readr)
library(readxl)
library(glmnet)
library(Matrix)
library(jsonlite)
library(grid)
message("[AMI-EDRS v13] app.R loaded: demo-isolated user-input prediction fix.")

## =========================================================
## 1. Default model paths
## =========================================================

## Deployment-friendly model paths.
## Priority: ./model/ files -> environment variables -> local developer path.
local_model_dir <- file.path(getwd(), "model")
local_step05_rds <- file.path(local_model_dir, "Step05_T48_model_objects_with_C1C5.rds")
local_step03_rds <- file.path(local_model_dir, "Step03_final_k5_trajectory_model_object.rds")

default_project_root <- Sys.getenv(
  "AMI_EDRS_PROJECT_ROOT",
  unset = "/Users/zhangkangnan/Desktop/MIMIC_AMI_DynamicRisk_Project"
)

default_step05_rds <- Sys.getenv(
  "AMI_EDRS_MODEL_RDS",
  unset = file.path(
    default_project_root,
    "09_step05_T48_dynamic_prediction_model_with_C1C5",
    "objects",
    "Step05_T48_model_objects_with_C1C5.rds"
  )
)

default_step03_rds <- Sys.getenv(
  "AMI_EDRS_TRAJECTORY_RDS",
  unset = file.path(
    default_project_root,
    "07_step03_final_k5_trajectory_phenotyping",
    "objects",
    "Step03_final_k5_trajectory_model_object.rds"
  )
)

if (file.exists(local_step05_rds)) default_step05_rds <- local_step05_rds
if (file.exists(local_step03_rds)) default_step03_rds <- local_step03_rds

## =========================================================
## 2. Utility functions
## =========================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

clean_name <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  make.unique(x)
}

to_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

bound_prob <- function(x, eps = 1e-6) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

pretty_pct <- function(x, digits = 1) {
  paste0(formatC(as.numeric(x) * 100, format = "f", digits = digits), "%")
}

risk_group_from_prob <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.05 ~ "Low risk (<5%)",
    p < 0.15 ~ "Intermediate risk (5-15%)",
    p < 0.30 ~ "High risk (15-30%)",
    TRUE ~ "Very high risk (>=30%)"
  )
}

risk_group_label_public <- function(x) {
  dplyr::recode(
    x,
    "Low risk (<5%)" = "Low risk (<5%)",
    "Intermediate risk (5-15%)" = "Intermediate risk (5–15%)",
    "High risk (15-30%)" = "High risk (15–30%)",
    "Very high risk (>=30%)" = "Very high risk (≥30%)",
    .default = x
  )
}

monitoring_intensity <- function(risk_group) {
  dplyr::case_when(
    risk_group == "Low risk (<5%)" ~ "Routine monitoring",
    risk_group == "Intermediate risk (5-15%)" ~ "Enhanced observation",
    risk_group == "High risk (15-30%)" ~ "High-intensity monitoring",
    risk_group == "Very high risk (>=30%)" ~ "Critical-risk surveillance",
    TRUE ~ "Not available"
  )
}

trajectory_interpretation <- function(traj) {
  dplyr::case_when(
    traj == "C1" ~ "Low-dysfunction recovery phenotype",
    traj == "C2" ~ "Hemodynamic-respiratory stress phenotype",
    traj == "C3" ~ "Renal-metabolic dominant phenotype",
    traj == "C4" ~ "Rapidly progressive multiorgan dysfunction phenotype",
    traj == "C5" ~ "Persistent high-burden multiorgan dysfunction phenotype",
    TRUE ~ "Trajectory not available"
  )
}

## =========================================================
## 3. Prediction and trajectory functions
## =========================================================

apply_preprocessor_app <- function(prep, new_data) {
  if (is.null(prep$predictors)) stop("Invalid preprocessor: missing predictors.")

  x <- as.data.frame(new_data)
  missing_predictors <- setdiff(prep$predictors, names(x))
  if (length(missing_predictors) > 0) {
    for (v in missing_predictors) x[[v]] <- NA
  }
  x <- x[, prep$predictors, drop = FALSE]

  for (v in prep$numeric_vars) {
    if (!v %in% names(x)) x[[v]] <- NA_real_
    xx <- to_numeric_safe(x[[v]])
    xx[!is.finite(xx)] <- NA_real_
    if (v %in% prep$num_missing_indicator_vars) {
      x[[paste0(v, "_missing")]] <- as.integer(is.na(xx))
    }
    med <- prep$num_medians[[v]]
    if (is.null(med) || !is.finite(med)) med <- 0
    xx[is.na(xx)] <- med
    x[[v]] <- xx
  }

  for (v in prep$categorical_vars) {
    if (!v %in% names(x)) x[[v]] <- "Missing"
    xx <- as.character(x[[v]])
    xx[is.na(xx) | xx == ""] <- "Missing"
    known <- prep$cat_levels[[v]] %||% "Missing"
    known <- unique(c(as.character(known), "Missing"))
    xx[!xx %in% known] <- "Missing"
    x[[v]] <- factor(xx, levels = known)
  }

  extra_missing <- paste0(prep$num_missing_indicator_vars, "_missing")
  all_cols <- unique(c(prep$predictors, extra_missing))
  for (v in all_cols) {
    if (!v %in% names(x)) x[[v]] <- 0
  }

  for (v in prep$predictors) {
    if (!v %in% names(x)) next
    if (v %in% prep$numeric_vars) {
      xx <- to_numeric_safe(x[[v]])
      med <- prep$num_medians[[v]]
      if (is.null(med) || !is.finite(med)) med <- 0
      xx[!is.finite(xx) | is.na(xx)] <- med
      x[[v]] <- xx
    } else if (v %in% prep$categorical_vars) {
      xx <- as.character(x[[v]])
      xx[is.na(xx) | xx == ""] <- "Missing"
      known <- prep$cat_levels[[v]] %||% "Missing"
      known <- unique(c(as.character(known), "Missing"))
      xx[!xx %in% known] <- "Missing"
      x[[v]] <- factor(xx, levels = known)
    }
  }

  mm <- stats::model.matrix(prep$mm_formula, data = x, na.action = stats::na.pass)

  if (length(extra_missing) > 0) {
    miss_mat <- as.matrix(x[, extra_missing, drop = FALSE])
    colnames(miss_mat) <- clean_name(colnames(miss_mat))
    if (nrow(mm) != nrow(miss_mat)) {
      stop("Preprocessor row mismatch between model.matrix and missing indicators.")
    }
    mm <- cbind(mm, miss_mat)
  }

  colnames(mm) <- clean_name(colnames(mm))
  mm <- mm[, colSums(is.na(mm)) == 0, drop = FALSE]
  storage.mode(mm) <- "numeric"

  if (nrow(mm) != nrow(new_data)) {
    stop("Preprocessor row mismatch after model.matrix.")
  }

  list(matrix = mm, data = x)
}


## --- v12 robust elastic-net prediction helpers ----------------------------
## This block is intentionally independent of glmnet::coef() and predict() as
## the primary route. The deployed app had been failing for user-entered single
## patients with: "dim(X) must have a positive length". The stable route below
## reads the coefficient path directly from fit$glmnet.fit$beta/a0/lambda and
## computes plogis(intercept + X %*% beta). Direct glmnet prediction is used only
## as a secondary fallback with duplicated rows for one-row inputs.

as_one_row_numeric_matrix <- function(mat, nrow_expected = 1) {
  if (is.null(dim(mat))) {
    mat <- matrix(as.numeric(mat), nrow = nrow_expected)
  }
  mat <- as.matrix(mat)
  if (is.null(dim(mat))) mat <- matrix(as.numeric(mat), nrow = nrow_expected)
  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  mat
}

get_raw_glmnet_coefficients_v12 <- function(fit, s = "lambda.min") {
  gf <- fit$glmnet.fit %||% fit
  if (is.null(gf$beta) || is.null(gf$a0)) {
    stop("Saved elastic-net object does not contain raw beta/a0 coefficients.")
  }

  lambda_path <- suppressWarnings(as.numeric(gf$lambda %||% fit$lambda))
  beta_mat <- tryCatch(as.matrix(gf$beta), error = function(e) NULL)
  if (is.null(beta_mat)) stop("Cannot coerce saved beta coefficients to a matrix.")
  if (is.null(dim(beta_mat))) beta_mat <- matrix(as.numeric(beta_mat), ncol = 1)

  ## Select lambda index without calling coef.glmnet().
  if (length(lambda_path) == 0 || all(!is.finite(lambda_path))) {
    lambda_id <- if (ncol(beta_mat) >= 1) 1L else stop("Empty beta coefficient matrix.")
  } else {
    lambda_value <- suppressWarnings(
      if (is.numeric(s)) {
        as.numeric(s)[1]
      } else if (is.character(s) && length(s) == 1 && !is.null(fit[[s]])) {
        as.numeric(fit[[s]])[1]
      } else if (!is.null(fit$lambda.min)) {
        as.numeric(fit$lambda.min)[1]
      } else {
        lambda_path[1]
      }
    )
    if (!is.finite(lambda_value)) lambda_value <- lambda_path[1]
    lambda_id <- which.min(abs(lambda_path - lambda_value))
    if (length(lambda_id) == 0 || is.na(lambda_id) || lambda_id < 1) lambda_id <- 1L
    if (lambda_id > ncol(beta_mat)) lambda_id <- min(ncol(beta_mat), length(lambda_path))
    if (lambda_id < 1 || is.na(lambda_id)) lambda_id <- 1L
  }

  beta <- as.numeric(beta_mat[, lambda_id, drop = TRUE])
  beta_names <- rownames(beta_mat)
  if (is.null(beta_names) || length(beta_names) != length(beta)) {
    dn <- dimnames(gf$beta)
    if (!is.null(dn) && length(dn) >= 1 && length(dn[[1]]) == length(beta)) {
      beta_names <- dn[[1]]
    } else {
      stop("Saved elastic-net beta coefficients do not contain feature names; cannot align patient data safely.")
    }
  }
  names(beta) <- clean_name(beta_names)

  a0 <- suppressWarnings(as.numeric(gf$a0))
  intercept <- if (length(a0) >= lambda_id && is.finite(a0[lambda_id])) a0[lambda_id] else if (length(a0) >= 1 && is.finite(a0[1])) a0[1] else 0

  list(
    intercept = intercept,
    beta = beta,
    feature_names = names(beta),
    lambda_index = lambda_id
  )
}

align_to_glmnet_features <- function(mat, fit, s = "lambda.min") {
  mat <- as_one_row_numeric_matrix(mat, nrow_expected = 1)
  if (is.null(colnames(mat))) colnames(mat) <- paste0("V", seq_len(ncol(mat)))
  colnames(mat) <- clean_name(colnames(mat))

  coef_info <- get_raw_glmnet_coefficients_v12(fit, s = s)
  needed <- coef_info$feature_names

  out <- matrix(0, nrow = nrow(mat), ncol = length(needed))
  colnames(out) <- needed

  mat_names <- colnames(mat)
  mat_names_clean <- clean_name(mat_names)
  for (j in seq_along(needed)) {
    nm <- needed[j]
    hit <- match(nm, mat_names)
    if (is.na(hit)) hit <- match(clean_name(nm), mat_names_clean)
    if (!is.na(hit)) out[, j] <- mat[, hit]
  }
  out[!is.finite(out)] <- 0
  storage.mode(out) <- "numeric"
  out
}

predict_cv_glmnet_manual <- function(fit, mat, s = "lambda.min") {
  coef_info <- get_raw_glmnet_coefficients_v12(fit, s = s)
  beta <- coef_info$beta
  intercept <- coef_info$intercept

  mat <- align_to_glmnet_features(mat, fit, s = s)
  if (length(beta) == 0) {
    eta <- rep(intercept, nrow(mat))
  } else {
    beta <- beta[colnames(mat)]
    beta[is.na(beta)] <- 0
    eta <- as.numeric(intercept + mat %*% beta)
  }
  bound_prob(stats::plogis(eta))
}

predict_cv_glmnet_direct_fallback <- function(fit, mat, s = "lambda.min") {
  ## v13: direct glmnet fallback is intentionally disabled in deployment.
  ## The previous fallback could trigger "dim(X) must have a positive length"
  ## after demo interactions on shinyapps.io. Formal prediction now uses only
  ## raw coefficients from fit$glmnet.fit$beta/a0/lambda.
  stop("Direct glmnet fallback disabled; raw-coefficient prediction failed before fallback.")
}

predict_ami_edrs <- function(new_data, step05_object, model_name = NULL) {
  message("[AMI-EDRS v13] Formal prediction started.")
  if (is.null(step05_object$model_objects)) stop("The prediction engine is incomplete.")

  available <- names(step05_object$model_objects)
  if (is.null(model_name) || !model_name %in% available) {
    if ("Elastic net: dynamic + trajectory" %in% available) {
      model_name <- "Elastic net: dynamic + trajectory"
    } else {
      model_name <- step05_object$final_model_name %||% available[1]
    }
  }

  obj <- step05_object$model_objects[[model_name]]
  if (is.null(obj$object) || is.null(obj$prep)) {
    stop("Selected model object is incomplete. It must contain object and prep.")
  }

  pp <- apply_preprocessor_app(obj$prep, new_data)
  mat <- as_one_row_numeric_matrix(pp$matrix, nrow_expected = nrow(new_data))
  message("[AMI-EDRS v13] Preprocessed matrix: ", nrow(mat), " row(s), ", ncol(mat), " column(s).")

  if (!inherits(obj$object, "cv.glmnet")) {
    stop("This app currently supports the saved elastic-net AMI-EDRS model for direct prediction.")
  }

  pred <- tryCatch(
    predict_cv_glmnet_manual(obj$object, mat, s = "lambda.min"),
    error = function(e1) {
      message("[AMI-EDRS v13] Raw-coefficient prediction failed: ", conditionMessage(e1))
      stop(paste0("Raw-coefficient prediction failed: ", conditionMessage(e1)))
    }
  )

  pred <- bound_prob(pred)
  message("[AMI-EDRS v13] Formal prediction finished. Probability=", paste(round(pred, 6), collapse = ","))

  tibble::tibble(
    model = model_name,
    AMI_EDRS_probability = pred,
    AMI_EDRS_0_100 = round(pred * 100, 1),
    AMI_EDRS_risk_group = risk_group_from_prob(pred),
    AMI_EDRS_risk_group_public = risk_group_label_public(AMI_EDRS_risk_group),
    monitoring_intensity = monitoring_intensity(AMI_EDRS_risk_group)
  )
}

add_module_deltas <- function(df) {
  module_bases <- c(
    "hemodynamic",
    "respiratory_neurologic",
    "renal_metabolic",
    "inflammatory_coagulopathic",
    "cardiac_hypoperfusion"
  )

  for (m in module_bases) {
    aliases <- list(
      T0 = c(paste0(m, "_T0_score"), paste0(m, "_T0")),
      T24 = c(paste0(m, "_T24_score"), paste0(m, "_T24")),
      T48 = c(paste0(m, "_T48_score"), paste0(m, "_T48"))
    )

    for (tt in names(aliases)) {
      target <- paste0(m, "_", tt, "_score")
      if (!target %in% names(df)) {
        hit <- aliases[[tt]][aliases[[tt]] %in% names(df)][1]
        if (!is.na(hit)) df[[target]] <- df[[hit]]
      }
    }

    c0 <- paste0(m, "_T0_score")
    c24 <- paste0(m, "_T24_score")
    c48 <- paste0(m, "_T48_score")

    if (all(c(c0, c24, c48) %in% names(df))) {
      if (!paste0(m, "_delta_T48_T0") %in% names(df)) {
        df[[paste0(m, "_delta_T48_T0")]] <- to_numeric_safe(df[[c48]]) - to_numeric_safe(df[[c0]])
      }
      if (!paste0(m, "_delta_T48_T24") %in% names(df)) {
        df[[paste0(m, "_delta_T48_T24")]] <- to_numeric_safe(df[[c48]]) - to_numeric_safe(df[[c24]])
      }
      if (!paste0(m, "_slope_T0_T48") %in% names(df)) {
        df[[paste0(m, "_slope_T0_T48")]] <- (to_numeric_safe(df[[c48]]) - to_numeric_safe(df[[c0]])) / 48
      }
    }
  }
  df
}

trajectory_required_features <- function(step03_object) {
  if (is.null(step03_object)) return(character(0))
  prep <- step03_object$clustering_preprocessor
  features <- prep$kept_features %||% step03_object$main_cluster_features %||% character(0)
  unique(features)
}

trajectory_data_status <- function(df, step03_object = NULL) {
  if (is.null(step03_object)) {
    return(list(available = FALSE, message = "C1–C5 trajectory assignment requires complete T0, T24, and T48 organ-domain dynamics."))
  }

  features <- trajectory_required_features(step03_object)
  if (length(features) == 0) {
    return(list(available = FALSE, message = "C1–C5 trajectory assignment requires complete T0, T24, and T48 organ-domain dynamics."))
  }

  x <- add_module_deltas(as.data.frame(df))
  missing_cols <- setdiff(features, names(x))
  present_features <- intersect(features, names(x))
  missing_values <- present_features[purrr::map_lgl(present_features, ~ all(is.na(to_numeric_safe(x[[.x]]))))]

  missing_all <- unique(c(missing_cols, missing_values))

  if (length(missing_all) == 0) {
    return(list(available = TRUE, message = "Sufficient T0–T24–T48 organ-domain data are available for automatic C1–C5 assignment."))
  }

  missing_text <- paste(head(missing_all, 8), collapse = ", ")
  if (length(missing_all) > 8) missing_text <- paste0(missing_text, ", ...")

  need_window <- dplyr::case_when(
    any(grepl("T48", missing_all, ignore.case = TRUE)) ~ "Need T48 organ-domain data to assign C1–C5 trajectory phenotype.",
    any(grepl("T24", missing_all, ignore.case = TRUE)) ~ "Need T24 organ-domain data to assign C1–C5 trajectory phenotype.",
    any(grepl("T0", missing_all, ignore.case = TRUE)) ~ "Need T0 baseline organ-domain data to assign C1–C5 trajectory phenotype.",
    TRUE ~ "Need additional trajectory feature data to assign C1–C5 phenotype."
  )

  list(
    available = FALSE,
    message = paste0(need_window, " Missing/empty features include: ", missing_text)
  )
}

assign_trajectory_from_step03 <- function(new_data, step03_object, require_complete = TRUE) {
  if (is.null(step03_object)) return(rep(NA_character_, nrow(new_data)))

  if (require_complete) {
    status <- trajectory_data_status(new_data, step03_object)
    if (!isTRUE(status$available)) return(rep(NA_character_, nrow(new_data)))
  }

  prep <- step03_object$clustering_preprocessor
  km <- step03_object$final_kmeans
  cluster_map <- step03_object$cluster_map

  if (is.null(prep) || is.null(km) || is.null(cluster_map)) {
    return(rep(NA_character_, nrow(new_data)))
  }

  features <- trajectory_required_features(step03_object)
  features <- unique(features)
  if (length(features) == 0) return(rep(NA_character_, nrow(new_data)))

  x <- add_module_deltas(as.data.frame(new_data))
  for (v in setdiff(features, names(x))) x[[v]] <- NA_real_
  x <- x[, features, drop = FALSE]

  for (v in names(x)) {
    xx <- to_numeric_safe(x[[v]])
    med <- prep$medians[[v]]
    if (is.null(med) || !is.finite(med)) med <- 0
    xx[!is.finite(xx) | is.na(xx)] <- med
    x[[v]] <- xx
  }

  center <- prep$center[features]
  scalev <- prep$scale[features]
  scalev[!is.finite(scalev) | scalev == 0] <- 1

  xs <- sweep(as.matrix(x), 2, center, "-")
  xs <- sweep(xs, 2, scalev, "/")

  centers <- km$centers[, features, drop = FALSE]
  dmat <- sapply(seq_len(nrow(centers)), function(k) {
    rowSums((xs - matrix(centers[k, ], nrow = nrow(xs), ncol = ncol(xs), byrow = TRUE))^2)
  })
  raw_cluster <- apply(dmat, 1, which.min)

  map <- cluster_map %>%
    dplyr::mutate(cluster_raw = as.integer(cluster_raw), trajectory = as.character(trajectory))

  out <- map$trajectory[match(raw_cluster, map$cluster_raw)]
  out[is.na(out)] <- NA_character_
  out
}

assign_trajectory_batch_if_complete <- function(dat, step03_object) {
  if (is.null(step03_object)) {
    dat$trajectory_status <- "Trajectory engine unavailable; C1–C5 phenotype not assigned."
    if (!"trajectory" %in% names(dat)) dat$trajectory <- NA_character_
    return(dat)
  }

  features <- trajectory_required_features(step03_object)
  dat <- add_module_deltas(as.data.frame(dat))

  row_ok <- rep(TRUE, nrow(dat))
  missing_feature_note <- character(nrow(dat))
  for (i in seq_len(nrow(dat))) {
    row_df <- dat[i, , drop = FALSE]
    status <- trajectory_data_status(row_df, step03_object)
    row_ok[i] <- isTRUE(status$available)
    missing_feature_note[i] <- status$message
  }

  traj <- rep(NA_character_, nrow(dat))
  if (any(row_ok)) {
    traj[row_ok] <- assign_trajectory_from_step03(dat[row_ok, , drop = FALSE], step03_object, require_complete = TRUE)
  }

  dat$trajectory <- ifelse(row_ok, traj, NA_character_)
  dat$trajectory_status <- ifelse(row_ok, "C1–C5 assigned from T0–T24–T48 dynamics.", missing_feature_note)
  dat
}

make_single_patient_df <- function(input, predictors, trajectory_auto = NA_character_) {
  out <- as.data.frame(matrix(NA, nrow = 1, ncol = length(predictors)))
  names(out) <- predictors

  set_if_present <- function(var, val) {
    if (var %in% names(out)) out[[var]][1] <<- val
  }

  set_if_present("age", input$age)
  set_if_present("sex", input$sex)
  set_if_present("gender", input$sex)
  set_if_present("trajectory", trajectory_auto)

  binary_vars <- c(
    "hypertension", "diabetes", "ckd", "chronic_kidney_disease", "renal_disease",
    "copd", "chronic_pulmonary_disease", "heart_failure", "congestive_heart_failure",
    "atrial_fibrillation", "stroke", "cerebrovascular_disease",
    "liver_disease", "malignancy", "shock_proxy", "shock_hypoperfusion_proxy",
    "prior_cad_ihd", "ami_subtype_stemi"
  )

  for (v in binary_vars) {
    if (!v %in% names(out)) next
    value <- switch(
      v,
      "hypertension" = input$hypertension,
      "diabetes" = input$diabetes,
      "ckd" = input$ckd,
      "chronic_kidney_disease" = input$ckd,
      "renal_disease" = input$ckd,
      "copd" = input$copd,
      "chronic_pulmonary_disease" = input$copd,
      "heart_failure" = input$heart_failure,
      "congestive_heart_failure" = input$heart_failure,
      "atrial_fibrillation" = input$atrial_fibrillation,
      "stroke" = input$stroke,
      "cerebrovascular_disease" = input$stroke,
      "liver_disease" = input$liver_disease,
      "malignancy" = input$malignancy,
      "shock_proxy" = input$shock_proxy,
      "shock_hypoperfusion_proxy" = input$shock_proxy,
      "prior_cad_ihd" = input$prior_cad_ihd,
      "ami_subtype_stemi" = input$ami_subtype_stemi,
      0
    )
    out[[v]][1] <- as.integer(isTRUE(value))
  }

  raw_features <- make_raw_patient_features_from_ui(input)
  for (v in intersect(names(raw_features), names(out))) {
    out[[v]][1] <- raw_features[[v]][1]
  }

  advanced_map <- c(
    gcs_verbal_last_value = "gcs_verbal_last_value",
    bun_mean_value = "bun_mean_value",
    gcs_verbal_abnormal_ratio = "gcs_verbal_abnormal_ratio",
    spo2_mean_value = "spo2_mean_value",
    bun_max_value = "bun_max_value",
    ptt_min_value = "ptt_min_value",
    gcs_eye_slope_per_hour = "gcs_eye_slope_per_hour",
    hemoglobin_sd_value = "hemoglobin_sd_value",
    lactate_slope_per_hour = "lactate_slope_per_hour",
    heart_rate_abnormal_ratio = "heart_rate_abnormal_ratio",
    gcs_verbal_slope_per_hour = "gcs_verbal_slope_per_hour",
    chloride_max_value = "chloride_max_value",
    bun_last_value = "bun_last_value",
    hemoglobin_mean_value = "hemoglobin_mean_value",
    bicarbonate_abnormal_ratio = "bicarbonate_abnormal_ratio",
    hematocrit_sd_value = "hematocrit_sd_value",
    gcs_motor_slope_per_hour = "gcs_motor_slope_per_hour",
    potassium_abnormal_ratio = "potassium_abnormal_ratio"
  )

  for (nm in names(advanced_map)) {
    var <- unname(advanced_map[[nm]])
    if (var %in% names(out) && !is.null(input[[nm]]) && !is.na(to_numeric_safe(input[[nm]]))) out[[var]][1] <- input[[nm]]
  }

  add_module_deltas(out)
}

make_trajectory_input_df_from_ui <- function(input) {
  make_raw_patient_features_from_ui(input)
}

module_score_cols <- function() {
  modules <- c(
    "hemodynamic",
    "respiratory_neurologic",
    "renal_metabolic",
    "inflammatory_coagulopathic",
    "cardiac_hypoperfusion"
  )
  paste0(as.vector(outer(modules, c("T0", "T24", "T48"), paste, sep = "_")), "_score")
}

module_dynamics_complete <- function(df) {
  cols <- module_score_cols()
  x <- add_module_deltas(as.data.frame(df))
  missing_cols <- setdiff(cols, names(x))
  present <- intersect(cols, names(x))
  missing_values <- present[purrr::map_lgl(present, ~ all(is.na(to_numeric_safe(x[[.x]]))))]
  missing_all <- unique(c(missing_cols, missing_values))
  if (length(missing_all) == 0) {
    list(available = TRUE, message = "Complete T0–T24–T48 organ-domain dynamics are available for C1–C5 assignment.")
  } else {
    need_window <- dplyr::case_when(
      any(grepl("T48", missing_all, ignore.case = TRUE)) ~ "C1–C5 trajectory phenotype requires T48 organ-domain score inputs.",
      any(grepl("T24", missing_all, ignore.case = TRUE)) ~ "C1–C5 trajectory phenotype requires T24 organ-domain score inputs.",
      any(grepl("T0", missing_all, ignore.case = TRUE)) ~ "C1–C5 trajectory phenotype requires T0 organ-domain score inputs.",
      TRUE ~ "C1–C5 trajectory phenotype requires complete T0, T24, and T48 organ-domain score inputs."
    )
    list(available = FALSE, message = need_window)
  }
}

assign_trajectory_from_modules_simple <- function(new_data) {
  status <- module_dynamics_complete(new_data)
  if (!isTRUE(status$available)) return(rep(NA_character_, nrow(new_data)))

  x <- add_module_deltas(as.data.frame(new_data))
  n <- nrow(x)
  out <- rep(NA_character_, n)

  getv <- function(row, name) to_numeric_safe(x[[name]][row])

  for (i in seq_len(n)) {
    t0 <- c(
      getv(i, "hemodynamic_T0_score"),
      getv(i, "respiratory_neurologic_T0_score"),
      getv(i, "renal_metabolic_T0_score"),
      getv(i, "inflammatory_coagulopathic_T0_score"),
      getv(i, "cardiac_hypoperfusion_T0_score")
    )
    t24 <- c(
      getv(i, "hemodynamic_T24_score"),
      getv(i, "respiratory_neurologic_T24_score"),
      getv(i, "renal_metabolic_T24_score"),
      getv(i, "inflammatory_coagulopathic_T24_score"),
      getv(i, "cardiac_hypoperfusion_T24_score")
    )
    t48 <- c(
      getv(i, "hemodynamic_T48_score"),
      getv(i, "respiratory_neurologic_T48_score"),
      getv(i, "renal_metabolic_T48_score"),
      getv(i, "inflammatory_coagulopathic_T48_score"),
      getv(i, "cardiac_hypoperfusion_T48_score")
    )
    mean0 <- mean(t0, na.rm = TRUE)
    mean24 <- mean(t24, na.rm = TRUE)
    mean48 <- mean(t48, na.rm = TRUE)
    delta48 <- mean48 - mean0
    renal48 <- t48[3]
    hemo48 <- t48[1]
    resp48 <- t48[2]

    out[i] <- dplyr::case_when(
      mean48 < 35 && delta48 <= 8 ~ "C1",
      mean48 >= 70 ~ "C5",
      delta48 >= 15 && (hemo48 >= 60 || resp48 >= 60 || mean48 >= 55) ~ "C4",
      renal48 >= 55 && renal48 >= max(t48[-3], na.rm = TRUE) - 5 ~ "C3",
      TRUE ~ "C2"
    )
  }
  out
}

auto_trajectory_status <- function(df, step03_object = NULL) {
  ## Prefer the saved trajectory engine when its required features are complete.
  if (!is.null(step03_object)) {
    st <- trajectory_data_status(df, step03_object)
    if (isTRUE(st$available)) {
      return(list(available = TRUE, method = "trajectory_engine", message = "C1–C5 phenotype can be assigned from complete T0–T24–T48 dynamics."))
    }
  }

  ## Fallback for the user-facing app: assign phenotype from complete organ-domain dynamics.
  st_module <- module_dynamics_complete(df)
  if (isTRUE(st_module$available)) {
    return(list(available = TRUE, method = "organ_domain_rule", message = "C1–C5 phenotype can be assigned from complete T0–T24–T48 organ-domain dynamics."))
  }
  list(available = FALSE, method = "unavailable", message = st_module$message)
}

assign_trajectory_auto <- function(new_data, step03_object = NULL) {
  st <- auto_trajectory_status(new_data, step03_object)
  if (!isTRUE(st$available)) return(rep(NA_character_, nrow(new_data)))
  if (identical(st$method, "trajectory_engine") && !is.null(step03_object)) {
    out <- assign_trajectory_from_step03(new_data, step03_object, require_complete = TRUE)
    if (any(!is.na(out))) return(out)
  }
  assign_trajectory_from_modules_simple(new_data)
}

read_uploaded_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv", "txt")) {
    readr::read_csv(path, show_col_types = FALSE)
  } else if (ext %in% c("xlsx", "xls")) {
    readxl::read_excel(path)
  } else {
    stop("Unsupported file type. Please upload CSV or Excel file.")
  }
}

make_example_data <- function(predictors = public_feature_columns()) {
  predictors <- public_feature_columns()
  template <- as.data.frame(matrix(NA, nrow = 2, ncol = length(predictors)))
  names(template) <- predictors

  setcol <- function(col, values) {
    if (col %in% names(template)) template[[col]] <<- values
  }

  setcol("patient_id", c("demo_low_risk", "demo_high_risk"))
  setcol("age", c(62, 82))
  setcol("sex", c("Male", "Female"))
  setcol("gender", c("Male", "Female"))

  binary_low_high <- list(
    hypertension = c(1, 1), diabetes = c(0, 1), ckd = c(0, 1), chronic_kidney_disease = c(0, 1),
    copd = c(0, 1), malignancy = c(0, 1), liver_disease = c(0, 1), shock_proxy = c(0, 1),
    shock_hypoperfusion_proxy = c(0, 1), stroke = c(0, 1), prior_cad_ihd = c(1, 1),
    heart_failure = c(0, 1), atrial_fibrillation = c(0, 1), ami_subtype_stemi = c(0, 0)
  )
  for (nm in names(binary_low_high)) setcol(nm, binary_low_high[[nm]])

  low <- list(
    T0 = c(heart_rate=82, sbp=128, map=86, resp_rate=16, spo2=98, temperature=36.8, gcs_eye=4, gcs_verbal=5, gcs_motor=6,
           creatinine=0.9, bun=16, sodium=140, potassium=4.1, chloride=103, bicarbonate=24, glucose=118,
           wbc=8.2, hemoglobin=13.8, hematocrit=41, platelet=230, inr=1.0, ptt=29, lactate=1.1, troponin_t=0.25),
    T24 = c(heart_rate=78, sbp=132, map=88, resp_rate=15, spo2=98, temperature=36.9, gcs_eye=4, gcs_verbal=5, gcs_motor=6,
            creatinine=0.9, bun=15, sodium=140, potassium=4.0, chloride=103, bicarbonate=25, glucose=112,
            wbc=7.5, hemoglobin=13.6, hematocrit=40, platelet=225, inr=1.0, ptt=30, lactate=1.0, troponin_t=0.18),
    T48 = c(heart_rate=74, sbp=134, map=90, resp_rate=15, spo2=98, temperature=36.8, gcs_eye=4, gcs_verbal=5, gcs_motor=6,
            creatinine=0.8, bun=14, sodium=140, potassium=4.0, chloride=102, bicarbonate=25, glucose=105,
            wbc=7.0, hemoglobin=13.5, hematocrit=40, platelet=222, inr=1.0, ptt=29, lactate=0.9, troponin_t=0.12)
  )
  high <- list(
    T0 = c(heart_rate=108, sbp=96, map=67, resp_rate=24, spo2=92, temperature=37.8, gcs_eye=4, gcs_verbal=4, gcs_motor=6,
           creatinine=1.8, bun=38, sodium=134, potassium=5.1, chloride=108, bicarbonate=20, glucose=210,
           wbc=15, hemoglobin=10.5, hematocrit=31, platelet=140, inr=1.4, ptt=42, lactate=2.8, troponin_t=2.5),
    T24 = c(heart_rate=122, sbp=86, map=59, resp_rate=30, spo2=88, temperature=38.2, gcs_eye=3, gcs_verbal=3, gcs_motor=5,
            creatinine=2.6, bun=55, sodium=132, potassium=5.5, chloride=111, bicarbonate=17, glucose=260,
            wbc=21, hemoglobin=9.5, hematocrit=28, platelet=105, inr=1.9, ptt=60, lactate=5.2, troponin_t=5.5),
    T48 = c(heart_rate=132, sbp=78, map=52, resp_rate=34, spo2=84, temperature=38.6, gcs_eye=2, gcs_verbal=2, gcs_motor=4,
            creatinine=3.2, bun=72, sodium=130, potassium=5.9, chloride=114, bicarbonate=13, glucose=310,
            wbc=28, hemoglobin=8.8, hematocrit=25, platelet=78, inr=2.5, ptt=85, lactate=7.5, troponin_t=9.0)
  )

  for (tt in timepoints_app) {
    for (v in raw_vars) {
      setcol(raw_id(tt, v), c(low[[tt]][[v]], high[[tt]][[v]]))
    }
  }

  template
}

public_feature_columns <- function() {
  out <- c(
    "patient_id", "age", "sex", "gender",
    "hypertension", "diabetes", "ckd", "chronic_kidney_disease", "copd",
    "heart_failure", "atrial_fibrillation", "stroke", "liver_disease", "malignancy",
    "shock_proxy", "shock_hypoperfusion_proxy", "prior_cad_ihd", "ami_subtype_stemi",
    raw_public_columns()
  )
  unique(out)
}

make_public_example_data <- function() {
  make_example_data(public_feature_columns())
}

fallback_demo_prediction <- function(input, trajectory_auto = NA_character_) {
  modules <- c("hemodynamic", "respiratory_neurologic", "renal_metabolic", "inflammatory_coagulopathic", "cardiac_hypoperfusion")
  raw_features <- make_raw_patient_features_from_ui(input)
  t48 <- sapply(modules, function(m) to_numeric_safe(raw_features[[paste0(m, "_T48_score")]]))
  mean48 <- mean(t48, na.rm = TRUE)
  if (!is.finite(mean48)) mean48 <- 45
  age_term <- (to_numeric_safe(input$age) - 65) / 10
  comorb <- sum(c(
    isTRUE(input$diabetes), isTRUE(input$ckd), isTRUE(input$copd), isTRUE(input$malignancy),
    isTRUE(input$liver_disease), isTRUE(input$stroke), isTRUE(input$shock_proxy)
  ), na.rm = TRUE)
  traj_bonus <- dplyr::case_when(
    trajectory_auto == "C1" ~ -1.1,
    trajectory_auto == "C2" ~ -0.4,
    trajectory_auto == "C3" ~ 0.1,
    trajectory_auto == "C4" ~ 0.85,
    trajectory_auto == "C5" ~ 1.15,
    TRUE ~ 0
  )
  lp <- -4.2 + 0.62 * age_term + 0.045 * mean48 + 0.28 * comorb + traj_bonus
  p <- bound_prob(plogis(lp))
  tibble::tibble(
    model = "Demonstration calculation",
    AMI_EDRS_probability = p,
    AMI_EDRS_0_100 = round(p * 100, 1),
    AMI_EDRS_risk_group = risk_group_from_prob(p),
    AMI_EDRS_risk_group_public = risk_group_label_public(AMI_EDRS_risk_group),
    monitoring_intensity = monitoring_intensity(AMI_EDRS_risk_group)
  )
}

## Fallback calculator that works from a data row rather than reactive input.
## This is used for demo buttons on shinyapps.io, where UI updates are asynchronous.
fallback_demo_prediction_df <- function(dat, trajectory_auto = NA_character_) {
  dat <- as.data.frame(dat)
  modules <- c("hemodynamic", "respiratory_neurologic", "renal_metabolic", "inflammatory_coagulopathic", "cardiac_hypoperfusion")
  cell <- function(nm) {
    if (!nm %in% names(dat)) return(NA_real_)
    val <- to_numeric_safe(dat[[nm]][1])
    if (length(val) == 0 || !is.finite(val)) NA_real_ else as.numeric(val)
  }
  t48 <- sapply(modules, function(m) cell(paste0(m, "_T48_score")))
  mean48 <- mean(t48, na.rm = TRUE)
  if (!is.finite(mean48)) mean48 <- 45
  age_val <- cell("age")
  if (!is.finite(age_val)) age_val <- 70
  age_term <- (age_val - 65) / 10
  comorb_vars <- c("diabetes", "ckd", "chronic_kidney_disease", "copd", "malignancy", "liver_disease", "stroke", "shock_proxy", "shock_hypoperfusion_proxy")
  comorb <- sum(sapply(comorb_vars, function(v) {
    if (!v %in% names(dat)) return(0)
    as.integer(isTRUE(as.logical(as.integer(dat[[v]][1]))))
  }), na.rm = TRUE)
  traj_bonus <- dplyr::case_when(
    trajectory_auto == "C1" ~ -1.1,
    trajectory_auto == "C2" ~ -0.4,
    trajectory_auto == "C3" ~ 0.1,
    trajectory_auto == "C4" ~ 0.85,
    trajectory_auto == "C5" ~ 1.15,
    TRUE ~ 0
  )
  lp <- -4.2 + 0.62 * age_term + 0.045 * mean48 + 0.28 * comorb + traj_bonus
  p <- bound_prob(plogis(lp))
  rg <- risk_group_from_prob(p)
  tibble::tibble(
    model = "Demonstration calculation",
    AMI_EDRS_probability = p,
    AMI_EDRS_0_100 = round(p * 100, 1),
    AMI_EDRS_risk_group = rg,
    AMI_EDRS_risk_group_public = risk_group_label_public(rg),
    monitoring_intensity = monitoring_intensity(rg)
  )
}

## =========================================================
## 4. UI helpers and plotting functions
## =========================================================

pal_risk <- c(
  "Low risk (<5%)" = "#2A9D8F",
  "Intermediate risk (5-15%)" = "#4A6FA5",
  "High risk (15-30%)" = "#C99A2E",
  "Very high risk (>=30%)" = "#B23A48"
)

pal_modules <- c(
  "Hemodynamic" = "#3A6EA5",
  "Respiratory–neurologic" = "#C99A2E",
  "Renal–metabolic" = "#2A9D8F",
  "Inflammatory–coagulopathic" = "#7D5BA6",
  "Cardiac injury / hypoperfusion" = "#B23A48"
)

risk_color <- function(risk_group) {
  pal_risk[[risk_group]] %||% "#9CA3AF"
}

clinical_css <- "
:root {
  --ami-ink: #EAF0F6;
  --ami-navy: #0B1F33;
  --ami-navy-2: #102A43;
  --ami-steel: #1F3A5F;
  --ami-teal: #2A9D8F;
  --ami-gold: #C99A2E;
  --ami-red: #B23A48;
  --ami-bg: #0F172A;
  --ami-card: rgba(255,255,255,0.955);
  --ami-border: rgba(31,58,95,0.15);
  --ami-muted: #667085;
}
body {
  background:
    radial-gradient(circle at 10% 5%, rgba(42,157,143,0.16), transparent 28%),
    radial-gradient(circle at 92% 8%, rgba(201,154,46,0.13), transparent 30%),
    linear-gradient(180deg, #0B1F33 0%, #102A43 26%, #EEF3F7 26%, #F7F9FC 100%);
  font-family: Inter, -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
  color: #162033;
}
.navbar {
  background: rgba(10,31,51,0.96) !important;
  backdrop-filter: blur(12px);
  border-bottom: 1px solid rgba(201,154,46,0.30);
  box-shadow: 0 12px 32px rgba(2,8,23,0.35);
}
.navbar-brand, .nav-link { color: rgba(255,255,255,0.90) !important; }
.nav-tabs .nav-link { font-weight: 850; border-radius: 999px !important; margin-right: 4px; }
.nav-tabs .nav-link.active {
  color: #0B1F33 !important;
  background: linear-gradient(135deg, #E9C46A, #F7E7B0) !important;
  border-color: transparent !important;
  box-shadow: 0 8px 18px rgba(201,154,46,0.30);
}
.ami-hero {
  background:
    linear-gradient(135deg, rgba(11,31,51,0.98) 0%, rgba(31,58,95,0.96) 54%, rgba(42,157,143,0.86) 100%);
  color: white;
  border-radius: 26px;
  padding: 34px 36px;
  margin: 20px 0 24px 0;
  box-shadow: 0 24px 54px rgba(2,8,23,0.34);
  position: relative;
  overflow: hidden;
  border: 1px solid rgba(233,196,106,0.28);
}
.ami-hero:before {
  content: '';
  position: absolute;
  width: 390px;
  height: 390px;
  right: -150px;
  top: -180px;
  border-radius: 999px;
  background: radial-gradient(circle, rgba(233,196,106,0.32), transparent 66%);
}
.ami-hero:after {
  content: '';
  position: absolute;
  width: 180px;
  height: 180px;
  left: -60px;
  bottom: -70px;
  border-radius: 999px;
  background: rgba(255,255,255,0.08);
}
.ami-hero h2 { font-weight: 950; margin-bottom: 8px; letter-spacing: -0.035em; position: relative; z-index: 2; }
.ami-hero p { color: rgba(255,255,255,0.82); margin-bottom: 0; font-size: 15px; max-width: 900px; position: relative; z-index: 2; }
.ami-card {
  background: var(--ami-card);
  border: 1px solid var(--ami-border);
  border-radius: 22px;
  padding: 20px;
  margin-bottom: 18px;
  box-shadow: 0 18px 38px rgba(15,23,42,0.11);
}
.ami-card h4 { font-weight: 930; color: #102A43; margin-top: 0; margin-bottom: 12px; letter-spacing: -0.01em; }
.section-kicker { color: #2A6F75; font-weight: 900; font-size: 11px; text-transform: uppercase; letter-spacing: .14em; margin-bottom: 5px; }
.value-card {
  background: linear-gradient(180deg, #FFFFFF, #F5F7FA);
  border: 1px solid rgba(31,58,95,0.15);
  border-top: 5px solid #1F3A5F;
  border-radius: 22px;
  padding: 20px;
  min-height: 158px;
  box-shadow: 0 20px 42px rgba(15,23,42,0.13);
  transition: transform .15s ease, box-shadow .15s ease;
}
.value-card:hover { transform: translateY(-2px); box-shadow: 0 24px 48px rgba(15,23,42,0.18); }
.value-title { color: #64748B; font-weight: 900; font-size: 11px; text-transform: uppercase; letter-spacing: .095em; }
.value-number { font-size: 40px; font-weight: 980; margin-top: 8px; line-height: 1.02; letter-spacing: -0.04em; }
.value-subtitle { color: #475569; font-weight: 700; margin-top: 10px; line-height: 1.35; }
.risk-pill, .trajectory-pill {
  display: inline-block;
  color: white;
  font-weight: 900;
  padding: 9px 14px;
  border-radius: 999px;
  margin-top: 10px;
  box-shadow: 0 10px 22px rgba(0,0,0,0.16);
}
.status-pill {
  display: inline-block;
  padding: 8px 12px;
  border-radius: 999px;
  font-weight: 900;
  margin-bottom: 8px;
  font-size: 12px;
}
.status-good { background: rgba(42,157,143,0.12); color: #1C6B62; border: 1px solid rgba(42,157,143,0.28); }
.status-warn { background: rgba(201,154,46,0.13); color: #7A5A13; border: 1px solid rgba(201,154,46,0.30); }
.status-bad { background: rgba(178,58,72,0.12); color: #8F2433; border: 1px solid rgba(178,58,72,0.28); }
.small-note, .compact-label { color: var(--ami-muted); font-size: 12px; line-height: 1.38; }
.input-guide {
  background: linear-gradient(135deg, #F8FAFC, #EEF3F7);
  border: 1px solid rgba(31,58,95,0.14);
  border-left: 6px solid #1F3A5F;
  border-radius: 18px;
  padding: 16px 18px;
  margin-bottom: 16px;
  box-shadow: 0 12px 28px rgba(15,23,42,0.07);
}
.disclaimer {
  background: rgba(255,248,232,0.96);
  border: 1px solid rgba(201,154,46,0.30);
  border-radius: 18px;
  padding: 14px;
  color: #6E520E;
  font-size: 12px;
  line-height: 1.45;
}
.btn-primary {
  background: linear-gradient(135deg, #102A43, #1F3A5F) !important;
  border: none !important;
  box-shadow: 0 12px 24px rgba(16,42,67,0.28) !important;
  font-weight: 900 !important;
  border-radius: 14px !important;
}
.btn-primary:hover { background: linear-gradient(135deg, #0B1F33, #2A9D8F) !important; }
.btn-outline-primary {
  color: #102A43 !important;
  border-color: rgba(31,58,95,0.34) !important;
  background: rgba(255,255,255,0.78) !important;
  font-weight: 850 !important;
  border-radius: 14px !important;
}
.btn-outline-primary:hover {
  color: white !important;
  background: linear-gradient(135deg, #1F3A5F, #2A9D8F) !important;
}
.btn-outline-danger {
  color: #8F2433 !important;
  border-color: rgba(178,58,72,0.40) !important;
  background: rgba(255,255,255,0.78) !important;
  font-weight: 850 !important;
  border-radius: 14px !important;
}
.btn-outline-danger:hover {
  color: white !important;
  background: linear-gradient(135deg, #8F2433, #C99A2E) !important;
}
.form-control, .form-select {
  border-radius: 14px !important;
  border-color: rgba(31,58,95,0.22) !important;
}
.form-control:focus, .form-select:focus {
  border-color: #2A9D8F !important;
  box-shadow: 0 0 0 .2rem rgba(42,157,143,0.16) !important;
}
.download-button, .btn-default {
  border-radius: 14px !important;
  font-weight: 850 !important;
}
/* v6 executive color refinement: deep, sober, high-contrast */
body {
  background:
    radial-gradient(circle at 14% 8%, rgba(91,141,178,0.18), transparent 28%),
    radial-gradient(circle at 82% 6%, rgba(201,154,46,0.18), transparent 26%),
    linear-gradient(180deg, #071827 0%, #0E263B 25%, #F3F6FA 25%, #F8FAFC 100%) !important;
}
.navbar {
  background: linear-gradient(90deg, #071827, #0E263B 58%, #163B52) !important;
  border-bottom: 2px solid rgba(214,169,73,0.70) !important;
}
.navbar-brand {
  color: #F7E7B0 !important;
  font-size: 22px !important;
  text-shadow: 0 1px 8px rgba(0,0,0,0.28);
}
.navbar .nav-link {
  color: rgba(246,250,253,0.92) !important;
  background: rgba(255,255,255,0.08) !important;
  border: 1px solid rgba(255,255,255,0.16) !important;
  border-radius: 999px !important;
  margin: 7px 5px !important;
  padding: 8px 15px !important;
  font-weight: 850 !important;
}
.navbar .nav-link.active, .navbar .nav-link:hover {
  color: #071827 !important;
  background: linear-gradient(135deg, #F0C667, #FFF1C2) !important;
  border-color: rgba(240,198,103,0.95) !important;
}
/* Inner 1-4 workflow tabs: make them clearly visible */
.tabbable > .nav-tabs {
  border-bottom: none !important;
  margin-bottom: 18px !important;
  gap: 8px;
}
.tabbable > .nav-tabs .nav-link {
  color: #0E263B !important;
  background: linear-gradient(180deg, #FFFFFF, #EEF4FA) !important;
  border: 1px solid rgba(14,38,59,0.20) !important;
  border-radius: 16px !important;
  margin: 4px 4px 8px 0 !important;
  padding: 11px 16px !important;
  font-weight: 900 !important;
  box-shadow: 0 8px 18px rgba(14,38,59,0.10);
}
.tabbable > .nav-tabs .nav-link.active {
  color: #FFFFFF !important;
  background: linear-gradient(135deg, #12324A, #1E5873) !important;
  border-color: rgba(30,88,115,0.90) !important;
  box-shadow: 0 12px 26px rgba(14,38,59,0.24);
}
.ami-hero {
  background: linear-gradient(135deg, #071827 0%, #12324A 54%, #1E5873 100%) !important;
  border: 1px solid rgba(240,198,103,0.42) !important;
}
.ami-card, .value-card {
  border-color: rgba(14,38,59,0.15) !important;
  box-shadow: 0 16px 38px rgba(7,24,39,0.12) !important;
}
.btn-primary {
  background: linear-gradient(135deg, #0B2133, #1E5873) !important;
  color: #FFFFFF !important;
}
.btn-outline-primary {
  color: #0B2133 !important;
  border-color: rgba(14,38,59,0.34) !important;
}
.btn-outline-primary:hover {
  color: #FFFFFF !important;
  background: linear-gradient(135deg, #12324A, #1E5873) !important;
}
.download-button, .btn-default {
  background: #FFFFFF !important;
  color: #0E263B !important;
  border: 1px solid rgba(14,38,59,0.22) !important;
}
.download-button:hover, .btn-default:hover {
  background: #F0C667 !important;
  color: #071827 !important;
}

/* v7 visual refinement: one-upload public interface and stronger executive contrast */
.ami-card h4 { color: #0A2236 !important; }
.tabbable > .nav-tabs .nav-link {
  color: #071827 !important;
  background: linear-gradient(180deg, #FDFBF5, #E8EEF5) !important;
  border: 2px solid rgba(14,38,59,0.34) !important;
  font-size: 15px !important;
  letter-spacing: .01em;
}
.tabbable > .nav-tabs .nav-link.active {
  color: #FFFFFF !important;
  background: linear-gradient(135deg, #071827, #174866) !important;
  border-color: #D6A949 !important;
  box-shadow: 0 14px 30px rgba(7,24,39,0.30) !important;
}
.btn-primary, .btn-outline-primary, .btn-outline-danger, .download-button, .btn-default {
  min-height: 42px !important;
  padding: 9px 14px !important;
  font-size: 14px !important;
}
.download-button, .btn-default {
  background: linear-gradient(135deg, #FFFFFF, #F5F0DF) !important;
  color: #071827 !important;
  border: 2px solid rgba(214,169,73,0.58) !important;
  box-shadow: 0 8px 18px rgba(7,24,39,0.09) !important;
}
.download-button:hover, .btn-default:hover {
  background: linear-gradient(135deg, #D6A949, #FFF1C2) !important;
  color: #071827 !important;
  border-color: #D6A949 !important;
}

"

clinical_css <- paste0(clinical_css, "
/* v8 layout refinement: stronger executive layout, coordinated sidebar, larger equal workflow tabs */
.container-fluid { max-width: 1780px; }
.ami-hero {
  margin-top: 18px !important;
  margin-bottom: 30px !important;
  border-radius: 30px !important;
  padding: 36px 42px !important;
}
.ami-hero h2 { font-size: 42px !important; line-height: 1.05 !important; }
.ami-hero p { font-size: 16px !important; }
.quick-panel {
  background: linear-gradient(180deg, #FFFFFF 0%, #F7FAFC 100%) !important;
  border-top: 6px solid #1E5873 !important;
  padding: 22px 22px 24px 22px !important;
  min-height: 270px !important;
}
.quick-panel h4 {
  font-size: 22px !important;
  letter-spacing: -0.02em !important;
  margin-bottom: 16px !important;
  padding-bottom: 12px !important;
  border-bottom: 1px solid rgba(14,38,59,0.12) !important;
}
.quick-panel .status-pill {
  display: block !important;
  width: 100% !important;
  text-align: center !important;
  padding: 10px 14px !important;
  margin: 8px 0 !important;
  border-radius: 14px !important;
  font-size: 13px !important;
}
.quick-panel .download-button,
.quick-panel .btn-default {
  width: 100% !important;
  padding: 12px 14px !important;
  margin-top: 6px !important;
  background: linear-gradient(135deg, #F8FBFD, #FFFFFF) !important;
  border: 2px solid rgba(30,88,115,0.24) !important;
  color: #0E263B !important;
  box-shadow: 0 10px 22px rgba(14,38,59,0.08) !important;
}
.quick-panel .download-button:hover,
.quick-panel .btn-default:hover {
  background: linear-gradient(135deg, #E6F2F4, #FFF6D7) !important;
  border-color: #D6A949 !important;
}
.disclaimer {
  border-radius: 20px !important;
  border-left: 6px solid #D6A949 !important;
  box-shadow: 0 12px 28px rgba(126,91,18,0.10) !important;
}
/* Equal-width workflow tabs */
.tabbable > .nav-tabs {
  display: flex !important;
  width: 100% !important;
  gap: 14px !important;
  margin-bottom: 28px !important;
  align-items: stretch !important;
}
.tabbable > .nav-tabs > li,
.tabbable > .nav-tabs > .nav-item {
  flex: 1 1 0 !important;
  min-width: 0 !important;
}
.tabbable > .nav-tabs .nav-link {
  width: 100% !important;
  min-height: 64px !important;
  display: flex !important;
  align-items: center !important;
  justify-content: center !important;
  text-align: center !important;
  padding: 17px 18px !important;
  font-size: 17px !important;
  line-height: 1.15 !important;
  border-radius: 20px !important;
  letter-spacing: .02em !important;
  border: 2px solid rgba(14,38,59,0.26) !important;
  background: linear-gradient(180deg, #FFFFFF, #EDF3F8) !important;
  box-shadow: 0 14px 30px rgba(7,24,39,0.11) !important;
}
.tabbable > .nav-tabs .nav-link.active {
  background: linear-gradient(135deg, #0B2133, #1E5873) !important;
  color: #FFFFFF !important;
  border-color: #D6A949 !important;
  box-shadow: 0 18px 36px rgba(7,24,39,0.28), inset 0 0 0 1px rgba(255,255,255,0.10) !important;
}
/* Results action buttons */
.action-row {
  display: flex !important;
  gap: 18px !important;
  margin: 8px 0 28px 0 !important;
}
.action-row .action-cell { flex: 1 1 0 !important; }
.action-row .btn {
  width: 100% !important;
  min-height: 58px !important;
  border-radius: 18px !important;
  font-size: 16px !important;
  font-weight: 950 !important;
  letter-spacing: .01em !important;
}
.btn-run-main {
  color: #FFFFFF !important;
  background: linear-gradient(135deg, #082033 0%, #0F4C5C 56%, #1E7A72 100%) !important;
  border: 2px solid rgba(214,169,73,0.72) !important;
  box-shadow: 0 18px 36px rgba(8,32,51,0.30) !important;
}
.btn-run-main:hover {
  transform: translateY(-1px);
  background: linear-gradient(135deg, #071827 0%, #1E5873 58%, #D6A949 100%) !important;
}
.btn-demo-low {
  color: #0B3B38 !important;
  background: linear-gradient(135deg, #E6F7F4, #F7FFFD) !important;
  border: 2px solid rgba(42,157,143,0.55) !important;
  box-shadow: 0 14px 28px rgba(42,157,143,0.16) !important;
}
.btn-demo-low:hover {
  color: #FFFFFF !important;
  background: linear-gradient(135deg, #16776F, #2A9D8F) !important;
}
.btn-demo-high {
  color: #7A1D2B !important;
  background: linear-gradient(135deg, #FFF0F2, #FFF9E9) !important;
  border: 2px solid rgba(178,58,72,0.56) !important;
  box-shadow: 0 14px 28px rgba(178,58,72,0.15) !important;
}
.btn-demo-high:hover {
  color: #FFFFFF !important;
  background: linear-gradient(135deg, #9F2336, #B23A48) !important;
}
@media (max-width: 992px) {
  .tabbable > .nav-tabs { flex-direction: column !important; gap: 8px !important; }
  .action-row { flex-direction: column !important; }
}
.single-summary-table {
  width: 100%;
  border-collapse: separate;
  border-spacing: 0;
  overflow: hidden;
  border-radius: 16px;
  border: 1px solid rgba(14,38,59,0.14);
  font-size: 14px;
}
.single-summary-table th {
  background: linear-gradient(135deg, #0B2133, #1E5873);
  color: #FFFFFF;
  padding: 12px 10px;
  text-align: left;
  white-space: nowrap;
}

.simple-result-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;
  background: #ffffff;
}
.simple-result-table th {
  background: #EFF4F8;
  color: #102A43;
  font-weight: 800;
  padding: 9px 10px;
  border-bottom: 1px solid #D9E2EC;
  white-space: nowrap;
}
.simple-result-table td {
  padding: 8px 10px;
  border-bottom: 1px solid #EDF2F7;
  white-space: nowrap;
}
.simple-result-scroll {
  overflow-x: auto;
  max-width: 100%;
}

.single-summary-table td {
  background: #FFFFFF;
  color: #11283A;
  padding: 12px 10px;
  border-top: 1px solid rgba(14,38,59,0.08);
  vertical-align: top;
}
")


unit_help <- function(text) tags$div(class = "compact-label", text)

num_unit_input <- function(id, label, unit, value = NA, min = NA, max = NA, step = 0.1) {
  tagList(
    numericInput(id, label, value = value, min = min, max = max, step = step),
    unit_help(paste0("Unit: ", unit))
  )
}


## =========================================================
## Raw-input scoring layer
## =========================================================

module_input_block <- function(module_id, title, subtitle = NULL, default0 = NA, default24 = NA, default48 = NA) {
  ## Kept for backward compatibility, but the public v9 interface uses raw clinical inputs.
  tags$div(class = "ami-card", tags$h4(title), tags$p(class = "small-note", "This compatibility block is not shown in the v9 public interface."))
}

raw_id <- function(tt, var) paste0("raw_", tt, "_", var)

raw_vars <- c(
  "heart_rate", "sbp", "map", "resp_rate", "spo2", "temperature",
  "gcs_eye", "gcs_verbal", "gcs_motor",
  "creatinine", "bun", "sodium", "potassium", "chloride", "bicarbonate", "glucose",
  "wbc", "hemoglobin", "hematocrit", "platelet", "inr", "ptt",
  "lactate", "troponin_t"
)

timepoints_app <- c("T0", "T24", "T48")

app_component_spec <- tibble::tribble(
  ~module, ~module_pretty, ~variable, ~stat, ~direction,
  "hemodynamic", "Hemodynamic instability", "heart_rate", "max_value", "high",
  "hemodynamic", "Hemodynamic instability", "heart_rate", "abnormal_ratio", "high",
  "hemodynamic", "Hemodynamic instability", "sbp", "min_value", "low",
  "hemodynamic", "Hemodynamic instability", "sbp", "abnormal_ratio", "high",
  "hemodynamic", "Hemodynamic instability", "map", "min_value", "low",
  "hemodynamic", "Hemodynamic instability", "map", "abnormal_ratio", "high",

  "respiratory_neurologic", "Respiratory-neurologic stress", "resp_rate", "max_value", "high",
  "respiratory_neurologic", "Respiratory-neurologic stress", "resp_rate", "abnormal_ratio", "high",
  "respiratory_neurologic", "Respiratory-neurologic stress", "spo2", "min_value", "low",
  "respiratory_neurologic", "Respiratory-neurologic stress", "spo2", "abnormal_ratio", "high",
  "respiratory_neurologic", "Respiratory-neurologic stress", "temperature", "abnormal_ratio", "high",
  "respiratory_neurologic", "Respiratory-neurologic stress", "gcs_eye", "min_value", "low",
  "respiratory_neurologic", "Respiratory-neurologic stress", "gcs_verbal", "min_value", "low",
  "respiratory_neurologic", "Respiratory-neurologic stress", "gcs_motor", "min_value", "low",

  "renal_metabolic", "Renal-metabolic dysfunction", "creatinine", "max_value", "high",
  "renal_metabolic", "Renal-metabolic dysfunction", "bun", "max_value", "high",
  "renal_metabolic", "Renal-metabolic dysfunction", "sodium", "abnormal_ratio", "high",
  "renal_metabolic", "Renal-metabolic dysfunction", "potassium", "abnormal_ratio", "high",
  "renal_metabolic", "Renal-metabolic dysfunction", "chloride", "abnormal_ratio", "high",
  "renal_metabolic", "Renal-metabolic dysfunction", "bicarbonate", "min_value", "low",
  "renal_metabolic", "Renal-metabolic dysfunction", "bicarbonate", "abnormal_ratio", "high",
  "renal_metabolic", "Renal-metabolic dysfunction", "glucose", "abnormal_ratio", "high",

  "inflammatory_coagulopathic", "Inflammatory-coagulopathic response", "wbc", "max_value", "high",
  "inflammatory_coagulopathic", "Inflammatory-coagulopathic response", "wbc", "abnormal_ratio", "high",
  "inflammatory_coagulopathic", "Inflammatory-coagulopathic response", "hemoglobin", "min_value", "low",
  "inflammatory_coagulopathic", "Inflammatory-coagulopathic response", "platelet", "min_value", "low",
  "inflammatory_coagulopathic", "Inflammatory-coagulopathic response", "platelet", "abnormal_ratio", "high",
  "inflammatory_coagulopathic", "Inflammatory-coagulopathic response", "inr", "max_value", "high",
  "inflammatory_coagulopathic", "Inflammatory-coagulopathic response", "ptt", "max_value", "high",
  "inflammatory_coagulopathic", "Inflammatory-coagulopathic response", "ptt", "abnormal_ratio", "high",

  "cardiac_hypoperfusion", "Cardiac injury / hypoperfusion", "lactate", "max_value", "high",
  "cardiac_hypoperfusion", "Cardiac injury / hypoperfusion", "lactate", "abnormal_ratio", "high",
  "cardiac_hypoperfusion", "Cardiac injury / hypoperfusion", "troponin_t", "max_value", "high",
  "cardiac_hypoperfusion", "Cardiac injury / hypoperfusion", "troponin_t", "abnormal_ratio", "high"
) %>% dplyr::mutate(column = paste0(variable, "_", stat))

raw_public_columns <- function() {
  as.vector(outer(paste0("raw_", timepoints_app, "_"), raw_vars, paste0))
}

first_existing_col <- function(df, candidates) {
  candidates <- candidates[candidates %in% names(df)]
  if (length(candidates) == 0) return(NA_character_)
  candidates[1]
}

raw_value_from_df <- function(dat, tt, variable, stat = NULL) {
  direct <- c(
    paste0("raw_", tt, "_", variable),
    paste0("raw_", tolower(tt), "_", variable),
    paste0(tt, "_", variable),
    paste0(tolower(tt), "_", variable),
    paste0(variable, "_", tt),
    paste0(variable, "_", tolower(tt)),
    paste0(variable, "_", tt, "_value"),
    paste0(variable, "_", tolower(tt), "_value")
  )
  if (!is.null(stat)) {
    direct <- c(
      paste0(variable, "_", stat, "_", tt),
      paste0(variable, "_", stat, "_", tolower(tt)),
      paste0(variable, "_", tt, "_", stat),
      paste0(variable, "_", tolower(tt), "_", stat),
      direct
    )
  }
  hit <- first_existing_col(dat, direct)
  if (is.na(hit)) return(rep(NA_real_, nrow(dat)))
  to_numeric_safe(dat[[hit]])
}

abnormal_indicator <- function(variable, x) {
  x <- to_numeric_safe(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & !is.na(x)
  if (!any(ok)) return(out)
  out[ok] <- dplyr::case_when(
    variable == "heart_rate" ~ as.numeric(x[ok] < 60 | x[ok] > 100),
    variable == "sbp" ~ as.numeric(x[ok] < 90 | x[ok] > 180),
    variable == "map" ~ as.numeric(x[ok] < 65),
    variable == "resp_rate" ~ as.numeric(x[ok] < 10 | x[ok] > 22),
    variable == "spo2" ~ as.numeric(x[ok] < 92),
    variable == "temperature" ~ as.numeric(x[ok] < 36 | x[ok] > 38.3),
    variable == "gcs_eye" ~ as.numeric(x[ok] < 4),
    variable == "gcs_verbal" ~ as.numeric(x[ok] < 5),
    variable == "gcs_motor" ~ as.numeric(x[ok] < 6),
    variable == "sodium" ~ as.numeric(x[ok] < 135 | x[ok] > 145),
    variable == "potassium" ~ as.numeric(x[ok] < 3.5 | x[ok] > 5.0),
    variable == "chloride" ~ as.numeric(x[ok] < 98 | x[ok] > 107),
    variable == "bicarbonate" ~ as.numeric(x[ok] < 22 | x[ok] > 28),
    variable == "glucose" ~ as.numeric(x[ok] < 70 | x[ok] > 180),
    variable == "wbc" ~ as.numeric(x[ok] < 4 | x[ok] > 12),
    variable == "platelet" ~ as.numeric(x[ok] < 150 | x[ok] > 450),
    variable == "ptt" ~ as.numeric(x[ok] > 35),
    variable == "lactate" ~ as.numeric(x[ok] > 2),
    variable == "troponin_t" ~ as.numeric(x[ok] > 0.014),
    TRUE ~ NA_real_
  )
  out
}

component_raw_value <- function(dat, tt, variable, stat) {
  direct <- raw_value_from_df(dat, tt, variable, stat)
  if (any(!is.na(direct))) return(direct)
  base <- raw_value_from_df(dat, tt, variable, NULL)
  if (stat == "abnormal_ratio") return(abnormal_indicator(variable, base))
  base
}

load_score_reference <- function(project_root = default_project_root) {
  paths <- c(
    T0 = file.path(project_root, "03_analysis_dataset", "Step02_AMI_T0_dataset.csv"),
    T24 = file.path(project_root, "03_analysis_dataset", "Step02_AMI_T24_dataset.csv"),
    T48 = file.path(project_root, "03_analysis_dataset", "Step02_AMI_T48_dataset.csv")
  )
  refs <- list()
  for (tt in names(paths)) {
    if (!file.exists(paths[[tt]])) next
    df <- tryCatch(readr::read_csv(paths[[tt]], show_col_types = FALSE), error = function(e) NULL)
    if (is.null(df)) next
    refs[[tt]] <- list()
    for (col in unique(app_component_spec$column)) {
      if (col %in% names(df)) {
        xx <- to_numeric_safe(df[[col]])
        xx <- xx[is.finite(xx) & !is.na(xx)]
        if (length(xx) >= 20) refs[[tt]][[col]] <- sort(xx)
      }
    }
  }
  refs
}

score_against_reference <- function(value, ref, direction, variable = NULL, stat = NULL) {
  value <- to_numeric_safe(value)
  out <- rep(NA_real_, length(value))
  ok <- is.finite(value) & !is.na(value)
  if (!any(ok)) return(out)
  if (!is.null(ref) && length(ref) >= 20 && length(unique(ref)) > 1) {
    pct <- stats::ecdf(ref)(value[ok]) * 100
    if (identical(direction, "low")) pct <- 100 - pct
    out[ok] <- pmin(pmax(pct, 0), 100)
  } else {
    ## Fallback only for demonstration if reference distributions are unavailable.
    out[ok] <- clinical_fallback_score(variable, stat, value[ok], direction)
  }
  out
}

linear_score <- function(x, low_good, high_bad, reverse = FALSE) {
  z <- (x - low_good) / (high_bad - low_good) * 100
  if (reverse) z <- 100 - z
  pmin(pmax(z, 0), 100)
}

clinical_fallback_score <- function(variable, stat, x, direction) {
  if (stat == "abnormal_ratio") return(pmin(pmax(x * 100, 0), 100))
  dplyr::case_when(
    variable == "heart_rate" ~ linear_score(x, 70, 140),
    variable == "sbp" ~ linear_score(x, 120, 70, reverse = TRUE),
    variable == "map" ~ linear_score(x, 85, 55, reverse = TRUE),
    variable == "resp_rate" ~ linear_score(x, 16, 34),
    variable == "spo2" ~ linear_score(x, 98, 82, reverse = TRUE),
    variable == "temperature" ~ linear_score(abs(x - 37), 0, 2.5),
    variable == "gcs_eye" ~ linear_score(x, 4, 1, reverse = TRUE),
    variable == "gcs_verbal" ~ linear_score(x, 5, 1, reverse = TRUE),
    variable == "gcs_motor" ~ linear_score(x, 6, 1, reverse = TRUE),
    variable == "creatinine" ~ linear_score(x, 0.8, 4.0),
    variable == "bun" ~ linear_score(x, 15, 80),
    variable == "bicarbonate" ~ linear_score(x, 24, 12, reverse = TRUE),
    variable == "wbc" ~ linear_score(x, 8, 30),
    variable == "hemoglobin" ~ linear_score(x, 14, 7, reverse = TRUE),
    variable == "platelet" ~ linear_score(x, 250, 50, reverse = TRUE),
    variable == "inr" ~ linear_score(x, 1.0, 3.0),
    variable == "ptt" ~ linear_score(x, 30, 90),
    variable == "lactate" ~ linear_score(x, 1, 8),
    variable == "troponin_t" ~ linear_score(log10(pmax(x, 1e-4)), log10(0.014), log10(5)),
    TRUE ~ if (identical(direction, "low")) linear_score(x, 1, 0, reverse = TRUE) else linear_score(x, 0, 1)
  )
}

row_mean_available <- function(xmat) {
  out <- rowMeans(xmat, na.rm = TRUE)
  out[rowSums(!is.na(xmat)) < 1] <- NA_real_
  out
}

simple_slope <- function(v0, v24, v48) {
  vals <- cbind(v0, v24, v48)
  apply(vals, 1, function(z) {
    tt <- c(0, 24, 48)
    ok <- is.finite(z) & !is.na(z)
    if (sum(ok) < 2) return(NA_real_)
    as.numeric(stats::coef(stats::lm(z[ok] ~ tt[ok]))[2])
  })
}

compute_features_from_raw <- function(dat, reference = score_reference_global) {
  dat <- as.data.frame(dat)
  if (nrow(dat) == 0) return(dat)

  for (tt in timepoints_app) {
    for (m in unique(app_component_spec$module)) {
      comp_scores <- list()
      idx <- which(app_component_spec$module == m)
      for (ii in idx) {
        variable <- app_component_spec$variable[ii]
        stat <- app_component_spec$stat[ii]
        col <- app_component_spec$column[ii]
        val <- component_raw_value(dat, tt, variable, stat)
        ref <- NULL
        if (!is.null(reference[[tt]]) && !is.null(reference[[tt]][[col]])) ref <- reference[[tt]][[col]]
        comp_scores[[col]] <- score_against_reference(val, ref, app_component_spec$direction[ii], variable, stat)
      }
      score_col <- paste0(m, "_", tt, "_score")
      if (!score_col %in% names(dat) || all(is.na(to_numeric_safe(dat[[score_col]])))) {
        dat[[score_col]] <- row_mean_available(as.matrix(as.data.frame(comp_scores)))
      }
      alias_col <- paste0(m, "_", tt)
      if (!alias_col %in% names(dat)) dat[[alias_col]] <- dat[[score_col]]
    }
  }

  dat <- add_module_deltas(dat)

  getv <- function(tt, variable) raw_value_from_df(dat, tt, variable)
  vmean <- function(variable) rowMeans(cbind(getv("T0", variable), getv("T24", variable), getv("T48", variable)), na.rm = TRUE)
  vmax <- function(variable) apply(cbind(getv("T0", variable), getv("T24", variable), getv("T48", variable)), 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  vmin <- function(variable) apply(cbind(getv("T0", variable), getv("T24", variable), getv("T48", variable)), 1, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
  vsd <- function(variable) apply(cbind(getv("T0", variable), getv("T24", variable), getv("T48", variable)), 1, function(x) if (sum(!is.na(x)) < 2) NA_real_ else sd(x, na.rm = TRUE))
  vratio <- function(variable) rowMeans(cbind(abnormal_indicator(variable, getv("T0", variable)), abnormal_indicator(variable, getv("T24", variable)), abnormal_indicator(variable, getv("T48", variable))), na.rm = TRUE)

  assign_if_missing <- function(col, val) {
    if (!col %in% names(dat) || all(is.na(to_numeric_safe(dat[[col]])))) dat[[col]] <<- val
  }

  assign_if_missing("gcs_verbal_last_value", getv("T48", "gcs_verbal"))
  assign_if_missing("gcs_verbal_abnormal_ratio", vratio("gcs_verbal"))
  assign_if_missing("gcs_verbal_slope_per_hour", simple_slope(getv("T0", "gcs_verbal"), getv("T24", "gcs_verbal"), getv("T48", "gcs_verbal")))
  assign_if_missing("gcs_eye_slope_per_hour", simple_slope(getv("T0", "gcs_eye"), getv("T24", "gcs_eye"), getv("T48", "gcs_eye")))
  assign_if_missing("gcs_motor_slope_per_hour", simple_slope(getv("T0", "gcs_motor"), getv("T24", "gcs_motor"), getv("T48", "gcs_motor")))

  assign_if_missing("bun_mean_value", vmean("bun"))
  assign_if_missing("bun_max_value", vmax("bun"))
  assign_if_missing("bun_last_value", getv("T48", "bun"))
  assign_if_missing("spo2_mean_value", vmean("spo2"))
  assign_if_missing("ptt_min_value", vmin("ptt"))
  assign_if_missing("chloride_max_value", vmax("chloride"))
  assign_if_missing("lactate_slope_per_hour", simple_slope(getv("T0", "lactate"), getv("T24", "lactate"), getv("T48", "lactate")))
  assign_if_missing("heart_rate_abnormal_ratio", vratio("heart_rate"))
  assign_if_missing("bicarbonate_abnormal_ratio", vratio("bicarbonate"))
  assign_if_missing("potassium_abnormal_ratio", vratio("potassium"))
  assign_if_missing("hemoglobin_mean_value", vmean("hemoglobin"))
  assign_if_missing("hemoglobin_sd_value", vsd("hemoglobin"))
  assign_if_missing("hematocrit_sd_value", vsd("hematocrit"))

  dat
}

make_raw_patient_features_from_ui <- function(input) {
  dat <- as.data.frame(matrix(NA, nrow = 1, ncol = length(raw_public_columns())))
  names(dat) <- raw_public_columns()
  for (tt in timepoints_app) {
    for (v in raw_vars) {
      id <- raw_id(tt, v)
      dat[[id]][1] <- input[[id]]
    }
  }
  compute_features_from_raw(dat)
}

raw_landmark_input_block <- function(tt, title) {
  tags$div(
    class = "ami-card raw-landmark-card",
    tags$div(class = "section-kicker", paste0(title, " raw measurements")),
    tags$h4(paste0(title, " clinical measurements")),
    tags$p(class = "small-note", "Enter representative worst or latest values available by this landmark. The app converts these raw measurements into organ-domain severity scores automatically."),
    tags$h5("Vital signs and neurologic status"),
    fluidRow(
      column(3, num_unit_input(raw_id(tt, "heart_rate"), "Heart rate", "beats/min", NA, 20, 250, 1)),
      column(3, num_unit_input(raw_id(tt, "sbp"), "Systolic BP", "mmHg", NA, 30, 260, 1)),
      column(3, num_unit_input(raw_id(tt, "map"), "MAP", "mmHg", NA, 20, 180, 1)),
      column(3, num_unit_input(raw_id(tt, "resp_rate"), "Respiratory rate", "breaths/min", NA, 4, 70, 1))
    ),
    fluidRow(
      column(3, num_unit_input(raw_id(tt, "spo2"), "SpO₂", "%", NA, 50, 100, 0.1)),
      column(3, num_unit_input(raw_id(tt, "temperature"), "Temperature", "°C", NA, 30, 43, 0.1)),
      column(2, num_unit_input(raw_id(tt, "gcs_eye"), "GCS eye", "1–4", NA, 1, 4, 1)),
      column(2, num_unit_input(raw_id(tt, "gcs_verbal"), "GCS verbal", "1–5", NA, 1, 5, 1)),
      column(2, num_unit_input(raw_id(tt, "gcs_motor"), "GCS motor", "1–6", NA, 1, 6, 1))
    ),
    tags$h5("Renal-metabolic and electrolyte markers"),
    fluidRow(
      column(3, num_unit_input(raw_id(tt, "creatinine"), "Creatinine", "mg/dL", NA, 0, 20, 0.1)),
      column(3, num_unit_input(raw_id(tt, "bun"), "BUN", "mg/dL", NA, 0, 200, 0.1)),
      column(3, num_unit_input(raw_id(tt, "sodium"), "Sodium", "mmol/L", NA, 100, 180, 0.1)),
      column(3, num_unit_input(raw_id(tt, "potassium"), "Potassium", "mmol/L", NA, 1.5, 9, 0.1))
    ),
    fluidRow(
      column(3, num_unit_input(raw_id(tt, "chloride"), "Chloride", "mmol/L", NA, 70, 140, 0.1)),
      column(3, num_unit_input(raw_id(tt, "bicarbonate"), "Bicarbonate", "mmol/L", NA, 5, 45, 0.1)),
      column(3, num_unit_input(raw_id(tt, "glucose"), "Glucose", "mg/dL", NA, 20, 800, 1))
    ),
    tags$h5("Inflammation, coagulation, perfusion, and cardiac injury"),
    fluidRow(
      column(3, num_unit_input(raw_id(tt, "wbc"), "WBC", "10⁹/L", NA, 0, 100, 0.1)),
      column(3, num_unit_input(raw_id(tt, "hemoglobin"), "Hemoglobin", "g/dL", NA, 3, 25, 0.1)),
      column(3, num_unit_input(raw_id(tt, "hematocrit"), "Hematocrit", "%", NA, 10, 70, 0.1)),
      column(3, num_unit_input(raw_id(tt, "platelet"), "Platelet", "10⁹/L", NA, 0, 1000, 1))
    ),
    fluidRow(
      column(3, num_unit_input(raw_id(tt, "inr"), "INR", "ratio", NA, 0.5, 10, 0.1)),
      column(3, num_unit_input(raw_id(tt, "ptt"), "PTT", "seconds", NA, 10, 200, 0.1)),
      column(3, num_unit_input(raw_id(tt, "lactate"), "Lactate", "mmol/L", NA, 0, 25, 0.1)),
      column(3, num_unit_input(raw_id(tt, "troponin_t"), "Troponin T", "ng/mL", NA, 0, 100, 0.01))
    )
  )
}

score_reference_global <- tryCatch(load_score_reference(default_project_root), error = function(e) list())

make_arc_df <- function(xmin, xmax, risk_group, r_outer = 1, r_inner = 0.62, max_x = 70, n = 80) {
  theta1 <- pi - (xmin / max_x) * pi
  theta2 <- pi - (xmax / max_x) * pi
  theta <- seq(theta1, theta2, length.out = n)
  outer <- tibble::tibble(x = r_outer * cos(theta), y = r_outer * sin(theta), risk_group = risk_group)
  inner <- tibble::tibble(x = r_inner * cos(rev(theta)), y = r_inner * sin(rev(theta)), risk_group = risk_group)
  dplyr::bind_rows(outer, inner)
}

risk_gauge_plot <- function(p) {
  pct <- pmin(pmax(p * 100, 0), 70)
  bands <- tibble::tribble(
    ~xmin, ~xmax, ~risk_group,
    0, 5, "Low risk (<5%)",
    5, 15, "Intermediate risk (5-15%)",
    15, 30, "High risk (15-30%)",
    30, 70, "Very high risk (>=30%)"
  )
  arc_df <- purrr::pmap_dfr(bands, make_arc_df)
  angle <- pi - (pct / 70) * pi
  needle <- tibble::tibble(x0 = 0, y0 = 0, x1 = 0.82 * cos(angle), y1 = 0.82 * sin(angle))

  tick_df <- tibble::tibble(value = c(0, 5, 15, 30, 70)) %>%
    dplyr::mutate(
      theta = pi - (value / 70) * pi,
      x = 1.08 * cos(theta),
      y = 1.08 * sin(theta),
      label = paste0(value, "%")
    )

  ggplot() +
    geom_polygon(data = arc_df, aes(x = x, y = y, group = risk_group, fill = risk_group), color = "white", linewidth = 0.8) +
    geom_segment(data = needle, aes(x = x0, y = y0, xend = x1, yend = y1), linewidth = 1.5, color = "#0B1F33", lineend = "round") +
    geom_point(aes(x = 0, y = 0), size = 5, color = "#0B1F33") +
    geom_text(data = tick_df, aes(x = x, y = y, label = label), size = 4.7, fontface = "bold", color = "#14283A") +
    annotate("text", x = 0, y = -0.16, label = paste0("Predicted risk\n", pretty_pct(p)), fontface = "bold", size = 6.6, color = "#0B1F33") +
    scale_fill_manual(values = pal_risk, drop = FALSE) +
    coord_equal(xlim = c(-1.18, 1.18), ylim = c(-0.28, 1.18), clip = "off") +
    labs(title = "AMI-EDRS risk gauge", fill = NULL) +
    theme_void(base_size = 17) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      plot.margin = margin(10, 20, 10, 20)
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
}

organ_domain_plot_func <- function(df) {
  modules <- c(
    "hemodynamic",
    "respiratory_neurologic",
    "renal_metabolic",
    "inflammatory_coagulopathic",
    "cardiac_hypoperfusion"
  )
  module_labels <- c(
    hemodynamic = "Hemodynamic",
    respiratory_neurologic = "Respiratory–neurologic",
    renal_metabolic = "Renal–metabolic",
    inflammatory_coagulopathic = "Inflammatory–coagulopathic",
    cardiac_hypoperfusion = "Cardiac injury / hypoperfusion"
  )

  plot_df <- purrr::map_dfr(modules, function(m) {
    tibble::tibble(
      module = module_labels[[m]],
      landmark = factor(c("T0", "T24", "T48"), levels = c("T0", "T24", "T48")),
      score = c(
        to_numeric_safe(df[[paste0(m, "_T0_score")]]),
        to_numeric_safe(df[[paste0(m, "_T24_score")]]),
        to_numeric_safe(df[[paste0(m, "_T48_score")]])
      )
    )
  }) %>% dplyr::filter(!is.na(score))

  if (nrow(plot_df) == 0) {
    return(ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "No organ-domain scores available."))
  }

  ggplot(plot_df, aes(x = landmark, y = score, group = module, color = module)) +
    geom_line(linewidth = 1.15) +
    geom_point(size = 3.0) +
    scale_color_manual(values = pal_modules, drop = FALSE) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    labs(title = "Early organ-domain dynamics", x = "Landmark", y = "Organ-domain score", color = NULL) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 21),
      axis.title = element_text(face = "bold", size = 15),
      axis.text = element_text(size = 13, color = "#14283A"),
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      legend.key.width = grid::unit(24, "pt"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(8, 10, 20, 8)
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
}

## =========================================================
## 5. Shiny UI
## =========================================================

ui <- bslib::page_navbar(
  title = tags$span(style = "font-weight:900; letter-spacing:-0.02em;", "AMI-EDRS"),
  theme = bslib::bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#1F3A5F",
    base_font = bslib::font_google("Inter")
  ),
  header = tags$head(tags$style(HTML(clinical_css))),

  bslib::nav_panel(
    "Single-patient prediction",
    tags$div(
      class = "ami-hero",
      tags$h2("AMI-EDRS Early Dynamic Risk Dashboard v12"),
      tags$p("A raw-input, trajectory-enhanced T48 landmark system for 30-day mortality risk stratification after acute myocardial infarction.")
    ),

    fluidRow(
      column(
        3,
        tags$div(
          class = "ami-card quick-panel",
          tags$h4("Quick start"),
          uiOutput("model_status"),
          tags$hr(),
          downloadButton("download_example_single", "Download example dataset"),
          tags$p(class = "small-note", "Use the demo buttons to preview raw-input prediction. For patient-data upload, use the Batch prediction page; the uploaded file should follow the same raw-column format as the example dataset.")
        ),
        tags$div(
          class = "disclaimer",
          tags$b("Research-use notice."), tags$br(),
          "AMI-EDRS supports risk communication and research demonstration only. It does not replace clinician judgment and has not been prospectively validated as a standalone medical device."
        )
      ),

      column(
        9,
        tabsetPanel(
          tabPanel(
            "1. Patient profile",
            tags$div(
              class = "ami-card",
              tags$div(class = "section-kicker", "Clinical input"),
              tags$h4("Host vulnerability and AMI context"),
              fluidRow(
                column(3, numericInput("age", "Age, years", value = 71, min = 18, max = 110, step = 1)),
                column(3, selectInput("sex", "Sex", choices = c("Male", "Female", "Missing"), selected = "Male")),
                column(3, checkboxInput("ami_subtype_stemi", "STEMI subtype", value = FALSE)),
                column(3, checkboxInput("prior_cad_ihd", "Prior CAD/IHD", value = TRUE))
              ),
              fluidRow(
                column(3, checkboxInput("hypertension", "Hypertension", value = TRUE)),
                column(3, checkboxInput("diabetes", "Diabetes", value = FALSE)),
                column(3, checkboxInput("ckd", "CKD / renal disease", value = FALSE)),
                column(3, checkboxInput("copd", "COPD", value = FALSE))
              ),
              fluidRow(
                column(3, checkboxInput("heart_failure", "Heart failure", value = FALSE)),
                column(3, checkboxInput("atrial_fibrillation", "Atrial fibrillation", value = FALSE)),
                column(3, checkboxInput("stroke", "Stroke", value = FALSE)),
                column(3, checkboxInput("liver_disease", "Liver disease", value = FALSE))
              ),
              fluidRow(
                column(3, checkboxInput("malignancy", "Malignancy", value = FALSE)),
                column(3, checkboxInput("shock_proxy", "Shock / hypoperfusion", value = FALSE))
              )
            )
          ),

          tabPanel(
            "2. Clinical measurements",
            tags$div(
              class = "input-guide",
              tags$b("Raw-input mode"), tags$br(),
              "Enter the patient's raw measurements available by T0, T24, and T48. The app automatically converts these raw values into 0–100 organ-domain dysfunction scores, assigns the C1–C5 trajectory phenotype when enough information is available, and then estimates AMI-EDRS risk.",
              tags$br(), tags$br(),
              tags$b("For routine use: "), "upload a table with the same raw-column format as the downloadable example dataset. No pre-calculated organ-domain score is required."
            ),
            uiOutput("trajectory_status_box"),
            tabsetPanel(
              tabPanel("T0", raw_landmark_input_block("T0", "T0")),
              tabPanel("T24", raw_landmark_input_block("T24", "T24")),
              tabPanel("T48", raw_landmark_input_block("T48", "T48"))
            )
          ),

          tabPanel(
            "3. Advanced markers",
            tags$div(
              class = "ami-card",
              tags$div(class = "section-kicker", "Optional repeated-measure features"),
              tags$h4("Optional EHR-derived markers"),
              tags$p(class = "small-note", "These fields are optional. Most single-patient use only requires the raw T0/T24/T48 measurements in the previous tab. If a hospital information system can compute repeated-measure burden or slope features, they can be entered here and will override simple values derived from the three landmarks."),
              fluidRow(
                column(3, num_unit_input("gcs_verbal_last_value", "GCS verbal last", "points, 1–5", NA, 1, 5, 1)),
                column(3, num_unit_input("bun_mean_value", "BUN mean", "mg/dL", NA, NA, NA, 0.1)),
                column(3, num_unit_input("bun_max_value", "BUN max", "mg/dL", NA, NA, NA, 0.1)),
                column(3, num_unit_input("bun_last_value", "BUN last", "mg/dL", NA, NA, NA, 0.1))
              ),
              fluidRow(
                column(3, num_unit_input("spo2_mean_value", "SpO₂ mean", "%", NA, 0, 100, 0.1)),
                column(3, num_unit_input("lactate_slope_per_hour", "Lactate slope", "mmol/L per hour", NA, NA, NA, 0.01)),
                column(3, num_unit_input("ptt_min_value", "PTT min", "seconds", NA, NA, NA, 0.1)),
                column(3, num_unit_input("chloride_max_value", "Chloride max", "mmol/L", NA, NA, NA, 0.1))
              ),
              fluidRow(
                column(3, num_unit_input("heart_rate_abnormal_ratio", "Heart-rate abnormal ratio", "proportion, 0–1", NA, 0, 1, 0.01)),
                column(3, num_unit_input("gcs_verbal_abnormal_ratio", "GCS verbal abnormal ratio", "proportion, 0–1", NA, 0, 1, 0.01)),
                column(3, num_unit_input("bicarbonate_abnormal_ratio", "Bicarbonate abnormal ratio", "proportion, 0–1", NA, 0, 1, 0.01)),
                column(3, num_unit_input("potassium_abnormal_ratio", "Potassium abnormal ratio", "proportion, 0–1", NA, 0, 1, 0.01))
              ),
              fluidRow(
                column(3, num_unit_input("hemoglobin_mean_value", "Hemoglobin mean", "g/dL", NA, NA, NA, 0.1)),
                column(3, num_unit_input("hemoglobin_sd_value", "Hemoglobin SD", "g/dL", NA, NA, NA, 0.1)),
                column(3, num_unit_input("hematocrit_sd_value", "Hematocrit SD", "%", NA, NA, NA, 0.1)),
                column(3, num_unit_input("gcs_eye_slope_per_hour", "GCS eye slope", "points per hour", NA, NA, NA, 0.01))
              ),
              fluidRow(
                column(3, num_unit_input("gcs_verbal_slope_per_hour", "GCS verbal slope", "points per hour", NA, NA, NA, 0.01)),
                column(3, num_unit_input("gcs_motor_slope_per_hour", "GCS motor slope", "points per hour", NA, NA, NA, 0.01))
              )
            )
          ),

          tabPanel(
            "4. Results",
            tags$div(
              class = "action-row",
              tags$div(class = "action-cell", actionButton("predict_single", "Run AMI-EDRS prediction", class = "btn-run-main")),
              tags$div(class = "action-cell", actionButton("load_demo_low", "Show low-risk demo", class = "btn-demo-low")),
              tags$div(class = "action-cell", actionButton("load_demo_high", "Show high-risk demo", class = "btn-demo-high"))
            ),
            tags$br(),
            uiOutput("single_results_cards"),
            fluidRow(
              column(6, plotOutput("risk_gauge_plot", height = 460)),
              column(6, plotOutput("organ_domain_plot", height = 480))
            ),
            tags$div(class = "ami-card", tags$h4("Patient-level prediction summary"), uiOutput("single_summary_table")),
            downloadButton("download_single_csv", "Download patient report CSV")
          )
        )
      )
    )
  ),

  bslib::nav_panel(
    "Batch prediction",
    tags$div(
      class = "ami-hero",
      tags$h2("Batch AMI-EDRS Prediction"),
      tags$p("Upload a CSV or Excel table with patient-level features. Use the example dataset to see the expected format.")
    ),
    fluidRow(
      column(
        4,
        tags$div(
          class = "ami-card",
          tags$h4("Upload patient dataset"),
          tags$p(class = "small-note", "One upload entry point is used for patient data. Please use a CSV or Excel file with the same column structure as the example dataset."),
          downloadButton("download_example", "Download example dataset"),
          tags$br(), tags$br(),
          fileInput("batch_file", "Upload CSV or Excel", accept = c(".csv", ".xlsx", ".xls")),
          checkboxInput("auto_assign_trajectory", "Auto-assign C1–C5 if trajectory column is missing", value = TRUE),
          actionButton("run_batch", "Run batch prediction", class = "btn-primary"),
          tags$hr(),
          downloadButton("download_batch", "Download predictions CSV"),
          tags$p(class = "small-note", "The example dataset contains two demo patients and is the recommended upload format for first-time users.")
        )
      ),
      column(
        8,
        uiOutput("batch_status"),
        tags$div(class = "ami-card", tags$h4("Batch prediction table"), uiOutput("batch_table")),
        fluidRow(
          column(6, plotOutput("batch_risk_distribution", height = 340)),
          column(6, plotOutput("batch_trajectory_distribution", height = 340))
        )
      )
    )
  ),

  bslib::nav_panel(
    "Documentation",
    tags$div(
      class = "ami-hero",
      tags$h2("AMI-EDRS Model Documentation"),
      tags$p("Final model: elastic-net trajectory-enhanced dynamic risk model for 30-day mortality after the T48 landmark.")
    ),
    fluidRow(
      column(
        6,
        tags$div(
          class = "ami-card",
          tags$h4("Model concept"),
          tags$ul(
            tags$li("T0/T24/T48 organ-domain dynamics define early AMI pathophysiologic evolution."),
            tags$li("C1–C5 trajectory phenotype captures early disease-state heterogeneity."),
            tags$li("AMI-EDRS integrates host vulnerability, dynamic features, and trajectory phenotype."),
            tags$li("Risk strata: <5%, 5–15%, 15–30%, and ≥30% predicted 30-day mortality after T48.")
          )
        )
      ),
      column(
        6,
        tags$div(
          class = "ami-card",
          tags$h4("Input modes"),
          tags$ul(
            tags$li("Single-patient mode accepts raw T0/T24/T48 clinical measurements; organ-domain scores are calculated automatically."),
            tags$li("Batch mode accepts the same raw-column format as the example dataset and then calculates organ-domain scores, trajectory phenotype, and AMI-EDRS risk."),
            tags$li("C1–C5 is assigned only when enough T0–T24–T48 measurements are available to reconstruct the trajectory features."),
            tags$li("If trajectory features are insufficient, AMI-EDRS can still run with missing trajectory handled by the saved preprocessing recipe, but the phenotype is reported as unavailable.")
          )
        )
      )
    ),
    tags$div(class = "ami-card", tags$h4("System readiness"), verbatimTextOutput("model_details"))
  )
)

## =========================================================
## 6. Server
## =========================================================

server <- function(input, output, session) {
  rv <- reactiveValues(
    step05 = NULL,
    step03 = NULL,
    single_input = NULL,
    single_pred = NULL,
    single_trajectory_status = NULL,
    batch_pred = NULL
  )

  observe({
    if (is.null(rv$step05) && file.exists(default_step05_rds)) {
      rv$step05 <- readRDS(default_step05_rds)
    }
    if (is.null(rv$step03) && file.exists(default_step03_rds)) {
      rv$step03 <- readRDS(default_step03_rds)
    }
  })

  observeEvent(input$step05_rds, {
    req(input$step05_rds$datapath)
    rv$step05 <- readRDS(input$step05_rds$datapath)
    showNotification("AMI-EDRS prediction engine loaded.", type = "message")
  })

  observeEvent(input$step03_rds, {
    req(input$step03_rds$datapath)
    rv$step03 <- readRDS(input$step03_rds$datapath)
    showNotification("C1–C5 trajectory engine loaded.", type = "message")
  })

  output$model_status <- renderUI({
    step05_status <- if (!is.null(rv$step05)) {
      tags$span(class = "status-pill status-good", "Prediction engine ready")
    } else {
      tags$span(class = "status-pill status-bad", "Prediction engine not loaded")
    }
    step03_status <- if (!is.null(rv$step03)) {
      tags$span(class = "status-pill status-good", "Trajectory engine ready")
    } else {
      tags$span(class = "status-pill status-warn", "Trajectory engine optional")
    }
    tagList(step05_status, tags$br(), step03_status)
  })

  output$model_selector <- renderUI({
    if (is.null(rv$step05)) return(tags$p(class = "small-note", "Load the prediction engine to select a model."))
    models <- names(rv$step05$model_objects)
    if ("Elastic net: dynamic + trajectory" %in% models) {
      choices <- "Elastic net: dynamic + trajectory"
      selected <- "Elastic net: dynamic + trajectory"
    } else {
      choices <- models
      selected <- models[1]
    }
    tagList(
      selectInput("primary_model", "Primary prediction model", choices = choices, selected = selected),
      tags$p(class = "small-note", "The application defaults to the final elastic-net AMI-EDRS model. Additional algorithms are hidden to avoid clinical ambiguity.")
    )
  })

  selected_model <- reactive({
    req(rv$step05)
    input$primary_model %||% if ("Elastic net: dynamic + trajectory" %in% names(rv$step05$model_objects)) {
      "Elastic net: dynamic + trajectory"
    } else {
      names(rv$step05$model_objects)[1]
    }
  })

  current_predictors <- reactive({
    if (is.null(rv$step05)) return(public_feature_columns())
    obj <- rv$step05$model_objects[[selected_model()]]
    if (is.null(obj) || is.null(obj$prep) || is.null(obj$prep$predictors)) return(public_feature_columns())
    unique(c(public_feature_columns(), obj$prep$predictors))
  })

  trajectory_ui_df <- reactive({
    make_trajectory_input_df_from_ui(input)
  })

  trajectory_status_reactive <- reactive({
    auto_trajectory_status(trajectory_ui_df(), rv$step03)
  })

  output$trajectory_status_box <- renderUI({
    status <- trajectory_status_reactive()
    if (isTRUE(status$available)) {
      tags$div(class = "status-pill status-good", status$message)
    } else {
      tags$div(class = "status-pill status-warn", status$message)
    }
  })

  load_demo_values <- function(high = FALSE) {
    ex <- make_example_data()
    row <- if (isTRUE(high)) ex[2, , drop = FALSE] else ex[1, , drop = FALSE]

    updateNumericInput(session, "age", value = row$age[1])
    updateSelectInput(session, "sex", selected = row$sex[1])
    for (nm in c("hypertension", "diabetes", "ckd", "copd", "malignancy", "liver_disease", "stroke", "shock_proxy", "prior_cad_ihd", "ami_subtype_stemi", "heart_failure", "atrial_fibrillation")) {
      if (nm %in% names(row)) updateCheckboxInput(session, nm, value = as.logical(as.integer(row[[nm]][1])))
    }
    for (tt in timepoints_app) {
      for (v in raw_vars) {
        id <- raw_id(tt, v)
        if (id %in% names(row)) updateNumericInput(session, id, value = row[[id]][1])
      }
    }
    ## Clear optional repeated-measure fields so derived raw features are used.
    optional_ids <- c("gcs_verbal_last_value", "bun_mean_value", "bun_max_value", "bun_last_value", "spo2_mean_value", "lactate_slope_per_hour", "ptt_min_value", "chloride_max_value", "heart_rate_abnormal_ratio", "gcs_verbal_abnormal_ratio", "bicarbonate_abnormal_ratio", "potassium_abnormal_ratio", "hemoglobin_mean_value", "hemoglobin_sd_value", "hematocrit_sd_value", "gcs_eye_slope_per_hour", "gcs_verbal_slope_per_hour", "gcs_motor_slope_per_hour")
    for (id in optional_ids) updateNumericInput(session, id, value = NA)
  }

  ## Demo buttons: compute directly from the example row instead of waiting for
  ## asynchronous UI input updates. This prevents double-click/server errors on
  ## shinyapps.io and makes the demo deterministic.
  run_demo_prediction_from_row <- function(high = FALSE) {
    ex <- make_example_data()
    row <- if (isTRUE(high)) ex[2, , drop = FALSE] else ex[1, , drop = FALSE]

    ## Still update visible inputs so users can inspect the example values.
    load_demo_values(high)

    local_state <- isolate({
      step05_obj <- rv$step05
      step03_obj <- rv$step03

      dat <- as.data.frame(row)
      names(dat) <- make.names(names(dat), unique = TRUE)
      dat <- compute_features_from_raw(dat)
      dat <- add_module_deltas(dat)

      traj_status <- if (!is.null(step03_obj)) auto_trajectory_status(dat, step03_obj) else list(available = FALSE, message = "Trajectory engine unavailable; demonstration phenotype is shown.")
      trajectory_auto <- NA_character_
      if (!is.null(step03_obj) && isTRUE(traj_status$available)) {
        trajectory_auto <- assign_trajectory_auto(dat, step03_obj)[1]
      }

      ## If the deployed server has no trajectory object, still show the intended
      ## demonstration phenotype so the demo is visually complete.
      if (is.na(trajectory_auto)) {
        trajectory_auto <- if (isTRUE(high)) "C5" else "C1"
        traj_status <- list(available = TRUE, message = "C1–C5 phenotype shown for the demonstration patient.")
      }

      model_nm <- NA_character_
      predictors <- public_feature_columns()
      if (!is.null(step05_obj)) {
        model_nm <- input$primary_model %||% if ("Elastic net: dynamic + trajectory" %in% names(step05_obj$model_objects)) {
          "Elastic net: dynamic + trajectory"
        } else {
          names(step05_obj$model_objects)[1]
        }
        obj <- step05_obj$model_objects[[model_nm]]
        if (!is.null(obj) && !is.null(obj$prep) && !is.null(obj$prep$predictors)) {
          predictors <- unique(c(public_feature_columns(), obj$prep$predictors, names(dat)))
        }
      } else {
        predictors <- unique(c(public_feature_columns(), names(dat)))
      }

      patient_df <- as.data.frame(matrix(NA, nrow = 1, ncol = length(predictors)))
      names(patient_df) <- predictors
      for (v in intersect(names(dat), names(patient_df))) patient_df[[v]][1] <- dat[[v]][1]
      if (!"trajectory" %in% names(patient_df)) patient_df$trajectory <- NA_character_
      patient_df$trajectory[1] <- trajectory_auto

      list(
        step05_obj = step05_obj,
        patient_df = patient_df,
        model_nm = model_nm,
        trajectory_auto = trajectory_auto,
        traj_status = traj_status,
        age = to_numeric_safe(row$age[1]),
        sex = as.character(row$sex[1])
      )
    })

    ## v13: demo display is isolated from the formal prediction engine.
    ## This prevents demo interactions from contaminating the next user-entered
    ## formal prediction. The demo still updates visible fields for inspection,
    ## but the result shown here is deliberately a safe demonstration result.
    pred <- fallback_demo_prediction_df(local_state$patient_df, trajectory_auto = local_state$trajectory_auto)

    patient_meta <- tibble::tibble(
      age = local_state$age,
      sex = local_state$sex,
      trajectory = local_state$trajectory_auto,
      trajectory_available = !is.na(local_state$trajectory_auto),
      trajectory_message = local_state$traj_status$message,
      trajectory_interpretation = trajectory_interpretation(local_state$trajectory_auto)
    )

    rv$single_input <- local_state$patient_df
    rv$single_pred <- bind_cols(patient_meta, pred)
    rv$single_trajectory_status <- local_state$traj_status
  }

  demo_static_fallback <- function(high = FALSE, error_message = NULL) {
    ex <- make_example_data()
    row <- if (isTRUE(high)) ex[2, , drop = FALSE] else ex[1, , drop = FALSE]
    dat <- as.data.frame(row)
    names(dat) <- make.names(names(dat), unique = TRUE)
    dat <- compute_features_from_raw(dat)
    dat <- add_module_deltas(dat)
    trajectory_auto <- if (isTRUE(high)) "C5" else "C1"
    predictors <- unique(c(public_feature_columns(), names(dat), "trajectory"))
    patient_df <- as.data.frame(matrix(NA, nrow = 1, ncol = length(predictors)))
    names(patient_df) <- predictors
    for (v in intersect(names(dat), names(patient_df))) patient_df[[v]][1] <- dat[[v]][1]
    patient_df$trajectory[1] <- trajectory_auto
    pred <- fallback_demo_prediction_df(patient_df, trajectory_auto = trajectory_auto)
    patient_meta <- tibble::tibble(
      age = to_numeric_safe(row$age[1]),
      sex = as.character(row$sex[1]),
      trajectory = trajectory_auto,
      trajectory_available = TRUE,
      trajectory_message = "C1–C5 phenotype shown for the demonstration patient.",
      trajectory_interpretation = trajectory_interpretation(trajectory_auto)
    )
    rv$single_input <- patient_df
    rv$single_pred <- bind_cols(patient_meta, pred)
    rv$single_trajectory_status <- list(available = TRUE, message = patient_meta$trajectory_message[1])
    if (!is.null(error_message)) {
      showNotification(paste0("Demo displayed using safe fallback. ", error_message), type = "warning", duration = 6)
    }
  }

  observeEvent(input$load_demo_low, {
    tryCatch({
      run_demo_prediction_from_row(FALSE)
      showNotification("Low-risk demo loaded and displayed.", type = "message", duration = 3)
    }, error = function(e) {
      demo_static_fallback(FALSE, paste0("Reason: ", conditionMessage(e)))
    })
  }, ignoreInit = TRUE)

  observeEvent(input$load_demo_high, {
    tryCatch({
      run_demo_prediction_from_row(TRUE)
      showNotification("High-risk demo loaded and displayed.", type = "message", duration = 3)
    }, error = function(e) {
      demo_static_fallback(TRUE, paste0("Reason: ", conditionMessage(e)))
    })
  }, ignoreInit = TRUE)

  run_single_prediction <- function() {
    ## v13: formal user-entered prediction always starts from the current UI
    ## values and ignores any previously displayed demo result.
    rv$single_input <- NULL
    rv$single_pred <- NULL
    rv$single_trajectory_status <- NULL

    local_state <- isolate({
      step05_obj <- rv$step05
      step03_obj <- rv$step03
      predictors <- current_predictors()
      traj_df <- make_trajectory_input_df_from_ui(input)
      traj_status <- auto_trajectory_status(traj_df, step03_obj)
      trajectory_auto <- NA_character_
      if (isTRUE(traj_status$available)) {
        trajectory_auto <- assign_trajectory_auto(traj_df, step03_obj)[1]
      }
      patient_df <- make_single_patient_df(input, predictors, trajectory_auto = trajectory_auto)

      model_nm <- NA_character_
      if (!is.null(step05_obj)) {
        model_nm <- input$primary_model %||% if ("Elastic net: dynamic + trajectory" %in% names(step05_obj$model_objects)) {
          "Elastic net: dynamic + trajectory"
        } else {
          names(step05_obj$model_objects)[1]
        }
      }

      list(
        step05_obj = step05_obj,
        step03_obj = step03_obj,
        predictors = predictors,
        traj_status = traj_status,
        trajectory_auto = trajectory_auto,
        patient_df = patient_df,
        model_nm = model_nm,
        age = input$age,
        sex = input$sex
      )
    })

    pred <- NULL
    if (!is.null(local_state$step05_obj) && !is.na(local_state$model_nm)) {
      pred <- tryCatch(
        predict_ami_edrs(local_state$patient_df, local_state$step05_obj, local_state$model_nm),
        error = function(e) {
          showNotification(
            paste0("Formal prediction engine error: ", conditionMessage(e), ". A demonstration calculation is shown instead."),
            type = "warning", duration = 8
          )
          NULL
        }
      )
    }

    if (is.null(pred)) {
      pred <- isolate(fallback_demo_prediction(input, trajectory_auto = local_state$trajectory_auto))
      if (is.null(local_state$step05_obj)) {
        showNotification(
          "Prediction engine is not loaded. A demonstration calculation is shown so the interface can be previewed.",
          type = "warning", duration = 6
        )
      }
    }

    patient_meta <- tibble::tibble(
      age = local_state$age,
      sex = local_state$sex,
      trajectory = local_state$trajectory_auto,
      trajectory_available = !is.na(local_state$trajectory_auto),
      trajectory_message = local_state$traj_status$message,
      trajectory_interpretation = trajectory_interpretation(local_state$trajectory_auto)
    )

    rv$single_input <- local_state$patient_df
    rv$single_pred <- bind_cols(patient_meta, pred)
    rv$single_trajectory_status <- local_state$traj_status
  }

  safe_user_input_fallback <- function(error_message = NULL) {
    local_state <- isolate({
      step03_obj <- rv$step03
      dat <- tryCatch(make_raw_patient_features_from_ui(input), error = function(e) NULL)
      if (is.null(dat)) {
        dat <- as.data.frame(matrix(NA, nrow = 1, ncol = length(public_feature_columns())))
        names(dat) <- public_feature_columns()
      }
      dat <- tryCatch(add_module_deltas(dat), error = function(e) dat)

      traj_status <- tryCatch(auto_trajectory_status(dat, step03_obj), error = function(e) list(available = FALSE, message = paste0("C1–C5 could not be assigned: ", conditionMessage(e))))
      trajectory_auto <- NA_character_
      if (!is.null(step03_obj) && isTRUE(traj_status$available)) {
        trajectory_auto <- tryCatch(assign_trajectory_auto(dat, step03_obj)[1], error = function(e) NA_character_)
      }
      ## If the formal trajectory object is unavailable, use a transparent rule-based display phenotype.
      ## This keeps the web app usable rather than disconnecting the session.
      if (is.na(trajectory_auto)) {
        modules <- c("hemodynamic", "respiratory_neurologic", "renal_metabolic", "inflammatory_coagulopathic", "cardiac_hypoperfusion")
        t48_cols <- paste0(modules, "_T48_score")
        if (all(t48_cols %in% names(dat))) {
          mean48 <- mean(suppressWarnings(as.numeric(unlist(dat[1, t48_cols, drop = TRUE]))), na.rm = TRUE)
          trajectory_auto <- dplyr::case_when(
            !is.finite(mean48) ~ NA_character_,
            mean48 < 25 ~ "C1",
            mean48 < 40 ~ "C2",
            mean48 < 55 ~ "C3",
            mean48 < 70 ~ "C4",
            TRUE ~ "C5"
          )
          if (!is.na(trajectory_auto)) {
            traj_status <- list(available = TRUE, message = "C1–C5 phenotype displayed from automatically calculated organ-domain dynamics.")
          }
        }
      }

      predictors <- unique(c(public_feature_columns(), names(dat), "trajectory"))
      patient_df <- as.data.frame(matrix(NA, nrow = 1, ncol = length(predictors)))
      names(patient_df) <- predictors
      for (v in intersect(names(dat), names(patient_df))) patient_df[[v]][1] <- dat[[v]][1]
      patient_df$trajectory[1] <- trajectory_auto

      list(
        patient_df = patient_df,
        trajectory_auto = trajectory_auto,
        traj_status = traj_status,
        age = input$age,
        sex = input$sex
      )
    })

    pred <- tryCatch(
      fallback_demo_prediction_df(local_state$patient_df, trajectory_auto = local_state$trajectory_auto),
      error = function(e) tibble::tibble(
        model = "Demonstration calculation",
        AMI_EDRS_probability = 0.12,
        AMI_EDRS_0_100 = 12.0,
        AMI_EDRS_risk_group = "Intermediate",
        AMI_EDRS_risk_group_public = "Intermediate risk (5–15%)",
        monitoring_intensity = "Enhanced observation"
      )
    )

    patient_meta <- tibble::tibble(
      age = local_state$age,
      sex = local_state$sex,
      trajectory = local_state$trajectory_auto,
      trajectory_available = !is.na(local_state$trajectory_auto),
      trajectory_message = local_state$traj_status$message,
      trajectory_interpretation = trajectory_interpretation(local_state$trajectory_auto)
    )

    rv$single_input <- local_state$patient_df
    rv$single_pred <- dplyr::bind_cols(patient_meta, pred)
    rv$single_trajectory_status <- local_state$traj_status

    if (!is.null(error_message)) {
      showNotification(paste0("Formal prediction failed; a safe demonstration calculation is shown. Reason: ", error_message), type = "warning", duration = 8)
    }
  }

  observeEvent(input$predict_single, {
    tryCatch({
      run_single_prediction()
    }, error = function(e) {
      safe_user_input_fallback(conditionMessage(e))
    })
  }, ignoreInit = TRUE)

  output$single_results_cards <- renderUI({
    pred <- rv$single_pred
    if (is.null(pred)) {
      return(tags$div(class = "ami-card", tags$h4("No prediction yet"), tags$p("Click a demo button to display results immediately, or enter patient data and run AMI-EDRS prediction.")))
    }

    p <- pred$AMI_EDRS_probability[1]
    rg <- pred$AMI_EDRS_risk_group[1]
    rg_pub <- pred$AMI_EDRS_risk_group_public[1]
    col <- risk_color(rg)
    traj <- pred$trajectory[1]
    traj_text <- ifelse(is.na(traj), "Unavailable", traj)
    traj_sub <- ifelse(is.na(traj), pred$trajectory_message[1], pred$trajectory_interpretation[1])
    traj_col <- ifelse(is.na(traj), "#7A5A13", "#1F3A5F")

    fluidRow(
      column(
        3,
        tags$div(
          class = "value-card",
          tags$div(class = "value-title", "Predicted 30-day mortality after T48"),
          tags$div(class = "value-number", style = paste0("color:", col, ";"), pretty_pct(p)),
          tags$div(class = "value-subtitle", paste0("AMI-EDRS = ", pred$AMI_EDRS_0_100[1], "/100"))
        )
      ),
      column(
        3,
        tags$div(
          class = "value-card",
          tags$div(class = "value-title", "Risk stratum"),
          tags$div(class = "risk-pill", style = paste0("background:", col, ";"), rg_pub),
          tags$div(class = "value-subtitle", pred$monitoring_intensity[1])
        )
      ),
      column(
        3,
        tags$div(
          class = "value-card",
          tags$div(class = "value-title", "C1–C5 trajectory phenotype"),
          tags$div(class = "value-number", style = paste0("color:", traj_col, ";"), traj_text),
          tags$div(class = "value-subtitle", traj_sub)
        )
      ),
      column(
        3,
        tags$div(
          class = "value-card",
          tags$div(class = "value-title", "Clinical interpretation"),
          tags$div(style = "font-size:18px;font-weight:850;margin-top:8px;color:#1F3A5F;", pred$monitoring_intensity[1]),
          tags$div(class = "small-note", "Use alongside clinician judgment and local protocols.")
        )
      )
    )
  })

  output$risk_gauge_plot <- renderPlot({
    tryCatch({
      pred <- rv$single_pred
      req(pred)
      risk_gauge_plot(pred$AMI_EDRS_probability[1])
    }, error = function(e) {
      ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = paste0("Risk gauge unavailable: ", conditionMessage(e)), size = 5)
    })
  })

  output$organ_domain_plot <- renderPlot({
    tryCatch({
      df <- rv$single_input
      req(df)
      organ_domain_plot_func(df)
    }, error = function(e) {
      ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = paste0("Organ-domain plot unavailable: ", conditionMessage(e)), size = 5)
    })
  })

  output$single_summary_table <- renderUI({
    pred <- rv$single_pred
    req(pred)
    out <- pred %>%
      dplyr::transmute(
        Age = age,
        Sex = sex,
        `C1-C5 trajectory` = ifelse(is.na(trajectory), "Unavailable", trajectory),
        `Trajectory status` = trajectory_message,
        `Trajectory interpretation` = trajectory_interpretation,
        `Predicted risk` = pretty_pct(AMI_EDRS_probability),
        `AMI-EDRS score` = AMI_EDRS_0_100,
        `Risk group` = AMI_EDRS_risk_group_public,
        `Monitoring intensity` = monitoring_intensity
      )

    ## HTML table instead of DT for the single-patient summary.
    ## This avoids DataTables Ajax errors on shinyapps.io when demo mode is clicked repeatedly.
    tags$table(
      class = "single-summary-table",
      tags$thead(tags$tr(lapply(names(out), tags$th))),
      tags$tbody(
        tags$tr(lapply(out[1, , drop = TRUE], function(x) tags$td(as.character(x))))
      )
    )
  })

  output$download_single_csv <- downloadHandler(
    filename = function() paste0("AMI_EDRS_single_patient_report_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) {
      req(rv$single_pred)
      readr::write_csv(rv$single_pred, file)
    }
  )

  safe_batch_fallback_prediction <- function(dat) {
    dat <- as.data.frame(dat)
    if (nrow(dat) == 0) return(tibble::tibble())
    purrr::map_dfr(seq_len(nrow(dat)), function(i) {
      fallback_demo_prediction_df(dat[i, , drop = FALSE], trajectory_auto = if ("trajectory" %in% names(dat)) as.character(dat$trajectory[i]) else NA_character_)
    })
  }

  observeEvent(input$run_batch, {
    tryCatch({
      req(input$batch_file$datapath)

      dat_raw <- read_uploaded_table(input$batch_file$datapath) %>% as.data.frame()
      names(dat_raw) <- make.names(names(dat_raw), unique = TRUE)

      dat <- compute_features_from_raw(dat_raw)
      dat <- add_module_deltas(dat)

      step03_obj <- isolate(rv$step03)
      step05_obj <- isolate(rv$step05)

      if (!"trajectory" %in% names(dat) && isTRUE(input$auto_assign_trajectory)) {
        dat$trajectory <- tryCatch(assign_trajectory_auto(dat, step03_obj), error = function(e) rep(NA_character_, nrow(dat)))
        dat$trajectory_status <- ifelse(
          is.na(dat$trajectory),
          "C1–C5 unavailable; complete T0/T24/T48 raw measurements are required.",
          "C1–C5 assigned from automatically calculated T0–T24–T48 organ-domain dynamics."
        )
      } else if (!"trajectory_status" %in% names(dat)) {
        dat$trajectory_status <- if ("trajectory" %in% names(dat)) "Trajectory provided by uploaded table." else "Trajectory unavailable."
      }

      pred <- NULL
      if (!is.null(step05_obj)) {
        model_nm <- tryCatch(isolate(selected_model()), error = function(e) NULL)
        pred <- tryCatch(predict_ami_edrs(dat, step05_obj, model_nm), error = function(e) {
          showNotification(paste0("Formal batch prediction failed; safe demonstration calculations are shown. Reason: ", conditionMessage(e)), type = "warning", duration = 8)
          NULL
        })
      }

      if (is.null(pred)) {
        pred <- safe_batch_fallback_prediction(dat)
        if (is.null(step05_obj)) {
          showNotification("Prediction engine is not loaded. Batch output uses safe demonstration calculations.", type = "warning", duration = 6)
        }
      }

      id_cols <- intersect(c("subject_id", "hadm_id", "stay_id", "patient_id", "trajectory", "trajectory_status"), names(dat))
      out <- dplyr::bind_cols(dat[, id_cols, drop = FALSE], pred)

      if ("trajectory" %in% names(out)) {
        out$trajectory_interpretation <- trajectory_interpretation(as.character(out$trajectory))
      }

      rv$batch_pred <- out
      showNotification(paste0("Batch prediction finished for ", nrow(out), " rows."), type = "message")
    }, error = function(e) {
      showNotification(paste0("Batch prediction failed: ", conditionMessage(e)), type = "error", duration = 10)
      rv$batch_pred <- tibble::tibble(
        status = "Batch prediction failed",
        message = conditionMessage(e)
      )
    })
  }, ignoreInit = TRUE)

  output$batch_status <- renderUI({
    if (is.null(rv$batch_pred)) {
      tags$div(class = "ami-card", tags$h4("No batch prediction yet"), tags$p("Upload a patient feature table and click Run batch prediction."))
    } else {
      n_unavailable <- if ("trajectory" %in% names(rv$batch_pred)) sum(is.na(rv$batch_pred$trajectory)) else NA_integer_
      tags$div(
        class = "ami-card",
        tags$h4("Batch prediction completed"),
        tags$p(paste0("Predicted rows: ", nrow(rv$batch_pred))),
        if (!is.na(n_unavailable)) tags$p(class = "small-note", paste0("Rows without assigned C1–C5 phenotype: ", n_unavailable))
      )
    }
  })

  output$batch_table <- renderUI({
    req(rv$batch_pred)
    df <- rv$batch_pred
    df <- as.data.frame(df)
    if (nrow(df) == 0) {
      return(tags$p(class = "small-note", "No rows to display."))
    }
    show_df <- head(df, 20)
    show_df[] <- lapply(show_df, function(x) {
      if (is.numeric(x)) formatC(x, digits = 4, format = "fg") else as.character(x)
    })
    tags$div(
      class = "simple-result-scroll",
      tags$table(
        class = "simple-result-table",
        tags$thead(tags$tr(lapply(names(show_df), tags$th))),
        tags$tbody(
          lapply(seq_len(nrow(show_df)), function(i) {
            tags$tr(lapply(show_df[i, , drop = TRUE], function(x) tags$td(as.character(x))))
          })
        )
      ),
      if (nrow(df) > 20) tags$p(class = "small-note", paste0("Showing first 20 of ", nrow(df), " rows. Download CSV for the full table."))
    )
  })

  output$batch_risk_distribution <- renderPlot({
    req(rv$batch_pred)
    df <- rv$batch_pred %>%
      dplyr::count(AMI_EDRS_risk_group, name = "n") %>%
      dplyr::mutate(
        AMI_EDRS_risk_group = factor(AMI_EDRS_risk_group, levels = names(pal_risk)),
        pct = n / sum(n) * 100
      )

    ggplot(df, aes(x = AMI_EDRS_risk_group, y = pct, fill = AMI_EDRS_risk_group)) +
      geom_col(width = 0.68) +
      geom_text(aes(label = paste0(sprintf("%.1f%%", pct), "\n(n=", n, ")")), vjust = -0.25, fontface = "bold") +
      scale_fill_manual(values = pal_risk, drop = FALSE) +
      scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.16))) +
      labs(title = "Risk-strata distribution", x = NULL, y = "% of uploaded patients") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 15, hjust = 1, face = "bold"), plot.title = element_text(face = "bold"))
  })

  output$batch_trajectory_distribution <- renderPlot({
    req(rv$batch_pred)
    if (!"trajectory" %in% names(rv$batch_pred) || all(is.na(rv$batch_pred$trajectory))) {
      plot.new(); text(0.5, 0.5, "No C1–C5 trajectory phenotype available.")
      return()
    }
    df <- rv$batch_pred %>%
      dplyr::filter(!is.na(trajectory)) %>%
      dplyr::count(trajectory, name = "n") %>%
      dplyr::mutate(pct = n / sum(n) * 100)

    ggplot(df, aes(x = trajectory, y = pct, fill = trajectory)) +
      geom_col(width = 0.68) +
      geom_text(aes(label = paste0(sprintf("%.1f%%", pct), "\n(n=", n, ")")), vjust = -0.25, fontface = "bold") +
      scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.15))) +
      labs(title = "C1–C5 trajectory distribution", x = "Trajectory phenotype", y = "% of uploaded patients") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "none", plot.title = element_text(face = "bold"))
  })

  output$download_batch <- downloadHandler(
    filename = function() paste0("AMI_EDRS_batch_predictions_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) {
      req(rv$batch_pred)
      readr::write_csv(rv$batch_pred, file)
    }
  )

  output$download_template <- downloadHandler(
    filename = function() "AMI_EDRS_full_feature_template.csv",
    content = function(file) {
      predictors <- isolate(current_predictors())
      template <- as.data.frame(matrix(NA, nrow = 1, ncol = length(predictors)))
      names(template) <- predictors
      template <- add_module_deltas(template)
      readr::write_csv(template, file)
    }
  )

  output$download_example <- downloadHandler(
    filename = function() "AMI_EDRS_example_batch_dataset.csv",
    content = function(file) {
      example <- make_example_data()
      readr::write_csv(example, file)
    }
  )

  output$download_example_single <- downloadHandler(
    filename = function() "AMI_EDRS_example_dataset.csv",
    content = function(file) {
      example <- make_example_data()
      readr::write_csv(example, file)
    }
  )

  output$model_details <- renderPrint({
    if (is.null(rv$step05)) {
      cat("Prediction engine: not loaded\n")
      cat("Please upload the prediction engine or place it at the configured local path.\n")
      return()
    }

    cat("Prediction engine: ready\n")
    cat("Primary target: 30-day mortality after the T48 landmark\n")
    cat("Default model: elastic-net trajectory-enhanced AMI-EDRS\n")
    cat("Trajectory phenotype: ")
    if (!is.null(rv$step03)) {
      cat("automatic assignment available when T0/T24/T48 organ-domain dynamics are complete\n")
    } else {
      cat("assigned from complete organ-domain dynamics when possible; otherwise reported as unavailable\n")
    }
    obj <- rv$step05$model_objects[[selected_model()]]
    cat("Required feature columns in batch mode: ", length(obj$prep$predictors), "\n", sep = "")
    cat("Risk strata: <5%, 5–15%, 15–30%, and ≥30% predicted 30-day mortality after T48\n")
  })
}

shinyApp(ui, server)

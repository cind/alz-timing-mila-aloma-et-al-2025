# R/biomarker_modelling.R

#====================#
# Biomarker trajectories modelling
#====================#

final_dataset_plasma <- read.csv("dataset_plasma.csv")
final_dataset_fbp <- read.csv("amyloid_dataset.csv")
tau_dataset <- read.csv("tau_dataset.csv")
csf_dataset <- read.csv("csf_dataset.csv")
CDR_dataset <- read.csv("CDR_dataset.csv")
MRI_dataset <- read.csv("mri_dataset.csv")

final_dataset_plasma$Centiloids_calc <- 300.66*final_dataset_plasma$SUVR_compositeRef  - 208.84
final_dataset_fbp$Centiloids_calc <- 300.66*final_dataset_fbp$SUVR_compositeRef - 208.84
tau_dataset$Centiloids_calc <- 300.66*tau_dataset$SUVR_compositeRef - 208.84
csf_dataset$Centiloids_calc <- 300.66*csf_dataset$SUVR_compositeRef - 208.84
CDR_dataset$Centiloids_calc <- 300.66*CDR_dataset$SUVR_compositeRef - 208.84
MRI_dataset$Centiloids_calc <- 300.66*MRI_dataset$SUVR_compositeRef - 208.84

#CDR_DX groups
final_dataset_plasma$CDR.x <-  as.factor(final_dataset_plasma$CDR.x)
final_dataset_plasma$DX_123 <-  as.factor(final_dataset_plasma$DX_123)

final_dataset_plasma <- final_dataset_plasma %>%
  mutate(CDR_DX_Group = case_when(
    CDR.x=="0" & DX_123=="1" ~ "CN_CDR0"    ,
    CDR.x!="0" & DX_123!="1" & nonAD_CI=="0" ~ "CI_AD_CDRgt0",
    CDR.x!="0" & DX_123!="1" & nonAD_CI=="1" ~ "CI_nonAD_CDRgt0",
    TRUE ~ NA_character_  
  ))

################################################
## Functions for derivatives of GAM(M) models ##

Deriv <- function(mod, n = 200, eps = 1e-7, newdata, term) {
  if(inherits(mod, "gamm"))
    mod <- mod$gam
  m.terms <- attr(terms(mod), "term.labels")
  if(missing(newdata)) {
    newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                   function(x) seq(min(x), max(x), length = n))
    names(newD) <- m.terms
  } else {
    newD <- newdata
  }
  X0 <- predict(mod, data.frame(newD), type = "lpmatrix")
  newD <- newD + eps
  X1 <- predict(mod, data.frame(newD), type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  Xp.r <- NROW(Xp)
  Xp.c <- NCOL(Xp)
  ## dims of bs
  bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
  ## number of smooth terms
  t.labs <- attr(mod$terms, "term.labels")
  ## match the term with the the terms in the model
  if(!missing(term)) {
    want <- grep(term, t.labs)
    if(!identical(length(want), length(term)))
      stop("One or more 'term's not found in model!")
    t.labs <- t.labs[want]
  }
  nt <- length(t.labs)
  ## list to hold the derivatives
  lD <- vector(mode = "list", length = nt)
  names(lD) <- t.labs
  for(i in seq_len(nt)) {
    Xi <- Xp * 0
    want <- grep(t.labs[i], colnames(X1))
    Xi[, want] <- Xp[, want]
    df <- Xi %*% coef(mod)
    df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
    lD[[i]] <- list(deriv = df, se.deriv = df.sd)
  }
  class(lD) <- "Deriv"
  lD$gamModel <- mod
  lD$eps <- eps
  lD$eval <- newD - eps
  lD ##return
}
Second_Deriv <- function(mod, n = 200, eps = 1e-7, newdata, term) {
  if(inherits(mod, "gamm"))
    mod <- mod$gam
  m.terms <- attr(terms(mod), "term.labels")
  if(missing(newdata)) {
    newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                   function(x) seq(min(x), max(x), length = n))
    names(newD) <- m.terms
  } else {
    newD <- newdata
  }
  X0 <- predict(mod, data.frame(newD), type = "lpmatrix")
  newD <- newD + eps
  X1 <- predict(mod, data.frame(newD), type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  Xp.r <- NROW(Xp)
  Xp.c <- NCOL(Xp)
  # second derivative
  newDFeps_m <- newD - 2*eps
  X_1 <- predict(mod, data.frame(newDFeps_m), type = 'lpmatrix')
  # design matrix for second derivative
  Xpp <- (X1 + X_1 - 2*X0)  / eps^2
  # second derivative
  fd_d2 <- Xpp %*% coef(mod)
  ## dims of bs
  bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
  ## number of smooth terms
  t.labs <- attr(mod$terms, "term.labels")
  ## match the term with the the terms in the model
  if(!missing(term)) {
    want <- grep(term, t.labs)
    if(!identical(length(want), length(term)))
      stop("One or more 'term's not found in model!")
    t.labs <- t.labs[want]
  }
  nt <- length(t.labs)
  ## list to hold the derivatives
  lD <- vector(mode = "list", length = nt)
  names(lD) <- t.labs
  for(i in seq_len(nt)) {
    Xi <- Xpp * 0
    want <- grep(t.labs[i], colnames(X1))
    Xi[, want] <- Xpp[, want]
    df <- Xpp %*% coef(mod)
    df.sd <- rowSums(Xpp %*% mod$Vp * Xpp)^.5
    lD[[i]] <- list(deriv = df, se.deriv = df.sd)
  }
  class(lD) <- "Deriv"
  lD$gamModel <- mod
  lD$eps <- eps
  lD$eval <- newD - eps
  lD ##return
}
confint.Deriv <- function(object, term, alpha = 0.05, ...) {
  l <- length(object) - 3
  term.labs <- names(object[seq_len(l)])
  if(missing(term)) {
    term <- term.labs
  } else { ## how many attempts to get this right!?!?
    ##term <- match(term, term.labs)
    ##term <- term[match(term, term.labs)]
    term <- term.labs[match(term, term.labs)]
  }
  if(any(miss <- is.na(term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  res <- vector(mode = "list", length = length(term))
  names(res) <- term
  residual.df <- df.residual(object$gamModel)
  tVal <- qt(1 - (alpha/2), residual.df)
  ##for(i in term.labs[term]) {
  for(i in term) {
    upr <- object[[i]]$deriv + tVal * object[[i]]$se.deriv
    lwr <- object[[i]]$deriv - tVal * object[[i]]$se.deriv
    res[[i]] <- list(upper = drop(upr), lower = drop(lwr))
  }
  res$alpha = alpha
  res
}

################################################
                   
# Calculation of mean and CI of biomarker values in the Reference group 
                   
###plasma biomarkers###
ref_group_subset <- subset(final_dataset_plasma , Ref_group==1)

stats_refgroup <- ref_group_subset %>%
  summarize(across(
    c(C2N_plasma_Abeta42_Abeta40_out, C2N_plasma_ptau217_out, C2N_plasma_ptau217_ratio, Fuji_plasma_ptau217_out,Fuji_plasma_Ab42_Ab40_out, AlzPath_plasma_ptau217_out,AlzPath_plasma_ptau217_out, Janssen_plasma_ptau217_out,Roche_plasma_GFAP_out,Roche_plasma_NfL_out, Roche_plasma_ptau181,Roche_plasma_Ab42_Ab40_out,QX_plasma_ptau181_out,QX_plasma_GFAP_out,QX_plasma_NfL_out,QX_plasma_NfL_out,QX_plasma_Ab42_Ab40_out),
    list(
      mean = ~mean(.x, na.rm = TRUE),
      sd = ~sd(.x, na.rm = TRUE),
      n = ~sum(!is.na(.x))
    )
  ))
alpha <- 0.05


stats_refgroup <- stats_refgroup %>%
  mutate(across(
    ends_with("_mean"),
    list(
      ci_lower = ~. - qt(1 - alpha / 2, df = stats_refgroup[[gsub("_mean", "_n", cur_column())]] - 1) * 
        stats_refgroup[[gsub("_mean", "_sd", cur_column())]] / 
        sqrt(stats_refgroup[[gsub("_mean", "_n", cur_column())]]),
      ci_upper = ~. + qt(1 - alpha / 2, df = stats_refgroup[[gsub("_mean", "_n", cur_column())]] - 1) * 
        stats_refgroup[[gsub("_mean", "_sd", cur_column())]] / 
        sqrt(stats_refgroup[[gsub("_mean", "_n", cur_column())]])
    ),
    .names = "{.col}_{.fn}"
  ))


final_dataset_plasma <- final_dataset_plasma %>%
  mutate(across(
    c(C2N_plasma_Abeta42_Abeta40_out, C2N_plasma_ptau217_out, C2N_plasma_ptau217_ratio, Fuji_plasma_ptau217_out,Fuji_plasma_Ab42_Ab40_out, AlzPath_plasma_ptau217_out,AlzPath_plasma_ptau217_out, Janssen_plasma_ptau217_out,Roche_plasma_GFAP_out,Roche_plasma_NfL_out, Roche_plasma_ptau181,Roche_plasma_Ab42_Ab40_out,QX_plasma_ptau181_out,QX_plasma_GFAP_out,QX_plasma_NfL_out,QX_plasma_NfL_out,QX_plasma_Ab42_Ab40_out),
    list(
      mean = ~stats_refgroup[[1, paste0(cur_column(), "_mean")]],
      sd = ~stats_refgroup[[1, paste0(cur_column(), "_sd")]],
      ci_lower = ~stats_refgroup[[1, paste0(cur_column(), "_mean_ci_lower")]],
      ci_upper = ~stats_refgroup[[1, paste0(cur_column(), "_mean_ci_upper")]]
    )
  ))


### reference biomarkers ##

add_stats_to_dataset <- function(data, variable_name, ref_group_column = "Ref_group", alpha = 0.05) {
  # Subset data based on Ref_group == 1
  ref_data <- subset(data, data[[ref_group_column]] == 1)
  
  # Extract the specific variable
  var_data <- ref_data[[variable_name]]
  
  # Calculate statistics
  var_mean <- mean(var_data, na.rm = TRUE)
  var_sd <- sd(var_data, na.rm = TRUE)
  n <- sum(!is.na(var_data))  # number of non-NA values
  se <- var_sd / sqrt(n)
  
  # Confidence interval calculation
  t_score <- qt(1 - alpha / 2, df = n - 1)
  ci_lower <- var_mean - t_score * se
  ci_upper <- var_mean + t_score * se
  
  # Add new columns to the dataset
  data <- data %>%
    mutate(!!paste0(variable_name, "_mean") := var_mean,
           !!paste0(variable_name, "_sd") := var_sd,
           !!paste0(variable_name, "_n") := n,
           !!paste0(variable_name, "_se") := se,
           !!paste0(variable_name, "_ci_lower") := ci_lower,
           !!paste0(variable_name, "_ci_upper") := ci_upper)
  
  return(data)
}

# Named list of datasets 
datasets <- list(
  final_dataset_fbp = final_dataset_fbp, 
  tau_dataset = tau_dataset,  
  csf_dataset = csf_dataset, 
  MRI_dataset = MRI_dataset, 
  CDR_dataset = CDR_dataset)

# Corresponding variables for each dataset
variables <- list(
  final_dataset_fbp = c("SUVR_compositeRef"),   
  tau_dataset = c("MesialTemporal", "TemporoParietal"),  # Multiple variables in tau_dataset
  csf_dataset = c("PTAU_over_ABETA42"),   
  MRI_dataset = c("metaROI.AgeAdj"),
  CDR_dataset = c("CDRSB", "MMSCORE") )

variables <- list(
  final_dataset_fbp = c("SUVR_compositeRef","MesialTemporal", "TemporoParietal","PTAU_over_ABETA42","metaROI.AgeAdj","CDRSB"))

# Loop over each dataset and its variables
for (dataset_name in names(datasets)) {
  dataset <- datasets[[dataset_name]]
  
  # Loop over each variable in the dataset
  for (variable in variables[[dataset_name]]) {
    
    # Add the stats to the dataset as new columns
    dataset <- add_stats_to_dataset(dataset, variable)
  }
  
  # Update the global variable directly with the modified dataset
  assign(dataset_name, dataset)
}
################################################
# =======GET ABNORMALITY TIMES VS REF GROUP==========


##---- code for Ab42/40 ----
get_tipping_point <- function(data, biomarker, z_value, time_variable) {
  formula <- as.formula(paste(biomarker, "~ s(", time_variable, ", k=3) + s(RID, bs = 're')"))
  model <- gam(formula, data = data, method = "REML")
  predictions <- predict(model, newdata = data, se.fit = TRUE)
  
  results <- data.frame(
    x = data[[time_variable]],
    fit = predictions$fit,
    se.fit = predictions$se.fit)
  results <- results %>% arrange(x) # Sort results by the time variable
  
  # Calculate confidence intervals
  results <- results %>%
    mutate(lower_ci = fit - z_value * se.fit,
           upper_ci = fit + z_value * se.fit)
  
  # Get the tipping point where the confidence intervals don't overlap
  non_overlap_points <- results$x[which(results$upper_ci < data[[paste0(biomarker, "_ci_lower")]])][1]
  
  return(non_overlap_points)
}

# Define parameters
n_bootstrap <- 1000
z_value <- qnorm(0.975)  # Z-value for 95% confidence interval

# Bootstrap process generalized for different biomarkers and time variables
run_bootstrap_for_biomarker <- function(data, biomarker, time_variable) {
  data <- data %>% filter(!is.na(.data[[time_variable]]))
  # Generate bootstrap samples
  set.seed(123)  # For reproducibility
  bootstrap_samples <- replicate(n_bootstrap, 
                                 data[sample(nrow(data), size = floor(0.8 * nrow(data)), replace = FALSE), ], 
                                 simplify = FALSE)
  
  # Apply the tipping point function to each bootstrap sample
  bootstrap_tipping_points <- map_dbl(bootstrap_samples, get_tipping_point, biomarker = biomarker, z_value = z_value, time_variable = time_variable)
  
  # Calculate median tipping point and confidence intervals
  bootstrap_tipping_points_median <- median(bootstrap_tipping_points, na.rm = TRUE)
  ci_lower <- quantile(bootstrap_tipping_points, 0.025, na.rm = TRUE)
  ci_upper <- quantile(bootstrap_tipping_points, 0.975, na.rm = TRUE)
  
  return(list(median = bootstrap_tipping_points_median, ci_lower = ci_lower, ci_upper = ci_upper))
}

biomarkers <- c("C2N_plasma_Abeta42_Abeta40_out", "Roche_plasma_Ab42_Ab40_out", "QX_plasma_Ab42_Ab40_out", "Fuji_plasma_Ab42_Ab40_out")
time_variables <- c("Centiloids_calc", "MesialTemporal")
#time_variables <- c("years_amy_onset", "years_tau_onset", "years_symp_onset", "years_tau_TP_onset")

results <- list()

# Apply the bootstrap process to each biomarker and time variable
for (biomarker in biomarkers) {
  for (time_variable in time_variables) {
    result_name <- paste(biomarker, time_variable, sep = "_")
    results[[result_name]] <- run_bootstrap_for_biomarker(final_dataset_plasma, biomarker, time_variable)
  }
}
results_df3 <- do.call(rbind, lapply(names(results), function(name) { # Convert the results list to a data frame
  # Extract biomarker and time variable from name
  parts <- unlist(strsplit(name, "_"))
  biomarker <- paste(parts[1:4], collapse = "_")  # Adjust based on your naming convention
  time_variable <- paste(parts[5:length(parts)], collapse = "_")
  
  res <- results[[name]]
  data.frame(
    Biomarker = biomarker,
    Time_Variable = time_variable,
    Median_Tipping_Point = res$median,
    CI_Lower = res$ci_lower,
    CI_Upper = res$ci_upper,
    stringsAsFactors = FALSE
  )
}))
write.csv(results_df3, file = "results_tipping_plasmaab1.csv", row.names = FALSE)

##---- code for rest of plasma biomarkers - increasing with pathology ----

get_tipping_point <- function(data, biomarker, z_value, time_variable) {
  formula <- as.formula(paste(biomarker, "~ s(", time_variable, ", k=3) + s(RID, bs = 're')"))
  model <- gam(formula, data = data, method = "REML")
  predictions <- predict(model, newdata = data, se.fit = TRUE)
  
  results <- data.frame(
    x = data[[time_variable]],
    fit = predictions$fit,
    se.fit = predictions$se.fit)
  results <- results %>% arrange(x) # Sort results by the time variable
  
  # Calculate confidence intervals
  results <- results %>%
    mutate(lower_ci = fit - z_value * se.fit,
           upper_ci = fit + z_value * se.fit)
  
  # Get the tipping point where the confidence intervals don't overlap
  non_overlap_points <- results$x[which(results$lower_ci > data[[paste0(biomarker, "_ci_upper")]])][1]
  
  return(non_overlap_points)
}

# Define parameters
n_bootstrap <- 1000
z_value <- qnorm(0.975)  # Z-value for 95% confidence interval

# Bootstrap process generalized for different biomarkers and time variables
run_bootstrap_for_biomarker <- function(data, biomarker, time_variable) {
  data <- data %>% filter(!is.na(.data[[time_variable]]))
  
  # Generate bootstrap samples
  set.seed(123)  # For reproducibility
  bootstrap_samples <- replicate(n_bootstrap, 
                                 data[sample(nrow(data), size = floor(0.8 * nrow(data)), replace = FALSE), ], 
                                 simplify = FALSE)
  
  # Apply the tipping point function to each bootstrap sample
  bootstrap_tipping_points <- map_dbl(bootstrap_samples, get_tipping_point, biomarker = biomarker, z_value = z_value, time_variable = time_variable)
  
  # Calculate median tipping point and confidence intervals
  bootstrap_tipping_points_median <- median(bootstrap_tipping_points, na.rm = TRUE)
  ci_lower <- quantile(bootstrap_tipping_points, 0.025, na.rm = TRUE)
  ci_upper <- quantile(bootstrap_tipping_points, 0.975, na.rm = TRUE)
  
  return(list(median = bootstrap_tipping_points_median, ci_lower = ci_lower, ci_upper = ci_upper))
}

biomarkers <- c("C2N_plasma_ptau217_out", "C2N_plasma_ptau217_ratio", "Fuji_plasma_ptau217_out", "AlzPath_plasma_ptau217_out", "Janssen_plasma_ptau217_out", 
                "Roche_plasma_GFAP_out", "Roche_plasma_NfL_out", "Roche_plasma_ptau181", "QX_plasma_ptau181_out","QX_plasma_GFAP_out", "QX_plasma_NfL_out"  )  

#time_variables <- c("years_amy_onset", "years_tau_onset", "years_symp_onset")
time_variables <- c("Centiloids_calc", "MesialTemporal")


results <- list()

# Apply the bootstrap process to each biomarker and time variable
for (biomarker in biomarkers) {
  for (time_variable in time_variables) {
    result_name <- paste(biomarker, time_variable, sep = "_")
    results[[result_name]] <- run_bootstrap_for_biomarker(final_dataset_plasma, biomarker, time_variable)
  }
}
results_df2 <- do.call(rbind, lapply(names(results), function(name) { # Convert the results list to a data frame
  # Extract biomarker and time variable from name
  parts <- unlist(strsplit(name, "_"))
  biomarker <- paste(parts[1:4], collapse = "_")  # Adjust based on your naming convention
  time_variable <- paste(parts[5:length(parts)], collapse = "_")
  
  res <- results[[name]]
  data.frame(
    Biomarker = biomarker,
    Time_Variable = time_variable,
    Median_Tipping_Point = res$median,
    CI_Lower = res$ci_lower,
    CI_Upper = res$ci_upper,
    stringsAsFactors = FALSE
  )
}))


##---- code for reference biomarkers that decrease with pathology ----
get_tipping_point <- function(data, biomarker, z_value, time_variable) {
  # Model formula with smoothing and random effect
  formula <- as.formula(paste(biomarker, "~ s(", time_variable, ", k=3) + s(RID, bs = 're')"))
  
  # Fit the GAM model
  model <- gam(formula, data = data, method = "REML")
  
  # Get model predictions and standard errors
  predictions <- predict(model, newdata = data, se.fit = TRUE)
  
  # Create results dataframe with x (time_variable), fit (predictions), and se.fit
  results <- data.frame(
    x = data[[time_variable]],
    fit = predictions$fit,
    se.fit = predictions$se.fit
  )
  
  # Sort results by the time variable (optional, for orderly inspection)
  results <- results %>% arrange(x)
  
  # Calculate confidence intervals
  results <- results %>%
    mutate(lower_ci = fit - z_value * se.fit,
           upper_ci = fit + z_value * se.fit)
  
  # Get the upper confidence interval of the observed biomarker values
  biomarker_ci_lower <- data[[paste0(biomarker, "_ci_lower")]]
  
  # Step 1: Identify all non-overlapping points
  non_overlap_points <- which(results$upper_ci < biomarker_ci_lower)
  
  # Step 2: Ensure continuous non-overlapping points
  # Find the first point where the confidence intervals no longer overlap,
  # and from which they continue to not overlap
  if (length(non_overlap_points) > 0) {
    for (i in seq_along(non_overlap_points)) {
      # Check if from this point forward, all subsequent points are non-overlapping
      if (all(non_overlap_points[i:length(non_overlap_points)] == seq(non_overlap_points[i], length.out = length(non_overlap_points) - i + 1))) {
        return(results$x[non_overlap_points[i]])  # Return the first "persistent" non-overlapping point
      }
    }
  }
  
  # If no continuous non-overlap point is found, return NA
  return(NA)
}

get_tipping_point <- function(data, biomarker, z_value, time_variable) {
  formula <- as.formula(paste(biomarker, "~ s(", time_variable, ", k=3) + s(RID, bs = 're')"))
  model <- gam(formula, data = data, method = "REML")
  predictions <- predict(model, newdata = data, se.fit = TRUE)
  
  results <- data.frame(
    x = data[[time_variable]],
    fit = predictions$fit,
    se.fit = predictions$se.fit)
  results <- results %>% arrange(x) # Sort results by the time variable
  
  # Calculate confidence intervals
  results <- results %>%
    mutate(lower_ci = fit - z_value * se.fit,
           upper_ci = fit + z_value * se.fit)
  
  # Get the tipping point where the confidence intervals don't overlap
  non_overlap_points <- results$x[which(results$upper_ci < data[[paste0(biomarker, "_ci_lower")]])][1]
  
  return(non_overlap_points)
}

# Define parameters
n_bootstrap <- 1000
z_value <- qnorm(0.975)  # Z-value for 95% confidence interval

# Define bootstrap process for different biomarkers and time variables
run_bootstrap_for_biomarker <- function(data, biomarker, time_variable) {
  # Remove rows with missing values in the specific time variable for this analysis
  data_filtered <- data %>% filter(!is.na(.data[[time_variable]]))
  
  # Remove rows with missing values in the biomarker (optional, depending on your needs)
  data_filtered <- data_filtered %>% filter(!is.na(.data[[biomarker]]))
  
  # Generate bootstrap samples
  set.seed(123)  # For reproducibility
  bootstrap_samples <- replicate(n_bootstrap, 
                                 data_filtered[sample(nrow(data_filtered), size = floor(0.8 * nrow(data_filtered)), replace = FALSE), ], 
                                 simplify = FALSE)
  
  # Apply the tipping point function to each bootstrap sample
  bootstrap_tipping_points <- map_dbl(bootstrap_samples, get_tipping_point, biomarker = biomarker, z_value = z_value, time_variable = time_variable)
  
  # Calculate median tipping point and confidence intervals
  bootstrap_tipping_points_median <- median(bootstrap_tipping_points, na.rm = TRUE)
  ci_lower <- quantile(bootstrap_tipping_points, 0.025, na.rm = TRUE)
  ci_upper <- quantile(bootstrap_tipping_points, 0.975, na.rm = TRUE)
  
  return(list(median = bootstrap_tipping_points_median, ci_lower = ci_lower, ci_upper = ci_upper))
}

# Define biomarkers, datasets, and time variables
datasets <- list(
  metaROI.AgeAdj = MRI_dataset,
  MMSCORE = CDR_dataset
)
datasets <- list(
  metaROI.AgeAdj = final_dataset_fbp
  
)

biomarkers <- c( "metaROI.AgeAdj")
time_variables <- c("years_amy_onset", "years_tau_onset", "years_symp_onset", "years_tau_TP_onset")
time_variables <- c("Centiloids_calc", "MesialTemporal")

results <- list()

# Apply the bootstrap process to each biomarker and time variable
for (biomarker in biomarkers) {
  dataset_name <- names(datasets)[which(biomarkers == biomarker)]
  data <- datasets[[dataset_name]]
  
  for (time_variable in time_variables) {
    result_name <- paste(biomarker, time_variable, sep = "_")
    results[[result_name]] <- run_bootstrap_for_biomarker(data, biomarker, time_variable)
  }
}

# Convert the results list to a data frame
results_dfref <- do.call(rbind, lapply(names(results), function(name) {
  # Extract biomarker and time variable from name
  parts <- unlist(strsplit(name, "_"))
  biomarker <- paste(parts[1:(length(parts) - 1)], collapse = "_")  # Adjust based on your naming convention
  time_variable <- parts[length(parts)]
  
  res <- results[[name]]
  data.frame(
    Biomarker = biomarker,
    Time_Variable = time_variable,
    Median_Tipping_Point = res$median,
    CI_Lower = res$ci_lower,
    CI_Upper = res$ci_upper,
    stringsAsFactors = FALSE
  )
}))

# Save the results to a CSV file
write.csv(results_dfref, file = "bootstrap_results_refdec1.csv", row.names = FALSE)



##---- code for reference biomarkers that increase with pathology ----
get_tipping_point <- function(data, biomarker, z_value, time_variable) {
  # Model formula with smoothing and random effect
  formula <- as.formula(paste(biomarker, "~ s(", time_variable, ", k=3) + s(RID, bs = 're')"))
  
  # Fit the GAM model
  model <- gam(formula, data = data, method = "REML")
  
  # Get model predictions and standard errors
  predictions <- predict(model, newdata = data, se.fit = TRUE)
  
  # Create results dataframe with x (time_variable), fit (predictions), and se.fit
  results <- data.frame(
    x = data[[time_variable]],
    fit = predictions$fit,
    se.fit = predictions$se.fit
  )
  
  # Sort results by the time variable (optional, for orderly inspection)
  results <- results %>% arrange(x)
  
  # Calculate confidence intervals
  results <- results %>%
    mutate(lower_ci = fit - z_value * se.fit,
           upper_ci = fit + z_value * se.fit)
  
  # Get the upper confidence interval of the observed biomarker values
  biomarker_ci_upper <- data[[paste0(biomarker, "_ci_upper")]]
  
  # Step 1: Identify all non-overlapping points
  non_overlap_points <- which(results$lower_ci > biomarker_ci_upper)
  
  # Step 2: Ensure continuous non-overlapping points
  # Find the first point where the confidence intervals no longer overlap,
  # and from which they continue to not overlap
  if (length(non_overlap_points) > 0) {
    for (i in seq_along(non_overlap_points)) {
      # Check if from this point forward, all subsequent points are non-overlapping
      if (all(non_overlap_points[i:length(non_overlap_points)] == seq(non_overlap_points[i], length.out = length(non_overlap_points) - i + 1))) {
        return(results$x[non_overlap_points[i]])  # Return the first "persistent" non-overlapping point
      }
    }
  }
  
  # If no continuous non-overlap point is found, return NA
  return(NA)
}

get_tipping_point <- function(data, biomarker, z_value, time_variable) {
  formula <- as.formula(paste(biomarker, "~ s(", time_variable, ", k=3) + s(RID, bs = 're')"))
  model <- gam(formula, data = data, method = "REML")
  predictions <- predict(model, newdata = data, se.fit = TRUE)
  
  results <- data.frame(
    x = data[[time_variable]],
    fit = predictions$fit,
    se.fit = predictions$se.fit)
  results <- results %>% arrange(x) # Sort results by the time variable
  
  # Calculate confidence intervals
  results <- results %>%
    mutate(lower_ci = fit - z_value * se.fit,
           upper_ci = fit + z_value * se.fit)
  
  # Get the tipping point where the confidence intervals don't overlap
  non_overlap_points <- results$x[which(results$lower_ci > data[[paste0(biomarker, "_ci_upper")]])][1]
  
  return(non_overlap_points)
}

# Define parameters
n_bootstrap <- 1000
z_value <- qnorm(0.975)  # Z-value for 95% confidence interval

# Define bootstrap process for different biomarkers and time variables
run_bootstrap_for_biomarker <- function(data, biomarker, time_variable) {
  # Remove rows with missing values in the specific time variable for this analysis
  data_filtered <- data %>% filter(!is.na(.data[[time_variable]]))
  
  # Remove rows with missing values in the biomarker (optional, depending on your needs)
  data_filtered <- data_filtered %>% filter(!is.na(.data[[biomarker]]))
  
  # Generate bootstrap samples
  set.seed(123)  # For reproducibility
  bootstrap_samples <- replicate(n_bootstrap, 
                                 data_filtered[sample(nrow(data_filtered), size = floor(0.8 * nrow(data_filtered)), replace = FALSE), ], 
                                 simplify = FALSE)
  
  # Apply the tipping point function to each bootstrap sample
  bootstrap_tipping_points <- map_dbl(bootstrap_samples, get_tipping_point, biomarker = biomarker, z_value = z_value, time_variable = time_variable)
  
  # Calculate median tipping point and confidence intervals
  bootstrap_tipping_points_median <- median(bootstrap_tipping_points, na.rm = TRUE)
  ci_lower <- quantile(bootstrap_tipping_points, 0.025, na.rm = TRUE)
  ci_upper <- quantile(bootstrap_tipping_points, 0.975, na.rm = TRUE)
  
  return(list(median = bootstrap_tipping_points_median, ci_lower = ci_lower, ci_upper = ci_upper))
}

# Define biomarkers, datasets, and time variables
datasets <- list(
SUVR_compositeRef = final_dataset_fbp,
  MesialTemporal = tau_dataset,
 TemporoParietal = tau_dataset,
 PTAU_over_ABETA42 = csf_dataset,
  CDRSB = CDR_dataset )

datasets <- list(
  SUVR_compositeRef = final_dataset_fbp,
  MesialTemporal = final_dataset_fbp,
  TemporoParietal = final_dataset_fbp,
  PTAU_over_ABETA42 = final_dataset_fbp,
  CDRSB = final_dataset_fbp )
biomarkers <- c("SUVR_compositeRef", "MesialTemporal", "TemporoParietal",  "PTAU_over_ABETA42", "CDRSB")
time_variables <- c("years_amy_onset", "years_tau_onset", "years_symp_onset")
time_variables <- c("Centiloids_calc", "MesialTemporal")


results <- list()

# Apply the bootstrap process to each biomarker and time variable
for (biomarker in biomarkers) {
  dataset_name <- names(datasets)[which(biomarkers == biomarker)]
  data <- datasets[[dataset_name]]
  
  for (time_variable in time_variables) {
    result_name <- paste(biomarker, time_variable, sep = "_")
    results[[result_name]] <- run_bootstrap_for_biomarker(data, biomarker, time_variable)
  }
}

# Convert the results list to a data frame
results_dfref2 <- do.call(rbind, lapply(names(results), function(name) {
  # Extract biomarker and time variable from name
  parts <- unlist(strsplit(name, "_"))
  biomarker <- paste(parts[1:(length(parts) - 1)], collapse = "_")  # Adjust based on your naming convention
  time_variable <- parts[length(parts)]
  
  res <- results[[name]]
  data.frame(
    Biomarker = biomarker,
    Time_Variable = time_variable,
    Median_Tipping_Point = res$median,
    CI_Lower = res$ci_lower,
    CI_Upper = res$ci_upper,
    stringsAsFactors = FALSE
  )
}))

# Save the results to a CSV file
write.csv(results_dfref2, file = "bootstrap_results_refinc1.csv", row.names = FALSE)

# Print the results data frame
print(results_df)

 ################################################
# =======IDENTIFY TIME PERIODS WITH SIGNIFICANT RATE OF CHANGE ==========

##---- plasma biomarkers that increase ----
plasma_biomarkers <- c("C2N_plasma_ptau217_out", "C2N_plasma_ptau217_ratio", "Fuji_plasma_ptau217_out", "AlzPath_plasma_ptau217_out", "Janssen_plasma_ptau217_out", 
                       "Roche_plasma_GFAP_out", "Roche_plasma_NfL_out", "Roche_plasma_ptau181", "QX_plasma_ptau181_out","QX_plasma_GFAP_out", "QX_plasma_NfL_out")  
time_variables <- c("years_amy_onset", "years_tau_onset", "years_symp_onset", "years_tau_TP_onset")

# Initialize a list to store the results
derivative_results <- list()

for (biomarker in plasma_biomarkers) {
  for (time_var in time_variables) {
    
    # Define the formula for the GAM model
    formula <- as.formula(paste(biomarker, "~ s(", time_var, ", k=3) + s(RID, bs = 're')"))
    
    # Fit the model
    model <- gam(formula, data = final_dataset_plasma, method = "REML")
    
    # Generate predictions
    predictions <- predict(model, newdata = final_dataset_plasma, se.fit = TRUE)
    
    results <- data.frame(
      x = final_dataset_plasma[[time_var]],
      fit = predictions$fit,
      se.fit = predictions$se.fit
    )
    results <- results %>% arrange(x)
    
    #### Derivatives to get the range of significant increase and acceleration ###
    first_derivative <- Deriv(model)
    
    # Use column index to access eval points
    col_index <- match(time_var, colnames(first_derivative$eval))
    if (!is.na(col_index)) {
      eval_points <- first_derivative$eval[, col_index]
      
      # Calculate confidence intervals for the derivative
      ci_deriv <- confint.Deriv(first_derivative)
      lower_bound <- ci_deriv[[time_var]]$lower
      upper_bound <- ci_deriv[[time_var]]$upper
      
      # Find the time points where the derivative is significantly positive (increase)
      significant_pos <- which(lower_bound > 0 & upper_bound > 0)
      
      if (length( significant_pos) > 0) {
        significant_timepoints <- eval_points[ significant_pos]
      } else {
        significant_timepoints <- NA
      }
      
      # Store the results for this biomarker and time variable
      derivative_results[[paste(biomarker, time_var, sep = "_")]] <- data.frame(
        Biomarker = biomarker,
        Time_Variable = time_var,
        Significant_Timepoints = I(list(significant_timepoints))
      )
    } else {
      warning(paste("Time variable", time_var, "not found in first_derivative$eval"))
    }
  }
}



# Combine all results into a single data frame
final_derivative_results <- do.call(rbind, derivative_results)
write.csv(final_derivative_results, file = "derivatives_plasma_biomarkers_inc_new.csv", row.names = FALSE)

##---- plasma biomarkers that decrease ----

plasma_biomarkers <- c("C2N_plasma_Abeta42_Abeta40_out", "Roche_plasma_Ab42_Ab40_out", "QX_plasma_Ab42_Ab40_out", "Fuji_plasma_Ab42_Ab40_out")
time_variables <- c("years_amy_onset", "years_tau_onset", "years_symp_onset", "years_tau_TP_onset")

derivative_results <- list()

for (biomarker in plasma_biomarkers) {
  for (time_var in time_variables) {
    
    # Define the formula for the GAM model
    formula <- as.formula(paste(biomarker, "~ s(", time_var, ", k=3) + s(RID, bs = 're')"))
    
    # Fit the model
    model <- gam(formula, data = final_dataset_plasma, method = "REML")
    
    # Generate predictions
    predictions <- predict(model, newdata = final_dataset_plasma, se.fit = TRUE)
    
    results <- data.frame(
      x = final_dataset_plasma[[time_var]],
      fit = predictions$fit,
      se.fit = predictions$se.fit
    )
    results <- results %>% arrange(x)
    
    #### Derivatives to get the range of significant increase and acceleration ###
    first_derivative <- Deriv(model)
    
    # Use column index to access eval points
    col_index <- match(time_var, colnames(first_derivative$eval))
    if (!is.na(col_index)) {
      eval_points <- first_derivative$eval[, col_index]
      
      # Calculate confidence intervals for the derivative
      ci_deriv <- confint.Deriv(first_derivative)
      lower_bound <- ci_deriv[[time_var]]$lower
      upper_bound <- ci_deriv[[time_var]]$upper
      
      # Find the time points where the derivative is significantly positive (increase)
      significant_neg <- which(lower_bound < 0 & upper_bound < 0)
      
      if (length(significant_neg) > 0) {
        significant_timepoints <- eval_points[significant_neg]
      } else {
        significant_timepoints <- NA
      }
      
      # Store the results for this biomarker and time variable
      derivative_results[[paste(biomarker, time_var, sep = "_")]] <- data.frame(
        Biomarker = biomarker,
        Time_Variable = time_var,
        Significant_Timepoints = I(list(significant_timepoints))
      )
    } else {
      warning(paste("Time variable", time_var, "not found in first_derivative$eval"))
    }
  }
}


# Combine all results into a single data frame
final_derivative_results <- do.call(rbind, derivative_results)
write.csv(final_derivative_results, file = "derivatives_plasma_biomarkers_dec_new.csv", row.names = FALSE)

##---- reference biomarkers that increase ----
datasets <- list(
  SUVR_compositeRef = final_dataset_fbp,
  MesialTemporal = tau_dataset,
  TemporoParietal = tau_dataset,
  PTAU_over_ABETA42 = csf_dataset,
  CDRSB = CDR_dataset)

time_variables <- c("years_amy_onset", "years_tau_onset", "years_symp_onset", "years_tau_TP_onset")

# Initialize a list to store the results
derivative_results <- list()

for (biomarker in names(datasets)) {
  data <- datasets[[biomarker]]
  
  for (time_var in time_variables) {
    
    # Define the formula for the GAM model
    formula <- as.formula(paste(biomarker, "~ s(", time_var, ", k=3) + s(RID, bs = 're')"))
    
    # Fit the model
    model <- gam(formula, data = data, method = "REML")
    
    # Generate predictions
    predictions <- predict(model, newdata = data, se.fit = TRUE)
    
    results <- data.frame(
      x = data[[time_var]],
      fit = predictions$fit,
      se.fit = predictions$se.fit
    )
    results <- results %>% arrange(x)
    
    #### Derivatives to get the range of significant increase and acceleration ###
    first_derivative <- Deriv(model)
    
    # Use column index to access eval points
    col_index <- match(time_var, colnames(first_derivative$eval))
    if (!is.na(col_index)) {
      eval_points <- first_derivative$eval[, col_index]
      
      # Calculate confidence intervals for the derivative
      ci_deriv <- confint.Deriv(first_derivative)
      lower_bound <- ci_deriv[[time_var]]$lower
      upper_bound <- ci_deriv[[time_var]]$upper
      
      # Find the time points where the derivative is significantly positive (increase)
      significant_pos <- which(lower_bound > 0 & upper_bound > 0)
      
      
      if (length(significant_pos) > 0) {
        significant_timepoints <- eval_points[significant_pos]
      } else {
        significant_timepoints <- NA
      }
      
      # Store the results for this biomarker and time variable
      derivative_results[[paste(biomarker, time_var, sep = "_")]] <- data.frame(
        Biomarker = biomarker,
        Time_Variable = time_var,
        Significant_Timepoints = I(list(significant_timepoints))
      )
    } else {
      warning(paste("Time variable", time_var, "not found in first_derivative$eval"))
    }
  }
}


# Combine all results into a single data frame
final_derivative_results <- do.call(rbind, derivative_results)

# View the final results
print(final_derivative_results)

# Save the results to a CSV file
write.csv(final_derivative_results, file = "derivatives_biomarkers_all_times_rangesref_inc_new.csv", row.names = FALSE)
##---- reference biomarkers that decrease----
datasets <- list(
  metaROI.AgeAdj = MRI_dataset,
  MMSCORE = CDR_dataset)

time_variables <- c("years_amy_onset", "years_tau_onset", "years_symp_onset", "years_tau_TP_onset")

# Initialize a list to store the results
derivative_results <- list()

for (biomarker in names(datasets)) {
  data <- datasets[[biomarker]]
  
  for (time_var in time_variables) {
    
    # Define the formula for the GAM model
    formula <- as.formula(paste(biomarker, "~ s(", time_var, ", k=3) + s(RID, bs = 're')"))
    
    # Fit the model
    model <- gam(formula, data = data, method = "REML")
    
    # Generate predictions
    predictions <- predict(model, newdata = data, se.fit = TRUE)
    
    results <- data.frame(
      x = data[[time_var]],
      fit = predictions$fit,
      se.fit = predictions$se.fit
    )
    results <- results %>% arrange(x)
    
    #### Derivatives to get the range of significant increase and acceleration ###
    first_derivative <- Deriv(model)
    
    # Use column index to access eval points
    col_index <- match(time_var, colnames(first_derivative$eval))
    if (!is.na(col_index)) {
      eval_points <- first_derivative$eval[, col_index]
      
      # Calculate confidence intervals for the derivative
      ci_deriv <- confint.Deriv(first_derivative)
      lower_bound <- ci_deriv[[time_var]]$lower
      upper_bound <- ci_deriv[[time_var]]$upper
      
      # Find the time points where the derivative is significantly positive (increase)
      significant_neg <- which(lower_bound < 0 & upper_bound < 0)
      
      
      if (length(significant_neg) > 0) {
        significant_timepoints <- eval_points[significant_neg]
      } else {
        significant_timepoints <- NA
      }
      
      # Store the results for this biomarker and time variable
      derivative_results[[paste(biomarker, time_var, sep = "_")]] <- data.frame(
        Biomarker = biomarker,
        Time_Variable = time_var,
        Significant_Timepoints = I(list(significant_timepoints))
      )
    } else {
      warning(paste("Time variable", time_var, "not found in first_derivative$eval"))
    }
  }
}

# Combine all results into a single data frame
final_derivative_results <- do.call(rbind, derivative_results)

# View the final results
print(final_derivative_results)

# Save the results to a CSV file
write.csv(final_derivative_results, file = "derivatives_biomarkers_all_times_rangesref_decr2_new.csv", row.names = FALSE)

################################################
# =======TEST EFFECTS OF AGE, SEX, APOE in biomarkers in ref group at bl and correct values==========

##---- Test effect of covariates in the ref group ----
ref_group_subset  <- ref_group_subset  %>%
  arrange(RID,EXAMDATE)
ref_group_subset_bl <- distinct(ref_group_subset, RID, .keep_all = T)

#plasma bioamrkers

biomarkers <- c("C2N_plasma_Abeta42_Abeta40_out", "Roche_plasma_Ab42_Ab40_out", "QX_plasma_Ab42_Ab40_out", "Fuji_plasma_Ab42_Ab40_out", "C2N_plasma_ptau217_out", "C2N_plasma_ptau217_ratio", "Fuji_plasma_ptau217_out", "AlzPath_plasma_ptau217_out", "Janssen_plasma_ptau217_out", 
                "Roche_plasma_GFAP_out", "Roche_plasma_NfL_out", "Roche_plasma_ptau181", "QX_plasma_ptau181_out","QX_plasma_GFAP_out", "QX_plasma_NfL_out"  )  

# Initialize a list to store results for each biomarker
lm_results <- list()

# Loop over each biomarker to run the linear model
for (biomarker in biomarkers) {
  # Define the formula for the linear model (Age, Sex, and APOE as predictors)
  formula <- as.formula(paste(biomarker, "~ Age + Sex + APOE_binary")) 
  
  # Run the linear model
  model <- lm(formula, data = ref_group_subset_bl )
  
  # Store summary of the model (coefficients, p-values, R-squared, etc.)
  lm_summary <- summary(model)
  
  coefficients <- lm_summary$coefficients
  r_squared <- lm_summary$r.squared
  
  # Create a data frame to store results for this biomarker
  result_df <- data.frame(
    Biomarker = biomarker,
    Coefficients = rownames(coefficients),
    Estimate = coefficients[, "Estimate"],
    Std_Error = coefficients[, "Std. Error"],
    t_value = coefficients[, "t value"],
    p_value = coefficients[, "Pr(>|t|)"],
    R_squared = r_squared,
    stringsAsFactors = FALSE
  )
  
  # Add the result for this biomarker to the list
  lm_results[[biomarker]] <- result_df
}

# Combine all results into a single data frame
final_results_df <- do.call(rbind, lm_results)

# Save results to a CSV file
write.csv(final_results_df, file = "lm_covariates_plasma_biomarkers.csv", row.names = FALSE)

#reference biomarkers
ref_group_subset_amy <- subset(final_dataset_fbp, Ref_group==1)
ref_group_subset_amy  <- ref_group_subset_amy  %>%
  arrange(RID,EXAMDATE)
ref_group_subset_amy_bl <- distinct(ref_group_subset_amy, RID, .keep_all = T)

ref_group_subset_tau <- subset(tau_dataset, Ref_group==1)
ref_group_subset_tau  <- ref_group_subset_tau  %>%
  arrange(RID,SCANDATE)
ref_group_subset_tau_bl <- distinct(ref_group_subset_tau, RID, .keep_all = T)

ref_group_subset_csf <- subset(csf_dataset, Ref_group==1)
ref_group_subset_csf  <- ref_group_subset_csf  %>%
  arrange(RID,EXAMDATE)
ref_group_subset_csf_bl <- distinct(ref_group_subset_csf, RID, .keep_all = T)

ref_group_subset_CDR <- subset(CDR_dataset, Ref_group==1)
ref_group_subset_CDR  <- ref_group_subset_CDR  %>%
  arrange(RID,VISDATE)
ref_group_subset_CDR_bl <- distinct(ref_group_subset_CDR, RID, .keep_all = T)

biomarkers <- c("SUVR_compositeRef", "MesialTemporal","TemporoParietal","PTAU_over_ABETA42", "CDRSB", "MMSCORE")  

datasets <- list(
  SUVR_compositeRef = ref_group_subset_amy_bl,
  MesialTemporal = ref_group_subset_tau_bl,
  TemporoParietal = ref_group_subset_tau_bl,
  PTAU_over_ABETA42 = ref_group_subset_csf_bl,
  CDRSB = ref_group_subset_CDR_bl,
  MMSCORE = ref_group_subset_CDR_bl)
  
final_dataset_fbp <- final_dataset_fbp %>%
  rename(PTGENDER = Sex)  
# Initialize a list to store results for each biomarker
lm_results <- list()

for (biomarker in biomarkers) {
  # Determine the dataset associated with the current biomarker
  dataset <- datasets[[biomarker]]
  
  # Ensure the dataset is valid
  if (is.null(dataset)) {
    warning(paste("No dataset found for biomarker:", biomarker))
    next
  }
  
  # Define the formula for the linear model (Age, Sex, and APOE as predictors)
  formula <- as.formula(paste(biomarker, "~ Age + Sex + APOE_binary"))
  
  # Run the linear model
  model <- lm(formula, data = dataset)
  
  # Store summary of the model (coefficients, p-values, R-squared, etc.)
  lm_summary <- summary(model)
  
  coefficients <- lm_summary$coefficients
  r_squared <- lm_summary$r.squared
  
  # Create a data frame to store results for this biomarker
  result_df <- data.frame(
    Biomarker = biomarker,
    Coefficients = rownames(coefficients),
    Estimate = coefficients[, "Estimate"],
    Std_Error = coefficients[, "Std. Error"],
    t_value = coefficients[, "t value"],
    p_value = coefficients[, "Pr(>|t|)"],
    R_squared = r_squared,
    stringsAsFactors = FALSE
  )
  
  # Add the result for this biomarker to the list
  lm_results[[biomarker]] <- result_df
}

# Combine all results into a single data frame
final_results_df <- do.call(rbind, lm_results)
write.csv(final_results_df, file = "lm_covariates_refbiomarkers.csv", row.names = FALSE)

##---- correct by age affected biomarker values ----
#plasma biomarkers
biomarkers <- c( "Fuji_plasma_ptau217_out",  "AlzPath_plasma_ptau217_out", "Janssen_plasma_ptau217_out", 
                  "QX_plasma_ptau181_out","Roche_plasma_ptau181",
                  "Roche_plasma_GFAP_out", "QX_plasma_GFAP_out","Roche_plasma_NfL_out","QX_plasma_NfL_out")
# Calculate the mean age of the reference group
mean_AGEref <- mean(ref_group_subset_bl$Age, na.rm = TRUE)

age_effects <- data.frame(biomarker = biomarkers, age_effect = NA)

# Loop through each biomarker to fit the model and extract the age effect
for (biomarker in biomarkers) {
  fit <- lm(as.formula(paste(biomarker, "~ Age")), data = ref_group_subset_bl)
  age_effect <- coef(fit)["Age"]
  age_effects$age_effect[age_effects$biomarker == biomarker] <- age_effect
}

# Adjust the biomarker levels in the whole dataset
for (biomarker in biomarkers) {
  age_effect <- age_effects$age_effect[age_effects$biomarker == biomarker]
  adj_biomarker <- paste0(biomarker, "_adj")
  final_dataset_plasma[[adj_biomarker]] <- final_dataset_plasma[[biomarker]] -  age_effect * (final_dataset_plasma$Age - mean_AGEref)
}
#reference biomarkers
biomarker_dataset <- list("MesialTemporal" = list( "full"= "tau_dataset", "ref" = "ref_group_subset_tau_bl"),
                                  "TemporoParietal" = list("full"= "tau_dataset", "ref" = "ref_group_subset_tau_bl"))


# Create a data frame to store the age effects
age_effects <- data.frame(biomarker = names(biomarker_dataset), age_effect = NA)

# Loop through each biomarker, calculate the mean age for the reference group, and fit the model to get age effect
for (biomarker in names(biomarker_dataset)) {
  
  # Reference group dataset
  ref_group_dataset <- get(biomarker_dataset[[biomarker]]$ref)
  
  # Calculate the mean age of the reference group
  mean_AGEref <- mean(ref_group_dataset$Age, na.rm = TRUE)
  
  # Fit the model and extract age effect
  fit <- lm(as.formula(paste(biomarker, "~ Age")), data = ref_group_dataset)
  age_effect <- coef(fit)["Age"]
  
  # Store the age effect
  age_effects$age_effect[age_effects$biomarker == biomarker] <- age_effect
  
  # Full dataset
  full_dataset <- get(biomarker_dataset[[biomarker]]$full)
  
  # Adjust the biomarker levels in the full dataset
  adj_biomarker <- paste0(biomarker, "_adj")
  full_dataset[[adj_biomarker]] <- full_dataset[[biomarker]] - age_effect * (full_dataset$Age - mean_AGEref)
  
  # Update the dataset with adjusted biomarker levels
  assign(biomarker_dataset[[biomarker]]$full, full_dataset)
}

####### Calculation of mean and CI of age adjusted biomarker values in the Reference group #######
###plasma biomarkers###
ref_group_subset <- subset(final_dataset_plasma , Ref_group==1)
##without outliers##
stats_refgroup <- ref_group_subset %>%
  summarize(across(
    c(Fuji_plasma_ptau217_out_adj,AlzPath_plasma_ptau217_out_adj,Janssen_plasma_ptau217_out_adj,Roche_plasma_GFAP_out_adj,Roche_plasma_NfL_out_adj, Roche_plasma_ptau181_adj,QX_plasma_ptau181_out_adj,QX_plasma_GFAP_out_adj,QX_plasma_NfL_out_adj),
    list(
      mean = ~mean(.x, na.rm = TRUE),
      sd = ~sd(.x, na.rm = TRUE),
      n = ~sum(!is.na(.x))
    )
  ))
alpha <- 0.05


stats_refgroup <- stats_refgroup %>%
  mutate(across(
    ends_with("_mean"),
    list(
      ci_lower = ~. - qt(1 - alpha / 2, df = stats_refgroup[[gsub("_mean", "_n", cur_column())]] - 1) * 
        stats_refgroup[[gsub("_mean", "_sd", cur_column())]] / 
        sqrt(stats_refgroup[[gsub("_mean", "_n", cur_column())]]),
      ci_upper = ~. + qt(1 - alpha / 2, df = stats_refgroup[[gsub("_mean", "_n", cur_column())]] - 1) * 
        stats_refgroup[[gsub("_mean", "_sd", cur_column())]] / 
        sqrt(stats_refgroup[[gsub("_mean", "_n", cur_column())]])
    ),
    .names = "{.col}_{.fn}"
  ))


final_dataset_plasma <- final_dataset_plasma %>%
  mutate(across(
    c(Fuji_plasma_ptau217_out_adj,AlzPath_plasma_ptau217_out_adj, Janssen_plasma_ptau217_out_adj,Roche_plasma_GFAP_out_adj,Roche_plasma_NfL_out_adj, Roche_plasma_ptau181_adj,QX_plasma_ptau181_out_adj,QX_plasma_GFAP_out_adj,QX_plasma_NfL_out_adj),
    list(
      mean = ~stats_refgroup[[1, paste0(cur_column(), "_mean")]],
      sd = ~stats_refgroup[[1, paste0(cur_column(), "_sd")]],
      ci_lower = ~stats_refgroup[[1, paste0(cur_column(), "_mean_ci_lower")]],
      ci_upper = ~stats_refgroup[[1, paste0(cur_column(), "_mean_ci_upper")]]
    )
  ))


### reference biomarkers ##

add_stats_to_dataset <- function(data, variable_name, ref_group_column = "Ref_group", alpha = 0.05) {
  # Subset data based on Ref_group == 1
  ref_data <- subset(data, data[[ref_group_column]] == 1)
  
  # Extract the specific variable
  var_data <- ref_data[[variable_name]]
  
  # Calculate statistics
  var_mean <- mean(var_data, na.rm = TRUE)
  var_sd <- sd(var_data, na.rm = TRUE)
  n <- sum(!is.na(var_data))  # number of non-NA values
  se <- var_sd / sqrt(n)
  
  # Confidence interval calculation
  t_score <- qt(1 - alpha / 2, df = n - 1)
  ci_lower <- var_mean - t_score * se
  ci_upper <- var_mean + t_score * se
  
  # Add new columns to the dataset
  data <- data %>%
    mutate(!!paste0(variable_name, "_mean") := var_mean,
           !!paste0(variable_name, "_sd") := var_sd,
           !!paste0(variable_name, "_n") := n,
           !!paste0(variable_name, "_se") := se,
           !!paste0(variable_name, "_ci_lower") := ci_lower,
           !!paste0(variable_name, "_ci_upper") := ci_upper)
  
  return(data)
}

# Named list of datasets (using variables directly)
datasets <- list(
  
  tau_dataset = tau_dataset)

# Corresponding variables for each dataset
variables <- list(tau_dataset = c("MesialTemporal_adj", "TemporoParietal_adj"))  # Multiple variables in tau_dataset)

# Loop over each dataset and its variables
for (dataset_name in names(datasets)) {
  dataset <- datasets[[dataset_name]]
  
  # Loop over each variable in the dataset
  for (variable in variables[[dataset_name]]) {
    
    # Add the stats to the dataset as new columns
    dataset <- add_stats_to_dataset(dataset, variable)
  }
  
  # Update the global variable directly with the modified dataset
  assign(dataset_name, dataset)
}

################################################
##----Tipping point for age adjusted plasma biomarkers (increasing)----
biomarkers <- c( "Fuji_plasma_ptau217_out_adj",  "AlzPath_plasma_ptau217_out_adj", "Janssen_plasma_ptau217_out_adj", 
                 "QX_plasma_ptau181_out_adj","Roche_plasma_ptau181_adj",
                 "Roche_plasma_GFAP_out_adj", "QX_plasma_GFAP_out_adj","Roche_plasma_NfL_out_adj","QX_plasma_NfL_out_adj")


get_tipping_point <- function(data, biomarker, z_value, time_variable) {
  formula <- as.formula(paste(biomarker, "~ s(", time_variable, ", k=3) + s(RID, bs = 're')"))
  model <- gam(formula, data = data, method = "REML")
  predictions <- predict(model, newdata = data, se.fit = TRUE)
  
  results <- data.frame(
    x = data[[time_variable]],
    fit = predictions$fit,
    se.fit = predictions$se.fit)
  results <- results %>% arrange(x) # Sort results by the time variable
  
  # Calculate confidence intervals
  results <- results %>%
    mutate(lower_ci = fit - z_value * se.fit,
           upper_ci = fit + z_value * se.fit)
  
  # Get the tipping point where the confidence intervals don't overlap
  non_overlap_points <- results$x[which(results$lower_ci > data[[paste0(biomarker, "_ci_upper")]])][1]
  
  return(non_overlap_points)
}

# Define parameters
n_bootstrap <- 1000
z_value <- qnorm(0.975)  # Z-value for 95% confidence interval

# Bootstrap process generalized for different biomarkers and time variables
run_bootstrap_for_biomarker <- function(data, biomarker, time_variable) {
  data <- data %>% filter(!is.na(.data[[time_variable]]))
  
  # Generate bootstrap samples
  set.seed(123)  # For reproducibility
  bootstrap_samples <- replicate(n_bootstrap, 
                                 data[sample(nrow(data), size = floor(0.8 * nrow(data)), replace = FALSE), ], 
                                 simplify = FALSE)
  
  # Apply the tipping point function to each bootstrap sample
  bootstrap_tipping_points <- map_dbl(bootstrap_samples, get_tipping_point, biomarker = biomarker, z_value = z_value, time_variable = time_variable)
  
  # Calculate median tipping point and confidence intervals
  bootstrap_tipping_points_median <- median(bootstrap_tipping_points, na.rm = TRUE)
  ci_lower <- quantile(bootstrap_tipping_points, 0.025, na.rm = TRUE)
  ci_upper <- quantile(bootstrap_tipping_points, 0.975, na.rm = TRUE)
  
  return(list(median = bootstrap_tipping_points_median, ci_lower = ci_lower, ci_upper = ci_upper))
}


time_variables <- c("years_amy_onset", "years_tau_onset", "years_symp_onset")

results <- list()

# Apply the bootstrap process to each biomarker and time variable
for (biomarker in biomarkers) {
  for (time_variable in time_variables) {
    result_name <- paste(biomarker, time_variable, sep = "_")
    results[[result_name]] <- run_bootstrap_for_biomarker(final_dataset_plasma, biomarker, time_variable)
  }
}
results_df2 <- do.call(rbind, lapply(names(results), function(name) { # Convert the results list to a data frame
  # Extract biomarker and time variable from name
  parts <- unlist(strsplit(name, "_"))
  biomarker <- paste(parts[1:4], collapse = "_")  # Adjust based on your naming convention
  time_variable <- paste(parts[5:length(parts)], collapse = "_")
  
  res <- results[[name]]
  data.frame(
    Biomarker = biomarker,
    Time_Variable = time_variable,
    Median_Tipping_Point = res$median,
    CI_Lower = res$ci_lower,
    CI_Upper = res$ci_upper,
    stringsAsFactors = FALSE
  )
}))
write.csv(results_df2, file = "results_tipping_plasma_adj.csv", row.names = FALSE)

#for Nfl (continious non overlap)

get_tipping_point <- function(data, biomarker, z_value, time_variable) {
  # Model formula with smoothing and random effect
  formula <- as.formula(paste(biomarker, "~ s(", time_variable, ", k=3) + s(RID, bs = 're')"))
  
  # Fit the GAM model
  model <- gam(formula, data = data, method = "REML")
  
  # Get model predictions and standard errors
  predictions <- predict(model, newdata = data, se.fit = TRUE)
  
  # Create results dataframe with x (time_variable), fit (predictions), and se.fit
  results <- data.frame(
    x = data[[time_variable]],
    fit = predictions$fit,
    se.fit = predictions$se.fit
  )
  
  # Sort results by the time variable (optional, for orderly inspection)
  results <- results %>% arrange(x)
  
  # Calculate confidence intervals
  results <- results %>%
    mutate(lower_ci = fit - z_value * se.fit,
           upper_ci = fit + z_value * se.fit)
  
  # Get the upper confidence interval of the observed biomarker values
  biomarker_ci_upper <- data[[paste0(biomarker, "_ci_upper")]]
  
  # Step 1: Identify all non-overlapping points
  non_overlap_points <- which(results$lower_ci > biomarker_ci_upper)
  
  # Step 2: Ensure continuous non-overlapping points
  # Find the first point where the confidence intervals no longer overlap,
  # and from which they continue to not overlap
  if (length(non_overlap_points) > 0) {
    for (i in seq_along(non_overlap_points)) {
      # Check if from this point forward, all subsequent points are non-overlapping
      if (all(non_overlap_points[i:length(non_overlap_points)] == seq(non_overlap_points[i], length.out = length(non_overlap_points) - i + 1))) {
        return(results$x[non_overlap_points[i]])  # Return the first "persistent" non-overlapping point
      }
    }
  }
  
  # If no continuous non-overlap point is found, return NA
  return(NA)
}
# Define parameters
n_bootstrap <- 1000
z_value <- qnorm(0.975)  # Z-value for 95% confidence interval

# Bootstrap process generalized for different biomarkers and time variables
run_bootstrap_for_biomarker <- function(data, biomarker, time_variable) {
  data <- data %>% filter(!is.na(.data[[time_variable]]))
  
  # Generate bootstrap samples
  set.seed(123)  # For reproducibility
  bootstrap_samples <- replicate(n_bootstrap, 
                                 data[sample(nrow(data), size = floor(0.8 * nrow(data)), replace = FALSE), ], 
                                 simplify = FALSE)
  
  # Apply the tipping point function to each bootstrap sample
  bootstrap_tipping_points <- map_dbl(bootstrap_samples, get_tipping_point, biomarker = biomarker, z_value = z_value, time_variable = time_variable)
  
  # Calculate median tipping point and confidence intervals
  bootstrap_tipping_points_median <- median(bootstrap_tipping_points, na.rm = TRUE)
  ci_lower <- quantile(bootstrap_tipping_points, 0.025, na.rm = TRUE)
  ci_upper <- quantile(bootstrap_tipping_points, 0.975, na.rm = TRUE)
  
  return(list(median = bootstrap_tipping_points_median, ci_lower = ci_lower, ci_upper = ci_upper))
}

time_variables <- c("years_amy_onset", "years_tau_onset", "years_symp_onset")

results <- list()

# Apply the bootstrap process to each biomarker and time variable
for (biomarker in biomarkers) {
  for (time_variable in time_variables) {
    result_name <- paste(biomarker, time_variable, sep = "_")
    results[[result_name]] <- run_bootstrap_for_biomarker(final_dataset_plasma, biomarker, time_variable)
  }
}
results_df2 <- do.call(rbind, lapply(names(results), function(name) { # Convert the results list to a data frame
  # Extract biomarker and time variable from name
  parts <- unlist(strsplit(name, "_"))
  biomarker <- paste(parts[1:4], collapse = "_")  # Adjust based on your naming convention
  time_variable <- paste(parts[5:length(parts)], collapse = "_")
  
  res <- results[[name]]
  data.frame(
    Biomarker = biomarker,
    Time_Variable = time_variable,
    Median_Tipping_Point = res$median,
    CI_Lower = res$ci_lower,
    CI_Upper = res$ci_upper,
    stringsAsFactors = FALSE
  )
}))
write.csv(results_df2, file = "results_tipping_plasma_adj.csv", row.names = FALSE)

##----Tipping point for age adjusted ref biomarkers (increasing)----
biomarkers <- c( "MesialTemporal_adj", "TemporoParietal_adj")


get_tipping_point <- function(data, biomarker, z_value, time_variable) {
  formula <- as.formula(paste(biomarker, "~ s(", time_variable, ", k=3) + s(RID, bs = 're')"))
  model <- gam(formula, data = data, method = "REML")
  predictions <- predict(model, newdata = data, se.fit = TRUE)
  
  results <- data.frame(
    x = data[[time_variable]],
    fit = predictions$fit,
    se.fit = predictions$se.fit)
  results <- results %>% arrange(x) # Sort results by the time variable
  
  # Calculate confidence intervals
  results <- results %>%
    mutate(lower_ci = fit - z_value * se.fit,
           upper_ci = fit + z_value * se.fit)
  
  # Get the tipping point where the confidence intervals don't overlap
  non_overlap_points <- results$x[which(results$lower_ci > data[[paste0(biomarker, "_ci_upper")]])][1]
  
  return(non_overlap_points)
}

# Define parameters
n_bootstrap <- 1000
z_value <- qnorm(0.975)  # Z-value for 95% confidence interval

# Bootstrap process generalized for different biomarkers and time variables
run_bootstrap_for_biomarker <- function(data, biomarker, time_variable) {
  data <- data %>% filter(!is.na(.data[[time_variable]]))
  
  # Generate bootstrap samples
  set.seed(123)  # For reproducibility
  bootstrap_samples <- replicate(n_bootstrap, 
                                 data[sample(nrow(data), size = floor(0.8 * nrow(data)), replace = FALSE), ], 
                                 simplify = FALSE)
  
  # Apply the tipping point function to each bootstrap sample
  bootstrap_tipping_points <- map_dbl(bootstrap_samples, get_tipping_point, biomarker = biomarker, z_value = z_value, time_variable = time_variable)
  
  # Calculate median tipping point and confidence intervals
  bootstrap_tipping_points_median <- median(bootstrap_tipping_points, na.rm = TRUE)
  ci_lower <- quantile(bootstrap_tipping_points, 0.025, na.rm = TRUE)
  ci_upper <- quantile(bootstrap_tipping_points, 0.975, na.rm = TRUE)
  
  return(list(median = bootstrap_tipping_points_median, ci_lower = ci_lower, ci_upper = ci_upper))
}

time_variables <- c("years_amy_onset", "years_tau_onset", "years_symp_onset")

results <- list()

# Apply the bootstrap process to each biomarker and time variable
for (biomarker in biomarkers) {
  for (time_variable in time_variables) {
    result_name <- paste(biomarker, time_variable, sep = "_")
    results[[result_name]] <- run_bootstrap_for_biomarker(tau_dataset, biomarker, time_variable)
  }
}
results_df3 <- do.call(rbind, lapply(names(results), function(name) { # Convert the results list to a data frame
  # Extract biomarker and time variable from name
  parts <- unlist(strsplit(name, "_"))
  biomarker <- paste(parts[1:4], collapse = "_")  # Adjust based on your naming convention
  time_variable <- paste(parts[5:length(parts)], collapse = "_")
  
  res <- results[[name]]
  data.frame(
    Biomarker = biomarker,
    Time_Variable = time_variable,
    Median_Tipping_Point = res$median,
    CI_Lower = res$ci_lower,
    CI_Upper = res$ci_upper,
    stringsAsFactors = FALSE
  )
}))
write.csv(results_df3, file = "results_tipping_ref_adj.csv", row.names = FALSE)

### periods of significant change for age adjusted plasma biomarkers
plasma_biomarkers <- c("Fuji_plasma_ptau217_out_adj", "AlzPath_plasma_ptau217_out_adj", "Janssen_plasma_ptau217_out_adj", 
                       "Roche_plasma_GFAP_out_adj", "Roche_plasma_NfL_out_adj", "Roche_plasma_ptau181_adj", "QX_plasma_ptau181_out_adj","QX_plasma_GFAP_out_adj", "QX_plasma_NfL_out_adj")  
time_variables <- c("years_amy_onset", "years_tau_onset", "years_symp_onset")

# Initialize a list to store the results
derivative_results <- list()

for (biomarker in plasma_biomarkers) {
  for (time_var in time_variables) {
    
    # Define the formula for the GAM model
    formula <- as.formula(paste(biomarker, "~ s(", time_var, ", k=3) + s(RID, bs = 're')"))
    
    # Fit the model
    model <- gam(formula, data = final_dataset_plasma, method = "REML")
    
    # Generate predictions
    predictions <- predict(model, newdata = final_dataset_plasma, se.fit = TRUE)
    
    results <- data.frame(
      x = final_dataset_plasma[[time_var]],
      fit = predictions$fit,
      se.fit = predictions$se.fit
    )
    results <- results %>% arrange(x)
    
    #### Derivatives to get the range of significant increase and acceleration ###
    first_derivative <- Deriv(model)
    
    # Use column index to access eval points
    col_index <- match(time_var, colnames(first_derivative$eval))
    if (!is.na(col_index)) {
      eval_points <- first_derivative$eval[, col_index]
      
      # Calculate confidence intervals for the derivative
      ci_deriv <- confint.Deriv(first_derivative)
      lower_bound <- ci_deriv[[time_var]]$lower
      upper_bound <- ci_deriv[[time_var]]$upper
      
      # Find the time points where the derivative is significantly positive (increase)
      significant_pos <- which(lower_bound > 0 & upper_bound > 0)
      
      if (length( significant_pos) > 0) {
        significant_timepoints <- eval_points[ significant_pos]
      } else {
        significant_timepoints <- NA
      }
      
      # Store the results for this biomarker and time variable
      derivative_results[[paste(biomarker, time_var, sep = "_")]] <- data.frame(
        Biomarker = biomarker,
        Time_Variable = time_var,
        Significant_Timepoints = I(list(significant_timepoints))
      )
    } else {
      warning(paste("Time variable", time_var, "not found in first_derivative$eval"))
    }
  }
}



# Combine all results into a single data frame
final_derivative_results <- do.call(rbind, derivative_results)
write.csv(final_derivative_results, file = "derivatives_adj.csv", row.names = FALSE)


## ref biomarkers

plasma_biomarkers <- c("MesialTemporal_adj", "TemporoParietal_adj")  
time_variables <- c("years_amy_onset", "years_tau_onset", "years_symp_onset")

# Initialize a list to store the results
derivative_results <- list()

for (biomarker in plasma_biomarkers) {
  for (time_var in time_variables) {
    
    # Define the formula for the GAM model
    formula <- as.formula(paste(biomarker, "~ s(", time_var, ", k=3) + s(RID, bs = 're')"))
    
    # Fit the model
    model <- gam(formula, data = tau_dataset, method = "REML")
    
    # Generate predictions
    predictions <- predict(model, newdata = tau_dataset, se.fit = TRUE)
    
    results <- data.frame(
      x = tau_dataset[[time_var]],
      fit = predictions$fit,
      se.fit = predictions$se.fit
    )
    results <- results %>% arrange(x)
    
    #### Derivatives to get the range of significant increase and acceleration ###
    first_derivative <- Deriv(model)
    
    # Use column index to access eval points
    col_index <- match(time_var, colnames(first_derivative$eval))
    if (!is.na(col_index)) {
      eval_points <- first_derivative$eval[, col_index]
      
      # Calculate confidence intervals for the derivative
      ci_deriv <- confint.Deriv(first_derivative)
      lower_bound <- ci_deriv[[time_var]]$lower
      upper_bound <- ci_deriv[[time_var]]$upper
      
      # Find the time points where the derivative is significantly positive (increase)
      significant_pos <- which(lower_bound > 0 & upper_bound > 0)
      
      if (length( significant_pos) > 0) {
        significant_timepoints <- eval_points[ significant_pos]
      } else {
        significant_timepoints <- NA
      }
      
      # Store the results for this biomarker and time variable
      derivative_results[[paste(biomarker, time_var, sep = "_")]] <- data.frame(
        Biomarker = biomarker,
        Time_Variable = time_var,
        Significant_Timepoints = I(list(significant_timepoints))
      )
    } else {
      warning(paste("Time variable", time_var, "not found in first_derivative$eval"))
    }
  }
}



# Combine all results into a single data frame
final_derivative_results <- do.call(rbind, derivative_results)
write.csv(final_derivative_results, file = "derivatives_ref_adj.csv", row.names = FALSE)

################################################
########### sensitivity analyses with all plasma data available subsample n=331 #################

qx_subset <- final_dataset_plasma  %>%
 filter(!is.na(QX_plasma_Ab42_Ab40) & !is.na(Janssen_plasma_ptau217) & !is.na(C2N_plasma_ptau217_ratio) & !is.na(C2N_plasma_Abeta42_Abeta40) & !is.na(AlzPath_plasma_ptau217))
qx_subset_bl <- final_dataset_plasma_bl  %>%
  filter(!is.na(QX_plasma_Ab42_Ab40) & !is.na(Janssen_plasma_ptau217) & !is.na(C2N_plasma_ptau217_ratio) & !is.na(C2N_plasma_Abeta42_Abeta40) & !is.na(AlzPath_plasma_ptau217))

###plasma biomarkers###
ref_group_subset_all <- subset(qx_subset , Ref_group==1)

stats_refgroup_all <- ref_group_subset_all %>%
  summarize(across(
    c(C2N_plasma_Abeta42_Abeta40_out, C2N_plasma_ptau217_out, C2N_plasma_ptau217_ratio, Fuji_plasma_ptau217_out,Fuji_plasma_Ab42_Ab40_out, AlzPath_plasma_ptau217_out,AlzPath_plasma_ptau217_out, Janssen_plasma_ptau217_out,Roche_plasma_GFAP_out,Roche_plasma_NfL_out, Roche_plasma_ptau181,Roche_plasma_Ab42_Ab40_out,QX_plasma_ptau181_out,QX_plasma_GFAP_out,QX_plasma_NfL_out,QX_plasma_NfL_out,QX_plasma_Ab42_Ab40_out, PTAU_over_ABETA42,SUVR_compositeRef,MesialTemporal,atrophy,CDR_SOB),
    list(
      mean = ~mean(.x, na.rm = TRUE),
      sd = ~sd(.x, na.rm = TRUE),
      n = ~sum(!is.na(.x))
    )
  ))
alpha <- 0.05


stats_refgroup_all <- stats_refgroup_all %>%
  mutate(across(
    ends_with("_mean"),
    list(
      ci_lower = ~. - qt(1 - alpha / 2, df = stats_refgroup_all[[gsub("_mean", "_n", cur_column())]] - 1) * 
        stats_refgroup_all[[gsub("_mean", "_sd", cur_column())]] / 
        sqrt(stats_refgroup_all[[gsub("_mean", "_n", cur_column())]]),
      ci_upper = ~. + qt(1 - alpha / 2, df = stats_refgroup_all[[gsub("_mean", "_n", cur_column())]] - 1) * 
        stats_refgroup_all[[gsub("_mean", "_sd", cur_column())]] / 
        sqrt(stats_refgroup_all[[gsub("_mean", "_n", cur_column())]])
    ),
    .names = "{.col}_{.fn}"
  ))


qx_subset <- qx_subset %>%
  mutate(across(
    c(C2N_plasma_Abeta42_Abeta40_out, C2N_plasma_ptau217_out, C2N_plasma_ptau217_ratio, Fuji_plasma_ptau217_out,Fuji_plasma_Ab42_Ab40_out, AlzPath_plasma_ptau217_out,AlzPath_plasma_ptau217_out, Janssen_plasma_ptau217_out,Roche_plasma_GFAP_out,Roche_plasma_NfL_out, Roche_plasma_ptau181,Roche_plasma_Ab42_Ab40_out,QX_plasma_ptau181_out,QX_plasma_GFAP_out,QX_plasma_NfL_out,QX_plasma_NfL_out,QX_plasma_Ab42_Ab40_out, PTAU_over_ABETA42,SUVR_compositeRef,MesialTemporal,atrophy,CDR_SOB),
    list(
      mean = ~stats_refgroup_all[[1, paste0(cur_column(), "_mean")]],
      sd = ~stats_refgroup_all[[1, paste0(cur_column(), "_sd")]],
      ci_lower = ~stats_refgroup_all[[1, paste0(cur_column(), "_mean_ci_lower")]],
      ci_upper = ~stats_refgroup_all[[1, paste0(cur_column(), "_mean_ci_upper")]]
    )
  ))

##tipping point for biomarkers that decrease
get_tipping_point <- function(data, biomarker, z_value, time_variable) {
  formula <- as.formula(paste(biomarker, "~ s(", time_variable, ", k=3) + s(RID, bs = 're')"))
  model <- gam(formula, data = data, method = "REML")
  predictions <- predict(model, newdata = data, se.fit = TRUE)
  
  results <- data.frame(
    x = data[[time_variable]],
    fit = predictions$fit,
    se.fit = predictions$se.fit)
  results <- results %>% arrange(x) # Sort results by the time variable
  
  # Calculate confidence intervals
  results <- results %>%
    mutate(lower_ci = fit - z_value * se.fit,
           upper_ci = fit + z_value * se.fit)
  
  # Get the tipping point where the confidence intervals don't overlap
  non_overlap_points <- results$x[which(results$upper_ci < data[[paste0(biomarker, "_ci_lower")]])][1]
  
  return(non_overlap_points)
}

# Define parameters
n_bootstrap <- 1000
z_value <- qnorm(0.975)  # Z-value for 95% confidence interval

# Bootstrap process generalized for different biomarkers and time variables
run_bootstrap_for_biomarker <- function(data, biomarker, time_variable) {
  data <- data %>% filter(!is.na(.data[[time_variable]]))
  # Generate bootstrap samples
  set.seed(123)  # For reproducibility
  bootstrap_samples <- replicate(n_bootstrap, 
                                 data[sample(nrow(data), size = floor(0.8 * nrow(data)), replace = FALSE), ], 
                                 simplify = FALSE)
  
  # Apply the tipping point function to each bootstrap sample
  bootstrap_tipping_points <- map_dbl(bootstrap_samples, get_tipping_point, biomarker = biomarker, z_value = z_value, time_variable = time_variable)
  
  # Calculate median tipping point and confidence intervals
  bootstrap_tipping_points_median <- median(bootstrap_tipping_points, na.rm = TRUE)
  ci_lower <- quantile(bootstrap_tipping_points, 0.025, na.rm = TRUE)
  ci_upper <- quantile(bootstrap_tipping_points, 0.975, na.rm = TRUE)
  
  return(list(median = bootstrap_tipping_points_median, ci_lower = ci_lower, ci_upper = ci_upper))
}

biomarkers <- c("C2N_plasma_Abeta42_Abeta40_out", "Roche_plasma_Ab42_Ab40_out", "QX_plasma_Ab42_Ab40_out", "Fuji_plasma_Ab42_Ab40_out", "atrophy")
time_variables <- c("years_amy_onset", "years_tau_onset", "years_symp_onset")

results <- list()

# Apply the bootstrap process to each biomarker and time variable
for (biomarker in biomarkers) {
  for (time_variable in time_variables) {
    result_name <- paste(biomarker, time_variable, sep = "_")
    results[[result_name]] <- run_bootstrap_for_biomarker(qx_subset, biomarker, time_variable)
  }
}
results_df3 <- do.call(rbind, lapply(names(results), function(name) { # Convert the results list to a data frame
  # Extract biomarker and time variable from name
  parts <- unlist(strsplit(name, "_"))
  biomarker <- paste(parts[1:4], collapse = "_")  # Adjust based on your naming convention
  time_variable <- paste(parts[5:length(parts)], collapse = "_")
  
  res <- results[[name]]
  data.frame(
    Biomarker = biomarker,
    Time_Variable = time_variable,
    Median_Tipping_Point = res$median,
    CI_Lower = res$ci_lower,
    CI_Upper = res$ci_upper,
    stringsAsFactors = FALSE
  )
}))
write.csv(results_df3, file = "results_tipping_plasma_sensitivity.csv", row.names = FALSE)

##plasma biomarkers that increase
get_tipping_point <- function(data, biomarker, z_value, time_variable) {
  formula <- as.formula(paste(biomarker, "~ s(", time_variable, ", k=3) + s(RID, bs = 're')"))
  model <- gam(formula, data = data, method = "REML")
  predictions <- predict(model, newdata = data, se.fit = TRUE)
  
  results <- data.frame(
    x = data[[time_variable]],
    fit = predictions$fit,
    se.fit = predictions$se.fit)
  results <- results %>% arrange(x) # Sort results by the time variable
  
  # Calculate confidence intervals
  results <- results %>%
    mutate(lower_ci = fit - z_value * se.fit,
           upper_ci = fit + z_value * se.fit)
  
  # Get the tipping point where the confidence intervals don't overlap
  non_overlap_points <- results$x[which(results$lower_ci > data[[paste0(biomarker, "_ci_upper")]])][1]
  
  return(non_overlap_points)
}

# Define parameters
n_bootstrap <- 1000
z_value <- qnorm(0.975)  # Z-value for 95% confidence interval

# Bootstrap process generalized for different biomarkers and time variables
run_bootstrap_for_biomarker <- function(data, biomarker, time_variable) {
  data <- data %>% filter(!is.na(.data[[time_variable]]))
  
  # Generate bootstrap samples
  set.seed(123)  # For reproducibility
  bootstrap_samples <- replicate(n_bootstrap, 
                                 data[sample(nrow(data), size = floor(0.8 * nrow(data)), replace = FALSE), ], 
                                 simplify = FALSE)
  
  # Apply the tipping point function to each bootstrap sample
  bootstrap_tipping_points <- map_dbl(bootstrap_samples, get_tipping_point, biomarker = biomarker, z_value = z_value, time_variable = time_variable)
  
  # Calculate median tipping point and confidence intervals
  bootstrap_tipping_points_median <- median(bootstrap_tipping_points, na.rm = TRUE)
  ci_lower <- quantile(bootstrap_tipping_points, 0.025, na.rm = TRUE)
  ci_upper <- quantile(bootstrap_tipping_points, 0.975, na.rm = TRUE)
  
  return(list(median = bootstrap_tipping_points_median, ci_lower = ci_lower, ci_upper = ci_upper))
}


biomarkers <- c("C2N_plasma_ptau217_out", "C2N_plasma_ptau217_ratio", "Fuji_plasma_ptau217_out", "AlzPath_plasma_ptau217_out", "Janssen_plasma_ptau217_out", 
                "Roche_plasma_GFAP_out", "Roche_plasma_NfL_out", "Roche_plasma_ptau181", "QX_plasma_ptau181_out","QX_plasma_GFAP_out", "QX_plasma_NfL_out" ,"PTAU_over_ABETA42","SUVR_compositeRef","MesialTemporal","TemporoParietal", "CDR_SOB" )  
time_variables <- c("years_amy_onset", "years_tau_onset", "years_symp_onset")

results <- list()

# Apply the bootstrap process to each biomarker and time variable
for (biomarker in biomarkers) {
  for (time_variable in time_variables) {
    result_name <- paste(biomarker, time_variable, sep = "_")
    results[[result_name]] <- run_bootstrap_for_biomarker(qx_subset, biomarker, time_variable)
  }
}
results_df2 <- do.call(rbind, lapply(names(results), function(name) { # Convert the results list to a data frame
  # Extract biomarker and time variable from name
  parts <- unlist(strsplit(name, "_"))
  biomarker <- paste(parts[1:4], collapse = "_")  # Adjust based on your naming convention
  time_variable <- paste(parts[5:length(parts)], collapse = "_")
  
  res <- results[[name]]
  data.frame(
    Biomarker = biomarker,
    Time_Variable = time_variable,
    Median_Tipping_Point = res$median,
    CI_Lower = res$ci_lower,
    CI_Upper = res$ci_upper,
    stringsAsFactors = FALSE
  )
}))

write.csv(results_df2, file = "results_tipping_plasma_incsensitivity.csv", row.names = FALSE)


# =======IDENTIFY TIME PERIODS WITH SIGNIFICANT RATE OF CHANGE ==========
##---- biomarkers that increase ----
plasma_biomarkers <- c("C2N_plasma_ptau217_out", "C2N_plasma_ptau217_ratio", "Fuji_plasma_ptau217_out", "AlzPath_plasma_ptau217_out", "Janssen_plasma_ptau217_out", 
                       "Roche_plasma_GFAP_out", "Roche_plasma_NfL_out", "Roche_plasma_ptau181", "QX_plasma_ptau181_out","QX_plasma_GFAP_out", "QX_plasma_NfL_out","PTAU_over_ABETA42","SUVR_compositeRef","MesialTemporal","TemporoParietal", "CDR_SOB")  
time_variables <- c("years_amy_onset", "years_tau_onset", "years_symp_onset", "years_tau_TP_onset")

# Initialize a list to store the results
derivative_results <- list()

for (biomarker in plasma_biomarkers) {
  for (time_var in time_variables) {
    
    # Define the formula for the GAM model
    formula <- as.formula(paste(biomarker, "~ s(", time_var, ", k=3) + s(RID, bs = 're')"))
    
    # Fit the model
    model <- gam(formula, data = qx_subset, method = "REML")
    
    # Generate predictions
    predictions <- predict(model, newdata = qx_subset, se.fit = TRUE)
    
    results <- data.frame(
      x = qx_subset[[time_var]],
      fit = predictions$fit,
      se.fit = predictions$se.fit
    )
    results <- results %>% arrange(x)
    
    #### Derivatives to get the range of significant increase and acceleration ###
    first_derivative <- Deriv(model)
    
    # Use column index to access eval points
    col_index <- match(time_var, colnames(first_derivative$eval))
    if (!is.na(col_index)) {
      eval_points <- first_derivative$eval[, col_index]
      
      # Calculate confidence intervals for the derivative
      ci_deriv <- confint.Deriv(first_derivative)
      lower_bound <- ci_deriv[[time_var]]$lower
      upper_bound <- ci_deriv[[time_var]]$upper
      
      # Find the time points where the derivative is significantly positive (increase)
      significant_pos <- which(lower_bound > 0 & upper_bound > 0)
      
      if (length( significant_pos) > 0) {
        significant_timepoints <- eval_points[ significant_pos]
      } else {
        significant_timepoints <- NA
      }
      
      # Store the results for this biomarker and time variable
      derivative_results[[paste(biomarker, time_var, sep = "_")]] <- data.frame(
        Biomarker = biomarker,
        Time_Variable = time_var,
        Significant_Timepoints = I(list(significant_timepoints))
      )
    } else {
      warning(paste("Time variable", time_var, "not found in first_derivative$eval"))
    }
  }
}



# Combine all results into a single data frame
final_derivative_results <- do.call(rbind, derivative_results)
write.csv(final_derivative_results, file = "~/Documents/FNIH_paper_Amyclockplasma/Annals_Neurology_submission/Review/derivatives_plasma_biomarkers_inc_sensitiv2.csv", row.names = FALSE)

##---- biomarkers that decrease ----

plasma_biomarkers <- c("C2N_plasma_Abeta42_Abeta40_out", "Roche_plasma_Ab42_Ab40_out", "QX_plasma_Ab42_Ab40_out", "Fuji_plasma_Ab42_Ab40_out", "atrophy")
time_variables <- c("years_amy_onset", "years_tau_onset", "years_symp_onset")

derivative_results <- list()

for (biomarker in plasma_biomarkers) {
  for (time_var in time_variables) {
    
    # Define the formula for the GAM model
    formula <- as.formula(paste(biomarker, "~ s(", time_var, ", k=3) + s(RID, bs = 're')"))
    
    # Fit the model
    model <- gam(formula, data = qx_subset, method = "REML")
    
    # Generate predictions
    predictions <- predict(model, newdata = qx_subset, se.fit = TRUE)
    
    results <- data.frame(
      x = qx_subset[[time_var]],
      fit = predictions$fit,
      se.fit = predictions$se.fit
    )
    results <- results %>% arrange(x)
    
    #### Derivatives to get the range of significant increase and acceleration ###
    first_derivative <- Deriv(model)
    
    # Use column index to access eval points
    col_index <- match(time_var, colnames(first_derivative$eval))
    if (!is.na(col_index)) {
      eval_points <- first_derivative$eval[, col_index]
      
      # Calculate confidence intervals for the derivative
      ci_deriv <- confint.Deriv(first_derivative)
      lower_bound <- ci_deriv[[time_var]]$lower
      upper_bound <- ci_deriv[[time_var]]$upper
      
      # Find the time points where the derivative is significantly positive (increase)
      significant_neg <- which(lower_bound < 0 & upper_bound < 0)
      
      if (length(significant_neg) > 0) {
        significant_timepoints <- eval_points[significant_neg]
      } else {
        significant_timepoints <- NA
      }
      
      # Store the results for this biomarker and time variable
      derivative_results[[paste(biomarker, time_var, sep = "_")]] <- data.frame(
        Biomarker = biomarker,
        Time_Variable = time_var,
        Significant_Timepoints = I(list(significant_timepoints))
      )
    } else {
      warning(paste("Time variable", time_var, "not found in first_derivative$eval"))
    }
  }
}


# Combine all results into a single data frame
final_derivative_results <- do.call(rbind, derivative_results)
write.csv(final_derivative_results, file = "derivatives_plasma_biomarkers_dec_sensit2.csv", row.names = FALSE)




##DESCRIPTIVE TABLES##
install.packages("tableone")
library(tableone)
#######Descriptives table AMYLOID AND TAU LONG DATASETS#########


mean(data_suvr_interval$Age)
sd(data_suvr_interval$Age)
table(data_suvr_interval$Race)
#tau dataset
tau_dataset_bl  <- tau_long_subset %>%
  distinct(RID, .keep_all = TRUE)

mean(data_tausuvr_interval$Age)
table(data_tausuvr_interval$PTRACCAT)
#########Descriptives table non plasma cohorts#########

CDR_dataset <- CDR_dataset %>%
  arrange(RID, VISDATE)
CDR_dataset_bl  <- CDR_dataset %>%
  distinct(RID, .keep_all = TRUE)


#define groups
CDR_dataset_bl <- CDR_dataset_bl %>%
  mutate(Ref_amy_group = case_when(
    !is.na(`Amyloid_age_mean`) ~ 1,   # Set to 1 when "mean amyloid age" is not NA
    Ref_group == 1 ~ 0,               # Set to 0 when ref.group is 1
    TRUE ~ NA_real_                  # Keep as NA for other cases (optional, remove if not needed)
  ))

CDR_dataset_bl$Ref_amy_group <-factor(CDR_dataset_bl$Ref_amy_group, levels = c(0,1), labels = c("reference", "amy_accum"))

CDR_dataset_bl <- CDR_dataset_bl %>%
  mutate(Ref_tau_group = case_when(
    !is.na(`Tau_age_mean`) ~ 1,   # Set to 1 when "mean amyloid age" is not NA
    Ref_group == 1 ~ 0,               # Set to 0 when ref.group is 1
    TRUE ~ NA_real_                  # Keep as NA for other cases (optional, remove if not needed)
  ))
CDR_dataset_bl$Ref_tau_group <-factor(CDR_dataset_bl$Ref_tau_group, levels = c(0,1), labels = c("reference", "tau_accum"))
dataset_combined <- filter(CDR_dataset_bl,Ref_amy_group =="reference" | Ref_amy_group == "amy_accum" | Ref_tau_group == "tau_accum" )

CDR_dataset_bl$CDR <- factor(CDR_dataset_bl$CDR, levels=c(0,0.5,1, 2,3))
CDR_dataset_bl$DIAGNOSIS<- factor(CDR_dataset_bl$DIAGNOSIS, levels=c(1, 2,3))
CDR_dataset_bl$PTGENDER<- factor(CDR_dataset_bl$PTGENDER, levels=c(1, 2))
CDR_dataset_bl <- CDR_dataset_bl %>%
  mutate(PTRACCAT = as.character(PTRACCAT),
         PTRACCAT = case_when(
           PTRACCAT == 4 ~ "Black",
           PTRACCAT == 5 ~ "White",
    TRUE ~ "Other"  # Groups all other levels as "Other"
  ))%>%
  mutate(PTRACCAT  = factor(PTRACCAT ))
table <- CreateTableOne(vars = c("Age", "PTGENDER", "APOE_binary", "PTEDUCAT", "PTRACCAT", "CDR", "DIAGNOSIS", "MMSCORE","CDRSB"), strata = "Ref_tau_group", addOverall = T, data = CDR_dataset_bl, test=T, testNonNormal = TRUE)
df_table <- as.data.frame(print(table, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE))
df_table <- df_table %>%
  rownames_to_column(var = "Variable")
##Compute p values for continious variables
continuous_vars <- c("Age", "PTEDUCAT","CDRSB", "MMSCORE" )  
cat_vars <- c("PTGENDER", "PTRACCAT", "APOE_binary","CDR","DIAGNOSIS")  

p_values <- data.frame(Variable = character(), P_value = numeric(), stringsAsFactors = FALSE)

# Loop through each variable
for (var in continuous_vars ) {
  # Filter out rows with missing values for the current variable
  complete_cases <- !is.na(CDR_dataset_bl[[var]]) & !is.na(CDR_dataset_bl$Ref_amy_group)
  
  # Perform the t-test
  test_result <- t.test(CDR_dataset_bl[[var]][complete_cases] ~ CDR_dataset_bl$Ref_amy_group[complete_cases])
  
  # Store the p-value
  p_values <- rbind(p_values, data.frame(Variable = var, P_value = test_result$p.value, stringsAsFactors = FALSE))
}

print(p_values)
##chi sq test for categorical
for (var in cat_vars) {
  # Convert variables to factors if they are not already
  CDR_dataset_bl[[var]] <- as.factor(CDR_dataset_bl[[var]])
  CDR_dataset_bl$Ref_tau_group <- as.factor(CDR_dataset_bl$Ref_tau_group)
  
  # Filter out rows with missing values for the current variable
  complete_cases <- !is.na(CDR_dataset_bl[[var]]) & !is.na(CDR_dataset_bl$Ref_tau_group)
  
  # Create a contingency table
  contingency_table <- table(CDR_dataset_bl[[var]][complete_cases], CDR_dataset_bl$Ref_tau_group[complete_cases])
  
  # Perform the Chi-Square test
  test_result <- chisq.test(contingency_table)
  
  # Store the p-value
  p_values <- rbind(p_values, data.frame(Variable = var, P_value = test_result$p.value, stringsAsFactors = FALSE))
}

#########Descriptives table plasma dataset##########   
final_dataset_plasma$EXAMDATE <- as.Date(final_dataset_plasma$EXAMDATE, format = "%m/%d/%y")
final_dataset_plasma$SCANDATE_AMY <- as.Date(final_dataset_plasma$SCANDATE_AMY, format = "%m/%d/%y")
final_dataset_plasma$EXAMDATE_CSF.x <- as.Date(final_dataset_plasma$EXAMDATE_CSF.x, format = "%m/%d/%y")
final_dataset_plasma$CDR <- factor(final_dataset_plasma$CDR.x, levels=c(0,0.5,1, 2,3))
final_dataset_plasma$DX<- factor(final_dataset_plasma$DX_123, levels=c(1, 2,3))
final_dataset_plasma$binary_scan1 <- factor(final_dataset_plasma$binary_scan1 , levels=c("Ab-Negative","Ab-Positive"), labels=c("Ab-Negative", "Ab-Positive"))
final_dataset_plasma$Race <- factor(final_dataset_plasma$Race)

final_dataset_plasma$ABETA42 <- as.numeric(final_dataset_plasma$ABETA42)
final_dataset_plasma$PTAU <- as.numeric(final_dataset_plasma$PTAU)
final_dataset_plasma$TAU <- as.numeric(final_dataset_plasma$TAU)

final_dataset_plasma <- final_dataset_plasma %>%
  arrange(RID,EXAMDATE)
final_dataset_plasma_bl <- distinct(final_dataset_plasma, RID, .keep_all = T)

##final table paper#

final_dataset_plasma_bl <- final_dataset_plasma_bl %>%
  mutate(Ref_amy_group = case_when(
    !is.na(`Amyloid_age_mean`) ~ 1,   # Set to 1 when "mean amyloid age" is not NA
    Ref_group == 1 ~ 0,               # Set to 0 when ref.group is 1
    TRUE ~ NA_real_                  # Keep as NA for other cases (optional, remove if not needed)
  ))

final_dataset_plasma_bl$Ref_amy_group <-factor(final_dataset_plasma_bl$Ref_amy_group, levels = c(0,1), labels = c("reference", "amy_accum"))

final_dataset_plasma_bl <- final_dataset_plasma_bl %>%
  mutate(Ref_tau_group = case_when(
    !is.na(`Tau_age_mean`) ~ 1,   # Set to 1 when "mean amyloid age" is not NA
    Ref_group == 1 ~ 0,               # Set to 0 when ref.group is 1
    TRUE ~ NA_real_                  # Keep as NA for other cases (optional, remove if not needed)
  ))

table <- CreateTableOne(vars = c("Age", "Sex", "Race", "APOE_binary", "Years_of_Education","CDR", "DX",  "SUVR_compositeRef","mean_Amyloid_age1", "mean_tau_ageMestemp", "EYO1_tau", "binary_scan1", "C2N_plasma_Abeta42_Abeta40", "C2N_plasma_ptau217", "C2N_plasma_ptau217_ratio", "Fuji_plasma_ptau217", "Fuji_plasma_Ab42_Ab40","AlzPath_plasma_ptau217","Janssen_plasma_ptau217", "Roche_plasma_Ab42_Ab40" , "Roche_plasma_GFAP" ,"Roche_plasma_NfL", "Roche_plasma_ptau181", "QX_plasma_Ab42_Ab40" , "QX_plasma_ptau181" ,"QX_plasma_GFAP","QX_plasma_NfL", "atrophy", "PTAU_over_ABETA42", "MMSCORE", "CDRSB", "MesialTemporal","TemporoParietal"), strata = "Ref_tau_group", addOverall = T, data = final_dataset_plasma_bl, test=T, testNonNormal = TRUE)
df_table <- as.data.frame(print(table, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE))
df_table <- df_table %>%
  rownames_to_column(var = "Variable")


##Compute p values for continious variables
continuous_vars <- c("Age", "Years_of_Education","SUVR_compositeRef","C2N_plasma_Abeta42_Abeta40", "C2N_plasma_ptau217", "C2N_plasma_ptau217_ratio", "Fuji_plasma_ptau217", "Fuji_plasma_Ab42_Ab40","AlzPath_plasma_ptau217","Janssen_plasma_ptau217", "Roche_plasma_Ab42_Ab40" , "Roche_plasma_GFAP" ,"Roche_plasma_NfL", "Roche_plasma_ptau181", "QX_plasma_Ab42_Ab40" , "QX_plasma_ptau181" ,"QX_plasma_GFAP","QX_plasma_NfL", "atrophy", "PTAU_over_ABETA42", "MMSCORE", "CDRSB", "MesialTemporal","TemporoParietal" )  
cat_vars <- c("Sex", "APOE_binary","CDR","DX", "binary_scan1" )  

p_values <- data.frame(Variable = character(), P_value = numeric(), stringsAsFactors = FALSE)

# Loop through each variable
for (var in continuous_vars) {
  # Filter out rows with missing values for the current variable
  complete_cases <- !is.na(final_dataset_plasma_bl[[var]]) & !is.na(final_dataset_plasma_bl$Ref_amy_group)
  
  # Perform the t-test
  test_result <- t.test(final_dataset_plasma_bl[[var]][complete_cases] ~ final_dataset_plasma_bl$Ref_amy_group[complete_cases])
  
  # Store the p-value
  p_values <- rbind(p_values, data.frame(Variable = var, P_value = test_result$p.value, stringsAsFactors = FALSE))
}

print(p_values)

##chi sq test for categorical
for (var in cat_vars) {
  # Convert variables to factors if they are not already
  final_dataset_plasma_bl[[var]] <- as.factor(final_dataset_plasma_bl[[var]])
  final_dataset_plasma_bl$Ref_amy_group <- as.factor(final_dataset_plasma_bl$Ref_amy_group)
  
  # Filter out rows with missing values for the current variable
  complete_cases <- !is.na(final_dataset_plasma_bl[[var]]) & !is.na(final_dataset_plasma_bl$Ref_amy_group)
  
  # Create a contingency table
  contingency_table <- table(final_dataset_plasma_bl[[var]][complete_cases], final_dataset_plasma_bl$Ref_amy_group[complete_cases])
  
  # Perform the Chi-Square test
  test_result <- chisq.test(contingency_table)
  
  # Store the p-value
  p_values <- rbind(p_values, data.frame(Variable = var, P_value = test_result$p.value, stringsAsFactors = FALSE))
}


###calculate slopes for biomarkers to report longitudinal change
final_dataset_plasma <- final_dataset_plasma %>%
  group_by(RID) %>%
  mutate(first_plasma = min(EXAMDATE, na.rm = TRUE)) %>%
  ungroup()
final_dataset_plasma  <- final_dataset_plasma  %>%
  mutate(time_since_firstplasma = as.numeric(difftime(EXAMDATE, first_plasma, units = "days")) / 365.25)

biomarkers <- c("C2N_plasma_Abeta42_Abeta40_out", "Roche_plasma_Ab42_Ab40_out", "Fuji_plasma_Ab42_Ab40_out","QX_plasma_Ab42_Ab40_out",  "C2N_plasma_ptau217_ratio", "C2N_plasma_ptau217_out","AlzPath_plasma_ptau217_out", "Janssen_plasma_ptau217_out","Fuji_plasma_ptau217_out",  
                "Roche_plasma_ptau181", "QX_plasma_ptau181_out", "Roche_plasma_GFAP_out", "QX_plasma_GFAP_out","Roche_plasma_NfL_out",  "QX_plasma_NfL_out" ,"PTAU_over_ABETA42", "SUVR_compositeRef",  "atrophy", "MesialTemporal", "TemporoParietal", "CDRSB")


results <- data.frame(Biomarker = character(),
                      Slope = numeric(),
                      SE = numeric(),
                      p_value = numeric(),
                      stringsAsFactors = FALSE)

# Loop through each biomarker and fit the LME model
for (bio in biomarkers) {
  
  # Fit LME model
  model <- lme(as.formula(paste(bio, "~ time_since_firstplasma")), 
               random = ~ 1 | RID, 
               data = subset(final_dataset_plasma, Ref_amy_group=="reference"), 
               method = "REML", na.action = na.omit)
  
  # Extract fixed effect estimate (slope) for 'time'
  summary_model <- summary(model)
  slope <- summary_model$tTable["time_since_firstplasma", "Value"]
  se <- summary_model$tTable["time_since_firstplasma", "Std.Error"]
  p_value <- summary_model$tTable["time_since_firstplasma", "p-value"]
  print(summary_model)
  # Append results to data frame
  results <- rbind(results, data.frame(Biomarker = bio,  Slope = slope, SE = se, p_value = p_value))
}

# Save results to CSV
write.csv(results, "biomarker_change_combined.csv", row.names = FALSE)



final_dataset_plasma <- final_dataset_plasma %>%
  mutate(Ref_amy_group = case_when(
    !is.na(`Amyloid_age_mean`) ~ 1,   # Set to 1 when "mean amyloid age" is not NA
    Ref_group == 1 ~ 0,               # Set to 0 when ref.group is 1
    TRUE ~ NA_real_                  # Keep as NA for other cases (optional, remove if not needed)
  ))

final_dataset_plasma$Ref_amy_group <-factor(final_dataset_plasma$Ref_amy_group, levels = c(0,1), labels = c("reference", "amy_accum"))

final_dataset_plasma <- final_dataset_plasma %>%
  mutate(Ref_tau_group = case_when(
    !is.na(`Tau_age_mean`) ~ 1,   # Set to 1 when "mean amyloid age" is not NA
    Ref_group == 1 ~ 0,               # Set to 0 when ref.group is 1
    TRUE ~ NA_real_                  # Keep as NA for other cases (optional, remove if not needed)
  ))
final_dataset_plasma$Ref_tau_group <-factor(final_dataset_plasma$Ref_tau_group, levels = c(0,1), labels = c("reference", "tau_accum"))

dataset_combined <- filter(final_dataset_plasma,Ref_amy_group =="reference" | Ref_amy_group == "amy_accum" | Ref_tau_group == "tau_accum" )

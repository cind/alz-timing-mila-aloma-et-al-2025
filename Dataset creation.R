# R/dataset_creation.R

#====================#
# DATASETS CREATION
#====================#


# Load required packages
library(dplyr)
library(data.table)
library(lubridate)
library(readr)
library(tibble)

####### AMYLOID PET dataset #####
final_dataset_fbp <- subset(final_dataset_all, Tracer=="FBP") 
final_dataset_fbp <- merge(x=final_dataset_fbp, ages_dataset[, c("RID", "Amyloid_age_mean", "Tau_age_mean", "Tau_TP_age_mean","est_conversion_age_amy", "est_conversion_age_tau")], by=c("RID"), all.x = T)

final_dataset_fbp$years_amy_onset <- final_dataset_fbp$Age - final_dataset_fbp$Amyloid_age_mean
final_dataset_fbp$years_tau_onset <- final_dataset_fbp$Age - final_dataset_fbp$Tau_age_mean
final_dataset_fbp$years_tau_TP_onset <- final_dataset_fbp$Age - final_dataset_fbp$Tau_TP_age_mean
final_dataset_fbp$years_symp_onset <- final_dataset_fbp$Age - final_dataset_fbp$est_conversion_age_tau
##add DX##
DX <- read_csv("CDR_DX_demog0806.csv")
setDT(final_dataset_fbp)
setDT(DX)
final_dataset_fbp$EXAMDATE <- as.Date(final_dataset_fbp$EXAMDATE, format = "%m/%d/%y")
closest_match <- DX[
  final_dataset_fbp, 
  on = .(RID, DX_EXAMDATE = EXAMDATE), 
  roll = "nearest", # This performs the nearest date join
  mult = "first"    # Return the first match if multiple nearest dates are the same
]
closest_match <- closest_match[, .(RID, EXAMDATE = DX_EXAMDATE, DIAGNOSIS)]
final_dataset_fbp <- final_dataset_fbp[
  closest_match, 
  on = .(RID, EXAMDATE), 
  nomatch = NA
]

###ref group (stable a-neg and cu)
final_dataset_fbp$CDR <-  as.factor(final_dataset_fbp$CDR)
final_dataset_fbp$DIAGNOSIS <-  as.factor(final_dataset_fbp$DIAGNOSIS)

final_dataset_fbp <- final_dataset_fbp %>%
  mutate(CDR_DX_Group = case_when(
    CDR=="0" & DIAGNOSIS=="1" ~ "CN_CDR0"    ,
    CDR!="0" & DIAGNOSIS!="1" ~ "CI_CDRgt0",
    TRUE ~ NA_character_  
  ))

ref_group <- final_dataset_fbp %>%
  group_by(RID) %>%
  summarize(Ref_group = all(binary_scan1 == "Negative" & CDR_DX_Group == "CN_CDR0" & is.na(mean_tau_ageMestemp) )) %>%
  mutate(Ref_group = as.integer(Ref_group)) # Convert logical to integer (1/0)

final_dataset_fbp <- final_dataset_fbp %>%
  left_join(ref_group %>% select(RID, Ref_group), by = "RID")

ref_group_uniq <-  final_dataset_fbp %>% distinct(RID, Ref_group)

##add tau PET##
setDT(final_dataset_fbp)
setDT(tau_dataset)
tau_dataset$SCANDATE <- as.Date(tau_dataset$SCANDATE, format = "%m/%d/%y")
closest_match <- tau_dataset[
  final_dataset_fbp, 
  on = .(RID),           # Join based on RID first
  allow.cartesian = TRUE # Allow multiple matches for each RID, if necessary
]

closest_match[, diff_days := abs(as.numeric(difftime(SCANDATE, EXAMDATE, units = "days")))]
closest_match_filtered <- closest_match[diff_days <= 365]
closest_match_smallest <- closest_match_filtered[, .SD[which.min(diff_days)], by = .(RID, EXAMDATE)]

closest_match_smallest <- closest_match_smallest[, .(RID, EXAMDATE, SCANDATE_tau=SCANDATE, MesialTemporal, TemporoParietal, tau_MesTemp_binary)]
final_dataset_fbp <- merge(final_dataset_fbp, closest_match_smallest, 
                           by = c("RID", "EXAMDATE"), 
                           all.x = TRUE)
#add mri
setDT(final_dataset_fbp)
setDT(MRI_dataset)
MRI_dataset$EXAMDATE_MRI <- as.Date(MRI_dataset$EXAMDATE, format = "%m/%d/%y")
MRI_dataset$EXAMDATE <- NULL
closest_match <- MRI_dataset[
  final_dataset_fbp, 
  on = .(RID),           # Join based on RID first
  allow.cartesian = TRUE # Allow multiple matches for each RID, if necessary
]

closest_match[, diff_days := abs(as.numeric(difftime(EXAMDATE_MRI, EXAMDATE, units = "days")))]
closest_match_filtered <- closest_match[diff_days <= 365]
closest_match_smallest <- closest_match_filtered[, .SD[which.min(diff_days)], by = .(RID, EXAMDATE)]

closest_match_smallest <- closest_match_smallest[, .(RID, EXAMDATE, EXAMDATE_MRI, metaROI.AgeAdj)]
final_dataset_fbp <- merge(final_dataset_fbp, closest_match_smallest, 
                           by = c("RID", "EXAMDATE"), 
                           all.x = TRUE)
#add mmse
setDT(final_dataset_fbp)
setDT(MMSE_dataset)

closest_match <- MMSE_dataset[
  final_dataset_fbp, 
  on = .(RID),           # Join based on RID first
  allow.cartesian = TRUE # Allow multiple matches for each RID, if necessary
]

closest_match[, diff_days := abs(as.numeric(difftime(VISDATE, EXAMDATE, units = "days")))]
closest_match_filtered <- closest_match[diff_days <= 365]
closest_match_smallest <- closest_match_filtered[, .SD[which.min(diff_days)], by = .(RID, EXAMDATE)]

closest_match_smallest <- closest_match_smallest[, .(RID, EXAMDATE, VISDATE_MMSE=VISDATE, MMSCORE)]
final_dataset_fbp <- merge(final_dataset_fbp, closest_match_smallest, 
                           by = c("RID", "EXAMDATE"), 
                           all.x = TRUE)

final_dataset_fbp <- subset(final_dataset_fbp, Tracer=="FBP")
write.csv(final_dataset_fbp, "amyloid_dataset.csv")
#get centiloid values for FBP SUVR (composite ref region)
final_dataset_fbp$Centiloids_fbp_comp <- 300.66*final_dataset_fbp$SUVR_compositeRef - 208.84
#equation [18F]FBP: CLcomposite =300.66 × SUVR_compositeRef − 208.84


####### TAU PET dataset #####
tau_dataset <- read_csv("tau_dataset2107.csv")
tau_dataset <- merge(x=tau_dataset, estimated_ages[, c("RID", "Amyloid_age_mean", "Tau_age_mean","Tau_TP_age_mean", "est_conversion_age_amy", "est_conversion_age_tau")], by=c("RID"), all.x = T)
tau_dataset <- merge(x=tau_dataset, estimated_ages[, c("RID", "Tau_TP_age_mean")], by=c("RID"), all.x = T)


tau_dataset$years_amy_onset <- tau_dataset$Age - tau_dataset$Amyloid_age_mean
tau_dataset$years_tau_onset <- tau_dataset$Age - tau_dataset$Tau_age_mean
tau_dataset$years_tau_TP_onset <- tau_dataset$Age - tau_dataset$Tau_TP_age_mean
tau_dataset$years_symp_onset <- tau_dataset$Age - tau_dataset$est_conversion_age_tau

tau_dataset <- merge(x=tau_dataset, ref_group_uniq, by=c("RID"), all.x = T)

setDT(final_dataset_all)
setDT(tau_dataset)
tau_dataset$SCANDATE <- as.Date(tau_dataset$SCANDATE, format = "%m/%d/%y")
matched <- final_dataset_all[
  tau_dataset, 
  on = .(RID),           # Join based on RID first
  allow.cartesian = TRUE # Allow multiple matches for each RID, if necessary
]

matched[, diff_days := abs(as.numeric(difftime(SCANDATE, EXAMDATE, units = "days")))]
closest_match_filtered <- matched[diff_days <= 365]
closest_match_smallest <- closest_match_filtered[, .SD[which.min(diff_days)], by = .(RID, SCANDATE)]

closest_match_smallest <- closest_match_smallest[, .(RID, EXAMDATE_amy=EXAMDATE,SCANDATE, SUVR_compositeRef, binary_scan1)]
tau_dataset <- merge(tau_dataset, closest_match_smallest, 
                     by = c("RID", "SCANDATE"), 
                     all.x = TRUE)
#add csf
setDT(tau_dataset)
setDT(csf_dataset)

closest_match <- csf_dataset[
  tau_dataset, 
  on = .(RID),           # Join based on RID first
  allow.cartesian = TRUE # Allow multiple matches for each RID, if necessary
]

closest_match[, diff_days := abs(as.numeric(difftime(EXAMDATE, SCANDATE, units = "days")))]
closest_match_filtered <- closest_match[diff_days <= 365]
closest_match_smallest <- closest_match_filtered[, .SD[which.min(diff_days)], by = .(RID, SCANDATE)]

closest_match_smallest <- closest_match_smallest[, .(RID, SCANDATE, EXAMDATE_csf =EXAMDATE, PTAU_over_ABETA42)]
tau_dataset <- merge(tau_dataset, closest_match_smallest, 
                     by = c("RID", "SCANDATE"), 
                     all.x = TRUE)

#add mri
setDT(tau_dataset)
setDT(MRI_dataset)

closest_match <- MRI_dataset[
  tau_dataset, 
  on = .(RID),           # Join based on RID first
  allow.cartesian = TRUE # Allow multiple matches for each RID, if necessary
]

closest_match[, diff_days := abs(as.numeric(difftime(EXAMDATE_MRI, SCANDATE, units = "days")))]
closest_match_filtered <- closest_match[diff_days <= 365]
closest_match_smallest <- closest_match_filtered[, .SD[which.min(diff_days)], by = .(RID, SCANDATE)]

closest_match_smallest <- closest_match_smallest[, .(RID, SCANDATE, EXAMDATE_MRI, metaROI.AgeAdj)]
tau_dataset <- merge(tau_dataset, closest_match_smallest, 
                     by = c("RID", "SCANDATE"), 
                     all.x = TRUE)
#add mmse
setDT(tau_dataset)
setDT(MMSE_dataset)

closest_match <- MMSE_dataset[
  tau_dataset, 
  on = .(RID),           # Join based on RID first
  allow.cartesian = TRUE # Allow multiple matches for each RID, if necessary
]

closest_match[, diff_days := abs(as.numeric(difftime(VISDATE, SCANDATE, units = "days")))]
closest_match_filtered <- closest_match[diff_days <= 365]
closest_match_smallest <- closest_match_filtered[, .SD[which.min(diff_days)], by = .(RID, SCANDATE)]

closest_match_smallest <- closest_match_smallest[, .(RID, SCANDATE, VISDATE_MMSE=VISDATE, MMSCORE)]
tau_dataset <- merge(tau_dataset, closest_match_smallest, 
                     by = c("RID", "SCANDATE"), 
                     all.x = TRUE)

#### remove outliers- mesialtemporal and temporoparietal smaller than 4 !
tau_dataset <- subset(tau_dataset, MesialTemporal< 4 & TemporoParietal<4)
write.csv(tau_dataset, "tau_dataset.csv")

####### CSF PTAU/AB dataset #####
csf_dataset <- read_csv("UPENNBIOMK_ROCHE_ELECSYS_02Apr2024.csv")

PTDEMO <- fread("PTDEMOG.csv")
PTDEMO_unique <- PTDEMO[!duplicated(RID), .(RID, PTGENDER, PTDOBMM, PTDOBYY, PTEDUCAT)]
csf_dataset <- merge(csf_dataset, PTDEMO_unique, by = "RID", all.x = TRUE)
csf_dataset$EXAMDATE <- as.Date(csf_dataset$EXAMDATE, format = "%m/%d/%y")
csf_dataset$AGE_csf <-round(as.numeric(lubridate::year(csf_dataset$EXAMDATE)) -
                              csf_dataset$PTDOBYY +
                              ((as.numeric(lubridate::month(csf_dataset$EXAMDATE)) -
                                  csf_dataset$PTDOBMM) / 12), digits=1)

adnimerge <- fread("ADNIMERGE.csv")
adnimerge_unique <- adnimerge[!duplicated(RID), .(RID, APOE4)]
csf_dataset <- merge(csf_dataset, adnimerge_unique, by = "RID", all.x = TRUE)
csf_dataset <- csf_dataset %>%
  mutate(APOE_binary  = factor(ifelse(APOE4 == 0, "no carrier", "carrier")))
CDR <- read_csv("CDR_22May2024.csv")
CDR$VISDATE <- as.Date(CDR$VISDATE, format = "%m/%d/%y")

setDT(csf_dataset)
setDT(CDR)

closest_match <- CDR[
  csf_dataset, 
  on = .(RID, VISDATE = EXAMDATE), 
  roll = "nearest", # This performs the nearest date join
  mult = "first"    # Return the first match if multiple nearest dates are the same
]
closest_match <- closest_match[, .(RID, EXAMDATE = VISDATE, CDGLOBAL)]
csf_dataset <- csf_dataset[
  closest_match, 
  on = .(RID, EXAMDATE), 
  nomatch = NA
]


csf_dataset <- merge(x=csf_dataset, estimated_ages[, c("RID", "Amyloid_age_mean", "Tau_age_mean", "Tau_TP_age_mean","est_conversion_age_amy", "est_conversion_age_tau")], by=c("RID"), all.x = T)

csf_dataset$years_amy_onset <- csf_dataset$Age - csf_dataset$Amyloid_age_mean
csf_dataset$years_tau_onset <- csf_dataset$Age -csf_dataset$Tau_age_mean
csf_dataset$years_tau_TP_onset <- csf_dataset$Age -csf_dataset$Tau_TP_age_mean
csf_dataset$years_symp_onset <- csf_dataset$Age - csf_dataset$est_conversion_age_tau

csf_dataset$PTAU_over_ABETA42 <- csf_dataset$PTAU/csf_dataset$ABETA42
csf_dataset <- merge(x=csf_dataset, ref_group_uniq, by=c("RID"), all.x = T)
write.csv(csf_dataset, "csf_dataset.csv")
#######MRI Cortical thicknes#####

MRI_dataset <- read_csv("ADNI_metaROI_thickness_harmonized.csv")

MRI_dataset <- merge(MRI_dataset, adnimerge_unique[, c("RID", "APOE4")], by = "RID", all.x = TRUE)
MRI_dataset <- MRI_dataset %>%
  mutate(APOE_binary  = factor(ifelse(APOE4 == 0, "no carrier", "carrier")))

MRI_dataset <- merge(x=MRI_dataset, estimated_ages[, c("RID", "Amyloid_age_mean", "Tau_age_mean", "Tau_age_mean","est_conversion_age_amy", "est_conversion_age_tau")], by=c("RID"), all.x = T)

MRI_dataset$years_amy_onset <- MRI_dataset$Age - MRI_dataset$Amyloid_age_mean
MRI_dataset$years_tau_onset <- MRI_dataset$Age -MRI_dataset$Tau_age_mean
MRI_dataset$years_tau_TP_onset <- MRI_dataset$Age -MRI_dataset$Tau_TP_age_mean
MRI_dataset$years_symp_onset <- MRI_dataset$Age - MRI_dataset$est_conversion_age_tau


MRI_dataset <- merge(x=MRI_dataset, ref_group_uniq, by=c("RID"), all.x = T)
write.csv(MRI_dataset, "MRI_dataset.csv")

####### CDR dataset #####
CDR_dataset <- read_csv("CDR_22May2024.csv")

CDR_dataset <- merge(CDR_dataset, PTDEMO_unique, by = "RID", all.x = TRUE)
CDR_dataset$VISDATE <- as.Date(CDR_dataset$VISDATE, format = "%m/%d/%y")
CDR_dataset$Age_CDR <-round(as.numeric(lubridate::year(CDR_dataset$VISDATE)) -
                              CDR_dataset$PTDOBYY +
                              ((as.numeric(lubridate::month(CDR_dataset$VISDATE)) -
                                  CDR_dataset$PTDOBMM) / 12), digits=1)

CDR_dataset<- merge(CDR_dataset, adnimerge_unique[, c("RID", "APOE4")], by = "RID", all.x = TRUE)
CDR_dataset <- CDR_dataset %>%
  mutate(APOE_binary  = factor(ifelse(APOE4 == 0, "no carrier", "carrier")))

#add dx
DX_dataset <- read_csv("DXSUM_PDXCONV_21May2024.csv")
setDT(CDR_dataset)
setDT(DX_dataset)
DX_dataset$EXAMDATE_DX <- as.Date(DX_dataset$EXAMDATE , format = "%m/%d/%y")
closest_match <- DX_dataset[
  CDR_dataset, 
  on = .(RID),           # Join based on RID first
  allow.cartesian = TRUE # Allow multiple matches for each RID, if necessary
]

closest_match[, diff_days := abs(as.numeric(difftime(EXAMDATE_DX, VISDATE, units = "days")))]
closest_match_filtered <- closest_match[diff_days <= 365]
closest_match_smallest <- closest_match_filtered[, .SD[which.min(diff_days)], by = .(RID, VISDATE)]

closest_match_smallest <- closest_match_smallest[, .(RID, EXAMDATE_DX, VISDATE, DIAGNOSIS)]
CDR_dataset <- merge(CDR_dataset, closest_match_smallest, 
                     by = c("RID", "VISDATE"), 
                     all.x = TRUE)
#add mmse
MMSE_dataset <- read_csv("MMSE_22Aug2024.csv")
setDT(CDR_dataset)
setDT(MMSE_dataset)
MMSE_dataset$VISDATE_mmse <- as.Date(MMSE_dataset$VISDATE, format = "%m/%d/%y")
closest_match <- MMSE_dataset[
  CDR_dataset, 
  on = .(RID),           # Join based on RID first
  allow.cartesian = TRUE # Allow multiple matches for each RID, if necessary
]

closest_match[, diff_days := abs(as.numeric(difftime(VISDATE_mmse, VISDATE, units = "days")))]
closest_match_filtered <- closest_match[diff_days <= 365]
closest_match_smallest <- closest_match_filtered[, .SD[which.min(diff_days)], by = .(RID, VISDATE)]

closest_match_smallest <- closest_match_smallest[, .(RID, VISDATE_mmse, VISDATE, MMSCORE)]
CDR_dataset <- merge(CDR_dataset, closest_match_smallest, 
                     by = c("RID", "VISDATE"), 
                     all.x = TRUE)

CDR_dataset <- merge(x=CDR_dataset, estimated_ages[, c("RID", "Amyloid_age_mean", "Tau_age_mean", "Tau_TP_age_mean", "est_conversion_age_amy", "est_conversion_age_tau")], by=c("RID"), all.x = T)

CDR_dataset$years_amy_onset <- CDR_dataset$Age - CDR_dataset$Amyloid_age_mean
CDR_dataset$years_tau_onset <- CDR_dataset$Age -CDR_dataset$Tau_age_mean
CDR_dataset$years_tau_TP_onset <- CDR_dataset$Age -CDR_dataset$Tau_TP_age_mean
CDR_dataset$years_symp_onset <- CDR_dataset$Age - CDR_dataset$est_conversion_age_tau

CDR_dataset <- merge(x=CDR_dataset, ref_group_uniq, by=c("RID"), all.x = T)
write.csv(CDR_dataset, "CDR_dataset.csv")

#######PLASMA DATASET#####

longitudinal_dataset_2024_04_15 <- read.csv("dataset.csv")
longitudinal_dataset_2024_04_15$SCANDATE_AMY <- as.Date(longitudinal_dataset_2024_04_15$SCANDATE_AMY, format = "%m/%d/%y")
longitudinal_dataset_2024_04_15$SCANDATE_TAU <- as.Date(longitudinal_dataset_2024_04_15$SCANDATE_TAU, format = "%m/%d/%y")
longitudinal_dataset_2024_04_15$EXAMDATE <- as.Date(longitudinal_dataset_2024_04_15$EXAMDATE, format = "%m/%d/%y")

#add Amyloid PET and clock variables #

final_dataset_all <- read_csv("final_dataset_all_0207.csv")
final_dataset_all <- final_dataset_all %>% rename(SCANDATE_AMY = EXAMDATE)
final_dataset_all$SCANDATE_AMY <- as.Date(final_dataset_all$SCANDATE_AMY, format = "%m/%d/%y")

final_dataset_all <- final_dataset_all %>%
  rename(RID = ID)

final_dataset_plasma <- merge(x=longitudinal_dataset_2024_04_15, final_dataset_all, by=c("RID", "SCANDATE_AMY"))

#add tau PET PVC #
tau_dataset <- tau_dataset %>%
  rename(c(MesialTemporal = MesialTemporal_pvc, MetaTemporal = Metatemporal_pvc))
tau_alll <- tau_alll %>%
  rename(RID = ID)

final_dataset_plasma <- final_dataset_plasma %>%
  rename(MesialTemporal_notpvc = MesialTemporal)

tau_alll <- tau_alll %>%
  rename(SCANDATE_TAU = SCANDATE)

tau_selected <- tau_alll %>% select(RID, SCANDATE_TAU, MesialTemporal_notpvc, MesialTemporal,TemporoParietal)
plasma_selected <- final_dataset_plasma %>% select(RID, SCANDATE_TAU, MesialTemporal_notpvc)

final_dataset_plasma <- merge(x=final_dataset_plasma, tau_alll[, c("RID", "SCANDATE_TAU", "MesialTemporal_notpvc", "MesialTemporal", "TemporoParietal")], by=c("RID", "SCANDATE_TAU", "MesialTemporal_notpvc"), all.x = T)
final_dataset_plasma <- final_dataset_plasma[!duplicated(final_dataset_plasma), ]
final_dataset_plasma <- final_dataset_plasma %>% rename(c(TemporoParietal_notpvc = TemporoParietal.x, TemporoParietal = TemporoParietal.y))

###CALCULATION OF years since amyloid onset, tau onset and symptom onset BASED ON AGE AT PLASMA MEASURE (for plasma bioamrkers)###
final_dataset_plasma <- merge(x=final_dataset_plasma, ages_dataset[, c("RID", "Amyloid_age_mean", "Tau_age_mean", "Tau_TP_age_mean" ,"est_conversion_age_amy", "est_conversion_age_tau")], by=c("RID"), all.x = T)

final_dataset_plasma$years_amy_onset <- final_dataset_plasma$Age - final_dataset_plasma$Amyloid_age_mean
final_dataset_plasma$years_tau_onset <- final_dataset_plasma$Age -final_dataset_plasma$Tau_age_mean
final_dataset_plasma$years_tau_TP_onset <- final_dataset_plasma$Age -final_dataset_plasma$Tau_TP_age_mean
final_dataset_plasma$years_symp_onset <- final_dataset_plasma$Age - final_dataset_plasma$est_conversion_age_tau

final_dataset_plasma <- merge(x=final_dataset_plasma, ref_group_uniq, by=c("RID"), all.x = T)
write.csv(final_dataset_plasma, "final_dataset_plasma.csv")

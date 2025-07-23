# R/tau_clock_main.R

#====================#
# Tau Clock 
#====================#
# Load required libraries
required_packages <- c("tidyverse", "nlme", "mgcv", "lme4", "data.table", "here")
lapply(required_packages, library, character.only = TRUE)

# Constants
SUVR_MIN <- 0.98
SUVR_MAX <- 2.04
SUVR_POSITIVE_THRESHOLD <- 1.41005
# Load dataset
data_path <- here("data", "tau_dataset.csv")
dataset <- read_csv(data_path)
# Ensure date formatting
dataset$SCANDATE <- as.Date(dataset$SCANDATE, format = "%m/%d/%y")
# Calculate conversion age in PET converters (used later for control checks and validation)
dataset <- dataset %>%
  arrange(RID,SCANDATE) %>%
  group_by(RID) %>% 
  mutate(Previous_tau_status = lag(tau_MesTemp_binary))

tau_conversion <- dataset  %>%
  filter(Previous_tau_status  == "Negative", tau_MesTemp_binary == "Positive") %>%
  group_by(RID) %>% distinct(RID)

tau_converters_subset <- dataset %>%
  semi_join(tau_conversion, by = "RID")

midpoint_age_tau <- tau_converters_subset %>%
  group_by(RID) %>%
  summarize(
    Last_0_Age_tau = last(AGE_tau[tau_MesTemp_binary == "Negative"]),  # Age at last 0
    First_1_Age_tau = first(AGE_tau[tau_MesTemp_binary == "Positive"]),  # Age at first 1
    Midpoint_Age_tau_MesTemp = (Last_0_Age_tau  +  First_1_Age_tau) / 2  # Midpoint age as average of ages at last 0 and first 1
  )
dataset <- merge(dataset, midpoint_age_tau, by = "RID", all.x = TRUE)
dataset$Tau_time_convMesTemp <- (dataset$AGE_tau - dataset$ Midpoint_Age_tau_MesTemp)

# Filter participants with follow-up
Multiple_observations_data  <- dataset %>% 
  group_by(RID) %>% 
  filter(n() > 1) %>% 
  ungroup() %>% 
  arrange(RID, SCANDATE) %>%
  group_by(RID) %>%
  mutate(first_PET_scan = min(SCANDATE),
         last_PET_scan = max(SCANDATE),
         total_fu_time = as.numeric(difftime(last_PET_scan, first_PET_scan, units = "days")) / 365.25,
         scan_int = as.numeric(difftime(SCANDATE, first_PET_scan, units = "days")) / 365.25)


#Get slopes of tau accumulation
multiple_observation_data <- multiple_observation_data  %>%
  filter(!is.na(MesialTemporal))
fit_MesialTemprate <- lme(MesialTemporal ~ scan_int,  random= ~1 +scan_int| RID, data=multiple_observation_data , na.action = na.omit, method = "ML")
su
tau_slopes_df <- coef(fit_MesialTemprate) %>%
  as.data.frame() %>%
  rownames_to_column("RID") %>%
  mutate(RID = as.integer(RID)) %>%
  select(RID, slope = scan_int)
# Add slopes to dataset
Multiple_observations_data <- left_join(multiple_observations_data, tau_ slopes_df, by = "RID")

# calculate SUVR midpoint
baseline_MesialTemporal <- multiple_observations_data %>% arrange(RID, SCANDATE) %>% group_by(RID) %>% summarise(Baseline_SUVR = first(MesialTemporal))
multiple_observations_data  <- left_join(multiple_observations_data, baseline_MesialTemporal, by = "RID") %>%
 mutate(MesialTemporal_midpoint = Baseline_SUVR + (total_fu_time / 2) * slope)

dataset_fu_uq <- multiple_observations_data %>% distinct(RID, .keep_all = TRUE)

# GAM of slope vs SUVR midpoint
fit_slope_gam <- gam(slope ~ s(MesialTemporal_midpoint, bs = "cr"), data = dataset_fu_uq)
summary(fit_slope_gam)
plot(fit_slope_gam)
#plot
ggplot(dataset_fu_uq , aes(x=MesialTemporal_midpoint, y=slope)) + 
  geom_point(aes(color= APOE_binary))+ scale_color_manual(values=c("darkblue", "forestgreen")) +
  geom_smooth(method="gam") + theme_classic()
# Variance check
dataset_fu_uq <- dataset_fu_uq[order(dataset_fu_uq$MesialTemporal_midpoint), ]
predictions <- predict(fit_slope_gam, newdata = dataset_fu_uq, se.fit = TRUE)
variance <- predictions$se.fit^2
mean_variance <- mean(variance, na.rm = TRUE)
sd_variance <- sd(variance, na.rm = TRUE)
cutpoint <- quantile(variance, 0.90) 
high_variance_points <- dataset_fu_uq$MesialTemporal_midpoint[variance > cutpoint]
plot(dataset_fu_uq$MesialTemporal_midpoint, variance, type = "l", lwd = 2,
     ylab = "Variance", xlab = "SUVR Midpoint",
     main = "Variance Across SUVR Midpoint")
abline(h = cutpoint, col = "red", lty = 2)
points(high_variance_points, variance[variance > cutpoint], col = "blue", pch = 19)

# Run GAM again with restricted interval 
data_suvr_interval <- subset(dataset_fu_uq,  between(MesialTemporal_midpoint, SUVR_MIN, SUVR_MAX))
fit_SUVRrate_1 <- gam(slope ~ s(MesialTemporal_midpoint, bs="cr"), data=
                      data_suvr_interval)
summary(fit_SUVRrate_1)
plot(fit_SUVRrate_1)

# Calculation of Tau time 

predicted <- predict(fit_SUVRrate_1, newdata = data_suvr_interval)
cor(data_suvr_interval$slope, predicted)
plot(fit_SUVRrate_1)
data_suvr_interval$estim_rate_tau <- predict(fit_SUVRrate_1, type = "response")  

#calculate time intervals between SUVR
##create a data frame with intervals 
MesialTemporal_midpoint <- seq(from= SUVR_MIN,to=SUVR_MAX,by=.0001)
MesialTemporal_midpoint <- MesialTemporal_midpoint + 0.00005
MesialTemporal_midpoint <- as.data.frame(MesialTemporal_midpoint)

#extract predicted slopes for all SUVR within interval
MesialTemporal_midpoint$estim_rate_i <- predict(fit_SUVRrate,newdata =MesialTemporal_midpoint, type="response")
MesialTemporal_midpoint$recip_rate_i <- 1/SUVR_midpoint$estim_rate_i
MesialTemporal_midpoint$Tau_TIME_INT_i  <- 0.0001*MesialTemporal_midpoint$recip_rate_i
MesialTemporal_midpoint$TimeSum_i <- cumsum(MesialTemporal_midpoint$Tau_TIME_INT_i) ###cumulative sum of time###
MesialTemporal_midpoint$Tau_time <- MesialTemporal_midpoint$TimeSum_i - MesialTemporal_midpoint$Time[MesialTemporal_midpoint$MesialTemporal_midpoint == SUVR_POSITIVE_THRESHOLD] # Calculate the time difference from the target SUVR for each observation

#Merge with original dataset and estimate tau time for the rest of scans ###
data_SUVR_int <- subset(dataset,  between(MesialTemporal, SUVR_MIN, SUVR_MAX))
data_SUVR_int  <- data_SUVR_int  %>%
rename(SUVR = MesialTemporal)
MesialTemporal_midpoint <- MesialTemporal_midpoint %>%
rename(SUVR = MesialTemporal_midpoint)

seq_dataset <- as.data.frame(MesialTemporal_midpoint)
library(data.table)
setDT(data_SUVR_int)
setDT(seq_dataset )
# Define a function to find the nearest value 
find_nearest <- function(x, y) {
  idx <- findInterval(x, y, all.inside = TRUE)
  y[ifelse(idx == 0, 1, ifelse(idx == length(y), length(y), idx))]
}

data_SUVR_int[, Nearest_Value := find_nearest(SUVR, seq_dataset$SUVR)] data_SUVR_int <- merge(data_SUVR_int, seq_dataset, by.x = "Nearest_Value", by.y = "SUVR", all.x = TRUE)


##estimates tau time for all pet scans
model_tau <- gam(Tau_time ~ s(SUVR, bs="cr"), data = data_SUVR_int)
summary(model_tau)
data_SUVR_int$Tau_time <- predict(model_tau ,newdata = data_SUVR_int, type="response")
# Estimate tau age
subset_taupos <- data_SUVR_int %>% 
  filter(SUVR > SUVR_POSITIVE_THRESHOLD) %>%
  mutate(Tau_age = Age - Tau_time) %>%
  group_by(RID) %>%
  mutate(Tau_age_mean = mean(Tau_age, na.rm = TRUE))

tau_ages <- subset_taupos %>%
  distinct(RID, Tau_age_mean)
tau_dataset <- left_join(tau_dataset, tau_ages, by = "RID")
tau_dataset <- left_join(tau_dataset, data_SUVR_int %>% select(RID, SCANDATE, Tau_time), 
                     by = c("RID", "SCANDATE"))
########################################################################################
# Control checks
#estimated vs actual time in converters

cor.test(data_SUVR_int$Tau_time_convMesTemp, data_SUVR_int$Tau_time, method = "spearman")
model_check <- lm(Tau_time_convMesTemp ~ Tau_time, data = data_SUVR_int)
summary(model_check)
#estimated tau age vs midpoint age for converters
subset_taupos_uq <- subset_taupos %>%
  distinct(RID, .keep_all = TRUE)
ggplot(subset_taupos, aes(x=Midpoint_Age_tau_MesTemp, y=Tau_age_mean)) + geom_point() +
  geom_smooth(method="lm") + theme_grey()
cor.test(subset_taupos_uq$Midpoint_Age_tau_MesTemp, subset_taupos_uq$Tau_age_mean, method = "spearman")
model_check2 <- lm(Tau_age_mean ~ Midpoint_Age_tau_MesTemp, data = subset_taupos_uq)
summary(model_check2)
#Tau time interval vs actual time interval
data_SUVR_int <- data_SUVR_int %>% arrange(RID, SCANDATE)

first_tau_times <- data_SUVR_int %>% # Calculate the difference from the first observation
  group_by(RID) %>%
  summarise(first_tau_time_new = first(Tau_time))

data_SUVR_int <- left_join(data_SUVR_int, first_tau_times, by = "RID")
data_SUVR_int <- data_SUVR_int %>%
  group_by(RID) %>%
  mutate(pre_time_int = Tau_time - first_tau_time_new)

subset_taupos_times <- subset(data_SUVR_int ,  binary_scan_tau==”Postive”)
cor.test(subset_taupos_times$scan_int, subset_taupos_times$pre_time_int, method = "spearman")
ggplot(subset_taupos_times, aes(x=scan_int, y=pre_time_int)) + geom_point(aes(color=binary_scan_tau))  + scale_color_manual(values=c("darkred","blue")) +
  geom_smooth(method="lm") + theme_classic()
model_check3 <- lm(pre_time_int ~ scan_int, data = subset_taupos_times)
summary(model_check3)
#MAE between estimated and actual tau ages
subset_taupos<- subset_taupos %>%
  group_by(RID) %>%
  mutate(Tau_age_mean = ifelse(n() == 1, Tau_age, mean(Tau_age)))

MAE_tauage <- mean(abs(subset_taupos$Midpoint_Age_tau_MesTemp - subset_taupos$Tau_age_mean), na.rm=T)
ggplot(subset_taupos, aes(x=Midpoint_Age_tau_MesTemp, y=Tau_age_mean)) + geom_point() + geom_smooth(method = "lm") +labs( x = "Conversion age",y = "Estimated age") 



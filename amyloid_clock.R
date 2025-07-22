# R/amyloid_clock_main.R

#====================#
# Amyloid Clock 
#====================#

# Load required libraries
required_packages <- c("tidyverse", "nlme", "mgcv", "lme4", "data.table", "here")
lapply(required_packages, library, character.only = TRUE)

# Constants
SUVR_MIN <- 0.62
SUVR_MAX <- 1.11
SUVR_POSITIVE_THRESHOLD <- 0.78005

# Load dataset
data_path <- here("data", "final_dataset_fbp0509.csv")
dataset <- read_csv(data_path)
# Ensure date formatting
dataset$EXAMDATE <- as.Date(dataset$EXAMDATE, format = "%m/%d/%y")

# Calculate conversion age in PET converters (used later for control checks and validation)
dataset  <- dataset %>%
  arrange(RID,EXAMDATE) %>%
  group_by(RID) %>% 
  mutate(Previous_Amyloidstatus_078 = lag(binary_scan078))
conversion_078 <- dataset %>%
  filter(Previous_Amyloidstatus_078  == "Negative", binary_scan078 == "Positive") %>%
  group_by(RID) %>% distinct(RID)
converters078_subset <- dataset %>%
  semi_join(conversion_078, by = "RID")
midpoint_ages078 <- converters078_subset %>%
  group_by(RID) %>%
  summarize(
    Last_0_Age078 = last(Age[binary_scan078 == "Negative"]),  # Age at last 0
    First_1_Age078 = first(Age[binary_scan078 == "Positive"]),  # Age at first 1
    Midpoint_Age1 = (Last_0_Age078 + First_1_Age078) / 2  # Midpoint age as average of ages at last 0 and first 1
  )

dataset <- merge(dataset, midpoint_ages078, by = "RID", all.x = TRUE)
dataset$Amy_time_conv1 <- (dataset$Age - dataset$Midpoint_Age1)

# Filter participants with follow-up
dataset_fu <- dataset %>% 
  group_by(RID) %>% 
  filter(n() > 1) %>% 
  ungroup() %>% 
  arrange(RID, EXAMDATE) %>%
  group_by(RID) %>%
  mutate(first_PET_scan = min(EXAMDATE),
         last_PET_scan = max(EXAMDATE),
         total_fu_time = as.numeric(difftime(last_PET_scan, first_PET_scan, units = "days")) / 365.25,
         scan_int = as.numeric(difftime(EXAMDATE, first_PET_scan, units = "days")) / 365.25)

# LME model for SUVR accumulation rate
fit_SUVRrate <- lme(SUVR_compositeRef ~ scan_int, random = ~1 + scan_int | RID, data = dataset_fu, method = "ML")
slopes_df <- coef(fit_SUVRrate) %>%
  as.data.frame() %>%
  rownames_to_column("RID") %>%
  mutate(RID = as.integer(RID)) %>%
  select(RID, slope = scan_int)

# Add slopes to dataset
dataset_fu <- left_join(dataset_fu, slopes_df, by = "RID")

# Calculate SUVR midpoint
baseline_suvr <- dataset_fu %>% arrange(RID, EXAMDATE) %>% group_by(RID) %>% summarise(Baseline_SUVR = first(SUVR_compositeRef))
dataset_fu <- left_join(dataset_fu, baseline_suvr, by = "RID") %>%
  mutate(SUVR_midpoint = Baseline_SUVR + (total_fu_time / 2) * slope)

dataset_fu_uq <- dataset_fu %>% distinct(RID, .keep_all = TRUE)

# GAM of slope vs SUVR midpoint
fit_slope_gam <- gam(slope ~ s(SUVR_midpoint, bs = "cr"), data = dataset_fu_uq)
summary(fit_slope_gam)
plot(fit_slope_gam)
#plot
ggplot(dataset_fu_uq , aes(x=SUVR_midpoint, y=slopes)) + 
  geom_point(aes(color= APOE_binary))+ scale_color_manual(values=c("darkblue", "forestgreen")) +
  geom_smooth(method="gam") + theme_classic()
# Variance check
dataset_fu_uq <- dataset_fu_uq[order(dataset_fu_uq$SUVR_midpoint), ]

predictions <- predict(fit_SUVRrate, newdata = dataset_fu_uq, se.fit = TRUE)
dataset_fu_uq$variance <- predictions$se.fit^2
mean_variance <- mean(variance, na.rm = TRUE)
sd_variance <- sd(variance, na.rm = TRUE)

cutpoint <- quantile(variance, 0.90) 
high_variance_points <- dataset_fu_uq$SUVR_midpoint[variance > cutpoint]
plot(dataset_fu_uq$SUVR_midpoint, variance, type = "l", lwd = 2,
     ylab = "Variance", xlab = "SUVR Midpoint",
     main = "Variance Across SUVR Midpoint")
abline(h = cutpoint, col = "red", lty = 2)
points(high_variance_points, variance[variance > cutpoint], col = "blue", pch = 19)
# Run GAM again with restricted interval (0.61, 1.11)
data_suvr_interval <- subset(dataset_fu_uq,  between(SUVR_midpoint, SUVR_MIN, SUVR_MAX))
fit_SUVRrate <- gam(slopes ~ s(SUVR_midpoint, bs="cr"), data=
                      data_suvr_interval)
summary(fit_SUVRrate)
plot(fit_SUVRrate)
# Calculate time intervals between SUVR#
##create a data frame with intervals from 0.62 to 1.11##
SUVR_midpoint <- seq(from= SUVR_MIN,to=SUVR_MAX,by=.0001)
SUVR_midpoint <- SUVR_midpoint + 0.00005
SUVR_midpoint <- as.data.frame(SUVR_midpoint)

###extract predicted slopes for all SUVR between 0.6 to 1.21###
SUVR_midpoint$estim_rate_i <- predict(fit_SUVRrate,newdata =SUVR_midpoint, type="response")
SUVR_midpoint$recip_rate_i <- 1/SUVR_midpoint$estim_rate_i
SUVR_midpoint$Amyloid_TIME_INT_i  <- 0.0001*SUVR_midpoint$recip_rate_i

SUVR_midpoint$TimeSum_i <- cumsum(SUVR_midpoint$Amyloid_TIME_INT_i) ###cumulative sum of time###
SUVR_midpoint$Amyloid_time <- SUVR_midpoint$TimeSum_i - SUVR_midpoint$Time[SUVR_midpoint$SUVR_midpoint == SUVR_POSITIVE_THRESHOLD] # Calculate the time difference from the target SUVR for each observation

### Merge with original dataset and estimate amyloid time for the rest of scans between 0.6 and 1.21###
data_SUVR_int <- subset(dataset,  between(SUVR_compositeRef, SUVR_MIN, SUVR_MAX))

data_SUVR_int  <- data_SUVR_int  %>%
rename(SUVR = SUVR_compositeRef)
SUVR_midpoint <- SUVR_midpoint %>%
rename(SUVR = SUVR_midpoint)

seq_dataset <- as.data.frame(SUVR_midpoint)
library(data.table)
setDT(data_SUVR_int)
setDT(seq_dataset )
### Define a function to find the nearest value in df2 for each value in df1###
find_nearest <- function(x, y) {
  idx <- findInterval(x, y, all.inside = TRUE)
  y[ifelse(idx == 0, 1, ifelse(idx == length(y), length(y), idx))]
}

data_SUVR_int[, Nearest_Value := find_nearest(SUVR, seq_dataset$SUVR)] # Add a new column to df1 with the nearest value from df2
data_SUVR_int <- merge(data_SUVR_int, seq_dataset, by.x = "Nearest_Value", by.y = "SUVR", all.x = TRUE)

#estimates amyloid time for all pet scans

model_amy <- gam(Amyloid_time ~ s(SUVR, bs="cr"), data = data_SUVR_int)
summary(model_amy)

data_SUVR_int$Amyloid_time <- predict(model_amy ,newdata = data_SUVR_int, type="response")
# Estimate amyloid age
subset_pos <- data_SUVR_int %>% filter(SUVR > SUVR_POSITIVE_THRESHOLD) %>%
  mutate(Amyloid_age = Age - Amyloid_time) %>%
  group_by(RID) %>%
  mutate(Amyloid_age_mean = mean(Amyloid_age, na.rm = TRUE))

########################################################################################
# Control checks

#estimated vs actual time in converters

cor.test(data_SUVR_int$Amy_time_conv1, data_SUVR_int$Amyloid_time, method = "spearman")
model_check <- lm(Amy_time_conv1 ~ Amyloid_time, data = data_SUVR_int)
summary(model_check)

#estimated amyloid age vs midpoint age for converters
subset_pos_uq <- subset_pos %>%
  distinct(RID, .keep_all = TRUE)
ggplot(subset_pos, aes(x=Midpoint_Age1, y=Amyloid_age_mean)) + geom_point() +
  geom_smooth(method="lm") + theme_grey()
cor.test(subset_pos_uq$Midpoint_Age1, subset_pos_uq$Amyloid_age_mean, method = "spearman")
model_check2 <- lm(Amyloid_age_mean ~ Midpoint_Age1, data = subset_pos_uq)
summary(model_check2)

#amyloid time interval vs actual time interval
data_SUVR_int <- data_SUVR_int %>% arrange(RID, EXAMDATE)

first_amy_times <- data_SUVR_int %>% # Calculate the difference from the first observation
  group_by(RID) %>%
  summarise(first_amy_time_new = first(Amyloid_time))

data_SUVR_int <- left_join(data_SUVR_int, first_amy_times, by = "RID")
data_SUVR_int <- data_SUVR_int %>%
  group_by(RID) %>%
  mutate(pre_time_int = Amyloid_time - first_amy_time_new)

subset_pos_times <- subset(data_SUVR_int , binary_scan1=="Positive")
cor.test(subset_pos_times$scan_int, subset_pos_times$pre_time_int, method = "spearman")
ggplot(subset_pos_times, aes(x=scan_int, y=pre_time_int)) + geom_point(aes(color=binary_scan1))  + scale_color_manual(values=c("darkred","blue")) +
  geom_smooth(method="lm") + theme_classic()
model_check3 <- lm(pre_time_int ~ scan_int, data = subset_pos_times)
summary(model_check3)

#MAE between actual and estimated amyloid ages in converters
MAE_age <- median(abs(subset_pos$Midpoint_Age1 - subset_pos$Amyloid_age_mean), na.rm=T)
ggplot(subset_pos, aes(x=Midpoint_Age1, y=Amyloid_age_mean)) + geom_point() + geom_smooth(method = "lm") +labs( x = "Conversion age",y = "Estimated age")


########################################################################################
# Prediction of age symptom onset##

dataset$firstCDRgt0CI_date <- as.Date(dataset$firstCDRgt0CI_date, format = "%m/%d/%y")
dataset$last_cdr0_date <- as.Date(dataset$last_cdr0_date, format = "%m/%d/%y")

# Exclude converters with negative PET at time of conversion or later
dataset <- dataset %>%
  left_join(subset_pos_uq %>% select(RID, Amyloid_age_mean), by = "RID")

converters_subset <- subset(dataset, converter==TRUE)
ids_with_negative_scan_atconv <- converters_subset %>%
  filter(binary_scan1 == "Negative" & EXAMDATE >= firstCDRgt0CI_date) %>%
  pull(RID)

dataset_conv <- converters_subset %>%
  filter(!RID %in% ids_with_negative_scan_atconv )

dataset_conv_uq <- dataset_conv%>%
  distinct(RID, .keep_all = TRUE)

ggplot(dataset_conv_uq, aes(x=Amyloid_age_mean, y=conversion_age)) + geom_point()  +
  geom_smooth(method="lm") 
cor.test(dataset_conv_uq$Amyloid_age_mean, dataset_conv_uq$conversion_age, method = "spearman") 

EYO_model1<- lm(conversion_age ~ Amyloid_age_mean ,  data=dataset_conv_uq)
summary(EYO_model1)

# Predict estimated conversion age for the rest of participants with at least 1 positive scan (those that have amyloid age calculated
dataset_uq <- dataset %>%
  distinct(RID, .keep_all = TRUE)
dataset_uq$est_conversion_age_amy <- predict(EYO_model1, newdata = dataset_uq, type="response")
dataset <- dataset %>%
  left_join(dataset_uq %>% select(RID, est_conversion_age_amy), by = "RID")

dataset_conv_uq <- dataset_conv_uq %>%
  left_join(dataset_uq %>% select(RID, est_conversion_age_amy), by = "RID")

cor.test(dataset_conv_uq$est_conversion_age_amy, dataset_conv_uq$conversion_age)
ggplot(dataset_uq, aes(x=conversion_age, y=est_conversion_age_amy)) + geom_point()  +
  geom_smooth(method="lm")

model_check4 <- lm(est_conversion_age_amy ~ conversion_age, data = dataset_conv_uq)
summary(model_check4)




ages_dataset <- dataset_uq %>%
  select(RID, Amyloid_age_mean, est_conversion_age_amy) %>%
  full_join(dataset_tau_uq %>% select(RID, Tau_age_mean, est_conversion_age_tau), by = "RID")
write_csv(ages_dataset, "ages_dataset.csv")

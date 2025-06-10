####1. LME models to estimate slopes/rate of change in tau####
tau_dataset <- read_csv("tau_dataset.csv")
dataset <- tau_dataset

##subset with more than one obs##
id_counts <- dataset %>%
  group_by(RID) %>%
  summarise(num_observations = n())

# Filter IDs with more than one observation
ids_with_multiple_observations <- id_counts %>%
  filter(num_observations > 1)

# Filter the original dataset based on IDs with more than one observation
multiple_observation_data <- dataset %>%
  filter(RID %in% ids_with_multiple_observations$RID)

library(nlme)
library(lme4)
library(effects)  
library(performance)

##calculate time intervals
multiple_observation_data$SCANDATE <- as.Date(multiple_observation_data$SCANDATE, format = "%m/%d/%y")

multiple_observation_data <- multiple_observation_data %>%
  arrange(RID,SCANDATE)

first_taupetscan_dates <- multiple_observation_data %>% 
  group_by(RID) %>%
  summarise(first_tauPET_scan = min(SCANDATE))
multiple_observation_data <- left_join(multiple_observation_data, first_taupetscan_dates, by = "RID")

last_taupetscan_dates <- multiple_observation_data %>% ##find the last PET date for each ID and add to dataset
  group_by(RID) %>%
  summarise(last_tauPET_scan = max(SCANDATE))
multiple_observation_data <- left_join(multiple_observation_data, last_taupetscan_dates, by = "RID")

#calculate actual time#
multiple_observation_data$total_tau_fu_time <- ((multiple_observation_data$last_tauPET_scan - multiple_observation_data$first_tauPET_scan)/365.25)
multiple_observation_data$tau_scan_int <- ((multiple_observation_data$SCANDATE- multiple_observation_data$first_tauPET_scan)/365.25)


##Get slopes of tau accumulation##
multiple_observation_data <- multiple_observation_data  %>%
  filter(!is.na(MesialTemporal))

fit_MesialTemprate <- lme(MesialTemporal ~ tau_scan_int,  random= ~1 +tau_scan_int| RID, data=multiple_observation_data , na.action = na.omit, method = "ML")
summary(fit_MesialTemprate) 
plot(fit_MesialTemprate)
coefficients <- coef(fit_MesialTemprate) #exctract slopes#
Tauslopes <- coefficients$tau_scan_int

Tauslopes_df <- data.frame(RID = unique(multiple_observation_data$RID), Slopes = Tauslopes)

multiple_observation_data <- left_join(multiple_observation_data, Tauslopes_df, by = "RID")  #add to dataset#

############### calculate SUVR midpoint####################

multiple_observation_data <- multiple_observation_data %>%
  arrange(RID, SCANDATE) %>%
  group_by(RID) %>%
  mutate(Baseline_MesialTemporal = first(MesialTemporal)) %>%
  ungroup()

multiple_observation_data$MesialTemp_midpoint=((multiple_observation_data$total_tau_fu_time/2)*multiple_observation_data$Slopes )+ multiple_observation_data$Baseline_MesialTemporal
multiple_observation_data$MesialTemp_midpoint <-  as.numeric(multiple_observation_data$MesialTemp_midpoint)

############ 2. plot SUVR midpoinct as function of slopes ###########

multiple_observation_data_uq <- multiple_observation_data %>%
  distinct(RID, .keep_all = TRUE)


ggplot(multiple_observation_data_uq, aes(x=MesialTemp_midpoint, y=Slopes)) + geom_point()+
  geom_smooth(method="gam") +  geom_point(aes(color= APOE_binary))+ scale_color_manual(values=c("darkblue", "forestgreen")) + 
  theme_minimal() + xlab("Tau PET Mesial-Temporal Midpoint SUVR") + ylab("Tau PET Rate of change (SUVR/year)")

############ 3. Calculation of Tau time ###########

library(mgcv)
data_tausuvr_interval <- subset(multiple_observation_data_uq,  between(MesialTemp_midpoint,0.98, 2.04))

fit_TauSUVRrate_1 <- gam(TauSlopes ~ s(MesialTemp_midpoint, bs="cr"), data=data_tausuvr_interval)
summary(fit_TauSUVRrate_1)
plot(fit_TauSUVRrate_1, residuals = T)

predicted <- predict(fit_TauSUVRrate_1, newdata = data_tausuvr_interval)
cor(data_tausuvr_interval$TauSlopes, predicted)
plot(fit_TauSUVRrate_1)
data_tausuvr_interval$estim_rate_tau <- predict(fit_TauSUVRrate_1, type = "response")  

##variance check
multiple_observation_data_uq <- multiple_observation_data_uq[order(multiple_observation_data_uq$MesialTemp_midpoint), ]

predictions <- predict(fit_TauSUVRrate_1, newdata = multiple_observation_data_uq, se.fit = TRUE)
multiple_observation_data_uq$variance <- predictions$se.fit^2
cutpoint <- quantile(variance, 0.90)
high_variance_points <- multiple_observation_data_uq$MesialTemp_midpoint[variance > cutpoint]
plot(multiple_observation_data_uq$MesialTemp_midpoint, variance, type = "l", lwd = 2,
     ylab = "Variance", xlab = "MesialTemporal SUVR Midpoint",
     main = "Variance Across SUVR Midpoint")
abline(h = cutpoint, col = "red", lty = 2)
points(high_variance_points, variance[variance > cutpoint], col = "blue", pch = 19)
#######calculate time intervals between SUVR#######
###create a data frame with intervals from 0.80/ 1st perccentile) to 2.75###
MesialTemp_midpoint <- seq(from= 0.98,to=2.04,by=.0001)
MesialTemp_midpoint <- MesialTemp_midpoint + 0.00005
MesialTemp_midpoint <- as.data.frame(MesialTemp_midpoint)

###extract predicted slopes for all SUVR between 0.6 to 1.21###
MesialTemp_midpoint$estim_rate_i <- predict(fit_TauSUVRrate_1,newdata =MesialTemp_midpoint, type="response")
MesialTemp_midpoint$recip_rate_i <- 1/MesialTemp_midpoint$estim_rate_i
MesialTemp_midpoint$Tau_TIME_INT_i  <- 0.0001*MesialTemp_midpoint$recip_rate_i

MesialTemp_midpoint$TimeSum_i <- cumsum(MesialTemp_midpoint$Tau_TIME_INT_i) ###cumulative sum of time###
MesialTemp_midpoint$MesialTemp_midpoint <- round(MesialTemp_midpoint$MesialTemp_midpoint, 5)
target_SUVR <- 1.41005

MesialTemp_midpoint$Tau_time <- MesialTemp_midpoint$TimeSum_i - MesialTemp_midpoint$TimeSum_i[MesialTemp_midpoint$MesialTemp_midpoint == target_SUVR] # Calculate the time difference from the target SUVR for each observation
write_csv(MesialTemp_midpoint, "MesialTemp_midpoint.csv")
### Merge with original dataset and estimate tau time for the rest of scans ###

tau_dataset <- tau_dataset %>%
  rename(SUVR = MesialTemporal)
MesialTemp_midpoint <- MesialTemp_midpoint %>%
  rename(SUVR = MesialTemp_midpoint)

seq_dataset <- as.data.frame(MesialTemp_midpoint)
library(data.table)
setDT(tau_dataset)
setDT(seq_dataset )
### Define a function to find the nearest value in df2 for each value in df1###
find_nearest <- function(x, y) {
  idx <- findInterval(x, y, all.inside = TRUE)
  y[ifelse(idx == 0, 1, ifelse(idx == length(y), length(y), idx))]
}

tau_dataset[, Nearest_Value := find_nearest(SUVR, seq_dataset$SUVR)] # Add a new column to df1 with the nearest value from df2
tau_dataset <- merge(tau_dataset, seq_dataset, by.x = "Nearest_Value", by.y = "SUVR", all.x = TRUE)


##estimates tau time for all pet scans

tau_interval <- subset(tau_dataset, SUVR >0.98 & SUVR < 2.04)
model_tau <- gam(Tau_time ~ s(SUVR, bs="cr"), data = tau_interval)
summary(model_tau)

tau_interval$Tau_time <- predict(model_tau ,newdata = tau_interval, type="response")

subset_taupos <- tau_interval %>% 
  filter(SUVR > 1.41)

subset_taupos$Tau_age <- (subset_taupos$Age - subset_taupos$Tau_time)
#plot estimated ages vs estimated amyloid time
ggplot(subset_taupos, aes(x=Tau_time, y=Tau_age, group=RID, color=APOE_binary)) + geom_point() + geom_line() 
install.packages("lmerTest")  
model_bias <- lmer(Tau_age ~ Tau_time + (1 | RID), data = subset_taupos)
summary(model_bias)

#get individual slopes for amyloid age
slopes_by_id_tau <- subset_taupos %>%
  group_by(RID) %>%
  do({
    model <- lm(Tau_age ~ Tau_time, data = .)
    data.frame(slope = coef(model)[2])  # Extract the slope (amyloid_time coefficient)
  })
ggplot(slopes_by_id_tau, aes(x = slope)) +
  geom_histogram(aes(y = ..count../sum(..count..)),  # Calculate percentages
                 binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +  # Show y-axis as percentages
  labs(title = "Distribution of Random Slopes for Tau Age vs. Tau Time",
       x = "Slope",
       y = "Percentage of IDs") +
  theme_minimal() + xlim(-10, +10)
median(slopes_by_id_tau$slope, na.rm=T)
wilcox.test(slopes_by_id_tau$slope, mu = 0)

## get sd for estimated ages
tau_age_sd <- subset_taupos  %>%  group_by(RID) %>%   summarise(sd_age_tau = sd(Tau_age, na.rm=T))
ggplot(tau_age_sd , aes(x=sd_age_tau)) + geom_histogram(aes(y=(..count..) / sum(..count..)*100),fill="lightblue", color="black", bins=10) +
  labs(title="Percentage of people by SD of amyloid age", x= "SD of amyloid age", y="percentage of people") +
  theme_minimal()

median(tau_age_sd$sd_age_tau, na.rm=T)
wilcox.test(tau_age_sd$sd_age_tau, mu = 0)

#MAE between ages
subset_taupos<- subset_taupos %>%
  group_by(RID) %>%
  mutate(Tau_age_mean = ifelse(n() == 1, Tau_age, mean(Tau_age)))

MAE_tauage <- mean(abs(subset_taupos$Midpoint_Age_tau_MesTemp - subset_taupos$Tau_age_mean), na.rm=T)
ggplot(subset_taupos, aes(x=Midpoint_Age_tau_MesTemp, y=Tau_age_mean)) + geom_point() + geom_smooth(method = "lm") +labs( x = "Conversion age",y = "Estimated age") 


###control checks###
##plot estimated vs actual time in converters

cor.test(tau_interval$Tau_time_convMesTemp, tau_interval$Tau_time, method = "spearman")

model_check <- lm(Tau_time_convMesTemp ~ Tau_time, data = tau_interval)
summary(model_check)

###estimated amyloid age vs midpoint age for converters###
subset_taupos_uq <- subset_taupos %>%
  distinct(RID, .keep_all = TRUE)
ggplot(subset_taupos_uq, aes(x=Midpoint_Age_tau_MesTemp, y=Tau_age_mean)) + geom_point() +
  geom_smooth(method="lm") + theme_grey()
cor.test(subset_taupos_uq$Midpoint_Age_tau_MesTemp, subset_taupos_uq$Tau_age_mean, method = "spearman")
model_check2 <- lm(Tau_age_mean ~Midpoint_Age_tau_MesTemp, data = subset_taupos_uq)
summary(model_check2)
##amyloid time interval vs actual time interval
tau_interval$SCANDATE <- as.Date(tau_interval$SCANDATE, format = "%m/%d/%y")

tau_interval <- tau_interval %>% arrange(RID, SCANDATE)

first_tau_times <- tau_interval%>% # Calculate the difference from the first observation
  group_by(RID) %>%
  summarise(first_tau_time_new = first(Tau_time))

tau_interval <- left_join(tau_interval, first_tau_times, by = "RID")
tau_interval <- tau_interval %>%
  group_by(RID) %>%
  mutate(pre_time_int = Tau_time - first_tau_time_new)

subset_pos_tautimes <- subset(tau_interval , binary_scan_tau=="1")
cor.test(subset_pos_tautimes$Tauscan_int, subset_pos_tautimes$pre_time_int, method = "spearman")
ggplot(subset_pos_tautimes, aes(x=Tauscan_int, y=pre_time_int)) + geom_point(aes(color=binary_scan_tau))  + scale_color_manual(values=c("darkred","blue")) +
  geom_smooth(method="lm") + theme_classic()
model_check3 <- lm(pre_time_int ~ Tauscan_int, data = subset_pos_tautimes)
summary(model_check3)
##prediction of age symptom onset##

converters_all <- read_csv("/Users/martamilaaloma/Documents/Datasets/ADNI/FNIH project datasets/converters_final1106.csv")
converters_all_uq <- converters_all %>%
  distinct(RID, .keep_all = TRUE)
tau_dataset$SCANDATE <- as.Date(tau_dataset$SCANDATE, format = "%m/%d/%y")
tau_dataset <- tau_dataset %>%
  left_join(subset_taupos_uq %>% select(RID, Tau_age_mean), by = "RID")
tau_dataset <- tau_dataset %>%
  left_join(converters_all_uq %>% select(RID,conversion_age,firstCDRgt0CI_date, last_cdr0_date), by = "RID")

tau_dataset$firstCDRgt0CI_date <- as.Date(tau_dataset$firstCDRgt0CI_date, format = "%m/%d/%y")
tau_dataset$last_cdr0_date <- as.Date(tau_dataset$last_cdr0_date, format = "%m/%d/%y")
tau_dataset$EXAMDATE_amy <- as.Date(tau_dataset$EXAMDATE_amy, format = "%m/%d/%y")

tau_dataset <- tau_dataset %>%
  mutate(converter = !is.na(conversion_age))

###exclude converters with negative PET at time of conversion or later###
converters_tausubset <- subset(tau_dataset, converter==TRUE)
ids_with_negative_scan_atconvtau <- converters_tausubset %>%
  filter(binary_scan1 == "Negative" & EXAMDATE_amy >= firstCDRgt0CI_date) %>%
  pull(RID)

dataset_conv_tau <- converters_tausubset %>%
  filter(!RID %in% ids_with_negative_scan_atconvtau )

dataset_conv_tau_uq <- dataset_conv_tau%>%
  distinct(RID, .keep_all = TRUE)



ggplot(dataset_conv_tau_uq, aes(x=Tau_age_mean, y=conversion_age)) + geom_point()  +
  geom_smooth(method="lm") 
cor.test(dataset_conv_tau_uq$Tau_age_mean, dataset_conv_tau_uq$conversion_age) 

EYO_model1_tau<- lm(conversion_age ~ Tau_age_mean ,  data=dataset_conv_tau_uq)
summary(EYO_model1_tau)

###predict estimated conversion age for the rest of participants with at least 1 positive scan (those that have amyloid age calculated##
dataset_tau_uq <- tau_dataset %>%
  distinct(RID, .keep_all = TRUE)
dataset_tau_uq$est_conversion_age_tau <- predict(EYO_model1_tau, newdata = dataset_tau_uq, type="response")
tau_dataset <- tau_dataset %>%
  left_join(dataset_tau_uq %>% select(RID, est_conversion_age_tau), by = "RID")


dataset_conv_tau_uq <- dataset_conv_tau_uq %>%
  left_join(dataset_tau_uq %>% select(RID, est_conversion_age_tau), by = "RID")

cor.test(dataset_conv_tau_uq$est_conversion_age_tau, dataset_conv_tau_uq$conversion_age)
ggplot(dataset_conv_tau_uq, aes(x=conversion_age, y=est_conversion_age_tau)) + geom_point()  +
  geom_smooth(method="lm")

model_check4 <- lm(est_conversion_age_tau ~ conversion_age, data = dataset_conv_tau_uq)
summary(model_check4)


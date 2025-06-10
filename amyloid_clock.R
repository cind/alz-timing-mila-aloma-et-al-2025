###################### AMYLOID CLOCLK #################
setwd("/Users/martamilaaloma/Documents/Datasets/ADNI/FNIH project datasets/last")
final_dataset_fbp <- read_csv("final_dataset_fbp0509.csv")
dataset <- final_dataset_fbp


###### LME to get rates of change ####
#get rows with multiple obs
id_counts <- dataset %>%
  group_by(RID) %>%
  summarise(num_observations = n())

ids_with_multiple_observations <- id_counts %>%
  filter(num_observations > 1)

dataset_fu <- dataset %>%
  filter(RID %in% ids_with_multiple_observations$RID)

dataset_fu$EXAMDATE <- as.Date(dataset_fu$EXAMDATE, format = "%m/%d/%y")


##calculate time intervals
dataset_fu <- dataset_fu %>%
  arrange(RID,EXAMDATE)

dataset_fu <- dataset_fu %>% 
  group_by(RID) %>%
  mutate(first_PET_scan = min(EXAMDATE))

dataset_fu <- dataset_fu %>% 
  group_by(RID) %>%
  mutate(last_PET_scan = max(EXAMDATE))

#calculate actual time#
dataset_fu$total_fu_time <- ((dataset_fu$last_PET_scan - dataset_fu$first_PET_scan)/365.25)
dataset_fu$scan_int <- ((dataset_fu$EXAMDATE- dataset_fu$first_PET_scan)/365.25)

#get rates of change#
library(nlme)

fit_SUVRrate <- lme(SUVR_compositeRef ~ scan_int,  random= ~1 +scan_int| RID, data=dataset_fu, na.action = na.omit, method = "ML")
summary(fit_SUVRrate)
plot(fit_SUVRrate)

#extract slope 
coefficients <- coef(fit_SUVRrate) 
slopes <- coefficients$scan_int

slopes_df <- data.frame(RID = unique(dataset_fu$RID), slopes = slopes)
slopes_df <- slopes_df %>% rename(slopes = slopes)

dataset_fu <- left_join(dataset_fu, slopes_df, by = "RID")  


########calculate SUVR midpoint######

dataset_fu <- dataset_fu %>%
  arrange(RID, EXAMDATE) %>%
  group_by(RID) %>%
  mutate(Baseline_SUVR = first(SUVR_compositeRef)) %>%
  ungroup()
dataset_fu$SUVR_midpoint=((dataset_fu$total_fu_time/2)*dataset_fu$slopes)+ dataset_fu$Baseline_SUVR
dataset_fu$SUVR_midpoint <- as.numeric(dataset_fu$SUVR_midpoint)
########rates of change as function of SUVR midpoint ######
dataset_fu_uq <- dataset_fu %>%
distinct(RID, .keep_all = TRUE)

ggplot(dataset_fu_uq , aes(x=SUVR_midpoint, y=slopes)) + 
  geom_point(aes(color= APOE_binary))+ scale_color_manual(values=c("darkblue", "forestgreen")) +
  geom_smooth(method="gam") + theme_classic()

fit_SUVRrate <- gam(slopes ~ s(SUVR_midpoint, bs="cr"), data=
                      dataset_fu_uq)
summary(fit_SUVRrate)
plot(fit_SUVRrate)

##variance check
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

###run GAM again with restricted interval
data_suvr_interval <- subset(dataset_fu_uq,  between(SUVR_midpoint,0.62, 1.11))
fit_SUVRrate <- gam(slopes ~ s(SUVR_midpoint, bs="cr"), data=
                      data_suvr_interval)
summary(fit_SUVRrate)
plot(fit_SUVRrate)

#######calculate time intervals between SUVR#######
###create a data frame with intervals from 0.60/ 1st percentile) to 1.21###

SUVR_midpoint <- seq(from= 0.62,to=1.11,by=.0001)
SUVR_midpoint <- SUVR_midpoint + 0.00005
SUVR_midpoint <- as.data.frame(SUVR_midpoint)

###extract predicted slopes for all SUVR between 0.6 to 1.21###
SUVR_midpoint$estim_rate_i <- predict(fit_SUVRrate,newdata =SUVR_midpoint, type="response")
SUVR_midpoint$recip_rate_i <- 1/SUVR_midpoint$estim_rate_i
SUVR_midpoint$Amyloid_TIME_INT_i  <- 0.0001*SUVR_midpoint$recip_rate_i

SUVR_midpoint$TimeSum_i <- cumsum(SUVR_midpoint$Amyloid_TIME_INT_i) ###cumulative sum of time###
target_SUVR <- 0.78005

SUVR_midpoint$Amyloid_time <- SUVR_midpoint$TimeSum_i - SUVR_midpoint$Time[SUVR_midpoint$SUVR_midpoint == target_SUVR] # Calculate the time difference from the target SUVR for each observation

### Merge with original dataset and estimate amyloid time for the rest of scans between 0.6 and 1.21###
data_SUVR_int <- subset(dataset,  between(SUVR_compositeRef, 0.62, 1.11))

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

#####estimates amyloid time for all pet scans#####

model_amy <- gam(Amyloid_time ~ s(SUVR, bs="cr"), data = data_SUVR_int)
summary(model_amy)

data_SUVR_int$Amyloid_time <- predict(model_amy ,newdata = data_SUVR_int, type="response")

#####calculate age at amyloid onset and test bias #####
subset_pos <- data_SUVR_int %>% 
    filter(SUVR > 0.78)

subset_pos$Amyloid_age <- (subset_pos$Age - subset_pos$Amyloid_time)

subset_pos<- subset_pos %>%
  group_by(RID) %>%
  mutate(Amyloid_age_mean = ifelse(n() == 1, Amyloid_age, mean(Amyloid_age)))

library(lme4)
model_bias <- lmer(Amyloid_age ~ Amyloid_time + (1 | RID), data = subset_pos)
summary(model_bias)
#get individual slopes for amyloid age
slopes_by_id <- subset_pos %>%
  group_by(RID) %>%
  do({
    model <- lm(Amyloid_age ~ Amyloid_time, data = .)
    data.frame(slope = coef(model)[2])  # Extract the slope (amyloid_time coefficient)
  })

subset_pos <- subset_pos %>%
  left_join(slopes_by_id, by = "RID")
ggplot(slopes_by_id, aes(x = slope)) +
  geom_histogram(aes(y = ..count../sum(..count..)),  # Calculate percentages
                 binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +  # Show y-axis as percentages
  labs(title = "Distribution of Random Slopes for Amyloid Age vs. Amyloid Time",
       x = "Slope",
       y = "Percentage of IDs") +
  theme_minimal() + xlim(-10, +10)
median(slopes_by_id$slope, na.rm=T)
wilcox.test(slopes_by_id$slope, mu = 0)

#plot slopes vs rates
dataset_fu_uq <- dataset_fu_uq %>%
  left_join(slopes_by_id, by = "RID")
ggplot(dataset_fu_uq, aes(x=slopes, y=slope)) + geom_point() + labs( x = "rates of accumulation",y = "individual slopes of estimated ages") + ylim(-20,20)
#plot slopes  vs suvr midpoint
ggplot(dataset_fu_uq, aes(x=SUVR_midpoint, y=slope)) + geom_point() + labs( x = "SUVR midpoint",y = "individual slopes of estimated ages") + ylim(-20,20)


## get sd for estimated ages
amy_age_sd <- subset_pos %>%  group_by(RID) %>%   summarise(sd_age_amy = sd(Amyloid_age, na.rm=T))
ggplot(amy_age_sd , aes(x=sd_age_amy)) + geom_histogram(aes(y=(..count..) / sum(..count..)*100),fill="lightblue", color="black", bins=10) +
labs(title="Percentage of people by SD of amyloid age", x= "SD of amyloid age", y="percentage of people") +
theme_minimal()

median(amy_age_sd$sd_age_amy, na.rm=T)
wilcox.test(amy_age_sd$sd_age_amy, mu = 0)
#plot slopes vs sd
dataset_fu_uq <- dataset_fu_uq %>%
  left_join(amy_age_sd, by = "RID")
ggplot(dataset_fu_uq, aes(x=slopes, y=sd_age_amy)) + geom_point() + labs( x = "rates of accumulation",y = "SD of estimated ages") 

#plot midpoint suvr vs sd
ggplot(dataset_fu_uq, aes(x=SUVR_midpoint, y=sd_age_amy)) + geom_point() + labs( x = "SUVR midpoint",y = "SD of estimated ages") 

#plot estimated ages vs estimated amyloid time
ggplot(subset_pos, aes(x=Amyloid_time, y=Amyloid_age, group=RID, color=APOE_binary)) + geom_point() + geom_line() 

#MAE between actual and estimated amyloid ages in converters
MAE_age <- median(abs(subset_pos$Midpoint_Age1 - subset_pos$Amyloid_age_mean), na.rm=T)
ggplot(subset_pos, aes(x=Midpoint_Age1, y=Amyloid_age_mean)) + geom_point() + geom_smooth(method = "lm") +labs( x = "Conversion age",y = "Estimated age") 



###control checks###
##plot estimated vs actual time in converters

cor.test(data_SUVR_int$Amy_time_conv1, data_SUVR_int$Amyloid_time, method = "spearman")

model_check <- lm(Amy_time_conv1 ~ Amyloid_time, data = data_SUVR_int)
summary(model_check)

###estimated amyloid age vs midpoint age for converters###
subset_pos_uq <- subset_pos %>%
  distinct(RID, .keep_all = TRUE)
ggplot(subset_pos, aes(x=Midpoint_Age1, y=Amyloid_age_mean)) + geom_point() +
  geom_smooth(method="lm") + theme_grey()
cor.test(subset_pos_uq$Midpoint_Age1, subset_pos_uq$Amyloid_age_mean, method = "spearman")
model_check2 <- lm(Amyloid_age_mean ~ Midpoint_Age1, data = subset_pos_uq)
summary(model_check2)
##amyloid time interval vs actual time interval
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


##prediction of age symptom onset##

dataset$firstCDRgt0CI_date <- as.Date(dataset$firstCDRgt0CI_date, format = "%m/%d/%y")
dataset$last_cdr0_date <- as.Date(dataset$last_cdr0_date, format = "%m/%d/%y")

###exclude converters with negative PET at time of conversion or later###
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

###predict estimated conversion age for the rest of participants with at least 1 positive scan (those that have amyloid age calculated##
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


write_csv(dataset_uq, "dataset_uq_amy.csv")
write_csv(dataset_tau_uq, "dataset_tau_uq.csv")

ages_dataset <- dataset_uq %>%
  select(RID, Amyloid_age_mean, est_conversion_age_amy) %>%
  full_join(dataset_tau_uq %>% select(RID, Tau_age_mean, est_conversion_age_tau), by = "RID")
write_csv(ages_dataset, "ages_dataset.csv")

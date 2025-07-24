# R/trajectory_plots.R

#====================#
# Biomarker trajectories plots
#====================#

# Load required libraries
library(tidyverse)
library(cowplot)


final_dataset_plasma  <- read_csv("dataset_plasma.csv")
final_dataset_plasma$EXAMDATE <- as.Date(final_dataset_plasma$EXAMDATE, format = "%m/%d/%y")
final_dataset_plasma$SCANDATE_AMY <- as.Date(final_dataset_plasma$SCANDATE_AMY, format = "%m/%d/%y")
final_dataset_plasma$CDR <- factor(final_dataset_plasma$CDR.x, levels=c(0,0.5,1, 2,3))
final_dataset_plasma$APOE_binary <- factor(final_dataset_plasma$APOE_binary, levels=c("carrier", "non-carrier"), labels = c(paste("APOE-\u03B54 carrier"), paste("APOE-\u03B54 non-carrier"))) ###
final_dataset_plasma$CDR_DX_Group <- factor(final_dataset_plasma$CDR_DX_Group , levels=c( "CN_CDR0", "CI_AD_CDRgt0", "CI_nonAD_CDRgt0"))

final_dataset_plasma<-final_dataset_plasma %>%
       mutate(CDR_3 = case_when(
            CDR == 0 ~ "0",
            CDR == 0.5 ~ "0.5",
             CDR %in% c("1", "2", "3") ~ ">=1"
         ))
final_dataset_plasma$CDR_3 <- factor(final_dataset_plasma$CDR_3, levels = c("0", "0.5", ">=1"), labels = c("CDR = 0", "CDR = 0.5", "CDR >= 1"))

final_dataset_plasma <- final_dataset_plasma %>%
 arrange(RID,EXAMDATE)


biomarker_labels <- c(
  C2N_plasma_Abeta42_Abeta40_out = expression(paste("C2N A",beta,"42/A",beta,"40")),
  Fuji_plasma_Ab42_Ab40_out = expression(paste("Fujirebio A",beta,"42/A",beta,"40")),
  Roche_plasma_Ab42_Ab40_out = expression(paste("Roche A",beta, "42/A", beta,"40")),
  QX_plasma_Ab42_Ab40_out = expression(paste("Quanterix A",beta, "42/A",beta,"40")),
  Roche_plasma_ptau181 = "Roche p-tau181 (pg/ml)",
  QX_plasma_ptau181_out = "Quanterix p-tau181 (pg/ml)",
  C2N_plasma_ptau217_out = "C2N p-tau217 (pg/ml)",
  Fuji_plasma_ptau217_out = "Fujirebio p-tau217 (pg/ml)",
  AlzPath_plasma_ptau217_out = "ALZPath p-tau217 (pg/ml)",
  Janssen_plasma_ptau217_out = "Janssen p-tau217 (pg/ml)",
  C2N_plasma_ptau217_ratio = "C2N %p-tau217 (%)",
  Roche_plasma_GFAP_out = "Roche GFAP (ng/ml)",
  QX_plasma_GFAP_out = "Quanterix GFAP (pg/ml)",
  Roche_plasma_NfL_out = "Roche NfL (pg/ml)",
  QX_plasma_NfL_out = "Quanterix NfL (pg/ml)",
  PTAU_over_ABETA42 = expression(paste("CSF p-tau181/A", beta, "42")),
  metaROI.AgeAdj= "Cortical thickness meta-ROI (mm)",
  CDRSB= "CDR-SB",
  MesialTemporal = expression(""^18*"F-flortaucipir mesial-temporal tau PET (SUVR)"),
  TemporoParietal = expression(""^18*"F-flortaucipir temporo-parietal tau PET (SUVR)"),
  SUVR_compositeRef = expression(""^18*"F-florbetapir amyloid PET (SUVR)"))

#######plot AGE ##########

create_spaghetti_plot_age <- function(data,Age,dep_var, CDR_3) {
  
  # Create segments from filtered data
  segments <- data%>%
    group_by(RID) %>%
    arrange(!!sym( Age)) %>%
    mutate(
      Age_next = lead(!!sym( Age)),
      plasma_bmks_next = lead(!!sym(dep_var)),
      group_next = lead(CDR_3)
    ) %>%
    filter(!is.na(Age), !is.na(group_next)) %>%
    ungroup()
  
  # Create two segments for each transition
  segments <- segments %>%
    mutate(
      Age_mid = (!!sym(Age) + Age_next) / 2,
      plasma_bmks_mid = (!!sym(dep_var) +  plasma_bmks_next) / 2
    )
  
  # Create first half segments
  first_half <- segments %>%
    transmute(
      Age = !!sym(Age),
      plasma_bmks = !!sym(dep_var),
      Age_end = Age_mid,
      plasma_bmks_end = plasma_bmks_mid,
      CDR_3 = CDR_3
    )
  
  # Create second half segments
  second_half <- segments %>%
    transmute(
      Age = Age_mid,
      plasma_bmks = plasma_bmks_mid,
      Age_end = Age_next,
      plasma_bmks_end = plasma_bmks_next,
      CDR_3 = group_next
    )
  
  # Combine both halves
  all_segments <- bind_rows(first_half, second_half)
  
  # Plot using geom_segment and geom_point
  plot <- ggplot() +
    geom_segment(data = all_segments, aes(
      x = Age, y = plasma_bmks,
      xend = Age_end, yend = plasma_bmks_end,
      color = CDR_3
    ), size = 0.3, show.legend =F) +
    geom_point(data = data, aes(x = !!sym(Age), y = !!sym(dep_var), color = CDR_3, shape = APOE_binary), size = 1, show.legend = F) +
    
    theme_classic() +
    labs(x = "Age (years)",
         y = biomarker_labels[[dep_var]],
         color = "Group", shape="APOE_binary") + theme( axis.title =  element_text(size = 12), legend.text = element_text(size=12, hjust = 0)) +
    scale_shape_manual(values = c(16, 17), labels = c(expression(italic("APOE-")~epsilon*"4 non-carrier"), expression(italic("APOE-")~epsilon*"4 carrier"))) +
    scale_color_manual(values = c("blue2","darkorange", "brown3" ), labels = c("CDR = 0","CDR = 0.5", "CDR >= 1"))+ xlim(50,105) + 
    guides(color = guide_legend(title = NULL, override.aes = list(size = 2), order = 1), shape = guide_legend(title = NULL,override.aes = list(size = 2), order = 2))
  
  return(plot)
}

# List of dependent variables
dep_var <- c( "C2N_plasma_ptau217_out", "Fuji_plasma_ptau217_out",  "AlzPath_plasma_ptau217_out", "Janssen_plasma_ptau217_out", 
               "C2N_plasma_ptau217_ratio", "QX_plasma_ptau181_out","Roche_plasma_ptau181",
               "C2N_plasma_Abeta42_Abeta40_out" ,"Fuji_plasma_Ab42_Ab40_out","Roche_plasma_Ab42_Ab40_out","QX_plasma_Ab42_Ab40_out", 
               "Roche_plasma_GFAP_out", "QX_plasma_GFAP_out","Roche_plasma_NfL_out","QX_plasma_NfL_out")

# Apply the function to each dependent variable
plot_list_age_plasma <- lapply(dep_var, function(dep_var) {
  create_spaghetti_plot_age(final_dataset_plasma, "Age", dep_var, "color = CDR_3")
})


###nonplasma biomarkers
dep_var <- "SUVR_compositeRef"
plot_list_age_amyPET <- lapply(dep_var, function(dep_var) {
  create_spaghetti_plot_age(final_dataset_fbp, "Age", dep_var, "color = CDR_3")
})

tau_dataset_cut <- subset(tau_dataset, MesialTemporal < 4)
dep_var <- "MesialTemporal"
plot_list_age_tauPET <- lapply(dep_var, function(dep_var) {
  create_spaghetti_plot_age(tau_dataset_cut, "Age", dep_var, "color = CDR_3")
})


dep_var <- "TemporoParietal"
plot_list_age_tauTPPET <- lapply(dep_var, function(dep_var) {
  create_spaghetti_plot_age(tau_dataset_cut, "Age", dep_var, "color = CDR_3")
})

dep_var <- "PTAU_over_ABETA42"
plot_list_age_CSF <- lapply(dep_var, function(dep_var) {
  create_spaghetti_plot_age(csf_dataset, "Age", dep_var, "color = CDR_3")
})
dep_var <- "metaROI.AgeAdj"
plot_list_age_MRI <- lapply(dep_var, function(dep_var) {
  create_spaghetti_plot_age(MRI_dataset, "Age", dep_var, "color = CDR_3")
})

dep_var <- "CDRSB"
plot_list_age_CDR <- lapply(dep_var, function(dep_var) {
  create_spaghetti_plot_age(CDR_dataset, "Age", dep_var, "color = CDR_3")
})



blank_plot <- ggplot() + theme_void()
row1 <- plot_grid(plot_list_age_plasma[[1]], plot_list_age_plasma[[2]], plot_list_age_plasma[[3]], plot_list_age_plasma[[4]], nrow = 1)
row2 <- plot_grid(plot_list_age_plasma[[5]], plot_list_age_plasma[[6]], plot_list_age_plasma[[7]], blank_plot, nrow = 1)
row3 <- plot_grid(plot_list_age_plasma[[8]], plot_list_age_plasma[[9]], plot_list_age_plasma[[10]], plot_list_age_plasma[[11]],  nrow = 1)
row4 <- plot_grid(plot_list_age_plasma[[12]],plot_list_age_plasma[[13]], plot_list_age_plasma[[14]], plot_list_age_plasma[[15]], nrow = 1)
row5 <- plot_grid(plot_list_age_CSF[[1]], plot_list_age_amyPET[[1]],plot_list_age_tauPET[[1]], plot_list_age_tauTPPET[[1]],nrow = 1)
row6 <- plot_grid(plot_list_age_MRI[[1]],plot_list_age_CDR[[1]], blank_plot, blank_plot,nrow = 1)

allplots_age <- plot_grid(row1, row2, row3, row4, row5, row6, nrow = 6)

get_legend <- function(my_plot) {
  g <- ggplotGrob(my_plot)
  legend <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")][[1]]
  legend
}
legend <- get_legend(plot_list_age[[1]])
legend_col <- plot_grid(legend, blank_plot, ncol = 1, rel_heights = c(0.20, 0.80))
final_plot <- plot_grid(allplots_age, legend_col, ncol = 2, rel_widths = c(0.8, 0.2))

ggsave("~/Documents/FNIH_paper_Amyclockplasma/Annals_Neurology_submission/Review/Supp_fig8.pdf", plot = allplots_age ,  width = 17, height = 20)



#######plot Ab PET SUVR ##########
create_spaghetti_plot_amypet <- function(data,SUVR_compositeRef,dep_var, CDR_3) {
  
  # Create segments from filtered data
  segments <- data%>%
    group_by(RID) %>%
    arrange(!!sym(SUVR_compositeRef)) %>%
    mutate(
      SUVR_compositeRef_next = lead(!!sym(SUVR_compositeRef)),
      plasma_bmks_next = lead(!!sym(dep_var)),
      group_next = lead(CDR_3)
    ) %>%
    filter(!is.na(SUVR_compositeRef), !is.na(group_next)) %>%
    ungroup()
  
  # Create two segments for each transition
  segments <- segments %>%
    mutate(
      SUVR_compositeRef_mid = (!!sym(SUVR_compositeRef) + SUVR_compositeRef_next) / 2,
      plasma_bmks_mid = (!!sym(dep_var) +  plasma_bmks_next) / 2
    )
  
  # Create first half segments
  first_half <- segments %>%
    transmute(
      SUVR_compositeRef = !!sym(SUVR_compositeRef),
      plasma_bmks = !!sym(dep_var),
      SUVR_compositeRef_end = SUVR_compositeRef_mid,
      plasma_bmks_end = plasma_bmks_mid,
      CDR_3 = CDR_3
    )
  
  # Create second half segments
  second_half <- segments %>%
    transmute(
      SUVR_compositeRef = SUVR_compositeRef_mid,
      plasma_bmks = plasma_bmks_mid,
      SUVR_compositeRef_end = SUVR_compositeRef_next,
      plasma_bmks_end = plasma_bmks_next,
      CDR_3 = group_next
    )
  
  # Combine both halves
  all_segments <- bind_rows(first_half, second_half)
  
  # Plot using geom_segment and geom_point
  plot <- ggplot() +
    geom_segment(data = all_segments, aes(
      x = SUVR_compositeRef, y = plasma_bmks,
      xend = SUVR_compositeRef_end, yend = plasma_bmks_end,
      color = CDR_3
    ), size = 0.3, show.legend =F ) +
    geom_point(data = data, aes(x = !!sym(SUVR_compositeRef), y = !!sym(dep_var), color = CDR_3, shape = APOE_binary), size = 1, show.legend = F) +
    
    theme_classic() +
    labs(x = expression(""^18*"F-florbetapir amyloid PET (SUVR)"),
         y = biomarker_labels[[dep_var]],
         color = "Group", shape="APOE_binary") + theme( axis.title =  element_text(size = 12), legend.text = element_text(size=12, hjust = 0)) +
    scale_shape_manual(values = c(16, 17), labels = c(expression(italic("APOE-")~epsilon*"4 non-carrier"), expression(italic("APOE-")~epsilon*"4 carrier"))) +
    scale_color_manual(values = c("blue2","darkorange", "brown3" ), labels = c("CDR = 0","CDR = 0.5", "CDR >= 1"))+
    guides(color = guide_legend(title = NULL, override.aes = list(size = 2), order = 1), shape = guide_legend(title = NULL,override.aes = list(size = 2), order = 2))
  
  return(plot)
}


# List of dependent variables
dep_var <- c( "C2N_plasma_ptau217_out", "Fuji_plasma_ptau217_out",  "AlzPath_plasma_ptau217_out", "Janssen_plasma_ptau217_out", 
              "C2N_plasma_ptau217_ratio", "QX_plasma_ptau181_out","Roche_plasma_ptau181",
              "C2N_plasma_Abeta42_Abeta40_out" ,"Fuji_plasma_Ab42_Ab40_out","Roche_plasma_Ab42_Ab40_out","QX_plasma_Ab42_Ab40_out", 
              "Roche_plasma_GFAP_out", "QX_plasma_GFAP_out","Roche_plasma_NfL_out","QX_plasma_NfL_out")

# Apply the function to each dependent variable
plot_list_amypet_plasma <- lapply(dep_var, function(dep_var) {
  create_spaghetti_plot_amypet(final_dataset_plasma, "SUVR_compositeRef", dep_var, "color = CDR_3")
})


###nonplasma biomarkers
final_dataset_fbp$MesialTemporal[final_dataset_fbp$MesialTemporal > 4] <- NA
final_dataset_fbp$TemporoParietal[final_dataset_fbp$TemporoParietal > 4] <- NA
dep_var <- c("PTAU_over_ABETA42", "SUVR_compositeRef", "MesialTemporal", "TemporoParietal", "metaROI.AgeAdj", "CDRSB")
plot_list_amypet_amyPET <- lapply(dep_var, function(dep_var) {
  create_spaghetti_plot_amypet(final_dataset_fbp, "SUVR_compositeRef", dep_var, "color = CDR_3")
})


blank_plot <- ggplot() + theme_void()
row1 <- plot_grid(plot_list_amypet_plasma[[1]], plot_list_amypet_plasma[[2]], plot_list_amypet_plasma[[3]], plot_list_amypet_plasma[[4]], nrow = 1)
row2 <- plot_grid(plot_list_amypet_plasma[[5]],plot_list_amypet_plasma[[6]], plot_list_amypet_plasma[[7]], blank_plot, nrow = 1)
row3 <- plot_grid(plot_list_amypet_plasma[[8]], plot_list_amypet_plasma[[9]], plot_list_amypet_plasma[[10]], plot_list_amypet_plasma[[11]],  nrow = 1)
row4 <- plot_grid(plot_list_amypet_plasma[[12]],plot_list_amypet_plasma[[13]], plot_list_amypet_plasma[[14]], plot_list_amypet_plasma[[15]], nrow = 1)
row5 <- plot_grid(plot_list_amypet_amyPET[[1]], plot_list_amypet_amyPET[[3]], plot_list_amypet_amyPET[[4]], blank_plot, nrow = 1)
row6 <- plot_grid(plot_list_amypet_amyPET[[5]],plot_list_amypet_amyPET[[6]], blank_plot, blank_plot, nrow = 1)

allplots_suvr <- plot_grid(row1, row2, row3, row4, row5, row6, nrow = 6)

get_legend <- function(my_plot) {
  g <- ggplotGrob(my_plot)
  legend <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")][[1]]
  legend
}
legend <- get_legend(plot_list_amypet_plasma[[1]])
legend_col <- plot_grid(legend, blank_plot, ncol = 1, rel_heights = c(0.20, 0.80))
final_plot <- plot_grid(allplots_suvr, legend_col, ncol = 2, rel_widths = c(0.8, 0.2))

ggsave("~/Documents/FNIH_paper_Amyclockplasma/Annals_Neurology_submission/Review/Supp_fig9.pdf", plot = allplots_suvr ,  width = 17, height = 20)


#######plot TAU PET SUVR ##########
create_spaghetti_plot5 <- function(final_dataset_plasma,MesialTemporal,dep_var, CDR_3) {
  
  # Create segments from filtered data
  segments <- final_dataset_plasma%>%
    group_by(RID) %>%
    arrange(!!sym(MesialTemporal)) %>%
    mutate(
      MesialTemporal_next = lead(!!sym(MesialTemporal)),
      plasma_bmks_next = lead(!!sym(dep_var)),
      group_next = lead(CDR_3)
    ) %>%
    filter(!is.na(MesialTemporal), !is.na(group_next)) %>%
    ungroup()
  
  # Create two segments for each transition
  segments <- segments %>%
    mutate(
      MesialTemporal_mid = (!!sym(MesialTemporal) + MesialTemporal_next) / 2,
      plasma_bmks_mid = (!!sym(dep_var) +  plasma_bmks_next) / 2
    )
  
  # Create first half segments
  first_half <- segments %>%
    transmute(
      MesialTemporal = !!sym(MesialTemporal),
      plasma_bmks = !!sym(dep_var),
      MesialTemporal_end = MesialTemporal_mid,
      plasma_bmks_end = plasma_bmks_mid,
      CDR_3 = CDR_3
    )
  
  # Create second half segments
  second_half <- segments %>%
    transmute(
      MesialTemporal = MesialTemporal_mid,
      plasma_bmks = plasma_bmks_mid,
      MesialTemporal_end = MesialTemporal_next,
      plasma_bmks_end = plasma_bmks_next,
      CDR_3 = group_next
    )
  
  # Combine both halves
  all_segments <- bind_rows(first_half, second_half)
  
  # Plot using geom_segment and geom_point
  plot <- ggplot() +
    geom_segment(data = all_segments, aes(
      x = MesialTemporal, y = plasma_bmks,
      xend = MesialTemporal_end, yend = plasma_bmks_end,
      color = CDR_3
    ), size = 0.3, show.legend =F ) +
    geom_point(data = final_dataset_plasma, aes(x = !!sym(MesialTemporal), y = !!sym(dep_var), color = CDR_3, shape = APOE_binary), size = 1, show.legend = F) +
    
    theme_classic() + 
    labs(x = expression(""^18*"F-flortaucipir mesial-temporal tau PET (SUVR)"),
         y = biomarker_labels[[dep_var]],
         color = "Group", shape="APOE_binary") + theme( axis.title =  element_text(size = 12), legend.text = element_text(size=12, hjust = 0)) +
    scale_shape_manual(values = c(16, 17), labels = c(expression(italic("APOE-")~epsilon*"4 non-carrier"), expression(italic("APOE-")~epsilon*"4 carrier"))) + xlim(0.8,4) +
    scale_color_manual(values = c("blue2","darkorange", "brown3" ), labels = c("CDR = 0","CDR = 0.5", "CDR >= 1"))+
    guides(color = guide_legend(title = NULL, override.aes = list(size = 2), order = 1), shape = guide_legend(title = NULL,override.aes = list(size = 2), order = 2))
  
  return(plot)
}


# Apply the function to each dependent variable
dep_var <- c( "C2N_plasma_ptau217_out", "Fuji_plasma_ptau217_out",  "AlzPath_plasma_ptau217_out", "Janssen_plasma_ptau217_out", 
              "C2N_plasma_ptau217_ratio", "QX_plasma_ptau181_out","Roche_plasma_ptau181",
              "C2N_plasma_Abeta42_Abeta40_out" ,"Fuji_plasma_Ab42_Ab40_out","Roche_plasma_Ab42_Ab40_out","QX_plasma_Ab42_Ab40_out", 
              "Roche_plasma_GFAP_out", "QX_plasma_GFAP_out","Roche_plasma_NfL_out","QX_plasma_NfL_out")

plot_list_tausuvr <- lapply(dep_var, function(dep_var) {
  create_spaghetti_plot5(final_dataset_plasma, "MesialTemporal", dep_var, "color = CDR_3")
})
###nonplasma biomarkers
dep_var <- c("PTAU_over_ABETA42", "SUVR_compositeRef",  "TemporoParietal", "metaROI.AgeAdj", "CDRSB")
plot_list_taupet_ref <- lapply(dep_var, function(dep_var) {
  create_spaghetti_plot5(tau_dataset, "MesialTemporal", dep_var, "color = CDR_3")
})


blank_plot <- ggplot() + theme_void()
row1 <- plot_grid(plot_list_tausuvr[[1]], plot_list_tausuvr[[2]], plot_list_tausuvr[[3]], plot_list_tausuvr[[4]], nrow = 1)
row2 <- plot_grid(plot_list_tausuvr[[5]], plot_list_tausuvr[[6]], plot_list_tausuvr[[7]], blank_plot, nrow = 1)
row3 <- plot_grid(plot_list_tausuvr[[8]], plot_list_tausuvr[[9]], plot_list_tausuvr[[10]], plot_list_tausuvr[[11]],  nrow = 1)
row4 <- plot_grid(plot_list_tausuvr[[12]],plot_list_tausuvr[[13]], plot_list_tausuvr[[14]], plot_list_tausuvr[[15]], nrow = 1)
row5 <- plot_grid(plot_list_taupet_ref[[1]],plot_list_taupet_ref[[2]], plot_list_taupet_ref[[3]], blank_plot, nrow = 1)
row6 <- plot_grid(plot_list_taupet_ref[[4]],plot_list_taupet_ref[[5]], blank_plot, blank_plot, nrow = 1)

allplots_tausuvr <- plot_grid(row1, row2, row3, row4, row5, row6, nrow = 6)

get_legend <- function(my_plot) {
  g <- ggplotGrob(my_plot)
  legend <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")][[1]]
  legend
}
legend <- get_legend(plot_list_tausuvr[[1]])
legend_col <- plot_grid(legend, blank_plot, ncol = 1, rel_heights = c(0.20, 0.80))
final_plot <- plot_grid(allplots_tausuvr, legend_col, ncol = 2, rel_widths = c(0.8, 0.2))

ggsave("~/Documents/FNIH_paper_Amyclockplasma/Annals_Neurology_submission/Review/Supp_fig10.pdf", plot = allplots_tausuvr  ,  width = 17, height = 20)

#######plot EYO TAU ##########

create_spaghetti_plot6 <- function(final_dataset_plasma,EYO1_tau ,dep_var, CDR_3) {
  
  # Create segments from filtered data
  segments <- final_dataset_plasma%>%
    group_by(RID) %>%
    arrange(!!sym(EYO1_tau)) %>%
    mutate(
      EYO1_tau_next = lead(!!sym(EYO1_tau)),
      plasma_bmks_next = lead(!!sym(dep_var)),
      group_next = lead(CDR_3)
    ) %>%
    filter(!is.na(EYO1_tau), !is.na(group_next)) %>%
    ungroup()
  
  # Create two segments for each transition
  segments <- segments %>%
    mutate(
      EYO1_tau_mid = (!!sym(EYO1_tau) + EYO1_tau_next) / 2,
      plasma_bmks_mid = (!!sym(dep_var) +  plasma_bmks_next) / 2
    )
  
  # Create first half segments
  first_half <- segments %>%
    transmute(
      EYO1_tau = !!sym(EYO1_tau),
      plasma_bmks = !!sym(dep_var),
      EYO1_tau_end = EYO1_tau_mid,
      plasma_bmks_end = plasma_bmks_mid,
      CDR_3 = CDR_3
    )
  
  # Create second half segments
  second_half <- segments %>%
    transmute(
      EYO1_tau = EYO1_tau_mid,
      plasma_bmks = plasma_bmks_mid,
      EYO1_tau_end = EYO1_tau_next,
      plasma_bmks_end = plasma_bmks_next,
      CDR_3 = group_next
    )
  
  # Combine both halves
  all_segments <- bind_rows(first_half, second_half)
  
  # Plot using geom_segment and geom_point
  plot <- ggplot() +
    geom_segment(data = all_segments, aes(
      x =EYO1_tau, y = plasma_bmks,
      xend = EYO1_tau_end, yend = plasma_bmks_end,
      color = CDR_3
    ), size = 0.3, show.legend =F) +
    geom_point(data = final_dataset_plasma, aes(x = !!sym(EYO1_tau), y = !!sym(dep_var), color = CDR_3, shape = APOE_binary), size = 1, show.legend = F) +
    
    theme_minimal() +
    labs(x = "EYO (years)",
         y = biomarker_labels[[dep_var]],
         color = "Group", shape="APOE_binary") + theme( axis.title =  element_text(size = 10), legend.text = element_text(size=10, hjust = 0)) +
    scale_shape_manual(values = c(16, 17), labels = c(expression(italic("APOE-")~epsilon*"4 non-carrier"), expression(italic("APOE-")~epsilon*"4 carrier"))) +
    scale_color_manual(values = c("blue2","brown3", "forestgreen" ), labels = c("CDR = 0","CDR = 0.5", "CDR >= 1"))+
    guides(color = guide_legend(title = NULL, override.aes = list(size = 2), order = 1), shape = guide_legend(title = NULL,override.aes = list(size = 2), order = 2))
  
  return(plot)
}


# Apply the function to each dependent variable
plot_list_eyo <- lapply(dep_vars, function(dep_var) {
  create_spaghetti_plot6(final_dataset_plasma, "EYO1_tau", dep_var, "color = CDR_3")
})


blank_plot <- ggplot() + theme_void()
row1 <- plot_grid(plot_list_eyo[[1]], plot_list_eyo[[2]], plot_list_eyo[[3]], plot_list_eyo[[4]], nrow = 1)
row2 <- plot_grid(plot_list_eyo[[5]], plot_list_eyo[[6]], plot_list_eyo[[7]], blank_plot, nrow = 1)
row3 <- plot_grid(plot_list_eyo[[8]], plot_list_eyo[[9]], plot_list_eyo[[10]], plot_list_eyo[[11]],  nrow = 1)
row4 <- plot_grid(plot_list_eyo[[12]],plot_list_eyo[[13]], plot_list_eyo[[14]], plot_list_eyo[[15]], nrow = 1)
row5 <- plot_grid(plot_list_eyo[[16]],plot_list_eyo[[17]], plot_list_eyo[[18]], plot_list_eyo[[19]], nrow = 1)
row6 <- plot_grid(plot_list_eyo[[20]],plot_list_eyo[[21]], plot_list_eyo[[22]], blank_plot, nrow = 1)

allplots_eyoTAU <- plot_grid(row1, row2, row3, row4, row5, row6, nrow = 6)

get_legend <- function(my_plot) {
  g <- ggplotGrob(my_plot)
  legend <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")][[1]]
  legend
}
legend <- get_legend(plot_list_eyo[[1]])
legend_col <- plot_grid(legend, blank_plot, ncol = 1, rel_heights = c(0.20, 0.80))
final_plot <- plot_grid(allplots_eyoTAU, legend_col, ncol = 2, rel_widths = c(0.8, 0.2))

ggsave("~/Documents/FNIH_paper_Amyclockplasma/Figures paper/Suppl figures/trajectories_EYO_TAU.pdf", plot = final_plot ,  width = 15, height = 15)




######## main figure aligment ###

row1 <- plot_grid(plot_list_age[[20]], plot_list_amy[[20]], plot_list_tau[[20]], plot_list_eyo[[20]], nrow = 1)
row2 <- plot_grid(plot_list_age[[21]], plot_list_amy[[21]], plot_list_tau[[21]], plot_list_eyo[[21]],  nrow = 1)
plot_main <- plot_grid(row1, row2, nrow = 2)
legend_col <- plot_grid(legend, blank_plot, ncol = 1, rel_heights = c(0.3, 0.7))
final_plot <- plot_grid(plot_main, legend_col, ncol = 2, rel_widths = c(0.8, 0.2))
ggsave("~/Documents/FNIH_paper_Amyclockplasma/Figures paper/Fig1_trajectories_main_n.pdf", plot = final_plot ,  width = 18, height = 8)


#############GAM MODELS PLOTS################

##define CDR factor
#amy dataset#
final_dataset_plasma$CDR <- factor(final_dataset_plasma$CDR, levels=c(0,0.5,1, 2,3))
final_dataset_plasma<- final_dataset_plasma %>%
  mutate(CDR_3 = case_when(
    CDR == 0 ~ "0",
    CDR == 0.5 ~ "0.5",
    CDR %in% c("1", "2", "3") ~ ">=1"
  ))
final_dataset_plasma$CDR_3 <- factor(final_dataset_plasma$CDR_3, levels = c("0", "0.5", ">=1"), labels = c("CDR = 0", "CDR = 0.5", "CDR >= 1"))

plasma_biomarkers <- c( "C2N_plasma_ptau217_out", "Fuji_plasma_ptau217_out",  "AlzPath_plasma_ptau217_out", "Janssen_plasma_ptau217_out", 
                        "C2N_plasma_ptau217_ratio", "QX_plasma_ptau181_out","Roche_plasma_ptau181",
                        "C2N_plasma_Abeta42_Abeta40_out" ,"Fuji_plasma_Ab42_Ab40_out","Roche_plasma_Ab42_Ab40_out","QX_plasma_Ab42_Ab40_out", 
                        "Roche_plasma_GFAP_out", "QX_plasma_GFAP_out","Roche_plasma_NfL_out","QX_plasma_NfL_out", "atrophy", "PTAU_over_ABETA42", "CDR_SOB", "MMSCORE", "SUVR_compositeRef", "MesialTemporal", "TemporoParietal")

colors <- c("forestgreen", "forestgreen","forestgreen", "forestgreen", "darkgreen", "green3", "green3","blue", "blue", "blue", "blue","red", "red",  "brown", "brown", "orange2","darkgrey", "purple", "purple", "coral", "darkblue", "darkblue")





segment_ranges <- list(
  Roche_plasma_Ab42_Ab40_out = c(-7.28, 14.34),
  C2N_plasma_Abeta42_Abeta40_out = c(-7.28, 13.72), 
  Fuji_plasma_Ab42_Ab40_out = c(-7.28, 13.52), 
    QX_plasma_Ab42_Ab40_out = c(-7.28, 12.91 ),
  C2N_plasma_ptau217_out = c(-7.28, 33.3 ),
  Fuji_plasma_ptau217_out = c(-7.28,33.3),
  
  Roche_plasma_ptau181 = c(-7.28, 33.3),
  QX_plasma_ptau181_out = c(-7.28,33.3),
  
  AlzPath_plasma_ptau217_out = c(-7.28,33.3 ),
  Janssen_plasma_ptau217_out = c(-7.28, 33.3 ),
  C2N_plasma_ptau217_ratio = c(-7.28,33.3),
  Roche_plasma_GFAP_out = c(-7.28, 20.25 ),
  QX_plasma_GFAP_out = c(-7.28 ,15.97),
  Roche_plasma_NfL_out = c(-7.28,13.72),
  QX_plasma_NfL_out = c(-7.28, 13.32 ),
  PTAU_over_ABETA42 = c(-7.28,30.7),
  atrophy= c(-7.28, 31.2 ),
  CDR_SOB= c(-7.28, 33.3),
  MMSCORE= c(-7.28, 33.3),
  MesialTemporal = c(-2.65,19.67 ),
  TemporoParietal = c(4.84,30.5),
  SUVR_compositeRef = c(-7.28,33.3)
)

tipping_points <-c(
  Roche_plasma_Ab42_Ab40_out = -7.28,
  C2N_plasma_Abeta42_Abeta40_out = -7.28, 
  Fuji_plasma_Ab42_Ab40_out = -4.52, 
  QX_plasma_Ab42_Ab40_out = -3.98,
  C2N_plasma_ptau217_out = -2.17 ,
  Fuji_plasma_ptau217_out = -1.24  ,
  
  Roche_plasma_ptau181 = -2.63 ,
  QX_plasma_ptau181_out = -2.08,
  
  AlzPath_plasma_ptau217_out = -3.09,
  Janssen_plasma_ptau217_out = -2.66,
  C2N_plasma_ptau217_ratio = -3.9,
  Roche_plasma_GFAP_out = -2.6,
  QX_plasma_GFAP_out = -3.45,
  Roche_plasma_NfL_out = 30.8,
  QX_plasma_NfL_out = 30.8,
  PTAU_over_ABETA42 = -3.54,
  atrophy= 12.69,
  CDR_SOB= 0.61,
  MMSCORE= 0.81,
  MesialTemporal = 2.53,
  TemporoParietal = 5.67,
  SUVR_compositeRef = -4.73
)
plot_list_amyloid <- list()

# Loop through each y variable and create a scatter plot
for (i in seq_along(plasma_biomarkers)) {
  y_var <- plasma_biomarkers[i]
  color <- colors[i]
  
  ref_mean <- subset_IDs_with_positive1_scans[[paste0(y_var, "_mean")]]
  
 if (y_var %in% names(segment_ranges)) {
   x_range <- segment_ranges[[y_var]]
    tipping_point <- tipping_points[[y_var]]
  # Create the scatter plot
  plot <- ggplot(final_dataset_plasma, aes(x = Amyloid_time, y = !!sym(y_var))) +
    geom_point(size=0.4, color=color, show.legend = F) + 
    geom_smooth( method="gam", color=color, show.legend = F, formula = y ~ s(x, k = 3), size=0.25) +
      geom_hline(yintercept = ref_mean, linetype = "dashed", color = "grey0", size= 0.65)+
    geom_vline(xintercept = tipping_point, linetype = "solid", color = "grey0", size= 0.65)+
    geom_ribbon(aes(ymin = !!sym(paste0(y_var, "_ci_lower")), ymax = !!sym(paste0(y_var, "_ci_upper"))), alpha = 0.2, fill = "black") +
    
    labs(x = "Amyloid time (years)", y = biomarker_labels[[y_var]]) +
    theme_minimal() + 
    theme(
      axis.title.x = element_text(size = 10),  # Change x-axis label size
      axis.title.y = element_text(size = 10))
  

  if (!is.null(x_range)){
    full_smooth <- ggplot_build(plot)$data[[2]]
    segment_smooth <- full_smooth %>% filter(x >= x_range[1] & x <= x_range[2])
   
  # Add the thicker line segment
 plot <- plot +
  geom_line(data = segment_smooth, aes(x = x, y = y), color = color, size = 2.5)
  }  
  
# Add the plot to the list
plot_list_amyloid[[y_var]] <- plot
}
 }
  
  
blank_plot <- ggplot() + theme_void()
row1 <- plot_grid( plot_list_amyloid[[1]],  plot_list_amyloid[[2]],  plot_list_amyloid[[3]], plot_list_amyloid[[4]], nrow = 1)
row2 <- plot_grid( plot_list_amyloid[[5]],  plot_list_amyloid[[6]],  plot_list_amyloid[[7]], blank_plot, nrow = 1)
row3 <- plot_grid( plot_list_amyloid[[8]],  plot_list_amyloid[[9]],  plot_list_amyloid[[10]],  plot_list_amyloid[[11]],  nrow = 1)
row4 <- plot_grid( plot_list_amyloid[[12]], plot_list_amyloid[[13]],  plot_list_amyloid[[14]],  plot_list_amyloid[[15]], nrow = 1)
row5 <- plot_grid( plot_list_amyloid[[16]], plot_list_amyloid[[17]],  plot_list_amyloid[[18]],  plot_list_amyloid[[19]], nrow=1)
row6 <- plot_grid( plot_list_amyloid[[20]], plot_list_amyloid[[21]],  plot_list_amyloid[[22]],blank_plot,nrow=1)

allplots_AT <- plot_grid(row1, row2, row3, row4, row5, row6,  nrow = 6)

##gam amyloid modified


plasma_biomarkers <- c( "C2N_plasma_ptau217_out", "Fuji_plasma_ptau217_out",  "AlzPath_plasma_ptau217_out", "Janssen_plasma_ptau217_out", 
                        "C2N_plasma_ptau217_ratio", "QX_plasma_ptau181_out","Roche_plasma_ptau181",
                        "C2N_plasma_Abeta42_Abeta40_out" ,"Fuji_plasma_Ab42_Ab40_out","Roche_plasma_Ab42_Ab40_out","QX_plasma_Ab42_Ab40_out", 
                        "Roche_plasma_GFAP_out", "QX_plasma_GFAP_out","Roche_plasma_NfL_out","QX_plasma_NfL_out")


biomarker_labels <- c(
  C2N_plasma_Abeta42_Abeta40_out = expression(paste("C2N A",beta,"42/A",beta,"40")),
  Fuji_plasma_Ab42_Ab40_out = expression(paste("Fujirebio A",beta,"42/A",beta,"40")),
  Roche_plasma_Ab42_Ab40_out = expression(paste("Roche A",beta, "42/A", beta,"40")),
  QX_plasma_Ab42_Ab40_out = expression(paste("Quanterix A",beta, "42/A",beta,"40")),
  Roche_plasma_ptau181 = "Roche p-tau181 (pg/ml)",
  QX_plasma_ptau181_out = "Quanterix p-tau181 (pg/ml)",
  C2N_plasma_ptau217_out = "C2N p-tau217 (pg/ml)",
  Fuji_plasma_ptau217_out = "Fujirebio p-tau217 (pg/ml)",
  AlzPath_plasma_ptau217_out = "ALZPath p-tau217 (pg/ml)",
  Janssen_plasma_ptau217_out = "Janssen p-tau217 (pg/ml)",
  C2N_plasma_ptau217_ratio = "C2N %p-tau217 (%)",
  Roche_plasma_GFAP_out = "Roche GFAP (ng/ml)",
  QX_plasma_GFAP_out = "Quanterix GFAP (pg/ml)",
  Roche_plasma_NfL_out = "Roche NfL (pg/ml)",
  QX_plasma_NfL_out = "Quanterix NfL (pg/ml)")


segment_ranges <- list(
  Roche_plasma_Ab42_Ab40_out = c(-7.94, 13),
  C2N_plasma_Abeta42_Abeta40_out = c(-7.94, 12.2), 
  Fuji_plasma_Ab42_Ab40_out = c(-7.94, 12.4), 
  QX_plasma_Ab42_Ab40_out = c(-7.09, 11.5 ),
  C2N_plasma_ptau217_out = c(-6.84, 30.3 ),
  Fuji_plasma_ptau217_out = c(-7.94,30.3),
  
  Roche_plasma_ptau181 = c(-7.94,30.3),
  QX_plasma_ptau181_out = c(-7.09,29.3),
  
  AlzPath_plasma_ptau217_out = c(-7.94,30.3),
  Janssen_plasma_ptau217_out = c(-7.94, 30.3),
  C2N_plasma_ptau217_ratio = c(-7.94, 30.3),
  Roche_plasma_GFAP_out = c(-7.94, 30.3),
  QX_plasma_GFAP_out = c(-7.09 ,29.3),
  Roche_plasma_NfL_out = c(-7.94,30.3),
  QX_plasma_NfL_out = c(-7.09, 29.3 )
)
  

tipping_points <-c(
  Roche_plasma_Ab42_Ab40_out = -7.94,
  C2N_plasma_Abeta42_Abeta40_out = -7.94, 
  Fuji_plasma_Ab42_Ab40_out = -5.81, 
  QX_plasma_Ab42_Ab40_out = -1.56,
  C2N_plasma_ptau217_out = -0.35 ,
  Fuji_plasma_ptau217_out = -1.09  ,
  
  Roche_plasma_ptau181 = -2.71 ,
  QX_plasma_ptau181_out = -0.77,
  
  AlzPath_plasma_ptau217_out = -2.63,
  Janssen_plasma_ptau217_out = -3.00,
  C2N_plasma_ptau217_ratio = -4.38,
  Roche_plasma_GFAP_out = 0.81,
  QX_plasma_GFAP_out = 2.37,
  Roche_plasma_NfL_out = 9.21,
  QX_plasma_NfL_out = 10.5
)
plot_list_amyloid <- list()

create_segments <- function(final_dataset_plasma,years_amy_onset,dep_var, CDR_3) {
  
  # Create segments from filtered data
  segments <- final_dataset_plasma%>%
    group_by(RID) %>%
    arrange(!!sym( years_amy_onset)) %>%
    mutate(
      years_amy_onset_next = lead(!!sym( years_amy_onset)),
      plasma_bmks_next = lead(!!sym(dep_var)),
      group_next = lead(CDR_3)
    ) %>%
    filter(!is.na(years_amy_onset_next), !is.na(group_next)) %>%
    ungroup()
  
  # Create two segments for each transition
  segments <- segments %>%
    mutate(
      years_amy_onset_mid = (!!sym(years_amy_onset) + years_amy_onset_next) / 2,
      plasma_bmks_mid = (!!sym(dep_var) +  plasma_bmks_next) / 2
    )
  
  # Create first half segments
  first_half <- segments %>%
    transmute(
      years_amy_onset = !!sym(years_amy_onset),
      plasma_bmks = !!sym(dep_var),
      years_amy_onset_end = years_amy_onset_mid,
      plasma_bmks_end = plasma_bmks_mid,
      CDR_3 = CDR_3
    )
  
  # Create second half segments
  second_half <- segments %>%
    transmute(
      years_amy_onset = years_amy_onset_mid,
      plasma_bmks = plasma_bmks_mid,
      years_amy_onset_end = years_amy_onset_next,
      plasma_bmks_end = plasma_bmks_next,
      CDR_3 = group_next
    )
  
  # Combine both halves
  all_segments <- bind_rows(first_half, second_half)
  return(all_segments)
}

# Loop through each y variable and create a scatter plot
for (i in seq_along(plasma_biomarkers)) {
  y_var <- plasma_biomarkers[i]
  
  
  ref_mean <- final_dataset_plasma[[paste0(y_var, "_mean")]]
  
  if (y_var %in% names(segment_ranges)) {
    x_range <- segment_ranges[[y_var]]
    tipping_point <- tipping_points[[y_var]]
# Create the segments for geom_segment
    all_segments <- create_segments(final_dataset_plasma, "years_amy_onset", y_var, "CDR_3")
    
    # Create the scatter plot
    plot <- ggplot(final_dataset_plasma, aes(x = years_amy_onset, y = !!sym(y_var))) +
      geom_point(size=1, show.legend = F,  aes(color=CDR_3, shape = APOE_binary)) + 
      geom_segment(data = all_segments, aes(x = years_amy_onset, y = plasma_bmks, xend = years_amy_onset_end, yend = plasma_bmks_end,color = CDR_3 ), size = 0.3, show.legend =F) +
      geom_smooth( method="gam", color="grey1", show.legend = F, formula = y ~ s(x, k = 3), size=0.25) +
      geom_hline(yintercept = ref_mean, linetype = "dashed", color = "grey0", size= 0.65)+
      geom_vline(xintercept = tipping_point, linetype = "solid", color = "darkgrey", size= 0.75)+
      geom_ribbon(aes(ymin = !!sym(paste0(y_var, "_ci_lower")), ymax = !!sym(paste0(y_var, "_ci_upper"))), alpha = 0.2, fill = "black") +
      scale_color_manual(values = c("blue2","orange2", "brown3") ) +
      
      labs(x = "Estimated years from amyloid PET positivity", y = biomarker_labels[[y_var]]) +
      theme_classic() + 
      theme(
        axis.title.x = element_text(size = 14),  # Change x-axis label size
        axis.title.y = element_text(size = 14))
    
    
    if (!is.null(x_range)){
      full_smooth <- ggplot_build(plot)$data[[3]]
      segment_smooth <- full_smooth %>% filter(x >= x_range[1] & x <= x_range[2])
      
      # Add the thicker line segment
      plot <- plot +
        geom_line(data = segment_smooth, aes(x = x, y = y), color = "grey1", size = 1.25)
    }  
    
    # Add the plot to the list
    plot_list_amyloid[[y_var]] <- plot
  }
}


blank_plot <- ggplot() + theme_void()
row1 <- plot_grid( plot_list_amyloid[[1]],  plot_list_amyloid[[2]],  plot_list_amyloid[[3]], plot_list_amyloid[[4]], nrow = 1)
row2 <- plot_grid( plot_list_amyloid[[5]],  plot_list_amyloid[[6]],  plot_list_amyloid[[7]], blank_plot, nrow = 1)
row3 <- plot_grid( plot_list_amyloid[[8]],  plot_list_amyloid[[9]],  plot_list_amyloid[[10]],  plot_list_amyloid[[11]],  nrow = 1)
row4 <- plot_grid( plot_list_amyloid[[12]], plot_list_amyloid[[13]],  plot_list_amyloid[[14]],  plot_list_amyloid[[15]], nrow = 1)

allplots_AT <- plot_grid(row1, row2, row3, row4,  nrow = 4)

##GAM for tau##
segment_ranges_tau <- list(
  Roche_plasma_Ab42_Ab40_out = c(-11.4, 17.9),
  C2N_plasma_Abeta42_Abeta40_out = c( -11.4, 17.9), 
  Fuji_plasma_Ab42_Ab40_out = c(-0.77, 4.37 ), 
  QX_plasma_Ab42_Ab40_out = c(-11.4, 5.99 ),
  C2N_plasma_ptau217_out = c( -0.99, 17.9 ),
  Fuji_plasma_ptau217_out = c(-11.4, 17.9),
  
  Roche_plasma_ptau181 = c(-11.4, 17.9 ),
  QX_plasma_ptau181_out = c(-11.4,3.64 ),
  
  AlzPath_plasma_ptau217_out = c(-11.4,17.9 ),
  Janssen_plasma_ptau217_out = c(-11.4,17.9),
  C2N_plasma_ptau217_ratio = c(-11.4,17.9),
  Roche_plasma_GFAP_out = c(-11.4,17.9 ),
  QX_plasma_GFAP_out = NULL,
  Roche_plasma_NfL_out = NULL,
  QX_plasma_NfL_out = NULL)

tipping_points_tau <-c(
  Roche_plasma_Ab42_Ab40_out = -7.1,
  C2N_plasma_Abeta42_Abeta40_out = -9.3, 
  Fuji_plasma_Ab42_Ab40_out = -7.5 , 
  QX_plasma_Ab42_Ab40_out = -3.26,
  C2N_plasma_ptau217_out = -4.72 ,
  Fuji_plasma_ptau217_out = -3.88,
  
  Roche_plasma_ptau181 = -6.7,
  QX_plasma_ptau181_out = -5.4 ,
  
  AlzPath_plasma_ptau217_out = -5.4 ,
  Janssen_plasma_ptau217_out = -5.2 ,
  C2N_plasma_ptau217_ratio = -4.86,
  Roche_plasma_GFAP_out = -11.4 ,
  QX_plasma_GFAP_out = -11.4 ,
  Roche_plasma_NfL_out = 5.28 ,
  QX_plasma_NfL_out = 15.3
)
plot_list_tau <- list()



create_segments <- function(final_dataset_plasma,years_tau_onset,dep_var, CDR_3) {
  
  # Create segments from filtered data
  segments <- final_dataset_plasma%>%
    group_by(RID) %>%
    arrange(!!sym( years_tau_onset)) %>%
    mutate(
      years_tau_onset_next = lead(!!sym( years_tau_onset)),
      plasma_bmks_next = lead(!!sym(dep_var)),
      group_next = lead(CDR_3)
    ) %>%
    filter(!is.na(years_tau_onset_next), !is.na(group_next)) %>%
    ungroup()
  
  # Create two segments for each transition
  segments <- segments %>%
    mutate(
      years_tau_onset_mid = (!!sym(years_tau_onset) + years_tau_onset_next) / 2,
      plasma_bmks_mid = (!!sym(dep_var) +  plasma_bmks_next) / 2
    )
  
  # Create first half segments
  first_half <- segments %>%
    transmute(
      years_tau_onset = !!sym(years_tau_onset),
      plasma_bmks = !!sym(dep_var),
      years_tau_onset_end = years_tau_onset_mid,
      plasma_bmks_end = plasma_bmks_mid,
      CDR_3 = CDR_3
    )
  
  # Create second half segments
  second_half <- segments %>%
    transmute(
      years_tau_onset = years_tau_onset_mid,
      plasma_bmks = plasma_bmks_mid,
      years_tau_onset_end = years_tau_onset_next,
      plasma_bmks_end = plasma_bmks_next,
      CDR_3 = group_next
    )
  
  # Combine both halves
  all_segments <- bind_rows(first_half, second_half)
  return(all_segments)
}

# Loop through each y variable and create a scatter plot
for (i in seq_along(plasma_biomarkers)) {
  y_var <- plasma_biomarkers[i]

  
  ref_mean <- final_dataset_plasma[[paste0(y_var, "_mean")]]
  
  if (y_var %in% names(segment_ranges_tau)) {
    x_range <- segment_ranges_tau[[y_var]]
    tipping_point <- tipping_points_tau [[y_var]]
    # Create the segments for geom_segment
    all_segments <- create_segments(final_dataset_plasma, "years_tau_onset", y_var, "CDR_3")
    
    # Create the scatter plot
    plot <- ggplot(final_dataset_plasma, aes(x = years_tau_onset, y = !!sym(y_var))) +
      geom_point(size=1, show.legend = F,  aes(color=CDR_3, shape = APOE_binary)) + 
      geom_segment(data = all_segments, aes(x = years_tau_onset, y = plasma_bmks, xend = years_tau_onset_end, yend = plasma_bmks_end,color = CDR_3 ), size = 0.3, show.legend =F) +
      geom_smooth( method="gam", color="grey1", show.legend = F, formula = y ~ s(x, k = 3), size=0.25) +
      geom_hline(yintercept = ref_mean, linetype = "dashed", color = "grey0", size= 0.65)+
      geom_vline(xintercept = tipping_point, linetype = "solid", color = "darkgrey", size= 0.75)+
      geom_ribbon(aes(ymin = !!sym(paste0(y_var, "_ci_lower")), ymax = !!sym(paste0(y_var, "_ci_upper"))), alpha = 0.2, fill = "black") +
      scale_color_manual(values = c("blue2","orange2", "brown3") ) +
      
      labs(x = "Estimated years from tau PET positivity", y = biomarker_labels[[y_var]]) +
      theme_classic() + 
      theme(
        axis.title.x = element_text(size = 14),  # Change x-axis label size
        axis.title.y = element_text(size = 14))
    
    
    if (!is.null(x_range)){
      full_smooth <- ggplot_build(plot)$data[[3]]
      segment_smooth <- full_smooth %>% filter(x >= x_range[1] & x <= x_range[2])
      
      # Add the thicker line segment
      plot <- plot +
        geom_line(data = segment_smooth, aes(x = x, y = y), color = "grey1", size = 1.25)
    }  
    
    # Add the plot to the list
    plot_list_tau[[y_var]] <- plot
  }
}



blank_plot <- ggplot() + theme_void()
row1 <- plot_grid(plot_list_tau[[1]], plot_list_tau[[2]], plot_list_tau[[3]],plot_list_tau[[4]], nrow = 1)
row2 <- plot_grid(plot_list_tau[[5]], plot_list_tau[[6]], plot_list_tau[[7]], blank_plot, nrow = 1)
row3 <- plot_grid(plot_list_tau[[8]], plot_list_tau[[9]], plot_list_tau[[10]], plot_list_tau[[11]],  nrow = 1)
row4 <- plot_grid(plot_list_tau[[12]],plot_list_tau[[13]], plot_list_tau[[14]], plot_list_tau[[15]], nrow = 1)

allplots_TT <- plot_grid(row1, row2, row3, row4,  nrow = 4)

##GAM for tau temporoparietal##
segment_ranges_tau <- list(
  Roche_plasma_Ab42_Ab40_out =  NULL,
  C2N_plasma_Abeta42_Abeta40_out =  NULL, 
  Fuji_plasma_Ab42_Ab40_out =  NULL, 
  QX_plasma_Ab42_Ab40_out =  NULL,
  C2N_plasma_ptau217_out = c( -10.2, 14.3 ),
  Fuji_plasma_ptau217_out = c(-10.2, 14.3),
  
  Roche_plasma_ptau181 = c(-10.2, 14.3 ),
  QX_plasma_ptau181_out = c(-10.2, 14.3 ),
  
  AlzPath_plasma_ptau217_out = c(-10.2, 14.3 ),
  Janssen_plasma_ptau217_out = c(-10.2, 14.3),
  C2N_plasma_ptau217_ratio = c(-10.2, 14.3),
  Roche_plasma_GFAP_out = NULL,
  QX_plasma_GFAP_out = c(1.7 , 14.3),
  Roche_plasma_NfL_out = c(-10.2, 14.3),
  QX_plasma_NfL_out = c(-10.2, 14.3)
)

tipping_points_tau <-c(
  Roche_plasma_Ab42_Ab40_out = -5.8 ,
  C2N_plasma_Abeta42_Abeta40_out = -5.8, 
  Fuji_plasma_Ab42_Ab40_out = -5.8 , 
  QX_plasma_Ab42_Ab40_out = -4.9,
  C2N_plasma_ptau217_out = -3.8  ,
  Fuji_plasma_ptau217_out = -3.8,
  
  Roche_plasma_ptau181 = -5.8 ,
  QX_plasma_ptau181_out = -1.8 ,
  
  AlzPath_plasma_ptau217_out = -5.8  ,
  Janssen_plasma_ptau217_out = -5.8  ,
  C2N_plasma_ptau217_ratio = -5.8 ,
  Roche_plasma_GFAP_out = -5.8  ,
  QX_plasma_GFAP_out = -5.8  ,
  Roche_plasma_NfL_out = -1.0 ,
  QX_plasma_NfL_out = -0.1
  )

plot_list_tau <- list()



create_segments <- function(final_dataset_plasma,years_tau_TP_onset,dep_var, CDR_3) {
  
  # Create segments from filtered data
  segments <- final_dataset_plasma%>%
    group_by(RID) %>%
    arrange(!!sym( years_tau_TP_onset)) %>%
    mutate(
      years_tau_TP_onset_next = lead(!!sym( years_tau_TP_onset)),
      plasma_bmks_next = lead(!!sym(dep_var)),
      group_next = lead(CDR_3)
    ) %>%
    filter(!is.na(years_tau_TP_onset_next), !is.na(group_next)) %>%
    ungroup()
  
  # Create two segments for each transition
  segments <- segments %>%
    mutate(
      years_tau_TP_onset_mid = (!!sym(years_tau_TP_onset) + years_tau_TP_onset_next) / 2,
      plasma_bmks_mid = (!!sym(dep_var) +  plasma_bmks_next) / 2
    )
  
  # Create first half segments
  first_half <- segments %>%
    transmute(
      years_tau_TP_onset = !!sym(years_tau_TP_onset),
      plasma_bmks = !!sym(dep_var),
      years_tau_TP_onset_end = years_tau_TP_onset_mid,
      plasma_bmks_end = plasma_bmks_mid,
      CDR_3 = CDR_3
    )
  
  # Create second half segments
  second_half <- segments %>%
    transmute(
      years_tau_TP_onset = years_tau_TP_onset_mid,
      plasma_bmks = plasma_bmks_mid,
      years_tau_TP_onset_end = years_tau_TP_onset_next,
      plasma_bmks_end = plasma_bmks_next,
      CDR_3 = group_next
    )
  
  # Combine both halves
  all_segments <- bind_rows(first_half, second_half)
  return(all_segments)
}

# Loop through each y variable and create a scatter plot
for (i in seq_along(plasma_biomarkers)) {
  y_var <- plasma_biomarkers[i]
  
  
  ref_mean <- final_dataset_plasma[[paste0(y_var, "_mean")]]
  
  if (y_var %in% names(segment_ranges_tau)) {
    x_range <- segment_ranges_tau[[y_var]]
    tipping_point <- tipping_points_tau [[y_var]]
    # Create the segments for geom_segment
    all_segments <- create_segments(final_dataset_plasma, "years_tau_TP_onset", y_var, "CDR_3")
    
    # Create the scatter plot
    plot <- ggplot(final_dataset_plasma, aes(x = years_tau_TP_onset, y = !!sym(y_var))) +
      geom_point(size=1, show.legend = F,  aes(color=CDR_3, shape = APOE_binary)) + 
      geom_segment(data = all_segments, aes(x = years_tau_TP_onset, y = plasma_bmks, xend = years_tau_TP_onset_end, yend = plasma_bmks_end,color = CDR_3 ), size = 0.3, show.legend =F) +
      geom_smooth( method="gam", color="grey1", show.legend = F, formula = y ~ s(x, k = 3), size=0.25) +
      geom_hline(yintercept = ref_mean, linetype = "dashed", color = "grey0", size= 0.65)+
      geom_vline(xintercept = tipping_point, linetype = "solid", color = "darkgrey", size= 0.75)+
      geom_ribbon(aes(ymin = !!sym(paste0(y_var, "_ci_lower")), ymax = !!sym(paste0(y_var, "_ci_upper"))), alpha = 0.2, fill = "black") +
      scale_color_manual(values = c("blue2","orange2", "brown3") ) +
      
      labs(x = "Estimated years from tau PET positivity", y = biomarker_labels[[y_var]]) +
      theme_classic() + 
      theme(
        axis.title.x = element_text(size = 14),  # Change x-axis label size
        axis.title.y = element_text(size = 14))
    
    
    if (!is.null(x_range)){
      full_smooth <- ggplot_build(plot)$data[[3]]
      segment_smooth <- full_smooth %>% filter(x >= x_range[1] & x <= x_range[2])
      
      # Add the thicker line segment
      plot <- plot +
        geom_line(data = segment_smooth, aes(x = x, y = y), color = "grey1", size = 1.25)
    }  
    
    # Add the plot to the list
    plot_list_tau[[y_var]] <- plot
  }
}



blank_plot <- ggplot() + theme_void()
row1 <- plot_grid(plot_list_tau[[1]], plot_list_tau[[2]], plot_list_tau[[3]],plot_list_tau[[4]], nrow = 1)
row2 <- plot_grid(plot_list_tau[[5]], plot_list_tau[[6]], plot_list_tau[[7]], blank_plot, nrow = 1)
row3 <- plot_grid(plot_list_tau[[8]], plot_list_tau[[9]], plot_list_tau[[10]], plot_list_tau[[11]],  nrow = 1)
row4 <- plot_grid(plot_list_tau[[12]],plot_list_tau[[13]], plot_list_tau[[14]], plot_list_tau[[15]], nrow = 1)

allplots_TT <- plot_grid(row1, row2, row3, row4,  nrow = 4)

##GAM for EYO ##
segment_ranges_eyo <- list(
  Roche_plasma_Ab42_Ab40_out = NULL,
  C2N_plasma_Abeta42_Abeta40_out = c(-11.4, -2.2 ), 
  Fuji_plasma_Ab42_Ab40_out = c(-16.7, -3.92 ), 
  QX_plasma_Ab42_Ab40_out = c(-16.7, -0.79 ),
  C2N_plasma_ptau217_out = c(-15.6, 7.28 ),
  Fuji_plasma_ptau217_out = c(-16.7, 7.28),
  
  Roche_plasma_ptau181 = c(-16.7, 2.22 ),
  QX_plasma_ptau181_out = c(-16.7, 7.28 ),
  
  AlzPath_plasma_ptau217_out = c(-16.7, 7.28 ),
  Janssen_plasma_ptau217_out = c(-16.7, 7.28 ),
  C2N_plasma_ptau217_ratio = c(-16.7, 7.28 ),
  Roche_plasma_GFAP_out = c(-16.7, 7.28 ),
  QX_plasma_GFAP_out = c(-16.7, 7.28),
  Roche_plasma_NfL_out = c(-16.7, 7.28 ),
  QX_plasma_NfL_out = c(-16.7, 7.28 ))

tipping_points_eyo <-c(
  Roche_plasma_Ab42_Ab40_out = -14.4,
  C2N_plasma_Abeta42_Abeta40_out = -14.4, 
  Fuji_plasma_Ab42_Ab40_out = -13.7 , 
  QX_plasma_Ab42_Ab40_out = -8.42,
  C2N_plasma_ptau217_out = -9.91 ,
  Fuji_plasma_ptau217_out = -9.91,
  
  Roche_plasma_ptau181 = -9.76,
  QX_plasma_ptau181_out = -9.44  ,
  
  AlzPath_plasma_ptau217_out = -10.2 ,
  Janssen_plasma_ptau217_out = -10.2 ,
  C2N_plasma_ptau217_ratio = -11.7,
  Roche_plasma_GFAP_out = -11.3 ,
  QX_plasma_GFAP_out = -10.2 ,
  Roche_plasma_NfL_out = -4.33,
  QX_plasma_NfL_out = -2.95)

plot_list_eyo <- list()


create_segments <- function(final_dataset_plasma,years_symp_onset,dep_var, CDR_3) {
  
  # Create segments from filtered data
  segments <- final_dataset_plasma%>%
    group_by(RID) %>%
    arrange(!!sym( years_symp_onset)) %>%
    mutate(
      years_symp_onset_next = lead(!!sym( years_symp_onset)),
      plasma_bmks_next = lead(!!sym(dep_var)),
      group_next = lead(CDR_3)
    ) %>%
    filter(!is.na(years_symp_onset_next), !is.na(group_next)) %>%
    ungroup()
  
  # Create two segments for each transition
  segments <- segments %>%
    mutate(
      years_symp_onset_mid = (!!sym(years_symp_onset) + years_symp_onset_next) / 2,
      plasma_bmks_mid = (!!sym(dep_var) +  plasma_bmks_next) / 2
    )
  
  # Create first half segments
  first_half <- segments %>%
    transmute(
      years_symp_onset = !!sym(years_symp_onset),
      plasma_bmks = !!sym(dep_var),
      years_symp_onset_end = years_symp_onset_mid,
      plasma_bmks_end = plasma_bmks_mid,
      CDR_3 = CDR_3
    )
  
  # Create second half segments
  second_half <- segments %>%
    transmute(
      years_symp_onset = years_symp_onset_mid,
      plasma_bmks = plasma_bmks_mid,
      years_symp_onset_end = years_symp_onset_next,
      plasma_bmks_end = plasma_bmks_next,
      CDR_3 = group_next
    )
  
  # Combine both halves
  all_segments <- bind_rows(first_half, second_half)
  return(all_segments)
}

# Loop through each y variable and create a scatter plot
for (i in seq_along(plasma_biomarkers)) {
  y_var <- plasma_biomarkers[i]
  
  
  ref_mean <- final_dataset_plasma[[paste0(y_var, "_mean")]]
  
  if (y_var %in% names(segment_ranges_eyo)) {
    x_range <- segment_ranges_eyo[[y_var]]
    tipping_point <- tipping_points_eyo [[y_var]]
    # Create the segments for geom_segment
    all_segments <- create_segments(final_dataset_plasma, "years_symp_onset", y_var, "CDR_3")
    
    # Create the scatter plot
    plot <- ggplot(final_dataset_plasma, aes(x = years_symp_onset, y = !!sym(y_var))) +
      geom_point(size=1, show.legend = F,  aes(color=CDR_3, shape = APOE_binary)) + 
      geom_segment(data = all_segments, aes(x = years_symp_onset, y = plasma_bmks, xend = years_symp_onset_end, yend = plasma_bmks_end,color = CDR_3 ), size = 0.3, show.legend =F) +
      geom_smooth( method="gam", color="grey1", show.legend = F, formula = y ~ s(x, k = 3), size=0.25) +
      geom_hline(yintercept = ref_mean, linetype = "dashed", color = "grey0", size= 0.65)+
      geom_vline(xintercept = tipping_point, linetype = "solid", color = "darkgrey", size= 0.75)+
      geom_ribbon(aes(ymin = !!sym(paste0(y_var, "_ci_lower")), ymax = !!sym(paste0(y_var, "_ci_upper"))), alpha = 0.2, fill = "black") +
      scale_color_manual(values = c("blue2","orange2", "brown3") ) +
      
      labs(x = "Estimated years from symptom onset", y = biomarker_labels[[y_var]]) +
      theme_classic() + 
      theme(
        axis.title.x = element_text(size = 14),  # Change x-axis label size
        axis.title.y = element_text(size = 14))
    
    
    if (!is.null(x_range)){
      full_smooth <- ggplot_build(plot)$data[[3]]
      segment_smooth <- full_smooth %>% filter(x >= x_range[1] & x <= x_range[2])
      
      # Add the thicker line segment
      plot <- plot +
        geom_line(data = segment_smooth, aes(x = x, y = y), color = "grey1", size = 1.25)
    }  
    
    # Add the plot to the list
    plot_list_eyo[[y_var]] <- plot
  }
}


blank_plot <- ggplot() + theme_void()
row1 <- plot_grid(plot_list_eyo[[1]], plot_list_eyo[[2]], plot_list_eyo[[3]],plot_list_eyo[[4]], nrow = 1)
row2 <- plot_grid(plot_list_eyo[[5]], plot_list_eyo[[6]], plot_list_eyo[[7]], blank_plot, nrow = 1)
row3 <- plot_grid(plot_list_eyo[[8]], plot_list_eyo[[9]], plot_list_eyo[[10]], plot_list_eyo[[11]],  nrow = 1)
row4 <- plot_grid(plot_list_eyo[[12]],plot_list_eyo[[13]], plot_list_eyo[[14]], plot_list_eyo[[15]], nrow = 1)

allplots_EYO <- plot_grid(row1, row2, row3, row4,  nrow = 4)


###ARRANGE FOR MAIN FIGURE
blank_plot <- ggplot() + theme_void()

row2 <- plot_grid(plot_list_amyloid[[5]], plot_list_tau[[5]], plot_list_eyo[[5]],  nrow = 1)
row3 <- plot_grid(plot_list_amyloid[[7]], plot_list_tau[[7]], plot_list_eyo[[7]],  nrow = 1)
row4 <- plot_grid(plot_list_amyloid[[10]],plot_list_tau[[10]], plot_list_eyo[[10]], nrow = 1)
row5 <- plot_grid(plot_list_amyloid[[12]],plot_list_tau[[12]], plot_list_eyo[[12]], nrow=1)
row6 <- plot_grid(plot_list_amyloid[[14]],plot_list_tau[[14]], plot_list_eyo[[14]],nrow=1)
allplots_main <- plot_grid(row2, row3, row4, row5, row6,  nrow = 5)

ggsave("~/Documents/FNIH_paper_Amyclockplasma/Figures paper/GAM_MAIN_median_last_presentation.pdf", plot = allplots_main , width = 15, height = 20)


###ARRANGE FOR supplementary FIGURE
blank_plot <- ggplot() + theme_void()
row01 <- plot_grid(plot_list_amyloid[[1]], plot_list_tau[[1]], plot_list_eyo[[1]], nrow = 1)
row1 <- plot_grid(plot_list_amyloid[[2]], plot_list_tau[[2]], plot_list_eyo[[2]], nrow = 1)
row2 <- plot_grid(plot_list_amyloid[[3]], plot_list_tau[[3]], plot_list_eyo[[3]],  nrow = 1)
row3 <- plot_grid(plot_list_amyloid[[4]], plot_list_tau[[4]], plot_list_eyo[[4]],  nrow = 1)
row4 <- plot_grid(plot_list_amyloid[[6]], plot_list_tau[[6]], plot_list_eyo[[6]],  nrow = 1)
allplots_supp1 <- plot_grid(row01, row1, row2, row3, row4,  nrow = 5)


row5 <- plot_grid(plot_list_amyloid[[8]], plot_list_tau[[8]], plot_list_eyo[[8]],  nrow = 1)
row6 <- plot_grid(plot_list_amyloid[[9]],plot_list_tau[[9]], plot_list_eyo[[9]], nrow = 1)
row7 <- plot_grid(plot_list_amyloid[[11]],plot_list_tau[[11]], plot_list_eyo[[11]], nrow = 1)
row8 <- plot_grid(plot_list_amyloid[[13]],plot_list_tau[[13]], plot_list_eyo[[13]], nrow=1)
row9 <- plot_grid(plot_list_amyloid[[15]],plot_list_tau[[15]], plot_list_eyo[[15]],nrow=1)
allplots_supp2 <- plot_grid(row5, row6,row7, row8, row9,  nrow = 5)

ggsave("~/Documents/FNIH_paper_Amyclockplasma/Figures paper/GAM_suppl_current2.pdf", plot = allplots_supp2 , width = 15, height = 20)

##########summary plots for selected biomarerks########
selected_biomarkers <- c( "C2N_plasma_ptau217_out",  
                        "C2N_plasma_ptau217_ratio", "Roche_plasma_ptau181",
                        "Roche_plasma_Ab42_Ab40_out", 
                        "Roche_plasma_GFAP_out", "Roche_plasma_NfL_out", "atrophy", 
                        "PTAU_over_ABETA42", "MMSCORE", "SUVR_compositeRef", 
                        "MesialTemporal", "TemporoParietal")
##calculate  z scores with mean and sd from ref group
final_dataset_plasma <- final_dataset_plasma %>%
  mutate(across(all_of(selected_biomarkers), ~ (. - get(paste0(cur_column(), "_mean"))) / get(paste0(cur_column(), "_sd")),.names = "{.col}_zscore"))

##calculate  z scores with mean from ref group and sd from pos group
sd_values <- subset_IDs_with_positive1_scans %>%
  summarize(across(all_of(selected_biomarkers), sd, na.rm = TRUE))

# Rename the columns to indicate they are standard deviations
sd_values <- sd_values %>%
  rename_with(~ paste0(.x, "_sdpos"), everything())
sd_df <- final_dataset_plasma %>%
  mutate(across(all_of(selected_biomarkers), ~ sd_values[[paste0(cur_column(), "_sdpos")]], .names = "{.col}_sdpos"))


final_dataset_plasma <- sd_df %>%
  mutate(across(all_of(selected_biomarkers), ~ (. - get(paste0(cur_column(), "_mean"))) / get(paste0(cur_column(), "_sdpos")),.names = "{.col}_zscore2"))


#invert z scores
final_dataset_plasma <- final_dataset_plasma %>%
  mutate(
    MMSCORE_zscore = -MMSCORE_zscore,
    atrophy_zscore = -atrophy_zscore,
    Roche_plasma_Ab42_Ab40_out_zscore = -Roche_plasma_Ab42_Ab40_out_zscore
  )



plot_biomarkers <-  c("C2N_plasma_ptau217_out_zscore2", "C2N_plasma_ptau217_ratio_zscore2", 
                      "Roche_plasma_ptau181_zscore2", "Roche_plasma_Ab42_Ab40_out_zscore2", 
                      "Roche_plasma_GFAP_out_zscore2", "Roche_plasma_NfL_out_zscore2", 
                      "PTAU_over_ABETA42_zscore2",  "atrophy_zscore2",
                      "MMSCORE_zscore2", "SUVR_compositeRef_zscore2", "MesialTemporal_zscore2", 
                      "TemporoParietal_zscore2")

long_dataset <- final_dataset_plasma %>%
  select(Amyloid_time, Tau_time, EYO1_tau, all_of(plot_biomarkers)) %>%
  pivot_longer(cols = all_of(plot_biomarkers), names_to = "Biomarker", values_to = "Z_Score")

long_dataset$Biomarker <- factor(long_dataset$Biomarker, levels = plot_biomarkers)

color_palette <- c("forestgreen", "darkgreen", "green", "blue",  "red","brown", "darkgrey","gold2",   "purple3", "coral", "cyan3", "black")


linetype_mapping <- c("solid", "solid", "solid", "solid", "solid", "solid", 
                      "dotted", "dotted", "dotted", "dotted", "dotted", "dotted")

plot_amyloid <-  ggplot(long_dataset, aes(x = Amyloid_time, y = Z_Score, color = Biomarker, linetype = Biomarker)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = FALSE,show.legend = F) +
  scale_color_manual(values = color_palette,  labels=c("C2N p-tau217 (pg/ml)", "C2N %p-tau217 (%)", "Roche p-tau181 (pg/ml)",
                                                       expression(paste("Roche A",beta,"42/A",beta,"40")),   "Roche GFAP (ng/ml)",
                                                       "Roche NfL (pg/ml)", expression(paste("CSF p-tau181/A", beta, "42")),"Cortical thickness meta-ROI", "MMSE", expression(paste("A",beta," PET SUVR")),
                                                       "Tau PET Mesial-Temporal SUVR", "Tau PET Temporo-Parietal SUVR" )) +
  labs(x = "Amyloid time (years)", y = "Z-Score", color = "Biomarker", linetype = "Biomarker") + ylim(-1,4.5) +
  theme_minimal() + scale_linetype_manual(values = linetype_mapping, labels=c("C2N p-tau217 (pg/ml)", "C2N %p-tau217 (%)", "Roche p-tau181 (pg/ml)",
                                                                       expression(paste("Roche A",beta,"42/A",beta,"40")), "Roche GFAP (ng/ml)",
                                                                       "Roche NfL (pg/ml)", expression(paste("CSF p-tau181/A", beta, "42")),
                                                                       "Cortical thickness meta-ROI", "MMSE", expression(paste("A",beta," PET SUVR")),
                                                                       "Tau PET Mesial-Temporal SUVR", "Tau PET Temporo-Parietal SUVR")) +
  theme(legend.position = "right", legend.title = element_blank())

plot_tau <-  ggplot(long_dataset, aes(x = Tau_time, y = Z_Score, color = Biomarker, linetype = Biomarker)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = FALSE,show.legend = F) +
  scale_color_manual(values = color_palette,  labels=c("C2N p-tau217 (pg/ml)", "C2N %p-tau217 (%)", "Roche p-tau181 (pg/ml)",
                                                       expression(paste("Roche A",beta,"42/A",beta,"40")),   "Roche GFAP (ng/ml)",
                                                       "Roche NfL (pg/ml)", expression(paste("CSF p-tau181/A", beta, "42")),"Cortical thickness meta-ROI", "MMSE", expression(paste("A",beta," PET SUVR")),
                                                       "Tau PET Mesial-Temporal SUVR", "Tau PET Temporo-Parietal SUVR" )) +
  scale_linetype_manual(values = linetype_mapping, labels=c("C2N p-tau217 (pg/ml)", "C2N %p-tau217 (%)", "Roche p-tau181 (pg/ml)",
                                                            expression(paste("Roche A",beta,"42/A",beta,"40")), "Roche GFAP (ng/ml)",
                                                            "Roche NfL (pg/ml)", expression(paste("CSF p-tau181/A", beta, "42")),
                                                            "Cortical thickness meta-ROI", "MMSE", expression(paste("A",beta," PET SUVR")),
                                                            "Tau PET Mesial-Temporal SUVR", "Tau PET Temporo-Parietal SUVR")) +
  
  labs(x = "Tau time (years)", y = "Z-Score", color = "Biomarker", linetype = "Biomarker") +ylim(-1,4.5) +
  theme_minimal() + 
  theme(legend.position = "right")

plot_eyo <-  ggplot(long_dataset, aes(x = EYO1_tau, y = Z_Score, color = Biomarker, linetype= Biomarker)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = FALSE, show.legend = F) +
  scale_color_manual(values = color_palette,  labels=c("C2N p-tau217 (pg/ml)", "C2N %p-tau217 (%)", "Roche p-tau181 (pg/ml)",
                                                       expression(paste("Roche A",beta,"42/A",beta,"40")),   "Roche GFAP (ng/ml)",
                                                       "Roche NfL (pg/ml)", expression(paste("CSF p-tau181/A", beta, "42")),"Cortical thickness meta-ROI", "MMSE", expression(paste("A",beta," PET SUVR")),
                                                       "Tau PET Mesial-Temporal SUVR", "Tau PET Temporo-Parietal SUVR" )) +
 
  scale_linetype_manual(values = linetype_mapping, labels=c("C2N p-tau217 (pg/ml)", "C2N %p-tau217 (%)", "Roche p-tau181 (pg/ml)",
                                                            expression(paste("Roche A",beta,"42/A",beta,"40")), "Roche GFAP (ng/ml)",
                                                            "Roche NfL (pg/ml)", expression(paste("CSF p-tau181/A", beta, "42")),
                                                            "Cortical thickness meta-ROI", "MMSE", expression(paste("A",beta," PET SUVR")),
                                                            "Tau PET Mesial-Temporal SUVR", "Tau PET Temporo-Parietal SUVR")) +
  
  
   labs(x = "EYO (years)", y = "Z-Score", color = "Biomarker", linetype = "Biomarker") + ylim(-1,4.5) +
  theme_minimal() +
  theme(legend.position = "right")


allplots_summary <- plot_grid(plot_amyloid, plot_tau, plot_eyo, nrow = 1)

get_legend <- function(my_plot) {
  g <- ggplotGrob(my_plot)
  legend <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")][[1]]
  legend
}
legend <- get_legend(plot_amyloid)
legend_col <- plot_grid(legend, blank_plot, ncol = 1, rel_heights = c(0.80, 0.20))
final_plot <- plot_grid(allplots_summary, legend_col, ncol = 2, rel_widths = c(0.8, 0.2))
ggsave("Summaryplots2.pdf", plot = final_plot , width = 15, height = 5)




#####Forest plots tipping points####
biomarker_labels <- c(
  C2N_plasma_Abeta42_Abeta40_out = expression(paste("C2N A",beta,"42/A",beta,"40")),
  Fuji_plasma_Ab42_Ab40_out = expression(paste("Fujirebio A",beta,"42/A",beta,"40")),
  Roche_plasma_Ab42_Ab40_out = expression(paste("Roche A",beta, "42/A", beta,"40")),
  QX_plasma_Ab42_Ab40_out = expression(paste("Quanterix A",beta, "42/A",beta,"40")),
  Roche_plasma_ptau181 = "Roche p-tau181",
  QX_plasma_ptau181_out = "Quanterix p-tau181",
  C2N_plasma_ptau217_out = "C2N p-tau217",
  Fuji_plasma_ptau217_out = "Fujirebio p-tau217",
  AlzPath_plasma_ptau217_out = "ALZPath p-tau217",
  Janssen_plasma_ptau217_out = "Janssen p-tau217",
  C2N_plasma_ptau217_ratio = "C2N %p-tau217",
  Roche_plasma_GFAP_out = "Roche GFAP",
  QX_plasma_GFAP_out = "Quanterix GFAP",
  Roche_plasma_NfL_out = "Roche NfL",
  QX_plasma_NfL_out = "Quanterix NfL",
  PTAU_over_ABETA42 = expression(paste("CSF p-tau181/A", beta, "42")),
  atrophy= "Cortical thickness meta-ROI",
  CDR_SOB= "CDR-SB",
  MesialTemporal = expression(""^18*"F-flortaucipir mesial-temporal tau PET"),
  TemporoParietal = expression(""^18*"F-flortaucipir temporo-parietal tau PET"),
  SUVR_compositeRef = expression(""^18*"F-florbetapir amyloid PET"))

##amyloid time plot##                               
times_AT <- read.csv("~/Documents/FNIH_paper_Amyclockplasma/times_AT_last.csv")
times_AT$biomarker <- factor(times_AT$biomarker, levels = unique(times_AT$biomarker))

color_palette <- c("blue", "blue","coral","blue","darkgrey", "darkgreen","red", "red","blue","forestgreen", "forestgreen", "green","forestgreen","green","darkblue","darkblue", "forestgreen", "purple3", "orange2", "brown", "brown")


amytime_forestplot_  <- ggplot(times_AT,
                              aes(y=biomarker, color= biomarker)) +
  geom_point(aes(x=time),show.legend = F) +
  geom_linerange(aes(xmin=ul, xmax=il),show.legend = FALSE) +
  scale_y_discrete(limits = rev(levels(times_AT$biomarker)), labels=biomarker_labels) + ggtitle("Amyloid timeline") + scale_color_manual(values=color_palette) +
  labs(y="", x="Estimated years from amyloid PET positivity") +  theme_base() +  theme(plot.title = element_text( size=14),axis.title.x = element_text(size = 12)) + xlim(-10, 35)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")

##tau time plot##
times_tau <- read.csv("~/Documents/FNIH_paper_Amyclockplasma/times_tau_last.csv")
times_tau$biomarker <- factor(times_tau$biomarker, levels = unique(times_tau$biomarker))
color_palette <- c("red", "red","blue","coral","blue","darkgrey","blue","green","purple3","green","forestgreen", "forestgreen", "darkgreen","darkblue","forestgreen","forestgreen", "blue","orange2","darkblue", "brown", "brown")

tautime_forestplot_  <- ggplot(times_tau,
                               aes(y=biomarker, color= biomarker)) +
  geom_point(aes(x=time),show.legend = F) +
  geom_linerange(aes(xmin=ul, xmax=il),show.legend = FALSE) +
  scale_y_discrete(limits = rev(levels(times_tau$biomarker)), labels=biomarker_labels) + ggtitle("Tau timeline") + scale_color_manual(values=color_palette) +
  labs(y="", x="Estimated years from tau PET positivity") +  theme_base() +  theme(plot.title = element_text( size=14),axis.title.x = element_text(size = 12)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") 

##EYO  plot##
times_eyo <- read.csv("~/Documents/FNIH_paper_Amyclockplasma/times_eyo_last.csv")
times_eyo$biomarker <- factor(times_eyo$biomarker, levels = unique(times_eyo$biomarker))
color_palette <- c("blue","coral","blue","darkgrey", "blue", "darkblue","darkgreen","darkblue","red","red","forestgreen", "forestgreen","forestgreen","forestgreen", "green","green","blue", "purple3",  "brown", "brown", "orange2" )

EYOtime_forestplot_  <- ggplot(times_eyo,
                               aes(y=biomarker, color= biomarker)) +
  geom_point(aes(x=time),show.legend = F) +
  geom_linerange(aes(xmin=ul, xmax=il),show.legend = FALSE) +
  scale_y_discrete(limits = rev(levels(times_eyo$biomarker)), labels=biomarker_labels) + ggtitle("Symptom timeline") + scale_color_manual(values=color_palette) +
  labs(y="", x="Estimated years from symptom onset") +  theme_base() +  theme(plot.title = element_text( size=14),axis.title.x = element_text(size = 12)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") 

allplots_forest <- plot_grid(amytime_forestplot_, tautime_forestplot_ , EYOtime_forestplot_, nrow = 1)
ggsave("~/Documents/FNIH_paper_Amyclockplasma/Annals_Neurology_submission/Fig4.pdf", plot = allplots_forest , width = 23, height = 8)


##centiloids plot
cl_AT <- read.csv("~/Documents/FNIH_paper_Amyclockplasma/Times_AT_Cent.csv")
cl_AT$biomarker <- factor(cl_AT$biomarker, levels = unique(cl_AT$biomarker))

color_palette <- c("blue", "blue","coral","blue", "darkgrey","darkgreen", "red", "red","blue", "forestgreen","forestgreen", "green", "forestgreen","green", "darkblue","darkblue", "forestgreen","purple3",  "orange2",  "brown", "brown")


cl_forestplot_  <- ggplot(cl_AT,
                          aes(y=biomarker, color= biomarker)) +
  geom_point(aes(x=time),show.legend = F) +
  geom_linerange(aes(xmin=ul, xmax=il),show.legend = FALSE) +
  scale_y_discrete(limits = rev(levels(cl_AT$biomarker)), labels=biomarker_labels) + ggtitle("Amyloid burden") + scale_color_manual(values=color_palette) +
  labs(y="", x="Centiloids") +  theme_base() +  theme(plot.title = element_text( size=14)) + scale_x_continuous(breaks = pretty_breaks(n = 10)) +
  geom_vline(xintercept = 25.7, linetype = "dashed", color = "black") 
ggsave("~/Documents/FNIH_paper_Amyclockplasma/Annals_Neurology_submission/review/forest_centiloids.pdf", plot = cl_forestplot_ , width = 7, height = 9)

##centiloids model plot
cl_AT <- read.csv("~/Documents/FNIH_paper_Amyclockplasma/Times_Cent_model.csv")
cl_AT$biomarker <- factor(cl_AT$biomarker, levels = unique(cl_AT$biomarker))

color_palette <- c("blue", "blue","coral","blue", "darkgrey","darkgreen", "red", "red","blue", "forestgreen","forestgreen", "green", "forestgreen","green", "darkblue","darkblue", "forestgreen","purple3",  "orange2",  "brown", "brown")


cl_mod_forestplot_  <- ggplot(cl_AT,
                          aes(y=biomarker, color= biomarker)) +
  geom_point(aes(x=time),show.legend = F) +
  geom_linerange(aes(xmin=ul, xmax=il),show.legend = FALSE) +
  scale_y_discrete(limits = rev(levels(cl_AT$biomarker)), labels=biomarker_labels) + ggtitle("Amyloid burden") + 
  labs(y="", x="Centiloids") +  theme_base() +  theme(plot.title = element_text( size=14)) + scale_x_continuous(breaks = pretty_breaks(n = 10)) +
  geom_vline(xintercept = 25.7, linetype = "dashed", color = "black") 
ggsave("~/Documents/FNIH_paper_Amyclockplasma/Annals_Neurology_submission/review/forest_centiloids.pdf", plot = cl_forestplot_ , width = 7, height = 9)

##tau SUVR plot
tau_suvr <- read.csv("~/Documents/FNIH_paper_Amyclockplasma/Times_tau_suvr.csv")
tau_suvr$biomarker <- factor(tau_suvr$biomarker, levels = unique(tau_suvr$biomarker))

color_palette <- c("red", "red","blue","coral","blue","darkgrey","blue","green","purple3","green","forestgreen", "forestgreen", "darkgreen","darkblue","forestgreen","forestgreen", "blue","orange2","darkblue", "brown", "brown")


tau_forestplot_  <- ggplot(tau_suvr,
                          aes(y=biomarker, color= biomarker)) +
  geom_point(aes(x=time),show.legend = F) +
  geom_linerange(aes(xmin=ul, xmax=il),show.legend = FALSE) +
  scale_y_discrete(limits = rev(levels(tau_suvr$biomarker)), labels=biomarker_labels) + ggtitle("Tau burden") + scale_color_manual(values=color_palette) +
  labs(y="", x=expression(""^18*"F-flortaucipir tau PET")) +  theme_base() +  theme(plot.title = element_text( size=14)) + scale_x_continuous(breaks = pretty_breaks(n = 10)) +
  geom_vline(xintercept = 1.41, linetype = "dashed", color = "black") 

ggsave("~/Documents/FNIH_paper_Amyclockplasma/Annals_Neurology_submission/review/forest_centiloids.pdf", plot = cl_forestplot_ , width = 7, height = 9)

suppl_forestplots <- plot_grid(cl_forestplot_ , tau_forestplot_, nrow = 1)
ggsave("~/Documents/FNIH_paper_Amyclockplasma/Annals_Neurology_submission/review/Suppl_forests.pdf", plot = suppl_forestplots , width = 15, height = 8)

blank_plot <- ggplot() + theme_void()
forestplots_all <- plot_grid(amytime_forestplot_, cl_forestplot_, blank_plot, tautime_forestplot_ ,tau_forestplot_, EYOtime_forestplot_,    nrow = 2)
ggsave("~/Documents/FNIH_paper_Amyclockplasma/Annals_Neurology_submission/review/forest_all.pdf", plot = forestplots_all , width = 23, height = 16)

##tau model SUVR plot
tau_suvr <- read.csv("~/Documents/FNIH_paper_Amyclockplasma/times_tau_suvr_model.csv")
tau_suvr$biomarker <- factor(tau_suvr$biomarker, levels = unique(tau_suvr$biomarker))

color_palette <- c("red", "red","blue","coral","blue","darkgrey","blue","green","purple3","green","forestgreen", "forestgreen", "darkgreen","darkblue","forestgreen","forestgreen", "blue","orange2","darkblue", "brown", "brown")


tau_forestplot_mod  <- ggplot(tau_suvr,
                           aes(y=biomarker, color= biomarker)) +
  geom_point(aes(x=time),show.legend = F) +
  geom_linerange(aes(xmin=ul, xmax=il),show.legend = FALSE) +
  scale_y_discrete(limits = rev(levels(tau_suvr$biomarker)), labels=biomarker_labels) + ggtitle("Tau burden") + 
  labs(y="", x=expression(""^18*"F-flortaucipir tau PET")) +  theme_base() +  theme(plot.title = element_text( size=14)) + scale_x_continuous(breaks = pretty_breaks(n = 10)) +
  geom_vline(xintercept = 1.41, linetype = "dashed", color = "black") 

forestplots_mods <- plot_grid(cl_mod_forestplot_,tau_forestplot_mod, nrow = 1)

#forest plots double axis
library(scales)
lookup_table <- read.csv("~/Documents/FNIH_paper_Amyclockplasma/suvr_midpoint.csv")
years_to_centiloid <- approxfun(lookup_table$Amyloid_time, lookup_table$Centiloid)

times_AT <- read.csv("~/Documents/FNIH_paper_Amyclockplasma/times_AT_last.csv")
times_AT$biomarker <- factor(times_AT$biomarker, levels = unique(times_AT$biomarker))

color_palette <- c("blue", "blue","coral","blue","darkgrey", "darkgreen","red", "red","blue","forestgreen", "forestgreen", "green","forestgreen","green","darkblue","darkblue", "forestgreen", "purple3", "orange2", "brown", "brown")


amytime_forestplot_xx  <- ggplot(times_AT,
                               aes(y=biomarker, color= biomarker)) +
  geom_point(aes(x=time),show.legend = F) +
  geom_linerange(aes(xmin=ul, xmax=il),show.legend = FALSE) +
  scale_y_discrete(limits = rev(levels(times_AT$biomarker)), labels=biomarker_labels) + ggtitle("Amyloid timeline") + scale_color_manual(values=color_palette) +
  labs(y="", x="Estimated years from amyloid PET positivity") +  theme_base() +  theme(plot.title = element_text( size=14),axis.title.x = element_text(size = 12)) + xlim(-10, 35)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + scale_x_continuous(
    name = "Estimated years from amyloid PET positivity",
    limits = c(-10, 20),
    sec.axis = sec_axis(~ years_to_centiloid(.), name = "Centiloids")
  )

lookup_table_tau <- read.csv("~/Documents/FNIH_paper_Amyclockplasma/MesialTemp_midpoint.csv")
years_to_suvr <- approxfun(lookup_table_tau$Tau_time, lookup_table_tau$MesialTemp_midpoint)

times_T <- read.csv("~/Documents/FNIH_paper_Amyclockplasma/times_tau_last.csv")
times_T$biomarker <- factor(times_T$biomarker, levels = unique(times_T$biomarker))

color_palette <- c("red", "red","blue","coral","blue","darkgrey","blue","green","purple3","green","forestgreen", "forestgreen", "darkgreen","darkblue","forestgreen","forestgreen", "blue","orange2","darkblue", "brown", "brown")


tautime_forestplot_xx  <- ggplot(times_T,
                                 aes(y=biomarker, color= biomarker)) +
  geom_point(aes(x=time),show.legend = F) +
  geom_linerange(aes(xmin=ul, xmax=il),show.legend = FALSE) +
  scale_y_discrete(limits = rev(levels(times_T$biomarker)), labels=biomarker_labels) + ggtitle("Tau timeline") + scale_color_manual(values=color_palette) +
  labs(y="", x="Estimated years from tau PET positivity") +  theme_base() +  theme(plot.title = element_text( size=14),axis.title.x = element_text(size = 12)) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + scale_x_continuous(
    name = "Estimated years from tau PET positivity",
    limits = c(-15, 15.7),
    sec.axis = sec_axis(~ years_to_suvr(.), name = "Mesial-temporal tau PET SUVR")
  )

forestplots_doublex_all <- plot_grid(amytime_forestplot_xx,tautime_forestplot_xx, EYOtime_forestplot_, nrow = 1)
ggsave("~/Documents/FNIH_paper_Amyclockplasma/Annals_Neurology_submission/review/forests_new.pdf", plot = forestplots_doublex_all , width = 23, height = 8)

#####3D plots #####
install.packages("scatterplot3d")
library(scatterplot3d)


plasma_biomarkers <- c( "C2N_plasma_ptau217_out", "Fuji_plasma_ptau217_out",  "AlzPath_plasma_ptau217_out", "Janssen_plasma_ptau217_out", 
                        "C2N_plasma_ptau217_ratio", "QX_plasma_ptau181_out","Roche_plasma_ptau181",
                        "C2N_plasma_Abeta42_Abeta40_out" ,"Fuji_plasma_Ab42_Ab40_out","Roche_plasma_Ab42_Ab40_out","QX_plasma_Ab42_Ab40_out", 
                        "Roche_plasma_GFAP_out", "QX_plasma_GFAP_out","Roche_plasma_NfL_out","QX_plasma_NfL_out", "atrophy", "PTAU_over_ABETA42", "CDR_SOB", "MMSCORE", "SUVR_compositeRef", "MesialTemporal", "TemporoParietal")

dir.create("3D_plots", showWarnings = FALSE)

color_map <- c("red", "blue")  # Red for one group, blue for the other
labels <- c("APOE-ε4 carrier", "APOE-ε4 non-carrier")
# Loop over each biomarker and create a 3D scatter plot
for (i in seq_along(plasma_biomarkers)) {
  biomarker <- plasma_biomarkers[i]  # Get the name of the biomarker
  
  colors <- color_map[as.factor(final_dataset_plasma$APOE_binary)]
  
  png(filename = sprintf("3D_plots/plot_%02d.png", i), width = 800, height = 600)
  
  # Create a 3D scatter plot for each biomarker
  scatterplot3d(final_dataset_plasma$Amyloid_time, 
                final_dataset_plasma[[biomarker]], 
                final_dataset_plasma$Tau_time, 
                main = biomarker,  # Set the title to the biomarker name
                xlab = "Estimated years from amyloid onset", 
                ylab = biomarker,  # Set y-axis label to the biomarker name
                zlab = "Estimated years from tau onset", 
                color = colors, pch = 19)
  legend("topright", legend = labels, col = color_map, pch = 19)
  dev.off()
  }

image_files <- list.files("3D_plots", pattern = "\\.png$", full.names = TRUE)
images <- lapply(image_files, readPNG)
plots <- lapply(images, rasterGrob)

combined_plot <- plot_grid(plotlist = plots, ncol = 3)  # Adjust ncol as needed
ggsave("combined_3D_plot_grid.pdf", plot = combined_plot, width = 6, height = 8)  # Save as PDF




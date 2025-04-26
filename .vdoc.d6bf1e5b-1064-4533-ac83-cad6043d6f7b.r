#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: setup
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, ggh4x, cowplot, fitdistrplus, MASS, goftest, forcats, nortest, fitdistrplus, ggh4x, PupillometryR, png, grid, remotes, ggthemes, bayestestR, HDInterval, DiagrammeR, magick)
#
#
#
#| label: cleandata
# Obtain the main df using "./R/1_data_process.R"
source(here("R", "data_process.R"))
#
#
#
#| label: sampleSize
# List with the sample sizes from the main database.
source(here("R", "func.R"))
#
hormone <- c("CORT", "Control")
temperature <- c("Cold", "Hot")
#
n_list <- list()
#
for(k in 1:length(hormone)){
  for(l in 1:length(temperature)){
    list_name <- paste0(hormone[k], "_", temperature[l])
    n_list[[list_name]] <- sample(df = clean_df, corti = hormone[k], therm = temperature[l])
  }
}
#
#
#
#
#| label: countclutches
# Count the number of clutches per species
#
clutches <- clean_df %>% 
  distinct(clutch) %>% 
  nrow()
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: fig-Methods
#| fig-cap: "Scheme of our experimental design. In panel A, we show the different stages of our experiment and the main manipulations. In panel B, we show the experimental device used to present the stimuli in the behavioural tests. In panel C, we show the experimental setup for the prey discrimination tests. In panel D, we show the relevant times from our behavioural tests."
#
knitr::include_graphics("./Others/Methods.png")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: models_mitochondrial
#
# Run models mitochondrial physiology (each region separately)
#
var_m <- c("mean_mitodensity", "mean_potential", "mean_ros", "mean_dnadamage", "mean_peroxidation")
regions <- c("OB", "OT")
for (p in var_m){
  for (h in regions){
    if (h == "OB"){
      df <- clean_df %>% filter(region == "OB")
      l <- "OB"
      if (p %in% c("mean_mitodensity", "mean_potential", "mean_ros")){
        formula <- paste0(p, "~ cort*temp + (1|clutch)")
      } else if (p == "mean_dnadamage"){
        formula <- paste0(p, "~ cort*temp + age_euthanasia + sex + (1|clutch)")
      } else {
        formula <- paste0(p, "~ cort*temp + age_euthanasia + (1|clutch)")
      }
    } else {
      df <- clean_df %>% filter(region == "OT")
      l <- "OT"
      if (p %in% c("mean_mitodensity", "mean_potential", "mean_ros")){
        formula <- paste0(p, "~ cort*temp + (1|clutch)")
      } else if (p == "mean_dnadamage"){
        formula <- paste0(p, "~ cort*temp + sex + (1|clutch)")
      } else {
        formula <- paste0(p, "~ cort*temp + age_euthanasia + (1|clutch)")
      }
    }

    pmodel_name <- paste0("m_def_", p, "_", h)
    assign(pmodel_name, fit_m(df = df,
                              cat = "def",
                              var = p,
                              formula = formula,
                              fam = gaussian(),
                              label = l,
                              refit = FALSE),
          envir = .GlobalEnv)  # Assign to the global environment
  }
} 
#
#
#
#| label: models_behaviour
# Fitting the model and extraction of posteriors for Detection Latency (log-normal).
source(here("R", "func.R"))
#
#
## Run model behaviour for both stimuli separately
#
beh_df <- clean_df
stimuli <- c("Chemical", "Visual")
for (s in stimuli){
  df_b <- beh_df %>% filter(stimulus == s)
  formula_t_D <- t_D ~ cort*temp + (1|clutch) + (1|lizard_id)

  pmodel_name <- paste0("m_def_t_D_", s)
  assign(pmodel_name, fit_m(df = df_b,
                             cat = "def",
                             var = "t_D",
                             formula = formula_t_D,
                             fam = gaussian(),
                             label = s,
                             refit = FALSE),
          envir = .GlobalEnv)  # Assign to the global environment
}
#
#
#
#| label: organise_posteriors
#
# Organising the posteriors of the previous models to fit the tables below
#
source(here("R", "func.R"))
#
# Databases for each region/stimulus
post_OB <- data.frame()
post_OT <- data.frame()
# Names posteriors:
region <- c("OB", "OT")
names_OB <- c("m_def_mean_mitodensity_OB", "m_def_mean_potential_OB", "m_def_mean_ros_OB", "m_def_mean_dnadamage_OB", "m_def_mean_peroxidation_OB", "m_def_t_D_Chemical")
names_OT <- c("m_def_mean_mitodensity_OT", "m_def_mean_potential_OT", "m_def_mean_ros_OT", "m_def_mean_dnadamage_OT", "m_def_mean_peroxidation_OT", "m_def_t_D_Visual")
#
# Organising the results
for (r in region) {
  model_select <- get(paste0("names_", r))   
  for (pos in model_select) {
    model <- get(pos)      # Get the model from the global environment
    post_result <- tidy_post(model)        # Apply tidy_post to each model
    # Add a new column to identify the region and model
    post_result$Region <- r
    post_result$Model <- pos
    # Append to the appropriate data frame
    if (r == "OB") {
      post_OB <- bind_rows(post_OB, post_result)
    } else {
      post_OT <- bind_rows(post_OT, post_result)
    }
  }
}
#
#
#
#| label: values_posteriors
#
# Extracting the posteriors for the models and the values of interest. Here, I am creating dfs for each of the variables with the values for all the prenatal conditions to make contrasts easier to write.
#
source(here("R", "func.R"))
#
# A) Olfactory bulbs/Chemical stimulus
#
MD_OB <- post_values(m_def_mean_mitodensity_OB, "none")
MP_OB <- post_values(m_def_mean_potential_OB, "none")
ROS_OB <- post_values(m_def_mean_ros_OB, "none")
DNA_OB <- post_values(m_def_mean_dnadamage_OB, "none")
LP_OB <- post_values(m_def_mean_peroxidation_OB, "none")
DET_OB <- post_values(m_def_t_D_Chemical, "none")
#
# B) Optic tecta/Visual stimulus
MD_OT <- post_values(m_def_mean_mitodensity_OT, "none")
MP_OT <- post_values(m_def_mean_potential_OT, "none")
ROS_OT <- post_values(m_def_mean_ros_OT, "none")
DNA_OT <- post_values(m_def_mean_dnadamage_OT, "sex")
LP_OT <- post_values(m_def_mean_peroxidation_OT, "none")
DET_OT <- post_values(m_def_t_D_Visual, "none")
#
#
#
#
#
#
#
#
#| label: fig-results_energy
#| fig-cap: "Estimates of mitochondrial density (A, C) and metabolic capacity (B, D) in the olfactory bulbs (A, B) and optic tecta (C, D) of L. delicata hatchlings as a function of the different prenatal conditions. Black dots indicate the posterior mean, and the bars represent the SD of the estimates. The y-axis represents the posterior estimates of the variable of interest, and the x-axis represents the different prenatal conditions. Lines with asterisks represent significant differences between groups based on pMCMC values (pMCMC < 0.05), no lines indicate no significant differences."
#| fig-name: "fig-results_energy"
#
source(here("R", "func.R"))
#
# A) Plotting the results for OB/Chemical stimulus using the df from before
plot_density_OB <- plotting(MD_OB, "Mit density")
plot_potential_OB <- plotting(MP_OB, "Mit potential")
#
fig_OB_energy <- plot_grid(plot_density_OB, NULL, plot_potential_OB, NULL,
                    nrow = 1, rel_widths = c(0.9, 0.1, 0.9, 0.3))
#
# B) Plotting the results for OT/Visual stimulus
plot_density_OT <- plotting(MD_OT, "Mit density")
plot_potential_OT <- plotting(MP_OT, "Mit potential")
#
fig_OT_energy <- plot_grid(plot_density_OT, NULL, plot_potential_OT, NULL,
                    nrow = 1, rel_widths = c(0.9, 0.1, 0.9, 0.3))
#
# Getting the legend (get_legend() does not work if the legend is on the top or bottom)
plot_legend_bottom <- plotting(MP_OB, "Mit potential") + theme(legend.position = "bottom", legend.title = element_blank())
gtable <- ggplot_gtable(ggplot_build(plot_legend_bottom))
legend_plot <- gtable$grobs[[which(sapply(gtable$grobs, function(x) x$name) == "guide-box")]]
# C) Merging plots
#
# Create the figure grid with extra space for images
fig_results_energy <- plot_grid(
  fig_OB_energy, fig_OT_energy, NULL,
  nrow = 3, rel_heights = c(0.9, 0.9, 0.15))
# Final composition: Merge everything, adding images and legend
final_plot_energy <- ggdraw(fig_results_energy) +
  # Insert OB image in the top right
  draw_image("./Others/OB.png", x = 0.87, y = 0.79, width = 0.12, height = 0.2) +
  # Insert OT image in the middle right
  draw_image("./Others/OT.png", x = 0.87, y = 0.34, width = 0.11, height = 0.18) +
  # Insert legend at the bottom-right
  draw_grob(legend_plot, x = 0.45, y = 0.05, width = 0.0001, height = 0.0001) +
  # Insert title for each plot
  annotate("text", x = 0.032, y = 0.97, label = "A", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.487, y = 0.97, label = "B", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.032, y = 0.515, label = "C", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.49, y = 0.515, label = "D", hjust = 1, vjust = 1, size = 7, fontface = "bold")
#
# Print final plot
ggsave(here("./output/figures/text/results_energy.png"), plot = final_plot_energy, width = 21, height = 14, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/text/results_energy.png")
#
#
#
#
#
#| label: fig-results_oxidative
#| fig-cap: "Estimates of ROS (A, D), DNA damage (B, E), and lipid peroxidation (C, F) in the olfactory bulbs (A - C) and optic tecta (D - F) of L. delicata hatchlings as a function of the different prenatal conditions. Black dots indicate the posterior mean, and the bars represent the SD of the estimates. The y-axis represents the posterior estimates of the variable of interest, and the x-axis represents the different prenatal conditions. Lines with asterisks represent significant differences between groups based on pMCMC values (pMCMC < 0.05), no lines indicate no significant differences."
#| fig-name: "fig-results_oxidative"
#
source(here("R", "func.R"))
#
# A) Plotting the results for OB/Chemical stimulus
plot_ros_OB <- plotting(ROS_OB, "ROS")
plot_dnadamage_OB <- plotting(DNA_OB, "DNA damage")
plot_peroxidation_OB <- plotting(LP_OB, "Lipid peroxidation")
#
fig_OB_oxidative_top <- plot_grid(plot_ros_OB, NULL, NULL,
                    nrow = 1, rel_widths = c(0.8, 0.1, 0.8))
fig_OB_oxidative_bottom <- plot_grid(plot_dnadamage_OB, NULL, plot_peroxidation_OB,
                    nrow = 1, rel_widths = c(0.8, 0.1, 0.8))
fig_OB_oxidative <- plot_grid(fig_OB_oxidative_top, fig_OB_oxidative_bottom,
                    nrow = 2)
#
# B) Plotting the results for OT/Visual stimulus
plot_ros_OT <- plotting(ROS_OT, "ROS")
plot_dnadamage_OT <- plotting(DNA_OT, "DNA damage")
plot_peroxidation_OT <- plotting(LP_OT, "Lipid peroxidation")
#
fig_OT_oxidative_top <- plot_grid(plot_ros_OT, NULL, NULL,
                    nrow = 1, rel_widths = c(0.8, 0.1, 0.8))
fig_OT_oxidative_bottom <- plot_grid(plot_dnadamage_OT, NULL, plot_peroxidation_OT,
                    nrow = 1, rel_widths = c(0.8, 0.1, 0.8))
fig_OT_oxidative <- plot_grid(fig_OT_oxidative_top, fig_OT_oxidative_bottom,
                    nrow = 2)
#
# Getting the legend (get_legend() does not work if the legend is on the top or bottom)
plot_legend_bottom <- plotting(ROS_OB, "ROS") + theme(legend.position = "bottom", legend.title = element_blank())
gtable <- ggplot_gtable(ggplot_build(plot_legend_bottom))
legend_plot <- gtable$grobs[[which(sapply(gtable$grobs, function(x) x$name) == "guide-box")]]
# C) Merging plots
#
# Create the figure grid with extra space for images
fig_results_oxidative <- plot_grid(
  fig_OB_oxidative, NULL, fig_OT_oxidative, NULL,
  nrow = 4, rel_heights = c(1, 0.05, 1, 0.1))
# Final composition: Merge everything, adding images and legend
final_plot_oxidative <- ggdraw(fig_results_oxidative) +
  # Insert OB image in the top right
  draw_image("./Others/OB.png", x = 0.56, y = 0.755, width = 0.2, height = 0.3) +
  # Insert OT image in the middle right
  draw_image("./Others/OT.png", x = 0.56, y = 0.28, width = 0.19, height = 0.28) +
  # Insert legend at the bottom-right
  draw_grob(legend_plot, x = 0.5, y = 0.025, width = 0.0001, height = 0.0001) +
  # Insert title for each plot
  annotate("text", x = 0.0365, y = 0.986, label = "A", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.0365, y = 0.755, label = "B", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.567, y = 0.755, label = "C", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.0365, y = 0.5, label = "D", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.0365, y = 0.268, label = "E", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.567, y = 0.268, label = "F", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("segment", x = 0.46, xend= 0.46, y = 0.21, yend = 0.26, size = 0.5, colour = "black") +
  annotate("segment", x = 0.45, xend= 0.46, y = 0.21, yend = 0.21, size = 0.5, colour = "black") +
  annotate("segment", x = 0.45, xend= 0.46, y = 0.26, yend = 0.26, size = 0.5, colour = "black") +
  annotate("text", x = 0.47, y = 0.228, label = "*", hjust = 0.5, vjust = 0.5, size = 7)
#
# Print final plot
ggsave(here("./output/figures/text/results_oxidative.png"), plot = final_plot_oxidative, width = 18, height = 24, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/text/results_oxidative.png")
#
#
#
#
#
#
#
#
#| label: fig-results_behaviour
#| fig-cap: "Estimates of detection latency of chemical (A) and visual (B) stimulus by L. delicata hatchlings as a function of the different prenatal conditions. Black dots indicate the posterior mean, and the bars represent the SD of the estimates. The y-axis represents the posterior estimates of the variable of interest, and the x-axis represents the different prenatal conditions. Lines with asterisks represent significant differences between groups based on pMCMC values (pMCMC < 0.05), no lines indicate no significant differences."
#| fig-name: "fig-results"
#
source(here("R", "func.R"))
# A) Plotting the results for OB/Chemical stimulus
plot_t_D_Chemical <- plotting(DET_OB, "Detection latency")
#
# B) Plotting the results for OT/Visual stimulus
plot_t_D_Visual <- plotting(DET_OT, "Detection latency")
#
# Getting the legend (get_legend() does not work if the legend is on the top or bottom)
plot_legend_bottom <- plotting(DET_OB, "Detection latency") + theme(legend.position = "bottom", legend.title = element_blank())
gtable <- ggplot_gtable(ggplot_build(plot_legend_bottom))
legend_plot <- gtable$grobs[[which(sapply(gtable$grobs, function(x) x$name) == "guide-box")]]
# C) Merging plots
#
# Create the figure grid with extra space for images
fig_results_det <- plot_grid(plot_t_D_Chemical, NULL, plot_t_D_Visual, NULL,
                    nrow = 1, rel_widths = c(1, 0.4, 1, 0.35))
fig_results_behaviour <- plot_grid(fig_results_det, NULL,
  nrow = 2, rel_heights = c(1, 0.15))
# Final composition: Merge everything, adding images and legend
final_plot_behaviour <- ggdraw(fig_results_behaviour) +
  # Insert CS image in the top middle
  draw_image("./Others/CS.png", x = 0.19, y = 0.36, width = 0.48, height = 0.6) +
  # Insert VS image in the top right
  draw_image("./Others/VS.png", x = 0.6, y = 0.35, width = 0.65, height = 0.6) +
  # Insert legend at the bottom-right
  draw_grob(legend_plot, x = 0.45, y = 0.07, width = 0.0001, height = 0.0001) +
  # Insert title for each plot
  annotate("text", x = 0.035, y = 0.96, label = "A", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.542, y = 0.96, label = "B", hjust = 1, vjust = 1, size = 7, fontface = "bold") +  
  annotate("segment", x = 0.34, xend= 0.34, y = 0.73, yend = 0.91, size = 0.5, colour = "black") +
  annotate("segment", x = 0.335, xend= 0.34, y = 0.73, yend = 0.73, size = 0.5, colour = "black") +
  annotate("segment", x = 0.335, xend= 0.34, y = 0.91, yend = 0.91, size = 0.5, colour = "black") +
  annotate("text", x = 0.35, y = 0.795, label = "*", hjust = 0.5, vjust = 0.5, size = 7) +
  annotate("segment", x = 0.292, xend= 0.292, y = 0.418, yend = 0.598, size = 0.5, colour = "black") +
  annotate("segment", x = 0.287, xend= 0.292, y = 0.418, yend = 0.418, size = 0.5, colour = "black") +
  annotate("segment", x = 0.287, xend= 0.292, y = 0.598, yend = 0.598, size = 0.5, colour = "black") +
  annotate("text", x = 0.302, y = 0.471, label = "*", hjust = 0.5, vjust = 0.5, size = 7)
#
# Print final plot
ggsave(here("./output/figures/text/results_behaviour.png"), plot = final_plot_behaviour, width = 21, height = 7, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/text/results_behaviour.png")
#
#
#
#
#
#
#| label: model_sem_OB
# Making the SEM model by using a multivariate brms. The aim is to test the relationships between mitochondrial physiology and detection latency.
# To simplify the random factors, and since there is no effect of experience with the prey in prey detection, latency was averaged across the two preys. Every other predictor added to the model was based on previous brms for each variable separated.
# All continue variables were standardized (var/2SD) before running the models (see data_process.R). 
#
source(here("R", "func.R"))
#
SEM_df <- clean_df %>%
  group_by(lizard_id, prey) %>%
  mutate(t_D = mean(t_D, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(prey == "Unknown")
#
#
data_sem_OB <- SEM_df %>%
  filter(region == "OB") %>%
  mutate(obs = as.integer(c(1:80))) %>%
  mutate(vec = rep(1, length(obs)))
# Create the models
#
refit <- FALSE
#
#
if(refit){
  m_OB <- brm(
    bf(t_D | se(vec, sigma = FALSE) ~ mean_mitodensity + mean_potential + mean_dnadamage + mean_peroxidation + (1|clutch) + (1|q|obs)) +
    bf(mean_dnadamage | se(vec, sigma = FALSE) ~ mean_ros + (1|clutch) + (1|p|obs)) +
    bf(mean_peroxidation | se(vec, sigma = FALSE) ~ mean_ros + (1|clutch)+ (1|p|obs)) +
    bf(mean_ros | se(vec, sigma = FALSE) ~ mean_mitodensity + mean_potential + (1|clutch) + (1|t|obs)) +
  set_rescor(FALSE),
  family = gaussian(),
  data = data_sem_OB,
  chains = 4, cores = 4, iter = 8000, warmup = 2000,
  control = list(adapt_delta = 0.99, max_treedepth = 11))
  # Save the model
  saveRDS(m_OB, file = here("output/m_SEM/m_OB.rds"))
} else {
  m_OB <- readRDS(here("output/m_SEM/m_OB.rds"))
}
#
#
#
#| label: sem_tidy_OB
source(here("R", "func.R"))
#
# I am extracting here all the values for getting the total effects of each of the variables in the model. I am using the posterior values for each of the variables to get the total effects assuming that:
## total effect = direct effect + indirect effect + residual correlation
# In other words:
## total effect = 
# Extract the posteriors for the SEM model
post_sem_OB <- as_draws_df(m_OB) 
#
#### A) Get the direct paths per each variable
# 
# A.1) Detection
OB_coeff_mitodensity_det <- post_sem_OB$b_tD_mean_mitodensity
OB_coeff_potential_det <- post_sem_OB$b_tD_mean_potential
OB_coeff_dna_det <- post_sem_OB$b_tD_mean_dnadamage
OB_coeff_perox_det <- post_sem_OB$b_tD_mean_peroxidation
# 
# A.2) DNA damage
OB_coeff_ros_dna <- post_sem_OB$b_meandnadamage_mean_ros
#
# A.3) Lipid peroxidation
OB_coeff_ros_perox <- post_sem_OB$b_meanperoxidation_mean_ros
#
# A.4) ROS
OB_coeff_mitodensity_ros <- post_sem_OB$b_meanros_mean_mitodensity
OB_coeff_potential_ros <- post_sem_OB$b_meanros_mean_potential
#
#
#### B) Get the indirect paths for each variable
#
# B.1) Detection
OB_undir_ros_det <- OB_coeff_ros_dna * OB_coeff_dna_det + OB_coeff_ros_perox * OB_coeff_perox_det
OB_undir_mitodensity_det <- OB_coeff_mitodensity_ros * OB_coeff_ros_dna * OB_coeff_dna_det + OB_coeff_mitodensity_ros * OB_coeff_ros_perox * OB_coeff_perox_det
OB_undir_potential_det <- OB_coeff_potential_ros * OB_coeff_ros_dna * OB_coeff_dna_det + OB_coeff_potential_ros * OB_coeff_ros_perox * OB_coeff_perox_det
#
# B.2) DNA damage
OB_undir_mitodensity_dna <- OB_coeff_mitodensity_ros * OB_coeff_ros_dna
OB_undir_potential_dna <- OB_coeff_potential_ros * OB_coeff_ros_dna
#
# B.3) Lipid peroxidation
OB_undir_mitodensity_perox <- OB_coeff_mitodensity_ros * OB_coeff_ros_perox
OB_undir_potential_perox <- OB_coeff_potential_ros * OB_coeff_ros_perox
#
#
#### C) Get the total effects for each variable
#
# C.1) Detection
OB_total_mitodensity_det <- OB_undir_mitodensity_det
OB_total_potential_det <- OB_undir_potential_det
OB_total_ros_det <- OB_undir_ros_det
OB_total_dna_det <- OB_coeff_dna_det
OB_total_perox_det <- OB_coeff_perox_det
#
# C.2) DNA damage
OB_total_ros_dna <- OB_coeff_ros_dna
OB_total_mitodensity_dna <- OB_undir_mitodensity_dna
OB_total_potential_dna <- OB_undir_potential_dna
#
# C.3) Lipid peroxidation
OB_total_ros_perox <- OB_coeff_ros_perox
OB_total_mitodensity_perox <- OB_undir_mitodensity_perox
OB_total_potential_perox <- OB_undir_potential_perox
#
# C.4) ROS
OB_total_mitodensity_ros <- OB_coeff_mitodensity_ros
OB_total_potential_ros <- OB_coeff_potential_ros
#
#
# D) Create a df with the values for each variable
#
# D.1) Detection (for example)
det_results_semOB <- data.frame(
  reference_variable = rep("t_D", 5),
  predictor_modulator = c("mean_mitodensity",
                          "mean_potential",
                          "ROS",
                          "mean_dnadamage",
                          "mean_peroxidation"),
  direct_effects = I(list(OB_coeff_mitodensity_det,
                          OB_coeff_potential_det,
                          NA,
                          OB_coeff_dna_det,
                          OB_coeff_perox_det)),
  indirect_effects = I(list(OB_undir_mitodensity_det,
                          OB_undir_potential_det,
                          OB_undir_ros_det,
                          NA,
                          NA)),
  total_effects = I(list(OB_total_mitodensity_det,
                        OB_total_potential_det,
                        OB_total_ros_det,
                        OB_total_dna_det,
                        OB_total_perox_det))
  ) %>% 
  group_by(predictor_modulator) %>%
  summarize(
    # Summarizing each list column (direct_effects, indirect_effects, etc.)
    mean_direct_effects = format_dec(mean(unlist(direct_effects), na.rm = TRUE), 3),
    mean_indirect_effects = format_dec(mean(unlist(indirect_effects), na.rm = TRUE), 3),
    mean_total_effects = format_dec(mean(unlist(total_effects), na.rm = TRUE), 3),
    q5_direct = format_dec(quantile(unlist(direct_effects), 0.05, na.rm = TRUE), 3),
    q95_direct = format_dec(quantile(unlist(direct_effects), 0.95, na.rm = TRUE), 3),
    q5_indirect = format_dec(quantile(unlist(indirect_effects), 0.05, na.rm = TRUE), 3),
    q95_indirect = format_dec(quantile(unlist(indirect_effects), 0.95, na.rm = TRUE), 3),
    q5_total = format_dec(quantile(unlist(total_effects), 0.05, na.rm = TRUE), 3),
    q95_total = format_dec(quantile(unlist(total_effects), 0.95, na.rm = TRUE), 3)) %>%
  mutate(source = "det_results_semOB")
#
# D.2) DNA damage
dna_results_semOB <- data.frame(
  reference_variable = rep("DNA damage", 3),
  predictor_modulator = c("mean_mitodensity",
                          "mean_potential",
                          "ROS"),
  direct_effects = I(list(NA,
                          NA,
                          OB_coeff_ros_dna)),
  indirect_effects = I(list(OB_undir_mitodensity_dna,
                          OB_undir_potential_dna,
                          NA)),
  total_effects = I(list(OB_total_mitodensity_dna,
                        OB_total_potential_dna,
                        OB_total_ros_dna))) %>% 
  group_by(predictor_modulator) %>%
  summarize(
    # Summarizing each list column (direct_effects, indirect_effects, etc.)
    mean_direct_effects = format_dec(mean(unlist(direct_effects), na.rm = TRUE), 3),
    mean_indirect_effects = format_dec(mean(unlist(indirect_effects), na.rm = TRUE), 3),
    mean_total_effects = format_dec(mean(unlist(total_effects), na.rm = TRUE), 3),
    q5_direct = format_dec(quantile(unlist(direct_effects), 0.05, na.rm = TRUE), 3),
    q95_direct = format_dec(quantile(unlist(direct_effects), 0.95, na.rm = TRUE), 3),
    q5_indirect = format_dec(quantile(unlist(indirect_effects), 0.05, na.rm = TRUE), 3),
    q95_indirect = format_dec(quantile(unlist(indirect_effects), 0.95, na.rm = TRUE), 3),
    q5_total = format_dec(quantile(unlist(total_effects), 0.05, na.rm = TRUE), 3),
    q95_total = format_dec(quantile(unlist(total_effects), 0.95, na.rm = TRUE), 3)) %>%
  mutate(source = "dna_results_semOB")
#
# D.3) Lipid peroxidation
perox_results_semOB <- data.frame(
  reference_variable = rep("lipid peroxidation", 3),
  predictor_modulator = c("mean_mitodensity",
                          "mean_potential",
                          "ROS"),
  direct_effects = I(list(NA,
                          NA,
                          OB_coeff_ros_perox)),
  indirect_effects = I(list(OB_undir_mitodensity_perox,
                          OB_undir_potential_perox,
                          NA)),
  total_effects = I(list(OB_total_mitodensity_perox,
                        OB_total_potential_perox,
                        OB_total_ros_perox))) %>% 
  group_by(predictor_modulator) %>%
  summarize(
    # Summarizing each list column (direct_effects, indirect_effects, etc.)
    mean_direct_effects = format_dec(mean(unlist(direct_effects), na.rm = TRUE), 3),
    mean_indirect_effects = format_dec(mean(unlist(indirect_effects), na.rm = TRUE), 3),
    mean_total_effects = format_dec(mean(unlist(total_effects), na.rm = TRUE), 3),
    q5_direct = format_dec(quantile(unlist(direct_effects), 0.05, na.rm = TRUE), 3),
    q95_direct = format_dec(quantile(unlist(direct_effects), 0.95, na.rm = TRUE), 3),
    q5_indirect = format_dec(quantile(unlist(indirect_effects), 0.05, na.rm = TRUE), 3),
    q95_indirect = format_dec(quantile(unlist(indirect_effects), 0.95, na.rm = TRUE), 3),
    q5_total = format_dec(quantile(unlist(total_effects), 0.05, na.rm = TRUE), 3),
    q95_total = format_dec(quantile(unlist(total_effects), 0.95, na.rm = TRUE), 3)) %>%
  mutate(source = "perox_results_semOB")
#
# D.4) ROS
ros_results_semOB <- data.frame(
  reference_variable = rep("ROS", 2),
  predictor_modulator = c("mean_mitodensity",
                          "mean_potential"),
  direct_effects = I(list(OB_coeff_mitodensity_ros,
                          OB_coeff_potential_ros)),
  indirect_effects = I(list(NA,
                          NA)),
  total_effects = I(list(OB_total_mitodensity_ros,
                        OB_total_potential_ros))) %>% 
  group_by(predictor_modulator) %>%
  summarize(
    # Summarizing each list column (direct_effects, indirect_effects, etc.)
    mean_direct_effects = format_dec(mean(unlist(direct_effects), na.rm = TRUE), 3),
    mean_indirect_effects = format_dec(mean(unlist(indirect_effects), na.rm = TRUE), 3),
    mean_total_effects = format_dec(mean(unlist(total_effects), na.rm = TRUE), 3),
    q5_direct = format_dec(quantile(unlist(direct_effects), 0.05, na.rm = TRUE), 3),
    q95_direct = format_dec(quantile(unlist(direct_effects), 0.95, na.rm = TRUE), 3),
    q5_indirect = format_dec(quantile(unlist(indirect_effects), 0.05, na.rm = TRUE), 3),
    q95_indirect = format_dec(quantile(unlist(indirect_effects), 0.95, na.rm = TRUE), 3),
    q5_total = format_dec(quantile(unlist(total_effects), 0.05, na.rm = TRUE), 3),
    q95_total = format_dec(quantile(unlist(total_effects), 0.95, na.rm = TRUE), 3)) %>%
  mutate(source = "ros_results_semOB")
#
#
# E) Merge everything into a single df
#
sem_results_OB <- bind_rows(det_results_semOB, dna_results_semOB, perox_results_semOB, ros_results_semOB)
#
#
#
#
#
#| label: fig-sem_results_OB
#| fig-cap: "Structural Equation Models for OB/Chemical stimulus"
#| fig-name: "fig-sem_results_OB"
#
# A) Getting all the direct coefficients
OB_density_det <- paste0(format_dec(mean(OB_coeff_mitodensity_det), 3),
                  " [", format_dec(quantile(OB_coeff_mitodensity_det, 0.05), 3),
                  ", ", format_dec(quantile(OB_coeff_mitodensity_det, 0.95), 3), "]")
OB_potential_det <- paste0(format_dec(mean(OB_coeff_potential_det), 3),
                  " [", format_dec(quantile(OB_coeff_potential_det, 0.05), 3),
                  ", ", format_dec(quantile(OB_coeff_potential_det, 0.95), 3), "]")
OB_dna_det <- paste0(format_dec(mean(OB_coeff_dna_det), 3),
                  " [", format_dec(quantile(OB_coeff_dna_det, 0.05), 3),
                  ", ", format_dec(quantile(OB_coeff_dna_det, 0.95), 3), "]")
OB_perox_det <- paste0(format_dec(mean(OB_coeff_perox_det), 3),
                  " [", format_dec(quantile(OB_coeff_perox_det, 0.05), 3),
                  ", ", format_dec(quantile(OB_coeff_perox_det, 0.95), 3), "]")
#
OB_ros_dna <- paste0(format_dec(mean(OB_coeff_ros_dna), 3),
                  " [", format_dec(quantile(OB_coeff_ros_dna, 0.05), 3),
                  ", ", format_dec(quantile(OB_coeff_ros_dna, 0.95), 3), "]")
#
OB_ros_perox <- paste0(format_dec(mean(OB_coeff_ros_perox), 3),
                  " [", format_dec(quantile(OB_coeff_ros_perox, 0.05), 3),
                  ", ", format_dec(quantile(OB_coeff_ros_perox, 0.95), 3), "]")
#
OB_density_ros <- paste0(format_dec(mean(OB_coeff_mitodensity_ros), 3),
                  " [", format_dec(quantile(OB_coeff_mitodensity_ros, 0.05), 3),
                  ", ", format_dec(quantile(OB_coeff_mitodensity_ros, 0.95), 3), "]")
OB_potential_ros <- paste0(format_dec(mean(OB_coeff_potential_ros), 3),
                  " [", format_dec(quantile(OB_coeff_potential_ros, 0.05), 3),
                  ", ", format_dec(quantile(OB_coeff_potential_ros, 0.95), 3), "]")
#
imgOB <- readPNG(here("Others", "SEM.png"))
plot_SEM_OB <- rasterGrob(imgOB, interpolate = TRUE)
#
fig_SEM_OB <- ggdraw(plot_SEM_OB) +
  draw_image("./Others/OB.png", x = 0.855, y = 0.67, width = 0.18, height = 0.3) +
  annotate("text", x = 0.52, y = 0.075, label = OB_density_det, hjust = 1, vjust = 1, size = 2.8, family = "Times") +
  annotate("text", x = 0.52, y = 0.94, label = OB_potential_det, hjust = 1, vjust = 1, size = 2.8, family = "Times") +
  annotate("text", x = 0.83, y = 0.715, label = OB_dna_det, hjust = 1, vjust = 1, size = 2.8, family = "Times") +
  annotate("text", x = 0.83, y = 0.335, label = OB_perox_det, hjust = 1, vjust = 1, size = 2.8, family = "Times") +
  annotate("text", x = 0.525, y = 0.805, label = OB_ros_dna, hjust = 1, vjust = 1, size = 2.8, family = "Times") +
  annotate("text", x = 0.515, y = 0.215, label = OB_ros_perox, hjust = 1, vjust = 1, size = 2.8, family = "Times") +
  annotate("text", x = 0.282, y = 0.2, label = OB_density_ros, hjust = 1, vjust = 1, size = 2.8, family = "Times") +
  annotate("text", x = 0.282, y = 0.82, label = OB_potential_ros, hjust = 1, vjust = 1, size = 2.8, family = "Times")
ggsave(here("./output/figures/text/SEM_OB.png"), plot = fig_SEM_OB, width = 21, height = 7.5, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/text/SEM_OB.png")
#
#
#
#| label: models_sem_OT
# Making the same models above but for OT/Visual stimulus
#
source(here("R", "func.R"))
#
#
data_sem_OT <- SEM_df %>%
  filter(region == "OT") %>%
  mutate(obs = as.integer(c(1:80))) %>%
  mutate(vec = rep(1, length(obs)))
#
# Create the models
#
refit <- FALSE
#
if(refit){
  m_OT <- brm(
    bf(t_D | se(vec, sigma = FALSE) ~ mean_mitodensity + mean_potential + mean_dnadamage + mean_peroxidation + (1|clutch) + (1|p|obs)) +
    bf(mean_dnadamage | se(vec, sigma = FALSE) ~ mean_ros + (1|clutch) + (1|q|obs)) +
    bf(mean_peroxidation | se(vec, sigma = FALSE) ~ mean_ros + (1|clutch) + (1|q|obs)) +
    bf(mean_ros | se(vec, sigma = FALSE) ~ mean_mitodensity + mean_potential + (1|clutch) + (1|t|obs)) +
  set_rescor(FALSE),
  family = gaussian(),
  data = data_sem_OT, 
  chains = 4, cores = 4, iter = 8000, warmup = 2000,
  control = list(adapt_delta = 0.99, max_treedepth = 11))
  # Save the model
  saveRDS(m_OT, file = here("output/m_SEM/m_OT.rds"))
} else {
  m_OT <- readRDS(here("output/m_SEM/m_OT.rds"))
}
#
#
#
#| label: sem_tidy_OT
source(here("R", "func.R"))
#
# I am extracting here all the values for getting the total effects of each of the variables in the model. I am using the posterior values for each of the variables to get the total effects assuming that:
## total effect = direct effect + indirect effect + residual correlation
# In other words:
## total effect = 
# Extract the posteriors for the SEM model
post_sem_OT <- as_draws_df(m_OT) 
#
#### A) Get the direct paths per each variable
# 
# A.1) Detection
OT_coeff_mitodensity_det <- post_sem_OT$b_tD_mean_mitodensity
OT_coeff_potential_det <- post_sem_OT$b_tD_mean_potential
OT_coeff_dna_det <- post_sem_OT$b_tD_mean_dnadamage
OT_coeff_perox_det <- post_sem_OT$b_tD_mean_peroxidation
#
# A.2) DNA damage
OT_coeff_ros_dna <- post_sem_OT$b_meandnadamage_mean_ros
#
# A.3) Lipid peroxidation
OT_coeff_ros_perox <- post_sem_OT$b_meanperoxidation_mean_ros
#
# A.4) ROS
OT_coeff_mitodensity_ros <- post_sem_OT$b_meanros_mean_mitodensity
OT_coeff_potential_ros <- post_sem_OT$b_meanros_mean_potential
#
#
#### B) Get the indirect paths for each variable
#
# B.1) Detection
OT_undir_ros_det <- OT_coeff_ros_dna * OT_coeff_dna_det + OT_coeff_ros_perox * OT_coeff_perox_det
OT_undir_mitodensity_det <- OT_coeff_mitodensity_ros * OT_coeff_ros_dna * OT_coeff_dna_det + OT_coeff_mitodensity_ros * OT_coeff_ros_perox * OT_coeff_perox_det
OT_undir_potential_det <- OT_coeff_potential_ros * OT_coeff_ros_dna * OT_coeff_dna_det + OT_coeff_potential_ros * OT_coeff_ros_perox * OT_coeff_perox_det
#
# B.2) DNA damage
OT_undir_mitodensity_dna <- OT_coeff_mitodensity_ros * OT_coeff_ros_dna
OT_undir_potential_dna <- OT_coeff_potential_ros * OT_coeff_ros_dna
#
# B.3) Lipid peroxidation
OT_undir_mitodensity_perox <- OT_coeff_mitodensity_ros * OT_coeff_ros_perox
OT_undir_potential_perox <- OT_coeff_potential_ros * OT_coeff_ros_perox
#
#
#### C) Get the total effects for each variable
#
# C.1) Detection
OT_total_mitodensity_det <- OT_undir_mitodensity_det
OT_total_potential_det <- OT_undir_potential_det
OT_total_ros_det <- OT_undir_ros_det 
OT_total_dna_det <- OT_coeff_dna_det 
OT_total_perox_det <- OT_coeff_perox_det 
#
# C.2) DNA damage
OT_total_ros_dna <- OT_coeff_ros_dna 
OT_total_mitodensity_dna <- OT_undir_mitodensity_dna
OT_total_potential_dna <- OT_undir_potential_dna
#
# C.3) Lipid peroxidation
OT_total_ros_perox <- OT_coeff_ros_perox 
OT_total_mitodensity_perox <- OT_undir_mitodensity_perox
OT_total_potential_perox <- OT_undir_potential_perox
#
# C.4) ROS
OT_total_mitodensity_ros <- OT_coeff_mitodensity_ros
OT_total_potential_ros <- OT_coeff_potential_ros
#
#
# D) Create a df with the values for each variable
#
# D.1) Detection
det_results_semOT <- data.frame(
  reference_variable = rep("t_D", 5),
  predictor_modulator = c("mean_mitodensity",
                          "mean_potential",
                          "ROS",
                          "mean_dnadamage",
                          "mean_peroxidation"),
  direct_effects = I(list(OT_coeff_mitodensity_det,
                          OT_coeff_potential_det,
                          NA,
                          OT_coeff_dna_det,
                          OT_coeff_perox_det)),
  indirect_effects = I(list(OT_undir_mitodensity_det,
                          OT_undir_potential_det,
                          OT_undir_ros_det,
                          NA,
                          NA)),
  total_effects = I(list(OT_total_mitodensity_det,
                        OT_total_potential_det,
                        OT_total_ros_det,
                        OT_total_dna_det,
                        OT_total_perox_det))
  ) %>% 
  group_by(predictor_modulator) %>%
  summarize(
    # Summarizing each list column (direct_effects, indirect_effects, etc.)
    mean_direct_effects = format_dec(mean(unlist(direct_effects), na.rm = TRUE), 3),
    mean_indirect_effects = format_dec(mean(unlist(indirect_effects), na.rm = TRUE), 3),
    mean_total_effects = format_dec(mean(unlist(total_effects), na.rm = TRUE), 3),
    q5_direct = format_dec(quantile(unlist(direct_effects), 0.05, na.rm = TRUE), 3),
    q95_direct = format_dec(quantile(unlist(direct_effects), 0.95, na.rm = TRUE), 3),
    q5_indirect = format_dec(quantile(unlist(indirect_effects), 0.05, na.rm = TRUE), 3),
    q95_indirect = format_dec(quantile(unlist(indirect_effects), 0.95, na.rm = TRUE), 3),
    q5_total = format_dec(quantile(unlist(total_effects), 0.05, na.rm = TRUE), 3),
    q95_total = format_dec(quantile(unlist(total_effects), 0.95, na.rm = TRUE), 3)
  )  %>%
  mutate(source = "det_results_semOT")
#
# D.2) DNA damage
dna_results_semOT <- data.frame(
  reference_variable = rep("DNA damage", 3),
  predictor_modulator = c("mean_mitodensity",
                          "mean_potential",
                          "ROS"),
  direct_effects = I(list(NA,
                          NA,
                          OT_coeff_ros_dna)),
  indirect_effects = I(list(OT_undir_mitodensity_dna,
                          OT_undir_potential_dna,
                          NA)),
  total_effects = I(list(OT_total_mitodensity_dna,
                        OT_total_potential_dna,
                        OT_total_ros_dna))
  ) %>% 
  group_by(predictor_modulator) %>%
  summarize(
    # Summarizing each list column (direct_effects, indirect_effects, etc.)
    mean_direct_effects = format_dec(mean(unlist(direct_effects), na.rm = TRUE), 3),
    mean_indirect_effects = format_dec(mean(unlist(indirect_effects), na.rm = TRUE), 3),
    mean_total_effects = format_dec(mean(unlist(total_effects), na.rm = TRUE), 3),
    q5_direct = format_dec(quantile(unlist(direct_effects), 0.05, na.rm = TRUE), 3),
    q95_direct = format_dec(quantile(unlist(direct_effects), 0.95, na.rm = TRUE), 3),
    q5_indirect = format_dec(quantile(unlist(indirect_effects), 0.05, na.rm = TRUE), 3),
    q95_indirect = format_dec(quantile(unlist(indirect_effects), 0.95, na.rm = TRUE), 3),
    q5_total = format_dec(quantile(unlist(total_effects), 0.05, na.rm = TRUE), 3),
    q95_total = format_dec(quantile(unlist(total_effects), 0.95, na.rm = TRUE), 3)
  )  %>%
  mutate(source = "dna_results_semOT")
#
# D.3) Lipid peroxidation
perox_results_semOT <- data.frame(
  reference_variable = rep("Lipid peroxidation", 3),
  predictor_modulator = c("mean_mitodensity",
                          "mean_potential",
                          "ROS"),
  direct_effects = I(list(NA,
                          NA,
                          OT_coeff_ros_perox)),
  indirect_effects = I(list(OT_undir_mitodensity_perox,
                          OT_undir_potential_perox,
                          NA)),
  total_effects = I(list(OT_total_mitodensity_perox,
                        OT_total_potential_perox,
                        OT_total_ros_perox))
  ) %>% 
  group_by(predictor_modulator) %>%
  summarize(
    # Summarizing each list column (direct_effects, indirect_effects, etc.)
    mean_direct_effects = format_dec(mean(unlist(direct_effects), na.rm = TRUE), 3),
    mean_indirect_effects = format_dec(mean(unlist(indirect_effects), na.rm = TRUE), 3),
    mean_total_effects = format_dec(mean(unlist(total_effects), na.rm = TRUE), 3),
    q5_direct = format_dec(quantile(unlist(direct_effects), 0.05, na.rm = TRUE), 3),
    q95_direct = format_dec(quantile(unlist(direct_effects), 0.95, na.rm = TRUE), 3),
    q5_indirect = format_dec(quantile(unlist(indirect_effects), 0.05, na.rm = TRUE), 3),
    q95_indirect = format_dec(quantile(unlist(indirect_effects), 0.95, na.rm = TRUE), 3),
    q5_total = format_dec(quantile(unlist(total_effects), 0.05, na.rm = TRUE), 3),
    q95_total = format_dec(quantile(unlist(total_effects), 0.95, na.rm = TRUE), 3)
  ) %>%
  mutate(source = "perox_results_semOT")
#
# D.4) ROS
ros_results_semOT <- data.frame(
  reference_variable = rep("ROS", 2),
  predictor_modulator = c("mean_mitodensity",
                          "mean_potential"),
  direct_effects = I(list(OT_coeff_mitodensity_ros,
                          OT_coeff_potential_ros)),
  indirect_effects = I(list(NA,
                          NA)),
  total_effects = I(list(OT_total_mitodensity_ros,
                        OT_total_potential_ros))
  ) %>% 
  group_by(predictor_modulator) %>%
  summarize(
    # Summarizing each list column (direct_effects, indirect_effects, etc.)
    mean_direct_effects = format_dec(mean(unlist(direct_effects), na.rm = TRUE), 3),
    mean_indirect_effects = format_dec(mean(unlist(indirect_effects), na.rm = TRUE), 3),
    mean_total_effects = format_dec(mean(unlist(total_effects), na.rm = TRUE), 3),
    q5_direct = format_dec(quantile(unlist(direct_effects), 0.05, na.rm = TRUE), 3),
    q95_direct = format_dec(quantile(unlist(direct_effects), 0.95, na.rm = TRUE), 3),
    q5_indirect = format_dec(quantile(unlist(indirect_effects), 0.05, na.rm = TRUE), 3),
    q95_indirect = format_dec(quantile(unlist(indirect_effects), 0.95, na.rm = TRUE), 3),
    q5_total = format_dec(quantile(unlist(total_effects), 0.05, na.rm = TRUE), 3),
    q95_total = format_dec(quantile(unlist(total_effects), 0.95, na.rm = TRUE), 3)
  )  %>%
  mutate(source = "ros_results_semOT")
#
# E) Merge everything into a single df
#
sem_results_OT <- bind_rows(det_results_semOT, dna_results_semOT, perox_results_semOT, ros_results_semOT)
#
#
#
#
#
#| label: fig-sem_results_OT
#| fig-cap: "Structural Equation Models for OT/Visual stimulus"
#| fig-name: "fig-sem_results_OT"
#
# A) Getting all the direct coefficients
OT_density_det <- paste0(format_dec(mean(OT_coeff_mitodensity_det), 3),
                  " [", format_dec(quantile(OT_coeff_mitodensity_det, 0.05), 3),
                  ", ", format_dec(quantile(OT_coeff_mitodensity_det, 0.95), 3), "]")
OT_potential_det <- paste0(format_dec(mean(OT_coeff_potential_det), 3),
                  " [", format_dec(quantile(OT_coeff_potential_det, 0.05), 3),
                  ", ", format_dec(quantile(OT_coeff_potential_det, 0.95), 3), "]")
OT_dna_det <- paste0(format_dec(mean(OT_coeff_dna_det), 3),
                  " [", format_dec(quantile(OT_coeff_dna_det, 0.05), 3),
                  ", ", format_dec(quantile(OT_coeff_dna_det, 0.95), 3), "]")
OT_perox_det <- paste0(format_dec(mean(OT_coeff_perox_det), 3),
                  " [", format_dec(quantile(OT_coeff_perox_det, 0.05), 3),
                  ", ", format_dec(quantile(OT_coeff_perox_det, 0.95), 3), "]")
#
OT_ros_dna <- paste0(format_dec(mean(OT_coeff_ros_dna), 3),
                  " [", format_dec(quantile(OT_coeff_ros_dna, 0.05), 3),
                  ", ", format_dec(quantile(OT_coeff_ros_dna, 0.95), 3), "]")
#
OT_ros_perox <- paste0(format_dec(mean(OT_coeff_ros_perox), 3),
                  " [", format_dec(quantile(OT_coeff_ros_perox, 0.05), 3),
                  ", ", format_dec(quantile(OT_coeff_ros_perox, 0.95), 3), "]")
#
OT_density_ros <- paste0(format_dec(mean(OT_coeff_mitodensity_ros), 3),
                  " [", format_dec(quantile(OT_coeff_mitodensity_ros, 0.05), 3),
                  ", ", format_dec(quantile(OT_coeff_mitodensity_ros, 0.95), 3), "]")
OT_potential_ros <- paste0(format_dec(mean(OT_coeff_potential_ros), 3),
                  " [", format_dec(quantile(OT_coeff_potential_ros, 0.05), 3),
                  ", ", format_dec(quantile(OT_coeff_potential_ros, 0.95), 3), "]")
#
#
imgOT <- readPNG(here("Others", "SEM.png"))
plot_SEM_OT <- rasterGrob(imgOT, interpolate = TRUE)
#
fig_SEM_OT <- ggdraw(plot_SEM_OT) +
  draw_image("./Others/OT.png", x = 0.86, y = 0.67, width = 0.16, height = 0.3) +
  annotate("text", x = 0.52, y = 0.075, label = OT_density_det, hjust = 1, vjust = 1, size = 2.8, family = "Times") +
  annotate("text", x = 0.52, y = 0.94, label = OT_potential_det, hjust = 1, vjust = 1, size = 2.8, family = "Times") +
  annotate("text", x = 0.83, y = 0.715, label = OT_dna_det, hjust = 1, vjust = 1, size = 2.8, family = "Times") +
  annotate("text", x = 0.83, y = 0.335, label = OT_perox_det, hjust = 1, vjust = 1, size = 2.8, family = "Times") +
  annotate("text", x = 0.525, y = 0.805, label = OT_ros_dna, hjust = 1, vjust = 1, size = 2.8, family = "Times") +
  annotate("text", x = 0.515, y = 0.215, label = OT_ros_perox, hjust = 1, vjust = 1, size = 2.8, family = "Times") +
  annotate("text", x = 0.285, y = 0.2, label = OT_density_ros, hjust = 1, vjust = 1, size = 2.8, family = "Times") +
  annotate("text", x = 0.285, y = 0.82, label = OT_potential_ros, hjust = 1, vjust = 1, size = 2.8, family = "Times")
ggsave(here("./output/figures/text/SEM_OT.png"), plot = fig_SEM_OT, width = 21, height = 7.5, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/text/SEM_OT.png")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#
#| label: table_bayesR2
#| tbl-cap: "BayesR2 values of the definite models for both regions and stimuli."
#
data_bayes <- data.frame(
  Model = character(0),
  Mean = numeric(0),
  Error = numeric(0),
  Q2_5 = numeric(0),
  Q97_5 = numeric(0)
)
models <- c("mean_mitodensity_def_OB",
            "mean_potential_def_OB",
            "mean_ros_def_OB",
            "mean_dnadamage_def_OB",
            "mean_peroxidation_def_OB",
            "t_D_def_Chemical",
            "mean_mitodensity_def_OT",
            "mean_potential_def_OT",
            "mean_ros_def_OT",
            "mean_dnadamage_def_OT",
            "mean_peroxidation_def_OT",
            "t_D_def_Visual")
#
for (m in models){
  mod <- readRDS(here("output/models/", paste0(m, ".rds")))
  bayes <- bayes_R2(mod)
  data_bayes <- rbind(data_bayes, data.frame(
    Model = m,
    Mean = format_dec(bayes[1], 3),
    Error = format_dec(bayes[2], 3),
    Q2_5 = format_dec(bayes[3], 3),
    Q97_5 = format_dec(bayes[4], 3)
  ))
}
#
bayes_table_df <- data_bayes %>%
  mutate(Region = gsub(".*_", "", Model)    # Extract everything after the last "_"
  ) %>%
  mutate(Model = factor(Model,
                        levels = c("mean_mitodensity_def_OB",
                                  "mean_potential_def_OB",
                                  "mean_ros_def_OB",
                                  "mean_dnadamage_def_OB",
                                  "mean_peroxidation_def_OB",
                                  "t_D_def_Chemical",
                                  "mean_mitodensity_def_OT",
                                  "mean_potential_def_OT",
                                  "mean_ros_def_OT",
                                  "mean_dnadamage_def_OT",
                                  "mean_peroxidation_def_OT",
                                  "t_D_def_Visual"),
                        labels = c("m_def_mean_mitodensity_OB" = "Mit density",
                                  "m_def_mean_potential_OB" = "Mit potential",
                                  "m_def_mean_ros_OB" = "ROS",
                                  "m_def_mean_dnadamage_OB" = "DNA damage",
                                  "m_def_mean_peroxidation_OB" = "Peroxidation",
                                  "m_def_t_D_Chemical" = "Detection lat",
                                  "m_def_mean_mitodensity_OT" = "Mit density",
                                  "m_def_mean_potential_OT" = "Mit potential",
                                  "m_def_mean_ros_OT" = "ROS",
                                  "m_def_mean_dnadamage_OT" = "DNA damage",
                                  "m_def_mean_peroxidation_OT" = "Peroxidation",
                                  "m_def_t_D_Visual" = "Detection lat")),
      Region = factor(Region, levels = c("OB", "Chemical", "OT", "Visual"),
                      labels = c("OB" = "Olfactory bulbs",
                                "Chemical" = "Chemical",
                                "OT" = "Optic tecta",
                                "Visual" = "Visual"))) %>%
  dplyr::select(Region, Model, Mean, Error, Q2_5, Q97_5) %>%
  arrange(Region, Model)
#
#
# Create the table
#
## Table format
set_flextable_defaults(
 font.family = "Times New Roman",
 font.size = 10)
#
bayes_table <- flextable(bayes_table_df) %>%
  set_header_labels(
    Region = "Region/Stimulus",
    Mean = "Mean",
    Error = "Error",
    Q2_5 = "2.5%",
    Q97_5 = "97.5%") %>%
  align(align = "center", j = c(3:5), part = "body") %>%
  align(align = "center", j = c(1:5), part = "header") %>%
  flextable::compose(i = c(2:5,8:11), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
  flextable::hline(i = 6, part = "body") %>% 
  autofit()
#
bayes_table
#
#
#
cat("\\newpage")
#
#
#
#
#
#| label: tbl-contrasts
#| tbl-cap: "Contrasts between prenatal conditions for mitochondrial physiology and learning."
#
source(here("R", "func.R"))
#
# A) Organise df 
var <- c("MD", "MP", "ROS", "DNA", "LP", "DET")
region <- c("OB", "OT")
#
data_table <- data.frame()
for(v in var){
  for(r in region){
    x <- paste0(v, "_", r)
    df <- get(x)
    Temperature <- format_dec(mean(c(df$CORT_Hot, df$Control_Hot)) - mean(c(df$CORT_Cold, df$Control_Cold)), 3)
    pMCMC_temp <- format_p(pmcmc(c(df$CORT_Hot, df$Control_Hot) - c(df$CORT_Cold, df$Control_Cold)), 3, equal = FALSE)
    CORT <- format_dec(mean(c(df$Control_Hot, df$Control_Cold)) - mean(c(df$CORT_Hot, df$CORT_Cold)), 3)
    pMCMC_cort <- format_p(pmcmc(c(df$Control_Hot, df$Control_Cold) - c(df$CORT_Hot, df$CORT_Cold)), 3, equal = FALSE)
    Interaction <- format_dec((mean(df$Control_Hot) - mean(df$CORT_Hot)) - (mean(df$Control_Cold) - mean(df$CORT_Cold)), 3)
    pMCMC_int <- format_p(pmcmc((df$Control_Hot - df$CORT_Hot) - (df$Control_Cold - df$CORT_Cold)), 3, equal = FALSE)
    data_temp <- data.frame(Variable = x,
                          Temperature = as.numeric(Temperature),
                          pMCMC_temp = pMCMC_temp,
                          CORT = as.numeric(CORT),
                          pMCMC_cort = pMCMC_cort,
                          Interaction = as.numeric(Interaction),
                          pMCMC_int = pMCMC_int)
    data_table <- dplyr::bind_rows(data_table, data_temp)
  }
}
# Modify the df
data_table_final <- data_table %>%
  pivot_longer(cols = c(Temperature, CORT, Interaction), 
               names_to = "Predictor", 
               values_to = "Contrast") %>%
  mutate(
    # Extract the pMCMC values from the corresponding columns
    `pMCMC contrast` = case_when(
      Predictor == "Temperature" ~ pMCMC_temp,
      Predictor == "CORT" ~ pMCMC_cort,
      Predictor == "Interaction" ~ pMCMC_int
    )
  ) %>%
  separate(col = Variable, into = c("Variable", "Region"), sep = "_") %>%
  mutate(
    Variable = case_when(
      Variable == "MD" ~ "Mit density",
      Variable == "MP" ~ "Mit potential",
      Variable == "ROS" ~ "ROS",
      Variable == "DNA" ~ "DNA damage",
      Variable == "LP" ~ "Lipid peroxidation",
      Variable == "DET" ~ "Detection latency",
      TRUE ~ Variable
      )
    ) %>%
  mutate(Variable = factor(Variable, levels = c("Mit density", "Mit potential", "ROS", "DNA damage", "Lipid peroxidation", "Detection latency"))) %>%
  mutate(
    Region = case_when(
      Region == "OB" ~ "Olfactory bulbs / Chemical",
      Region == "OT" ~ "Optic tecta / Visual",
      TRUE ~ Region
    )
  ) %>%
  dplyr::select(Region, Variable, Predictor, Contrast, `pMCMC contrast`) %>%
  arrange(Region, Variable)
#
# C) Make the contrasts table:
#
set_flextable_defaults(
 font.family = "Times New Roman",
 font.size = 10)
#
contrast_table <- flextable(data_table_final) %>%
  align(align = "center", j = c(4,5), part = "body") %>%
  align(align = "center", j = c(1:4), part = "header") %>%
  bold(~`pMCMC contrast` < 0.05, j = c("pMCMC contrast", "Contrast", "Predictor"),
       bold = TRUE) %>%  # Bold when PMCMC is "<0.05"
  bold(~`pMCMC contrast` <0.001, j = c("pMCMC contrast", "Contrast", "Predictor")) %>%  # Bold when PMCMC is "<0.001"
  flextable::compose(i = c(2:18,20:36), j = 1, value = as_paragraph(""), part = "body") %>%
  flextable::compose(i = c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30,32,33,35,36),
                    j = 2, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the second column
  flextable::hline(i = c(3,6,9,12,15,21,24,27,30,33,36), j = c(2:5), part = "body") %>% 
  flextable::hline(i = 18, part = "body") %>%
  autofit()
#
contrast_table
#
#
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#| label: results_OB_table
#| tbl-cap: "Results of the models testing for Olfactory Bulbs."
#| tbl-name: "results_OB"
#| tbl-label: "results_OB"
source(here("R", "func.R"))
# 
# A) Refining the df summarizing the posteriors for OB/Chemical stimulus (post_OB)
post_OB_refined <- refine_post(post_OB) %>%
  arrange(Variable, Predictors)
#
# B) Create table
## Table format
set_flextable_defaults(
 font.family = "Times New Roman",
 font.size = 10)
#
OB_table <- flextable(post_OB_refined) %>%
  align(align = "center", j = c(3:5), part = "body") %>%
  align(align = "center", j = c(1:5), part = "header") %>%
  bold(~`PMCMC` < 0.05, j = c("PMCMC", "Estimate Mean", "95% CI", "Predictors"),
       bold = TRUE) %>%  # Bold when PMCMC is "<0.05"
  bold(~`PMCMC` <0.001, j = c("PMCMC", "Estimate Mean", "95% CI", "Predictors")) %>%  # Bold when PMCMC is "<0.001"
  flextable::compose(i = c(2:4,6:8,10:12,14:18,20:23,25:27), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
  autofit()
#
OB_table
#
#
#
cat("\\newpage")
#
#
#
#
#| label: results_OT_table
#| tbl-cap: "Results of the models testing for Olfactory Bulbs."
#| tbl-name: "results_OT"
#| tbl-label: "results_OT"
source(here("R", "func.R"))
# 
# A) Refining the df summarizing the posteriors for OT/Visual stimulus (post_OT)
post_OT_refined <- refine_post(post_OT) %>%
  arrange(Variable, Predictors)
#
# B) Create table
## Table format
set_flextable_defaults(
 font.family = "Times New Roman",
 font.size = 10)
#
OT_table <- flextable(post_OT_refined) %>%
  align(align = "center", j = c(3:5), part = "body") %>%
  align(align = "center", j = c(1:5), part = "header") %>%
  bold(~`PMCMC` < 0.05, j = c("PMCMC", "Estimate Mean", "95% CI", "Predictors"),
       bold = TRUE) %>%  # Bold when PMCMC is "<0.05"
  bold(~`PMCMC` <0.001, j = c("PMCMC", "Estimate Mean", "95% CI", "Predictors")) %>%  # Bold when PMCMC is "<0.001"
  flextable::compose(i = c(2:4,6:8,10:12,14:17,19:22,24:26), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
  autofit()
#
OT_table
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#| label: table_sem_results_OB
#| fig-cap: ""
#
# Table created from the df sem_results_OB (see above)
source(here("R", "func.R"))
#
# Modify database
table_semOB_df <- sem_results_OB %>%
  mutate(`Direct effects` = paste0(mean_direct_effects, " [", q5_direct, ", ", q95_direct, "]"),
         `Indirect effects` = paste0(mean_indirect_effects, " [", q5_indirect, ", ", q95_indirect, "]"),
         `Total effects` = paste0(mean_total_effects, " [", q5_total, ", ", q95_total, "]")) %>%
  dplyr::select(source, predictor_modulator, `Direct effects`, `Indirect effects`, `Total effects`) %>%
  mutate(
      predictor_modulator = case_when(
        predictor_modulator == "mean_mitodensity" ~ "Mitochondrial density",
        predictor_modulator == "mean_potential" ~ "Metabolic capacity",
        predictor_modulator == "mean_dnadamage" ~ "DNA damage",
        predictor_modulator == "mean_peroxidation" ~ "Lipid peroxidation",
        predictor_modulator == "ROS" ~ "ROS",
        TRUE ~ predictor_modulator
      ),
      source = case_when(
        source == "det_results_semOB" ~ "Detection latency",
        source == "dna_results_semOB" ~ "DNA damage",
        source == "perox_results_semOB" ~ "Lipid peroxidation",
        source == "ros_results_semOB" ~ "ROS",
        TRUE ~ source
      ),
      `Direct effects` = case_when(`Direct effects` == "NaN [NA, NA]" ~ " - ", TRUE ~ `Direct effects`),
      `Indirect effects` = case_when(`Indirect effects` == "NaN [NA, NA]" ~ " - ", TRUE ~ `Indirect effects`),
      `Total effects` = case_when(`Total effects` == "NaN [NA, NA]" ~ " - ", TRUE ~ `Total effects`)) %>%
  rename(Predictor = predictor_modulator, Response = source) %>%
  mutate(Response = factor(Response, levels = c("Detection latency", "DNA damage", "Lipid peroxidation", "ROS")),
         Predictor = factor(Predictor, levels = c("Mitochondrial density",
                                                "Metabolic capacity",
                                                "ROS",
                                                "DNA damage",
                                                "Lipid peroxidation"))) %>%
  arrange(Response, Predictor)
#
# Table format
set_flextable_defaults(
 font.family = "Times New Roman",
 font.size = 10)
#
table_semOB <- flextable(table_semOB_df) %>%
  align(align = "center", part = "header") %>%
  flextable::compose(i = c(2:5,7,8,10,11,13), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
  flextable::hline(i = c(5,8,11), part = "body") %>% 
  autofit()
#
table_semOB
#
#
#
cat("\\newpage")
#
#
#
#
#
#| label: table_sem_results_OT
#| fig-cap: ""
#
# Table created from the df sem_results_OT (see above)
source(here("R", "func.R"))
#
# Modify database
table_semOT_df <- sem_results_OT %>%
  mutate(`Direct effects` = paste0(mean_direct_effects, " [", q5_direct, ", ", q95_direct, "]"),
         `Indirect effects` = paste0(mean_indirect_effects, " [", q5_indirect, ", ", q95_indirect, "]"),
         `Total effects` = paste0(mean_total_effects, " [", q5_total, ", ", q95_total, "]")) %>%
  dplyr::select(source, predictor_modulator, `Direct effects`, `Indirect effects`, `Total effects`) %>%
  mutate(
      predictor_modulator = case_when(
        predictor_modulator == "mean_mitodensity" ~ "Mitochondrial density",
        predictor_modulator == "mean_potential" ~ "Metabolic capacity",
        predictor_modulator == "mean_dnadamage" ~ "DNA damage",
        predictor_modulator == "mean_peroxidation" ~ "Lipid peroxidation",
        predictor_modulator == "ROS" ~ "ROS",
        TRUE ~ predictor_modulator
      ),
      source = case_when(
        source == "det_results_semOT" ~ "Detection latency",
        source == "dna_results_semOT" ~ "DNA damage",
        source == "perox_results_semOT" ~ "Lipid peroxidation",
        source == "ros_results_semOT" ~ "ROS",
        TRUE ~ source
      ),
      `Direct effects` = case_when(`Direct effects` == "NaN [NA, NA]" ~ " - ", TRUE ~ `Direct effects`),
      `Indirect effects` = case_when(`Indirect effects` == "NaN [NA, NA]" ~ " - ", TRUE ~ `Indirect effects`),
      `Total effects` = case_when(`Total effects` == "NaN [NA, NA]" ~ " - ", TRUE ~ `Total effects`)) %>%
  rename(Predictor = predictor_modulator, Response = source) %>%
  mutate(Response = factor(Response, levels = c("Detection latency", "DNA damage", "Lipid peroxidation", "ROS")),
         Predictor = factor(Predictor, levels = c("Mitochondrial density",
                                                "Metabolic capacity",
                                                "ROS",
                                                "DNA damage",
                                                "Lipid peroxidation"))) %>%
  arrange(Response, Predictor)
#
# Table format
set_flextable_defaults(
 font.family = "Times New Roman",
 font.size = 10)
#
table_semOT <- flextable(table_semOT_df) %>%
  align(align = "center", part = "header") %>%
  flextable::compose(i = c(2:5,7,8,10,11,13), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
  flextable::hline(i = c(5,8,11), part = "body") %>% 
  autofit()
#
table_semOT
#
#
#
#
#
#
#| label: models_preliminary
# Fitting intial models to see if sex and age are relevant for our models
source(here("R", "func.R"))
#
#
# Run models mitochondrial physiology (each region separately)
#
var_m <- c("mean_mitodensity", "mean_potential", "mean_ros", "mean_dnadamage", "mean_peroxidation")
regions <- c("OB", "OT")
formula_list_ <- list()
for (p in var_m){
  formula_list_[[p]] <- paste0(p, "~ cort*temp + age_euthanasia + sex + (1|clutch)")
  for (h in regions){
    if (h == "OB"){
      df <- clean_df %>% filter(region == "OB")
      l <- "OB"
    } else {
      df <- clean_df %>% filter(region == "OT")
      l <- "OT"
    }
  
  pmodel_name <- paste0("m_prel_", p, "_", h)
  assign(pmodel_name, fit_m(df = df,
                             cat = "prel",
                             var = p,
                             formula = formula_list_[[p]],
                             fam = gaussian(),
                             label = l,
                             refit = FALSE),
          envir = .GlobalEnv)  # Assign to the global environment
  }
}
#
#
# Run model behaviour (stimuli separated)
beh_df <- clean_df
stimulus <- c("Chemical", "Visual")
formula_t_D <- t_D ~ motivation + motivation:cort + cort*temp + age_trial + sex + prey +(1|clutch) + (1|lizard_id) 
for (k in stimulus){
  df <- beh_df %>% filter(stimulus == k)
  pmodel_name <- paste0("m_prel_t_D_", k)
  assign(pmodel_name, fit_m(df = df,
                            cat = "prel",
                            var = "t_D",
                            formula = formula_t_D,
                            fam = gaussian(),
                            label = k,
                            refit = FALSE),
        envir = .GlobalEnv)  # Assign to the global environment
}
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for Mitochondrial Density in OB"
#| label: results_preliminary_mitdensity
#
sum_m_mitdensity_OB_prel <- m_prel_mean_mitodensity_OB %>%
  dplyr::select(-starts_with("r_"), -sigma, -Intercept, -lprior, -lp__) %>%  # Removes random effect terms
  summarise_draws()  %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_m_mitdensity_OB_prel)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for Mitochondrial Density in OT"
#| label: results_preliminary_mitdensity_OT
#
sum_m_mitdensity_OT_prel <- m_prel_mean_mitodensity_OT %>%
  dplyr::select(-starts_with("r_"), -sigma, -Intercept, -lprior, -lp__) %>%  # Removes random effect terms
  summarise_draws() %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_m_mitdensity_OT_prel)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for metabolic capacity in OB"
#| label: results_preliminary_potential_OB
#
sum_m_potential_OB_prel <- m_prel_mean_potential_OB %>%
  dplyr::select(-starts_with("r_"), -sigma, -Intercept, -lprior, -lp__) %>%  # Removes random effect terms
  summarise_draws() %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_m_potential_OB_prel)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for Metabolic capacity in OT"
#| label: results_preliminary_potential_OT
#
sum_m_potential_OT_prel <- m_prel_mean_potential_OT %>%
  dplyr::select(-starts_with("r_"), -sigma, -Intercept, -lprior, -lp__) %>%  # Removes random effect terms
  summarise_draws() %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_m_potential_OT_prel)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for ROS Production in OB"
#| label: results_preliminary_ros_OB
#
sum_m_ros_OB_prel <- m_prel_mean_ros_OB %>%
  dplyr::select(-starts_with("r_"), -sigma, -Intercept, -lprior, -lp__) %>%  # Removes random effect terms
  summarise_draws() %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_m_ros_OB_prel)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for ROS Production in OT"
#| label: results_preliminary_ros_OT
#
sum_m_ros_OT_prel <- m_prel_mean_ros_OT %>%
  dplyr::select(-starts_with("r_"), -sigma, -Intercept, -lprior, -lp__) %>%  # Removes random effect terms
  summarise_draws() %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_m_ros_OT_prel)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for DNA Damage in OB"
#| label: results_preliminary_dnadamage_OB
#
sum_m_dnadamage_OB_prel <- m_prel_mean_dnadamage_OB %>%
  dplyr::select(-starts_with("r_"), -sigma, -Intercept, -lprior, -lp__) %>%  # Removes random effect terms
  summarise_draws() %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_m_dnadamage_OB_prel)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for DNA Damage in OT"
#| label: results_preliminary_dnadamage_OT
#
sum_m_dnadamage_OT_prel <- m_prel_mean_dnadamage_OT %>%
  dplyr::select(-starts_with("r_"), -sigma, -Intercept, -lprior, -lp__) %>%  # Removes random effect terms
  summarise_draws() %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_m_dnadamage_OT_prel)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for Lipid Peroxidation in OB"
#| label: results_preliminary_peroxidation_OB
#
sum_m_peroxidation_OB_prel <- m_prel_mean_peroxidation_OB %>%
  dplyr::select(-starts_with("r_"), -sigma, -Intercept, -lprior, -lp__) %>%  # Removes random effect terms
  summarise_draws() %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_m_peroxidation_OB_prel)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for Lipid Peroxidation in OT"
#| label: results_preliminary_peroxidation_OT
#
sum_m_peroxidation_OT_prel <- m_prel_mean_peroxidation_OT %>%
  dplyr::select(-starts_with("r_"), -sigma, -Intercept, -lprior, -lp__) %>%  # Removes random effect terms
  summarise_draws() %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_m_peroxidation_OT_prel)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for Detection Latency (t_D) of Chemical stimuli"
#| label: results_preliminary_t_D_Chemical
# 
sum_m_t_D_prel_Chem <- m_prel_t_D_Chemical %>%
  dplyr::select(-starts_with("r_"), -sigma, -Intercept, -lprior, -lp__) %>%  # Removes random effect terms
  summarise_draws() %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_m_t_D_prel_Chem)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#| tbl-cap: "Preliminary results of the models testing for Detection Latency (t_D) of Visual stimuli"
#| label: results_preliminary_t_D_Visual
#
sum_m_t_D_prel_Vis <- m_prel_t_D_Visual %>%
  dplyr::select(-starts_with("r_"), -sigma, -Intercept, -lprior, -lp__) %>%  # Removes random effect terms
  summarise_draws() %>%
  mutate(across(where(is.numeric), ~ as.numeric(format_dec(.x, 3))))
#
flextable(sum_m_t_D_prel_Vis)
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#
#| label: plotmod_mitdensity_OB
#| caption: "Posterior predictive checks for the model of Mitochondrial Density in Olfactory Bulbs."
#
mod <- readRDS(here("output/models/mean_mitodensity_def_OB.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX3.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX3.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#| label: plotmod_mitdensity_OT
#| caption: "Posterior predictive checks for the model of Mitochondrial Density in Optic Tecta."
#
mod <- readRDS(here("output/models/mean_mitodensity_def_OT.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX4.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX4.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#| label: plotmod_potential_OB
#| caption: "Posterior predictive checks for the model of Metabolic capacity in Olfactory Bulbs."
#
mod <- readRDS(here("output/models/mean_potential_def_OB.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX5.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX5.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#| label: plotmod_potential_OT
#| caption: "Posterior predictive checks for the model of Metabolic capacity in Optic Tecta."
#
mod <- readRDS(here("output/models/mean_potential_def_OT.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX6.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX6.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#| label: plotmod_ros_OB
#| caption: "Posterior predictive checks for the model of ROS Production in Olfactory Bulbs."
#
mod <- readRDS(here("output/models/mean_ros_def_OB.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX7.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX7.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#| label: plotmod_ros_OT
#| caption: "Posterior predictive checks for the model of ROS Production in Optic Tecta."
#
mod <- readRDS(here("output/models/mean_ros_def_OT.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX8.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX8.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#| label: plotmod_dnadamage_OB
#| caption: "Posterior predictive checks for the model of DNA Damage in Olfactory Bulbs."
#
mod <- readRDS(here("output/models/mean_dnadamage_def_OB.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX9.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX9.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#| label: plotmod_dnadamage_OT
#| caption: "Posterior predictive checks for the model of DNA Damage in Optic Tecta."
#
mod <- readRDS(here("output/models/mean_dnadamage_def_OT.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX10.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX10.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#| label: plotmod_peroxidation_OB
#| caption: "Posterior predictive checks for the model of Lipid Peroxidation in Olfactory Bulbs."
#
mod <- readRDS(here("output/models/mean_peroxidation_def_OB.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX11.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX11.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#| label: plotmod_peroxidation_OT
#| caption: "Posterior predictive checks for the model of Lipid Peroxidation in Optic Tecta."
#
mod <- readRDS(here("output/models/mean_peroxidation_def_OT.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX12.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX12.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#| label: plotmod_tD_Chemical
#| caption: "Posterior predictive checks for the model of Detection Latency (t_D) in Chemical trials."
#
mod <- readRDS(here("output/models/t_D_def_Chemical.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX1.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX1.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#| label: plotmod_tD_Visual
#| caption: "Posterior predictive checks for the model of Detection Latency (t_D) in Visual trials."
#
mod <- readRDS(here("output/models/t_D_def_Visual.rds"))
plot_mod <- plot(mod, plot = FALSE)
#
plot_mod_1 <- plot_mod[[1]]
plot_mod_2 <- plot_mod[[2]]
#
fig_mod <- plot_grid(plot_mod_1, plot_mod_2, ncol = 1)
ggsave("./output/figures/suppl/Figure_SX2.png", plot = fig_mod,
       width = 20, height = 28, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX2.png")
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#
#| label: data_mitochondrial_controlled
# 
mit_controlled <- mit_df %>%
    mutate(Id_flow = as.character(Id_flow)) %>%
    mutate(mean_dnadamage = log10(mean_dnadamage)) %>%
    mutate(mean_peroxidation = log10(mean_peroxidation)) %>%
    mutate(across(where(~ is.numeric(.)), ~.x / mean_mitodensity))
#
df_controlled <- merge(data_discr, mit_controlled, by = c("stimulus", "Id_flow")) %>%
  mutate(across(c(age_trial, age_euthanasia), as.numeric)) %>%
  mutate(across(where(~ is.double(.) & !is.integer(.)),
                ~ (. - mean(., na.rm = TRUE)) / (2 * sd(., na.rm = TRUE))))
#
#
#
#| label: models_mitochondrial_controlled
#
# Run models mitochondrial physiology from df_controlled (each region separately)
#
var_m_controlled <- c("mean_potential", "mean_ros", "mean_dnadamage", "mean_peroxidation")
regions <- c("OB", "OT")
for (p in var_m_controlled){
  for (h in regions){
    if (h == "OB"){
      df <- df_controlled %>% filter(region == "OB")
      l <- "OB"
      if (p %in% c("mean_mitodensity", "mean_potential", "mean_ros")){
        formula <- paste0(p, "~ cort*temp + (1|clutch)")
      } else if (p == "mean_dnadamage"){
        formula <- paste0(p, "~ cort*temp + age_euthanasia + sex + (1|clutch)")
      } else {
        formula <- paste0(p, "~ cort*temp + age_euthanasia + (1|clutch)")
      }
    } else {
      df <- df_controlled %>% filter(region == "OT")
      l <- "OT"
      if (p %in% c("mean_mitodensity", "mean_potential", "mean_ros")){
        formula <- paste0(p, "~ cort*temp + (1|clutch)")
      } else if (p == "mean_dnadamage"){
        formula <- paste0(p, "~ cort*temp + sex + (1|clutch)")
      } else {
        formula <- paste0(p, "~ cort*temp + age_euthanasia + (1|clutch)")
      }
    }

    pmodel_name <- paste0("m_controlled_", p, "_", h)
    assign(pmodel_name, fit_m(df = df,
                              cat = "def",
                              var = p,
                              formula = formula,
                              fam = gaussian(),
                              label = l,
                              refit = TRUE),
          envir = .GlobalEnv)  # Assign to the global environment
  }
} 
#
#
#
#| label: organise_posteriors
#
# Organising the posteriors of the previous models to fit the tables below
#
source(here("R", "func.R"))
#
# Databases for each region/stimulus
post_OB_cont <- data.frame()
post_OT_cont <- data.frame()
# Names posteriors:
region <- c("OB", "OT")
names_OB <- c("m_controlled_mean_potential_OB", "m_controlled_mean_ros_OB", "m_controlled_mean_dnadamage_OB", "m_controlled_mean_peroxidation_OB", "m_controlled_t_D_Chemical")
names_OT <- c("m_controlled_mean_potential_OT", "m_controlled_mean_ros_OT", "m_controlled_mean_dnadamage_OT", "m_controlled_mean_peroxidation_OT", "m_controlled_t_D_Visual")
#
# Organising the results
for (r in region) {
  model_select <- get(paste0("names_", r))   
  for (pos in model_select) {
    model <- get(pos)      # Get the model from the global environment
    post_result <- tidy_post(model)        # Apply tidy_post to each model
    # Add a new column to identify the region and model
    post_result$Region <- r
    post_result$Model <- pos
    # Append to the appropriate data frame
    if (r == "OB") {
      post_OB <- bind_rows(post_OB, post_result)
    } else {
      post_OT <- bind_rows(post_OT, post_result)
    }
  }
}
#
#
#
#| label: values_posteriors
#
# Extracting the posteriors for the models and the values of interest. Here, I am creating dfs for each of the variables with the values for all the prenatal conditions to make contrasts easier to write.
#
source(here("R", "func.R"))
#
# A) Olfactory bulbs/Chemical stimulus
#
MP_OB_con <- post_values(m_controlled_mean_potential_OB, "none")
ROS_OB_con <- post_values(m_controlled_mean_ros_OB, "none")
DNA_OB_con <- post_values(m_controlled_mean_dnadamage_OB, "none")
LP_OB_con <- post_values(m_controlled_mean_peroxidation_OB, "none")
DET_OB_con <- post_values(m_controlled_t_D_Chemical, "none")
#
# B) Optic tecta/Visual stimulus
MP_OT_con <- post_values(m_controlled_mean_potential_OT, "none")
ROS_OT_con <- post_values(m_controlled_mean_ros_OT, "none")
DNA_OT_con <- post_values(m_controlled_mean_dnadamage_OT, "sex")
LP_OT_con <- post_values(m_controlled_mean_peroxidation_OT, "none")
DET_OT_con <- post_values(m_controlled_t_D_Visual, "none")
#
#
#
#
#
#
#| label: table_bayesR2
#| tbl-cap: "BayesR2 values of the definite models for both regions and stimuli."
#
data_bayes <- data.frame(
  Model = character(0),
  Mean = numeric(0),
  Error = numeric(0),
  Q2_5 = numeric(0),
  Q97_5 = numeric(0)
)
models <- c("mean_mitodensity_def_OB",
            "mean_potential_def_OB",
            "mean_ros_def_OB",
            "mean_dnadamage_def_OB",
            "mean_peroxidation_def_OB",
            "t_D_def_Chemical",
            "mean_mitodensity_def_OT",
            "mean_potential_def_OT",
            "mean_ros_def_OT",
            "mean_dnadamage_def_OT",
            "mean_peroxidation_def_OT",
            "t_D_def_Visual")
#
for (m in models){
  mod <- readRDS(here("output/models/", paste0(m, ".rds")))
  bayes <- bayes_R2(mod)
  data_bayes <- rbind(data_bayes, data.frame(
    Model = m,
    Mean = format_dec(bayes[1], 3),
    Error = format_dec(bayes[2], 3),
    Q2_5 = format_dec(bayes[3], 3),
    Q97_5 = format_dec(bayes[4], 3)
  ))
}
#
bayes_table_df <- data_bayes %>%
  mutate(Region = gsub(".*_", "", Model)    # Extract everything after the last "_"
  ) %>%
  mutate(Model = factor(Model,
                        levels = c("mean_mitodensity_def_OB",
                                  "mean_potential_def_OB",
                                  "mean_ros_def_OB",
                                  "mean_dnadamage_def_OB",
                                  "mean_peroxidation_def_OB",
                                  "t_D_def_Chemical",
                                  "mean_mitodensity_def_OT",
                                  "mean_potential_def_OT",
                                  "mean_ros_def_OT",
                                  "mean_dnadamage_def_OT",
                                  "mean_peroxidation_def_OT",
                                  "t_D_def_Visual"),
                        labels = c("m_def_mean_mitodensity_OB" = "Mit density",
                                  "m_def_mean_potential_OB" = "Mit potential",
                                  "m_def_mean_ros_OB" = "ROS",
                                  "m_def_mean_dnadamage_OB" = "DNA damage",
                                  "m_def_mean_peroxidation_OB" = "Peroxidation",
                                  "m_def_t_D_Chemical" = "Detection lat",
                                  "m_def_mean_mitodensity_OT" = "Mit density",
                                  "m_def_mean_potential_OT" = "Mit potential",
                                  "m_def_mean_ros_OT" = "ROS",
                                  "m_def_mean_dnadamage_OT" = "DNA damage",
                                  "m_def_mean_peroxidation_OT" = "Peroxidation",
                                  "m_def_t_D_Visual" = "Detection lat")),
      Region = factor(Region, levels = c("OB", "Chemical", "OT", "Visual"),
                      labels = c("OB" = "Olfactory bulbs",
                                "Chemical" = "Chemical",
                                "OT" = "Optic tecta",
                                "Visual" = "Visual"))) %>%
  dplyr::select(Region, Model, Mean, Error, Q2_5, Q97_5) %>%
  arrange(Region, Model)
#
#
# Create the table
#
## Table format
set_flextable_defaults(
 font.family = "Times New Roman",
 font.size = 10)
#
bayes_table <- flextable(bayes_table_df) %>%
  set_header_labels(
    Region = "Region/Stimulus",
    Mean = "Mean",
    Error = "Error",
    Q2_5 = "2.5%",
    Q97_5 = "97.5%") %>%
  align(align = "center", j = c(3:5), part = "body") %>%
  align(align = "center", j = c(1:5), part = "header") %>%
  flextable::compose(i = c(2:5,8:11), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
  flextable::hline(i = 6, part = "body") %>% 
  autofit()
#
bayes_table
#
#
#
cat("\\newpage")
#
#
#
#
#
#| label: tbl-contrasts
#| tbl-cap: "Contrasts between prenatal conditions for mitochondrial physiology and learning."
#
source(here("R", "func.R"))
#
# A) Organise df 
var <- c("MD", "MP", "ROS", "DNA", "LP", "DET")
region <- c("OB", "OT")
#
data_table <- data.frame()
for(v in var){
  for(r in region){
    x <- paste0(v, "_", r)
    df <- get(x)
    Temperature <- format_dec(mean(c(df$CORT_Hot, df$Control_Hot)) - mean(c(df$CORT_Cold, df$Control_Cold)), 3)
    pMCMC_temp <- format_p(pmcmc(c(df$CORT_Hot, df$Control_Hot) - c(df$CORT_Cold, df$Control_Cold)), 3, equal = FALSE)
    CORT <- format_dec(mean(c(df$Control_Hot, df$Control_Cold)) - mean(c(df$CORT_Hot, df$CORT_Cold)), 3)
    pMCMC_cort <- format_p(pmcmc(c(df$Control_Hot, df$Control_Cold) - c(df$CORT_Hot, df$CORT_Cold)), 3, equal = FALSE)
    Interaction <- format_dec((mean(df$Control_Hot) - mean(df$CORT_Hot)) - (mean(df$Control_Cold) - mean(df$CORT_Cold)), 3)
    pMCMC_int <- format_p(pmcmc((df$Control_Hot - df$CORT_Hot) - (df$Control_Cold - df$CORT_Cold)), 3, equal = FALSE)
    data_temp <- data.frame(Variable = x,
                          Temperature = as.numeric(Temperature),
                          pMCMC_temp = pMCMC_temp,
                          CORT = as.numeric(CORT),
                          pMCMC_cort = pMCMC_cort,
                          Interaction = as.numeric(Interaction),
                          pMCMC_int = pMCMC_int)
    data_table <- dplyr::bind_rows(data_table, data_temp)
  }
}
# Modify the df
data_table_final <- data_table %>%
  pivot_longer(cols = c(Temperature, CORT, Interaction), 
               names_to = "Predictor", 
               values_to = "Contrast") %>%
  mutate(
    # Extract the pMCMC values from the corresponding columns
    `pMCMC contrast` = case_when(
      Predictor == "Temperature" ~ pMCMC_temp,
      Predictor == "CORT" ~ pMCMC_cort,
      Predictor == "Interaction" ~ pMCMC_int
    )
  ) %>%
  separate(col = Variable, into = c("Variable", "Region"), sep = "_") %>%
  mutate(
    Variable = case_when(
      Variable == "MD" ~ "Mit density",
      Variable == "MP" ~ "Mit potential",
      Variable == "ROS" ~ "ROS",
      Variable == "DNA" ~ "DNA damage",
      Variable == "LP" ~ "Lipid peroxidation",
      Variable == "DET" ~ "Detection latency",
      TRUE ~ Variable
      )
    ) %>%
  mutate(Variable = factor(Variable, levels = c("Mit density", "Mit potential", "ROS", "DNA damage", "Lipid peroxidation", "Detection latency"))) %>%
  mutate(
    Region = case_when(
      Region == "OB" ~ "Olfactory bulbs / Chemical",
      Region == "OT" ~ "Optic tecta / Visual",
      TRUE ~ Region
    )
  ) %>%
  dplyr::select(Region, Variable, Predictor, Contrast, `pMCMC contrast`) %>%
  arrange(Region, Variable)
#
# C) Make the contrasts table:
#
set_flextable_defaults(
 font.family = "Times New Roman",
 font.size = 10)
#
contrast_table <- flextable(data_table_final) %>%
  align(align = "center", j = c(4,5), part = "body") %>%
  align(align = "center", j = c(1:4), part = "header") %>%
  bold(~`pMCMC contrast` < 0.05, j = c("pMCMC contrast", "Contrast", "Predictor"),
       bold = TRUE) %>%  # Bold when PMCMC is "<0.05"
  bold(~`pMCMC contrast` <0.001, j = c("pMCMC contrast", "Contrast", "Predictor")) %>%  # Bold when PMCMC is "<0.001"
  flextable::compose(i = c(2:18,20:36), j = 1, value = as_paragraph(""), part = "body") %>%
  flextable::compose(i = c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30,32,33,35,36),
                    j = 2, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the second column
  flextable::hline(i = c(3,6,9,12,15,21,24,27,30,33,36), j = c(2:5), part = "body") %>% 
  flextable::hline(i = 18, part = "body") %>%
  autofit()
#
contrast_table
#
#
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#| label: results_OB_table
#| tbl-cap: "Results of the models testing for Olfactory Bulbs."
#| tbl-name: "results_OB"
#| tbl-label: "results_OB"
source(here("R", "func.R"))
# 
# A) Refining the df summarizing the posteriors for OB/Chemical stimulus (post_OB)
post_OB_refined <- refine_post(post_OB) %>%
  arrange(Variable, Predictors)
#
# B) Create table
## Table format
set_flextable_defaults(
 font.family = "Times New Roman",
 font.size = 10)
#
OB_table <- flextable(post_OB_refined) %>%
  align(align = "center", j = c(3:5), part = "body") %>%
  align(align = "center", j = c(1:5), part = "header") %>%
  bold(~`PMCMC` < 0.05, j = c("PMCMC", "Estimate Mean", "95% CI", "Predictors"),
       bold = TRUE) %>%  # Bold when PMCMC is "<0.05"
  bold(~`PMCMC` <0.001, j = c("PMCMC", "Estimate Mean", "95% CI", "Predictors")) %>%  # Bold when PMCMC is "<0.001"
  flextable::compose(i = c(2:4,6:8,10:12,14:18,20:23,25:27), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
  autofit()
#
OB_table
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#| label: results_OT_table
#| tbl-cap: "Results of the models testing for Olfactory Bulbs."
#| tbl-name: "results_OT"
#| tbl-label: "results_OT"
source(here("R", "func.R"))
# 
# A) Refining the df summarizing the posteriors for OT/Visual stimulus (post_OT)
post_OT_refined <- refine_post(post_OT) %>%
  arrange(Variable, Predictors)
#
# B) Create table
## Table format
set_flextable_defaults(
 font.family = "Times New Roman",
 font.size = 10)
#
OT_table <- flextable(post_OT_refined) %>%
  align(align = "center", j = c(3:5), part = "body") %>%
  align(align = "center", j = c(1:5), part = "header") %>%
  bold(~`PMCMC` < 0.05, j = c("PMCMC", "Estimate Mean", "95% CI", "Predictors"),
       bold = TRUE) %>%  # Bold when PMCMC is "<0.05"
  bold(~`PMCMC` <0.001, j = c("PMCMC", "Estimate Mean", "95% CI", "Predictors")) %>%  # Bold when PMCMC is "<0.001"
  flextable::compose(i = c(2:4,6:8,10:12,14:17,19:22,24:26), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
  autofit()
#
OT_table
#
#
#
#
#
#
#
#
#
#
#
#
#| label: fig-results_energy
#| fig-cap: "Estimates of mitochondrial density (A, C) and metabolic capacity (B, D) in the olfactory bulbs (A, B) and optic tecta (C, D) of L. delicata hatchlings as a function of the different prenatal conditions. Black dots indicate the posterior mean, and the bars represent the SD of the estimates. The y-axis represents the posterior estimates of the variable of interest, and the x-axis represents the different prenatal conditions. Lines with asterisks represent significant differences between groups based on pMCMC values (pMCMC < 0.05), no lines indicate no significant differences."
#| fig-name: "fig-results_energy"
#
source(here("R", "func.R"))
#
# A) Plotting the results for OB/Chemical stimulus using the df from before
plot_density_OB <- plotting(MD_OB, "Mit density")
plot_potential_OB <- plotting(MP_OB, "Mit potential")
#
fig_OB_energy <- plot_grid(plot_density_OB, NULL, plot_potential_OB, NULL,
                    nrow = 1, rel_widths = c(0.9, 0.1, 0.9, 0.3))
#
# B) Plotting the results for OT/Visual stimulus
plot_density_OT <- plotting(MD_OT, "Mit density")
plot_potential_OT <- plotting(MP_OT, "Mit potential")
#
fig_OT_energy <- plot_grid(plot_density_OT, NULL, plot_potential_OT, NULL,
                    nrow = 1, rel_widths = c(0.9, 0.1, 0.9, 0.3))
#
# Getting the legend (get_legend() does not work if the legend is on the top or bottom)
plot_legend_bottom <- plotting(MP_OB, "Mit potential") + theme(legend.position = "bottom", legend.title = element_blank())
gtable <- ggplot_gtable(ggplot_build(plot_legend_bottom))
legend_plot <- gtable$grobs[[which(sapply(gtable$grobs, function(x) x$name) == "guide-box")]]
# C) Merging plots
#
# Create the figure grid with extra space for images
fig_results_energy <- plot_grid(
  fig_OB_energy, fig_OT_energy, NULL,
  nrow = 3, rel_heights = c(0.9, 0.9, 0.15))
# Final composition: Merge everything, adding images and legend
final_plot_energy <- ggdraw(fig_results_energy) +
  # Insert OB image in the top right
  draw_image("./Others/OB.png", x = 0.87, y = 0.79, width = 0.12, height = 0.2) +
  # Insert OT image in the middle right
  draw_image("./Others/OT.png", x = 0.87, y = 0.34, width = 0.11, height = 0.18) +
  # Insert legend at the bottom-right
  draw_grob(legend_plot, x = 0.45, y = 0.05, width = 0.0001, height = 0.0001) +
  # Insert title for each plot
  annotate("text", x = 0.032, y = 0.97, label = "A", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.487, y = 0.97, label = "B", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.032, y = 0.515, label = "C", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.49, y = 0.515, label = "D", hjust = 1, vjust = 1, size = 7, fontface = "bold")
#
# Print final plot
ggsave(here("./output/figures/text/results_energy.png"), plot = final_plot_energy, width = 21, height = 14, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/text/results_energy.png")
#
#
#
#
#| label: fig-results_oxidative
#| fig-cap: "Estimates of ROS (A, D), DNA damage (B, E), and lipid peroxidation (C, F) in the olfactory bulbs (A - C) and optic tecta (D - F) of L. delicata hatchlings as a function of the different prenatal conditions. Black dots indicate the posterior mean, and the bars represent the SD of the estimates. The y-axis represents the posterior estimates of the variable of interest, and the x-axis represents the different prenatal conditions. Lines with asterisks represent significant differences between groups based on pMCMC values (pMCMC < 0.05), no lines indicate no significant differences."
#| fig-name: "fig-results_oxidative"
#
source(here("R", "func.R"))
#
# A) Plotting the results for OB/Chemical stimulus
plot_ros_OB <- plotting(ROS_OB, "ROS")
plot_dnadamage_OB <- plotting(DNA_OB, "DNA damage")
plot_peroxidation_OB <- plotting(LP_OB, "Lipid peroxidation")
#
fig_OB_oxidative_top <- plot_grid(plot_ros_OB, NULL, NULL,
                    nrow = 1, rel_widths = c(0.8, 0.1, 0.8))
fig_OB_oxidative_bottom <- plot_grid(plot_dnadamage_OB, NULL, plot_peroxidation_OB,
                    nrow = 1, rel_widths = c(0.8, 0.1, 0.8))
fig_OB_oxidative <- plot_grid(fig_OB_oxidative_top, fig_OB_oxidative_bottom,
                    nrow = 2)
#
# B) Plotting the results for OT/Visual stimulus
plot_ros_OT <- plotting(ROS_OT, "ROS")
plot_dnadamage_OT <- plotting(DNA_OT, "DNA damage")
plot_peroxidation_OT <- plotting(LP_OT, "Lipid peroxidation")
#
fig_OT_oxidative_top <- plot_grid(plot_ros_OT, NULL, NULL,
                    nrow = 1, rel_widths = c(0.8, 0.1, 0.8))
fig_OT_oxidative_bottom <- plot_grid(plot_dnadamage_OT, NULL, plot_peroxidation_OT,
                    nrow = 1, rel_widths = c(0.8, 0.1, 0.8))
fig_OT_oxidative <- plot_grid(fig_OT_oxidative_top, fig_OT_oxidative_bottom,
                    nrow = 2)
#
# Getting the legend (get_legend() does not work if the legend is on the top or bottom)
plot_legend_bottom <- plotting(ROS_OB, "ROS") + theme(legend.position = "bottom", legend.title = element_blank())
gtable <- ggplot_gtable(ggplot_build(plot_legend_bottom))
legend_plot <- gtable$grobs[[which(sapply(gtable$grobs, function(x) x$name) == "guide-box")]]
# C) Merging plots
#
# Create the figure grid with extra space for images
fig_results_oxidative <- plot_grid(
  fig_OB_oxidative, NULL, fig_OT_oxidative, NULL,
  nrow = 4, rel_heights = c(1, 0.05, 1, 0.1))
# Final composition: Merge everything, adding images and legend
final_plot_oxidative <- ggdraw(fig_results_oxidative) +
  # Insert OB image in the top right
  draw_image("./Others/OB.png", x = 0.56, y = 0.755, width = 0.2, height = 0.3) +
  # Insert OT image in the middle right
  draw_image("./Others/OT.png", x = 0.56, y = 0.28, width = 0.19, height = 0.28) +
  # Insert legend at the bottom-right
  draw_grob(legend_plot, x = 0.5, y = 0.025, width = 0.0001, height = 0.0001) +
  # Insert title for each plot
  annotate("text", x = 0.0365, y = 0.986, label = "A", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.0365, y = 0.755, label = "B", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.567, y = 0.755, label = "C", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.0365, y = 0.5, label = "D", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.0365, y = 0.268, label = "E", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("text", x = 0.567, y = 0.268, label = "F", hjust = 1, vjust = 1, size = 7, fontface = "bold") +
  annotate("segment", x = 0.46, xend= 0.46, y = 0.21, yend = 0.26, size = 0.5, colour = "black") +
  annotate("segment", x = 0.45, xend= 0.46, y = 0.21, yend = 0.21, size = 0.5, colour = "black") +
  annotate("segment", x = 0.45, xend= 0.46, y = 0.26, yend = 0.26, size = 0.5, colour = "black") +
  annotate("text", x = 0.47, y = 0.228, label = "*", hjust = 0.5, vjust = 0.5, size = 7)
#
# Print final plot
ggsave(here("./output/figures/text/results_oxidative.png"), plot = final_plot_oxidative, width = 18, height = 24, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/text/results_oxidative.png")
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#
#
#
#
#
#
#| label: neurons_size
# Getting the size of the neurons from the validated data
# Olfactory Bulbs
# A) Experiment cells
cells_exp_OB <- mit_df %>% filter(region == "OB")
mean_cells_size_OB <- format_dec(mean(cells_exp_OB$mean_size), 3)
min_cells_OB <- format_dec(quantile(cells_exp_OB$mean_size, 0.025), 3)
max_cells_OB <- format_dec(quantile(cells_exp_OB$mean_size, 0.975), 3)
# B) Neurons validation
neurons_OB <- read.csv(here("./data/size_neurons_OB.csv")) %>%
  filter(X == "Olfactory Bulb_Tube_001_001.fcs")
neurons_size_OB <- log10(neurons_OB$neurons...Geometric.Mean..FSC.A.)
##
##
# Optic Tecta
# A) Experiment cells
cells_exp_OT <- mit_df %>% filter(region == "OT")
mean_cells_size_OT <- format_dec(mean(cells_exp_OT$mean_size), 3)
min_cells_OT <- format_dec(quantile(cells_exp_OT$mean_size, 0.025), 3)
max_cells_OT <- format_dec(quantile(cells_exp_OT$mean_size, 0.975), 3)
# B) Neurons validation
neurons_OT <- read.csv(here("./data/size_neurons_OT.csv")) %>%
  filter(X == "Optical Tecta_Tube_002_002.fcs")
neurons_size_OT <- log10(neurons_OT$cells.neurons...Geometric.Mean..FSC.A.)
#
#
#
#
#
#| label: fig-S13
#| fig-cap: "Fluorescence microscopy image of neuron nuclei stained with NeuN-Alexa 488. The image was excited using a 488 nm wavelength, and green fluorescence corresponds to NeuN-positive cells. In panel A, NeuN-positive cells for the homogenate of the Olfactory Bulb (OB). In panel B, NeuN-positive cells for the homogenate of the Optic Tecta (OT)."
#
#
imgOB <- image_read(here("output/figures/microscope/OB_neurons_black.jpg"))  %>%
  image_resize("600x400!")
plot_OB <- rasterGrob(imgOB, interpolate = TRUE)
#
imgOT <- image_read(here("output/figures/microscope/OT_neurons_black.jpg")) %>%
  image_resize("600x400!")
plot_OT <- rasterGrob(imgOT, interpolate = TRUE)
#
plot_neun <- plot_grid(plot_OB, plot_OT, nrow = 1)
#
fig_neun <- ggdraw(plot_neun) +
  annotate("text", x = 0.05, y = 0.92, label = "A", hjust = 1, vjust = 1, size = 4, family = "Times", color = "white") +
  annotate("text", x = 0.55, y = 0.92, label = "B", hjust = 1, vjust = 1, size = 4, family = "Times", color = "white") 
ggsave(here("./output/figures/suppl/Figure_SX13.png"), plot = fig_neun, width = 15, height = 5, units = "cm", dpi = 600, bg = "white")
knitr::include_graphics("./output/figures/suppl/Figure_SX13.png")
#
#
#
#
#
#
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#

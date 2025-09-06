####################################
# Exploratory analyses
####################################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes,
ggh4x, cowplot, fitdistrplus, MASS, goftest)
#
source(here("R", "func.R"))
#
# 
data_expl <- df_merged %>%
       mutate(treatment = paste(cort, temp, sep = "-"))
#
# A) Make histograms for the main variables and figures factorizing main predictors
var <- c("t_D", "mean_mitodensity", "mean_potential", "mean_conpotential", "mean_ros", "mean_conros", "mean_dnadamage", "mean_peroxidation")

hist_list <- list()
for (i in var) {
  bin_width <- (max(data_expl[[i]], na.rm = TRUE) - min(data_expl[[i]], na.rm = TRUE)) / sqrt(length(data_expl[[i]]))
  # Labels
  if (i == "t_D"){
     lab <- "Detection latency (s)"
  } else if (i == "mean_mitodensity") {
     lab <- "Mitochondrial density"
  } else if (i == "mean_potential") {
     lab <- "Mitochondrial potential"
  } else if (i == "mean_ros") {
     lab <- "ROS production"
  } else if (i == "mean_dnadamage") {
     lab <- "DNA damage"
  } else if (i == "mean_peroxidation") {
     lab <- "Lipid peroxidation"
  } else if (i == "mean_conpotential") {
    lab <- "Mit potential/Mit density"
  } else if (i == "mean_conros"){
    lab <- "ROS/Mit density"
  }
  # Histograms
  hist_plot <- ggplot(data_expl, aes(x = .data[[i]])) +
    geom_histogram(binwidth = bin_width, fill = "#062d00", color = "black", alpha = 0.5) +
    theme_classic() +
    labs(x = lab, y = "Counts") +
    theme(axis.title = element_text(size = 12, family = "Times"),
          axis.text = element_text(size = 10, family = "Times"))
  hist_list[[i]] <- hist_plot
}
fig_hist <- plot_grid(plotlist = hist_list, ncol = 3)
#
for (i in var) {
  if(i == "t_D"){
    plot_age <- ggplot(data_expl, aes(x = age_trial, y = t_D, color = treatment)) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("CORT-Cold"="#00008B", "Control-Cold"="#68bde1", 
                                   "CORT-Hot"="#b50101", "Control-Hot"="#fa927d")) +
      facet_wrap(~ stimulus, ncol = 1) +
      theme_classic() +
      labs(y = "Detection latency (s)", x = "Age (days)") + 
      theme(axis.title = element_text(size = 12, family = "Times"),
            axis.text = element_text(size = 10, family = "Times"))    
    plot_sex <- ggplot(data_expl, aes(x = sex, y = t_D, fill = treatment)) +
      geom_violin(alpha = 0.5, color = "black") +
      scale_fill_manual(values = c("CORT-Cold"="#00008B", "Control-Cold"="#68bde1", 
                                  "CORT-Hot"="#b50101", "Control-Hot"="#fa927d")) +
      facet_wrap(~ stimulus, ncol = 1) +
      theme_classic() +
      labs(y = "Detection latency (s)", x = "") + 
      theme(axis.title = element_text(size = 12, family = "Times"),
            axis.text = element_text(size = 10, family = "Times"),
            legend.position = "none")
    
    plot_prey <- ggplot(data_expl, aes(x = prey, y = t_D, fill = treatment)) +
      geom_violin(alpha = 0.5, color = "black") +
      scale_fill_manual(values = c("CORT-Cold"="#00008B", "Control-Cold"="#68bde1", 
                                  "CORT-Hot"="#b50101", "Control-Hot"="#fa927d")) +
      facet_wrap(~ stimulus, ncol = 1) +
      theme_classic() +
      labs(y = "Detection latency (s)", x = "") + 
      theme(axis.title = element_text(size = 12, family = "Times"),
            axis.text = element_text(size = 10, family = "Times"))
    fig_expl_t_D <- plot_grid(plot_age, plot_sex, plot_prey, nrow = 1, rel_widths = c(0.34, 0.33, 0.33))
  } else {
  if (i == "t_D"){
     lab <- "Detection latency (s)"
  } else if (i == "mean_mitodensity") {
     lab <- "Mitochondrial density"
  } else if (i == "mean_potential") {
     lab <- "Mitochondrial potential"
  } else if (i == "mean_ros") {
     lab <- "ROS production"
  } else if (i == "mean_dnadamage") {
     lab <- "DNA damage"
  } else if (i == "mean_peroxidation") {
     lab <- "Lipid peroxidation"
  } else if (i == "mean_conpotential") {
    lab <- "Mit potential/Mit density"
  } else if (i == "mean_conros"){
    lab <- "ROS/Mit density"
  }
    plot_age <- ggplot(data_expl, aes(x = age_trial, y = .data[[i]], color = treatment)) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("CORT-Cold"="#00008B", "Control-Cold"="#68bde1", 
                                   "CORT-Hot"="#b50101", "Control-Hot"="#fa927d")) +
      facet_wrap(~ stimulus, ncol = 1) +
      theme_classic() +
      labs(y = lab, x = "Age (days)") + 
      theme(axis.title = element_text(size = 12, family = "Times"),
            axis.text = element_text(size = 10, family = "Times"))
    
    plot_sex <- ggplot(data_expl, aes(x = sex, y = .data[[i]], fill = treatment)) +
      geom_violin(alpha = 0.5, color = "black") +
      scale_fill_manual(values = c("CORT-Cold"="#00008B", "Control-Cold"="#68bde1", 
                                  "CORT-Hot"="#b50101", "Control-Hot"="#fa927d")) +
      facet_wrap(~ stimulus, ncol = 1) +
      theme_classic() +
      labs(y = lab, x = "") + 
      theme(axis.title = element_text(size = 12, family = "Times"),
            axis.text = element_text(size = 10, family = "Times"))
    
  assign(paste0("fig_expl_", i), plot_grid(plot_age, plot_sex, nrow = 1, rel_widths = c(0.5, 0.5)))
  }
  
}
#
ggsave(here("./output/figures/exploratory/fig_hist.png"), fig_hist, width = 20, height = 20, dpi = 300)
ggsave(here("./output/figures/exploratory/fig_expl_t_D.png"), fig_expl_t_D, width = 20, height = 10, dpi = 300)
ggsave(here("./output/figures/exploratory/fig_expl_mean_mitodensity.png"), fig_expl_mean_mitodensity, width = 20, height = 10, dpi = 300)
ggsave(here("./output/figures/exploratory/fig_expl_mean_potential.png"), fig_expl_mean_potential, width = 20, height = 10, dpi = 300)
ggsave(here("./output/figures/exploratory/fig_expl_mean_conpotential.png"), fig_expl_mean_conpotential, width = 20, height = 10, dpi = 300)
ggsave(here("./output/figures/exploratory/fig_expl_mean_ros.png"), fig_expl_mean_ros, width = 20, height = 10, dpi = 300)
ggsave(here("./output/figures/exploratory/fig_expl_mean_conros.png"), fig_expl_mean_conros, width = 20, height = 10, dpi = 300)
ggsave(here("./output/figures/exploratory/fig_expl_mean_dnadamage.png"), fig_expl_mean_dnadamage, width = 20, height = 10, dpi = 300)
ggsave(here("./output/figures/exploratory/fig_expl_mean_peroxidation.png"), fig_expl_mean_peroxidation, width = 20, height = 10, dpi = 300)
#
# B) Get the distribution of the main independent variables
source(here("R", "func.R"))
#
data_expl$t_D <- data_expl$t_D + 1

qq_plots_list <- lapply(var, function(v) {
  label <- if (v == "t_D") {
    "Detection latency (s)"
  } else if (v == "mean_mitodensity") {
    "Mitochondrial density"
  } else if (v == "mean_potential") {
    "Mitochondrial potential"
  } else if (v == "mean_ros") {
    "ROS production"
  } else if (v == "mean_dnadamage") {
    "DNA damage"
  } else if (v == "mean_peroxidation") {
    "Lipid peroxidation"
  } else if (i == "mean_conpotential") {
    "Mit potential/Mit density"
  } else if (i == "mean_conros"){
    "ROS/Mit density"
  }
  else {
    v  # Fallback to variable name if not found
  }
  qq_plots_single(data_expl, v, label)
})
#
final_qq_plot <- plot_grid(plotlist = qq_plots_list, ncol = 1)
ggsave(here("./output/figures/exploratory/fig_qqplots.png"), final_qq_plot, width = 20, height = 20, dpi = 300)
#
# Log-Normal: Latency, Mitochondrial Density, DNA Damage, Lipid Peroxidation
# Normal: Mitochondrial Potential, Mit potential/Mit density, ROS Production, ROS/Mit density
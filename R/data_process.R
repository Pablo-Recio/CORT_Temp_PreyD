####################################
# Data_process
####################################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2,
               lme4, zoo, lmerTest, broom, tidybayes, forcats)
#
source(here("R", "func.R"))
# A) Processing behaviour_df
data_discr <- read.csv(here("./data/data_main.csv")) %>%
  mutate(temp = gsub("[AB]_", "", trt),
         cort = gsub("_[2][38]", "", trt))  %>%
  mutate(temp = factor(temp,
                       levels = c("23", "28"),
                       labels = c("23" = "Cold", "28" = "Hot"))) %>%
  mutate(cort = factor(cort,
                       levels = c("A", "B"),
                       labels = c("A" = "Control", "B" = "CORT"))) %>%
  mutate(test = factor(test,
                       levels = c("CC", "CV", "TC", "TV"),
                       labels = c("CC" = "Chemical_Familiar",
                                  "CV" = "Visual_Familiar",
                                  "TC" = "Chemical_Unknown",
                                  "TV" = "Visual_Unknown"))) %>%
  mutate(stimulus = gsub("_Familiar|_Unknown", "", test),
         prey = gsub("Chemical_|Visual_", "", test)) %>%
  mutate(sex = factor(sex,
                      levels = c("m", "f"),
                      labels = c("m" = "Male", "f" = "Female"))) %>%
  mutate(t_D = as.numeric(t_D)) %>%
  dplyr::select(lizard_id, clutch, sex, temp, cort, age_trial, age_euthanasia,
                trial, stimulus, prey, t_D, motivation, Id_flow) %>%
  data.frame()
#
#
# B) Processing and merging 'mitochondrial activity' (data/mit_act) and
# 'mitochondrial damage' (data/mit_damage) dfs
#
mit_act <- read.csv(here("./data/mit_act.csv")) %>%
  group_by(Id_flow, plate, region) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop") %>%
  ungroup() %>%
  dplyr::select(region,
                Id_flow,
                mean_mitodensity,
                mean_potential,
                mean_ros,
                mean_size) %>%
  data.frame()
#
mit_damage <- read.csv(here("./data/mit_damage.csv")) %>%
  group_by(Id_flow, plate, region) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop") %>%
  ungroup() %>%
  dplyr::select(region, Id_flow, mean_dnadamage, mean_peroxidation) %>%
  data.frame()
#
mit_df <- merge(mit_act, mit_damage,
                by = c("Id_flow", "region"), all = TRUE) %>%
  mutate(stimulus = factor(region,
                           levels = c("OB", "OT"),
                           labels = c("OB" = "Chemical",
                                      "OT" = "Visual")))

write.csv(mit_df, "./output/checking/mit_df.csv")
# C) Merging all three dfs
#
df_merged <- merge(data_discr, mit_df, by = c("stimulus", "Id_flow"))
#
#
# Transform variables according to distribution (see below and exploratory analyses)
# 
# Log-Normal: Latency, Mitochondrial Density, DNA Damage, Lipid Peroxidation
# Normal: Mitochondrial Potential, ROS Production
#
log_var <- c("t_D", "mean_mitodensity", "mean_dnadamage", "mean_peroxidation")
df_merged_log <- df_merged %>%
  mutate(across(all_of(log_var), ~ log(. + 1)))
#
# Convert age variables to numeric and standardize all numerical variables
clean_df <- df_merged_log %>%
  mutate(across(c(age_trial, age_euthanasia), as.numeric)) %>%
  mutate(across(where(~ is.double(.) & !is.integer(.)),
                ~ (. - mean(., na.rm = TRUE)) / (2 * sd(., na.rm = TRUE))))

write.csv(clean_df, "./output/database/clean_df.csv", row.names = FALSE)

######## DFT RBD vs controls longitudinal analysis ########
library(markdown)
library(data.table)
library(tidyverse)
library(tidyr)
library(dplyr)
library(readr)
library(xlsx)
library(ggplot2)
library(pROC)
require(MASS)
library(lme4)
library(lmerTest)
library(readxl)
library(writexl)
library(e1071)
library(broom.mixed)
library(boot)


RBD_control_DFT_data <- read_csv("DFT/RBD_control_DFT_data.csv")

table(RBD_control_DFT_data$PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline)
#rename the 'yes' 
RBD_control_DFT_data <- RBD_control_DFT_data %>%
  mutate(PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline = recode(
    PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline,
    "Yes, Survey + SurveyBRAIN + SurveySmell + SurveyBRAINSmell + SurveyAll" = "Yes",
    "No" = "No"
  ))


## remove those with essential tremor and under 50
RBD_control_DFT_data <- RBD_control_DFT_data %>%
  filter(hx_essentialtremor != "Yes") %>% 
  filter(age_session_baseline >= 50)

RBD_control_DFT_data %>%
  filter(Code == "Control", rx_anti_psychotics_use_final_ever != "Unknown") %>%
  count(rx_anti_psychotics_use_final_ever)
# remove high risk controls
RBD_control_DFT_data <- RBD_control_DFT_data %>%
  filter(!(Code == "Control" & PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline == "Yes"))
# double check
RBD_control_DFT_data %>%
  filter(Code == "Control") %>%
  count(PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline)

#check followup data on AT
RBD_control_DFT_data %>%
  rowwise() %>%
  mutate(
    non_missing_sessions = sum(!is.na(c_across(c(
      onefingertaptest_worst_at_5,
      onefingertaptest_worst_at_7,
      onefingertaptest_worst_is_5,
      onefingertaptest_worst_is_7,
      onefingertaptest_worst_ks_5,
      onefingertaptest_worst_ks_7
    ))))
  ) %>%
  ungroup() %>%
  group_by(PID2, Code) %>%
  summarise(sessions_with_data = max(non_missing_sessions), .groups = "drop") %>%
  count(Code, sessions_with_data)

#remove transgender
RBD_control_DFT_data <- RBD_control_DFT_data %>%
  filter(gender != "Transgender/Nonbinary/Other")

#need to reformat df
RBD_control_DFT_data2 <- RBD_control_DFT_data %>%
  pivot_longer(
    cols = matches("onefingertaptest_worst_(at|is|ks)_(baseline|lastfollowup)$"),
    names_to = c("Measure", "Timepoint"),
    names_pattern = "onefingertaptest_worst_(.*)_(baseline|lastfollowup)",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Measure,
    values_from = Value,
    names_prefix = "onefingertaptest_worst_"
  ) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("baseline", "lastfollowup"))
  )

#get rid of '0'
RBD_control_DFT_data2 <- RBD_control_DFT_data2 %>%
  filter(
    abs(onefingertaptest_worst_at) > 1e-6,
    abs(onefingertaptest_worst_is) > 1e-6,
    abs(onefingertaptest_worst_ks) > 1e-6
  )

summary(RBD_control_DFT_data2)

valid_ids <- RBD_control_DFT_data2 %>%
  group_by(PID2, Timepoint) %>%
  summarise(
    has_at = all(!is.na(onefingertaptest_worst_at)),
    has_is = all(!is.na(onefingertaptest_worst_is)),
    has_ks = all(!is.na(onefingertaptest_worst_ks)),
    .groups = "drop"
  ) %>%
  filter(has_at & has_is & has_ks) %>%
  count(PID2) %>%
  filter(n == 2) %>%  # Must have all 3 parameters at both baseline and follow-up
  pull(PID2)

length(valid_ids)

# Filter the original data
RBD_control_DFT_data2_filtered <- RBD_control_DFT_data2 %>%
  filter(PID2 %in% valid_ids)

#check for NAs
RBD_control_DFT_data2_filtered %>%
  dplyr::select(onefingertaptest_worst_is, onefingertaptest_worst_at, onefingertaptest_worst_ks) %>% 
  summarise_all(~ sum(is.na(.)))
# no NAs


table(RBD_control_DFT_data2_filtered$Code, RBD_control_DFT_data2_filtered$hx_neuro_parkinsons) # 8 RBD participants with PD
### make subgroup for those with PD
RBD_and_PD_longitudinal_data <- RBD_control_DFT_data2_filtered %>%
  filter(hx_neuro_parkinsons == "Yes")
#remove those with PD
RBD_control_DFT_data2_filtered <- RBD_control_DFT_data2_filtered %>%
  anti_join(RBD_and_PD_longitudinal_data, by = c("PID2", "TokenFinal"))

write.csv(RBD_control_DFT_data2_filtered, "DFT/Pilot study dft/Longitudinal/RBD vs controls/RBD_con_long_BT.csv", row.names = FALSE)


# log transform IS
RBD_control_DFT_data2_filtered$onefingertaptest_worst_is <- log(RBD_control_DFT_data2_filtered$onefingertaptest_worst_is)
#remove NaN rows
RBD_control_DFT_data2_filtered <- RBD_control_DFT_data2_filtered %>%
  filter(!is.nan(onefingertaptest_worst_is))


# Winsorization function to replace outliers with the max/min within bounds
winsorize_with_max_min <- function(x) {
  qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- qnt[2] - qnt[1]
  lower <- qnt[1] - 1.5 * iqr
  upper <- qnt[2] + 1.5 * iqr
  
  # Find the maximum and minimum values within the bounds
  capped_min <- min(x[x >= lower], na.rm = TRUE)
  capped_max <- max(x[x <= upper], na.rm = TRUE)
  
  # Replace outliers with the max/min within the bounds
  x[x < lower] <- capped_min
  x[x > upper] <- capped_max
  
  return(x)
}

# Apply to the selected columns to controls
RBD_control_DFT_data2_filtered <- RBD_control_DFT_data2_filtered %>%
  mutate(across(
    c(onefingertaptest_worst_is, onefingertaptest_worst_at, onefingertaptest_worst_ks),
    ~ if_else(Code == "Control", winsorize_with_max_min(.), .)
  ))

summary(RBD_control_DFT_data2_filtered)


# Check for missing values in specific columns
sum(is.na(RBD_control_DFT_data2_filtered$gender))
sum(is.na(RBD_control_DFT_data2_filtered$age_session_baseline))
sum(is.na(RBD_control_DFT_data2_filtered$onefingertaptest_date_baseline))
sum(is.na(RBD_control_DFT_data2_filtered$onefingertaptest_date_lastfollowup))


## reverse direction of KS
RBD_control_DFT_data2_filtered$onefingertaptest_worst_ks_rev <- -RBD_control_DFT_data2_filtered$onefingertaptest_worst_ks
#check 
head(RBD_control_DFT_data2_filtered[, c("onefingertaptest_worst_ks", "onefingertaptest_worst_ks_rev")])


### check distribution of data - density plots
RBD_control_DFT_data2_filtered %>%
  filter(Timepoint == "baseline") %>%
  ggplot(aes(x = age_session_baseline, group = Code, fill = Code)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(45, 95)


DFT_baseline_male <- RBD_control_DFT_data2_filtered[RBD_control_DFT_data2_filtered$gender == "Male", ]

DFT_baseline_male %>%
  filter(Timepoint == "baseline") %>%
  ggplot(aes(x = age_session_baseline, group = Code, fill = Code)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(45, 95)


DFT_baseline_female <- RBD_control_DFT_data2_filtered[RBD_control_DFT_data2_filtered$gender == "Female", ]

DFT_baseline_female %>%
  filter(Timepoint == "baseline") %>%
  ggplot(aes(x = age_session_baseline, group = Code, fill = Code)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(45, 95)


# IS
ggplot(RBD_control_DFT_data2_filtered, aes(x = onefingertaptest_worst_is, group = interaction(Code, Timepoint), fill = Code)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Timepoint) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(4,12) +
  labs(x = "DFT Worst IS Adjusted", y = "Density", title = "Baseline vs Follow-up by Group")


# AT
ggplot(RBD_control_DFT_data2_filtered, aes(x = onefingertaptest_worst_at, group = interaction(Code, Timepoint), fill = Code)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Timepoint) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(0, 300) +
  labs(x = "DFT Worst AT Adjusted", y = "Density", title = "Baseline vs Follow-up by Group")

# KS
ggplot(RBD_control_DFT_data2_filtered, aes(x = onefingertaptest_worst_ks, group = interaction(Code, Timepoint), fill = Code)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Timepoint) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(0,150) +
  labs(x = "DFT Worst KS Adjusted", y = "Density", title = "Baseline vs Follow-up by Group")


df_filtered <- RBD_control_DFT_data2_filtered

#### demographics ####
table(df_filtered$Code, df_filtered$Timepoint)

t_test_age_result <- t.test(age_session_baseline ~ Code, data = df_filtered[df_filtered$Timepoint == "baseline", ])
print(t_test_age_result)

std_dev_ages <- df_filtered %>% filter(Timepoint == "baseline") %>% group_by(Code) %>% 
  summarize(std_dev_ages = sd(age_session_baseline, na.rm = TRUE))
print(std_dev_ages)

sex_percentage <- df_filtered %>%
  filter(Timepoint == "baseline") %>%
  count(Code, gender) %>%
  group_by(Code) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()
print(sex_percentage)

chi_test_result_sex <- df_filtered %>%
  filter(Timepoint == "baseline") %>%
  with(table(gender, Code)) %>%
  chisq.test()

print(chi_test_result_sex)

# Wilcoxon test for IS
df_filtered %>%
  group_by(Timepoint, Code) %>%
  summarise(median_bt = median(onefingertaptest_worst_is, na.rm = TRUE), .groups = "drop") %>%
  print() %>%
  invisible() %>%
  { df_filtered %>%
      group_split(Timepoint) %>%
      setNames(unique(df_filtered$Timepoint)) %>%
      lapply(function(x) wilcox.test(onefingertaptest_worst_is ~ Code, data = x)$p.value) }

t.test(onefingertaptest_worst_at ~ Code, data = subset(df_filtered, Timepoint == "baseline"))
t.test(onefingertaptest_worst_at ~ Code, data = subset(df_filtered, Timepoint == "lastfollowup"))

t.test(onefingertaptest_worst_ks ~ Code, data = subset(df_filtered, Timepoint == "baseline"))
t.test(onefingertaptest_worst_ks ~ Code, data = subset(df_filtered, Timepoint == "lastfollowup"))


##############   Linear Mixed Effect Regression    ##############
#assign cases and controls to data sheet
df_filtered2 <- df_filtered %>%
  mutate(Code = ifelse(Code == "iRBD", 1, 0))
df_filtered2$Code <- as.factor(df_filtered2$Code)
table(df_filtered2$Code, df_filtered2$Timepoint)


class(df_filtered2$onefingertaptest_date_baseline)
class(df_filtered2$onefingertaptest_date_lastfollowup)
df_filtered2 <- df_filtered2 %>%
  mutate(
    onefingertaptest_date_baseline = as.Date(onefingertaptest_date_baseline, format = "%d/%m/%Y"),
    onefingertaptest_date_lastfollowup = as.Date(onefingertaptest_date_lastfollowup, format = "%d/%m/%Y")
  )

#to avoid collinearity produce elapsed_time
df_filtered2 <- df_filtered2 %>%
  group_by(PID2) %>%
  mutate(
    elapsed_time = as.numeric(difftime(onefingertaptest_date_lastfollowup, min(onefingertaptest_date_baseline), units = "days")) / 365.25
  ) %>%
  ungroup()


# check mean followup time
df_filtered2 %>%
  group_by(PID2) %>%
  summarise(followup_years = max(elapsed_time, na.rm = TRUE), .groups = "drop") %>%
  summarise(
    mean_followup   = mean(followup_years, na.rm = TRUE),
    median_followup = median(followup_years, na.rm = TRUE),
    min_followup    = min(followup_years, na.rm = TRUE),
    max_followup    = max(followup_years, na.rm = TRUE)
  )


df_filtered2 %>%
  group_by(PID2, Code) %>%
  summarise(followup_years = max(elapsed_time, na.rm = TRUE), .groups = "drop") %>%
  group_by(Code) %>%
  summarise(
    mean_followup   = mean(followup_years, na.rm = TRUE),
    median_followup = median(followup_years, na.rm = TRUE),
    min_followup    = min(followup_years, na.rm = TRUE),
    max_followup    = max(followup_years, na.rm = TRUE),
    n_participants  = n()
  )


## need to add random slopes- to allow individual variation in slopes

# Incoordination score
IS_model <- lmer(
  onefingertaptest_worst_is ~ elapsed_time * Code + age_session_baseline + gender + (1 | PID2),
  data = df_filtered2
)
summary(IS_model)
# create a summary table with CIs
IS <- tidy(IS_model, effects = "fixed", conf.int = TRUE) %>%
  rename(
    Term = term,         
    Estimate = estimate, 
    Std_Error = std.error,
    t_value = statistic,
    DF = df,             
    P_value = p.value,   
    Lower_CI = conf.low, 
    Upper_CI = conf.high 
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(IS)


# Helper for median + IQR
median_iqr <- function(x) {
  data.frame(
    y = median(x, na.rm = TRUE),
    ymin = quantile(x, 0.25, na.rm = TRUE),
    ymax = quantile(x, 0.75, na.rm = TRUE)
  )
}

# Plot IS across years for cases vs controls (median + IQR)
ggplot(df_filtered2, aes(x = Timepoint, y = onefingertaptest_worst_is, color = Code, group = Code)) +
  stat_summary(fun = median, geom = "line", linewidth = 1.2) +
  stat_summary(fun = median, geom = "point", size = 3, shape = 21, fill = "white", stroke = 1) +
  stat_summary(fun.data = median_iqr, geom = "errorbar", width = 0.2, linewidth = 0.8) +
  labs(
    title = "DFT Incoordination Score at Baseline and Follow-up",
    subtitle = "PSG_RBD vs. Controls",
    x = "Timepoint",
    y = "Median Incoordination Score (with IQR)",
    color = "Group"
  ) +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2"),
                     labels = c("Controls", "Cases")) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12),         # keep subtitle readable
    axis.text = element_text(size = 11),            # axis tick labels
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )



# Kinesia score
KS_model <- lmer(
  onefingertaptest_worst_ks_rev ~ elapsed_time * Code + age_session_baseline + gender + (1 | PID2),
  data = df_filtered2
)
summary(KS_model)
# create a summary table with CIs
KS <- tidy(KS_model, effects = "fixed", conf.int = TRUE) %>%
  rename(
    Term = term,         
    Estimate = estimate, 
    Std_Error = std.error,
    t_value = statistic,
    DF = df,             
    P_value = p.value,   
    Lower_CI = conf.low, 
    Upper_CI = conf.high 
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(KS)


# Plot KS across years for cases vs controls
ggplot(df_filtered2, aes(x = Timepoint, y = onefingertaptest_worst_ks, color = Code, group = Code)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3, shape = 21, fill = "white", stroke = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.8) +  # Changed here
  labs(
    title = "DFT Kinesia Score at Baseline and Follow-up",
    subtitle = "PSG_RBD vs. Controls",
    x = "Timepoint",
    y = "Mean Kinesia Score",
    color = "Group"
  ) +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2"), labels = c("Controls", "Cases")) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )



# Akinesia Time
AT_model <- lmer(
  onefingertaptest_worst_at ~ elapsed_time * Code + age_session_baseline + gender + (1 | PID2),
  data = df_filtered2
)
summary(AT_model)
# create a summary table with CIs
AT <- tidy(AT_model, effects = "fixed", conf.int = TRUE) %>%
  rename(
    Term = term,         
    Estimate = estimate, 
    Std_Error = std.error,
    t_value = statistic,
    DF = df,             
    P_value = p.value,   
    Lower_CI = conf.low, 
    Upper_CI = conf.high 
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(AT)


# Plot AT across years for cases vs controls
ggplot(df_filtered2, aes(x = Timepoint, y = onefingertaptest_worst_at, color = Code, group = Code)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3, shape = 21, fill = "white", stroke = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.8) + 
  labs(
    title = "DFT Akinesia Time at Baseline and Follow-up",
    subtitle = "PSG_RBD vs. Controls",
    x = "Timepoint",
    y = "Mean Akinesia Time",
    color = "Group"
  ) +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2"), labels = c("Controls", "Cases")) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )


#### check linearity of age
#IS
df_filtered2$age_squared <- df_filtered2$age_session_baseline^2

#scale variables
df_filtered2 <- df_filtered2 %>%
  mutate(
    age_baseline_z = scale(age_session_baseline, center = TRUE, scale = TRUE),
    age_squared_z = scale(age_session_baseline^2, center = TRUE, scale = TRUE),
    elapsed_time_z = scale(elapsed_time, center = TRUE, scale = TRUE)
  )


IS_model2 <- lmer(
  onefingertaptest_worst_is ~ elapsed_time_z * Code + age_baseline_z + age_squared_z + gender + (1 | PID2),
  data = df_filtered2
)
summary(IS_model2)
# create a summary table with CIs
IS2 <- tidy(IS_model2, effects = "fixed", conf.int = TRUE) %>%
  rename(
    Term = term,         
    Estimate = estimate, 
    Std_Error = std.error,
    t_value = statistic,
    DF = df,             
    P_value = p.value,   
    Lower_CI = conf.low, 
    Upper_CI = conf.high 
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(IS2)
anova(IS_model, IS_model2)


#KS
KS_model2 <- lmer(
  onefingertaptest_worst_ks_rev ~ elapsed_time_z * Code + age_baseline_z + age_squared_z + gender + (1 | PID2),
  data = df_filtered2
)
summary(KS_model2)
# create a summary table with CIs
KS2 <- tidy(KS_model2, effects = "fixed", conf.int = TRUE) %>%
  rename(
    Term = term,         
    Estimate = estimate, 
    Std_Error = std.error,
    t_value = statistic,
    DF = df,             
    P_value = p.value,   
    Lower_CI = conf.low, 
    Upper_CI = conf.high 
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(KS2)
anova(KS_model, KS_model2)


#AT
AT_model2 <- lmer(
  onefingertaptest_worst_at ~ elapsed_time_z * Code + age_baseline_z + age_squared_z + gender + (1 | PID2),
  data = df_filtered2
)
summary(AT_model2)
# create a summary table with CIs
AT2 <- tidy(AT_model2, effects = "fixed", conf.int = TRUE) %>%
  rename(
    Term = term,         
    Estimate = estimate, 
    Std_Error = std.error,
    t_value = statistic,
    DF = df,             
    P_value = p.value,   
    Lower_CI = conf.low, 
    Upper_CI = conf.high 
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(AT2)
anova(AT_model, AT_model2)

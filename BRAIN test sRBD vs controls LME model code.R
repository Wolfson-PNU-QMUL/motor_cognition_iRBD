######## BRAIN test sRBD vs controls longitudinal analysis ########
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
library(ggdist)
library(emmeans)
library(nortest)

sRBD_control_BT_data <- read_csv("BRAIN test/sRBD_control_BT_data.csv")

table(sRBD_control_BT_data$PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline)
#rename the 'yes' 
sRBD_control_BT_data <- sRBD_control_BT_data %>%
  mutate(PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline = recode(
    PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline,
    "Yes, Survey + SurveyBRAIN + SurveySmell + SurveyBRAINSmell + SurveyAll" = "Yes",
    "No" = "No"
  ))


## remove those with essential tremor and under 50
sRBD_control_BT_data <- sRBD_control_BT_data %>%
  filter(hx_essentialtremor != "Yes") %>% 
  filter(age_session_baseline >= 50)

sRBD_control_BT_data %>%
  filter(Code == "Control", rx_anti_psychotics_use_final_ever != "Unknown") %>%
  count(rx_anti_psychotics_use_final_ever)
# remove high risk controls
sRBD_control_BT_data <- sRBD_control_BT_data %>%
  filter(!(Code == "Control" & PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline == "Yes"))
# double check
sRBD_control_BT_data %>%
  filter(Code == "Control") %>%
  count(PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline)


#check followup data on AT
sRBD_control_BT_data %>%
  rowwise() %>%
  mutate(
    non_missing_sessions = sum(!is.na(c_across(c(
      braintest_worst_at_0,
      braintest_worst_at_1,
      braintest_worst_at_2,
      braintest_worst_at_3,
      braintest_worst_at_adj_4,
      braintest_worst_at_adj_5,
      braintest_worst_at_adj_7,
      braintest_worst_at_adj_8
    ))))
  ) %>%
  ungroup() %>%
  group_by(PID2, Code) %>%
  summarise(sessions_with_data = max(non_missing_sessions), .groups = "drop") %>%
  count(Code, sessions_with_data)


#remove transgender
sRBD_control_BT_data <- sRBD_control_BT_data %>%
  filter(gender != "Transgender/Nonbinary/Other")


#need to reformat df so can run lme- year in separate rows rather than column
sRBD_control_BT_data2 <- sRBD_control_BT_data %>%
  pivot_longer(
    cols = c(
      braintest_worst_at_adj_4, braintest_worst_at_adj_5, braintest_worst_at_adj_7,
      braintest_worst_at_adj_8, braintest_worst_is_adj_4, braintest_worst_is_adj_5, 
      braintest_worst_is_adj_7, braintest_worst_is_adj_8,
      braintest_worst_ks_4, braintest_worst_ks_5, braintest_worst_ks_7, braintest_worst_ks_adj_8
    ),
    names_to = c("Measure", "Year"),
    names_pattern = "(.*)_(\\d+)$",  # <- updated to allow 2-digit years
    values_to = "Value"
  ) %>%
  mutate(
    Year = recode(Year, `4` = 1, `5` = 2, `7` = 3, `8` = 4),
    Year = as.numeric(Year)
  ) %>%
  pivot_wider(
    names_from = Measure,
    values_from = Value
  )


# combine KS non adjusted and adjusted columns
sRBD_control_BT_data2 <- sRBD_control_BT_data2 %>%
  mutate(braintest_worst_ks = if_else(!is.na(braintest_worst_ks_adj), braintest_worst_ks_adj, braintest_worst_ks)) %>%
  dplyr::select(-braintest_worst_ks_adj)

#get rid of '0'
sRBD_control_BT_data2 <- sRBD_control_BT_data2 %>%
  filter(
    abs(braintest_worst_at_adj) > 1e-6,
    abs(braintest_worst_is_adj) > 1e-6,
    abs(braintest_worst_ks) > 1e-6
  )

summary(sRBD_control_BT_data2)

## only take out baseline scores and most recent followup for each participant (2 timepoints)
# and remove participants with no follow up data

# Map Year 1–4 to Session4–Session8
sRBD_control_BT_data2 <- sRBD_control_BT_data2 %>%
  mutate(
    session_label = recode(as.character(Year),
                           `1` = "Session4",
                           `2` = "Session5",
                           `3` = "Session7",
                           `4` = "Session8")
  )

# Find IDs with both a baseline and follow-up session
valid_ids <- sRBD_control_BT_data2 %>%
  group_by(PID2) %>%
  summarise(
    has_baseline = any(session_label == braintest_longitudinal_baseline, na.rm = TRUE),
    has_followup = any(session_label == braintest_longitudinal_lastfollowup, na.rm = TRUE)
  ) %>%
  filter(has_baseline & has_followup) %>%
  pull(PID2)


# Now filter rows to just baseline and follow-up sessions *per participant*
sRBD_control_BT_data_clean <- sRBD_control_BT_data2 %>%
  filter(PID2 %in% valid_ids) %>%
  filter(session_label == braintest_longitudinal_baseline | session_label == braintest_longitudinal_lastfollowup)

#check for NAs
sRBD_control_BT_data_clean %>%
  dplyr::select(braintest_worst_at_adj, braintest_worst_is_adj, braintest_worst_ks) %>% 
  summarise_all(~ sum(is.na(.)))
# No NAs


## remove PDs
table(sRBD_control_BT_data_clean$Code, sRBD_control_BT_data_clean$hx_neuro_parkinsons) # 38 sRBD participants with PD
### make subgroup for those with PD
sRBD_and_PD_longitudinal_data <- sRBD_control_BT_data_clean %>%
  filter(hx_neuro_parkinsons == "Yes")
#remove those with PD
sRBD_control_BT_data_clean <- sRBD_control_BT_data_clean %>%
  anti_join(sRBD_and_PD_longitudinal_data, by = c("PID2", "TokenFinal"))


# log transform IS
sRBD_control_BT_data_clean$braintest_worst_is_adj <- log(sRBD_control_BT_data_clean$braintest_worst_is_adj)
#remove NaN rows
sRBD_control_BT_data_clean <- sRBD_control_BT_data_clean %>%
  filter(!is.nan(braintest_worst_is_adj))

#have outliers in control group but don't want to 'stack the deck' to get sig differences
#instead of excluding set the value at the cohort max/min
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
sRBD_control_BT_data_clean <- sRBD_control_BT_data_clean %>%
  group_by(Code) %>%
  mutate(across(c(braintest_worst_is_adj, braintest_worst_at_adj, braintest_worst_ks), ~ if (first(Code) == "Control") winsorize_with_max_min(.) else .)) %>%
  ungroup()

# check categories of data
str(sRBD_control_BT_data_clean)


# Check for missing values in specific columns
sum(is.na(sRBD_control_BT_data_clean$gender))
sum(is.na(sRBD_control_BT_data_clean$age_session_baseline))
sum(is.na(sRBD_control_BT_data_clean$braintest_date_baseline))
sum(is.na(sRBD_control_BT_data_clean$braintest_date_lastfollowup))


## reverse direction of KS
sRBD_control_BT_data_clean$bt_worst_ks_rev <- -sRBD_control_BT_data_clean$braintest_worst_ks
#check 
head(sRBD_control_BT_data_clean[, c("braintest_worst_ks", "bt_worst_ks_rev")])


### check distribution of data - density plots
sRBD_control_BT_data_clean %>%
  filter(session_label == braintest_longitudinal_baseline) %>%
  ggplot(aes(x = age_session_baseline, group = Code, fill = Code)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(45, 90)

BT_baseline_male <- sRBD_control_BT_data_clean[sRBD_control_BT_data_clean$gender == "Male", ]

BT_baseline_male %>%
  filter(session_label == braintest_longitudinal_baseline) %>%
  ggplot(aes(x = age_session_baseline, group = Code, fill = Code)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(45, 90)

BT_baseline_female <- sRBD_control_BT_data_clean[sRBD_control_BT_data_clean$gender == "Female", ]

BT_baseline_female %>%
  filter(session_label == braintest_longitudinal_baseline) %>%
  ggplot(aes(x = age_session_baseline, group = Code, fill = Code)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(45, 87)



df_filtered <- sRBD_control_BT_data_clean %>%
  filter(session_label %in% c(braintest_longitudinal_baseline, braintest_longitudinal_lastfollowup)) %>%
  mutate(Timepoint = case_when(
    session_label == braintest_longitudinal_baseline ~ "Baseline",
    session_label == braintest_longitudinal_lastfollowup ~ "Follow-up"
  ))


# IS
ggplot(df_filtered, aes(x = braintest_worst_is_adj, group = interaction(Code, Timepoint), fill = Code)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Timepoint) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(4,15) +
  labs(x = "Braintest Worst IS Adjusted", y = "Density", title = "Baseline vs Follow-up by Group")


# AT
ggplot(df_filtered, aes(x = braintest_worst_at_adj, group = interaction(Code, Timepoint), fill = Code)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Timepoint) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(-10, 300) +
  labs(x = "Braintest Worst AT Adjusted", y = "Density", title = "Baseline vs Follow-up by Group")

# KS
ggplot(df_filtered, aes(x = braintest_worst_ks, group = interaction(Code, Timepoint), fill = Code)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Timepoint) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(5,120) +
  labs(x = "Braintest Worst KS Adjusted", y = "Density", title = "Baseline vs Follow-up by Group")




#### demographics ####
table(df_filtered$Code, df_filtered$Timepoint)

t_test_age_result <- t.test(age_session_baseline ~ Code, data = df_filtered[df_filtered$Timepoint == "Baseline", ])
print(t_test_age_result)

std_dev_ages <- df_filtered %>% filter(Timepoint == "Baseline") %>% group_by(Code) %>% 
  summarize(std_dev_ages = sd(age_session_baseline, na.rm = TRUE))
print(std_dev_ages)

sex_percentage <- df_filtered %>%
  filter(Timepoint == "Baseline") %>%
  count(Code, gender) %>%
  group_by(Code) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()
print(sex_percentage)

chi_test_result_sex <- df_filtered %>%
  filter(Timepoint == "Baseline") %>%
  with(table(gender, Code)) %>%
  chisq.test()

print(chi_test_result_sex)

# Wilcoxon test for IS
df_filtered %>%
  group_by(Timepoint, Code) %>%
  summarise(median_bt = median(braintest_worst_is_adj, na.rm = TRUE), .groups = "drop") %>%
  print() %>%
  invisible() %>%
  { df_filtered %>%
      group_split(Timepoint) %>%
      setNames(unique(df_filtered$Timepoint)) %>%
      lapply(function(x) wilcox.test(braintest_worst_is_adj ~ Code, data = x)$p.value) }

t.test(braintest_worst_at_adj ~ Code, data = subset(df_filtered, Timepoint == "Baseline"))
t.test(braintest_worst_at_adj ~ Code, data = subset(df_filtered, Timepoint == "Follow-up"))

t.test(braintest_worst_ks ~ Code, data = subset(df_filtered, Timepoint == "Baseline"))
t.test(braintest_worst_ks ~ Code, data = subset(df_filtered, Timepoint == "Follow-up"))



##############   Linear Mixed Effect Regression    ##############
#assign cases and controls to data sheet
df_filtered2 <- df_filtered %>%
  mutate(Code = ifelse(Code == "sRBD", 1, 0))
df_filtered2$Code <- as.factor(df_filtered2$Code)
table(df_filtered_sRBD$Code, df_filtered_sRBD$Timepoint)


class(df_filtered2$braintest_date_baseline)
class(df_filtered2$braintest_date_lastfollowup)
df_filtered2 <- df_filtered2 %>%
  mutate(
    braintest_date_baseline = as.Date(braintest_date_baseline, format = "%d/%m/%Y"),
    braintest_date_lastfollowup = as.Date(braintest_date_lastfollowup, format = "%d/%m/%Y")
  )

#to avoid collinearity (similar data) produce elapsed_time- true time difference for each participant
df_filtered2 <- df_filtered2 %>%
  group_by(PID2) %>%
  mutate(
    elapsed_time = as.numeric(difftime(braintest_date_lastfollowup, min(braintest_date_baseline), units = "days")) / 365.25
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
# meaning participants can change at different rates over time! 


# Incoordination score
IS_model <- lmer(
  braintest_worst_is_adj ~ elapsed_time * Code + age_session_baseline + gender + (1 | PID2),
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

#bootstrapping due to small sample size
# Function to extract fixed effects (xtract the fixed-effect coefficients for each bootstrap sample)
fixef_fun <- function(fit) fixef(fit)

# Bootstrap 1000 times (resamples our dataset 1000 times with replacement, 
#refits the LME, and stores the fixed-effect estimates)
boot_results <- bootMer(IS_model, fixef_fun, nsim = 1000)

# 95% bootstrap CI
confint(boot_results, type = "perc")

# Helper for median + IQR
median_iqr <- function(x) {
  data.frame(
    y = median(x, na.rm = TRUE),
    ymin = quantile(x, 0.25, na.rm = TRUE),
    ymax = quantile(x, 0.75, na.rm = TRUE)
  )
}


# Kinesia score
KS_model <- lmer(
  bt_worst_ks_rev ~ elapsed_time * Code + age_session_baseline + gender + (1 | PID2),
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

#bootstrapping due to small sample size
# Function to extract fixed effects (xtract the fixed-effect coefficients for each bootstrap sample)
fixef_fun <- function(fit) fixef(fit)

# Bootstrap 1000 times (resamples our dataset 1000 times with replacement, 
#refits the LME, and stores the fixed-effect estimates)
boot_results <- bootMer(KS_model, fixef_fun, nsim = 1000)

# 95% bootstrap CI
confint(boot_results, type = "perc")


# Akinesia Time
AT_model <- lmer(
  braintest_worst_at_adj ~ elapsed_time * Code + age_session_baseline + gender + (1 | PID2),
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


#bootstrapping due to small sample size
# Function to extract fixed effects (xtract the fixed-effect coefficients for each bootstrap sample)
fixef_fun <- function(fit) fixef(fit)

# Bootstrap 1000 times (resamples our dataset 1000 times with replacement, 
#refits the LME, and stores the fixed-effect estimates)
boot_results <- bootMer(AT_model, fixef_fun, nsim = 1000)

# 95% bootstrap CI
confint(boot_results, type = "perc")



#### check linearity of age
#IS
ggplot(df_filtered2, aes(x = age_session_baseline, y = braintest_worst_is_adj, color = factor(Code))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE) +  # loess to explore shape
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), linetype = "dashed", se = FALSE) + # quadratic fit
  theme_minimal() +
  labs(
    x = "Baseline Age",
    y = "Worst IS Score",
    color = "Code"
  )

df_filtered2$age_squared <- df_filtered2$age_session_baseline^2

#scale variables
df_filtered2 <- df_filtered2 %>%
  mutate(
    age_baseline_z = scale(age_session_baseline, center = TRUE, scale = TRUE),
    age_squared_z = scale(age_session_baseline^2, center = TRUE, scale = TRUE),
    elapsed_time_z = scale(elapsed_time, center = TRUE, scale = TRUE)
  )


IS_model2 <- lmer(
  braintest_worst_is_adj ~ elapsed_time_z * Code + age_baseline_z + age_squared_z + gender + (1 | PID2),
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
ggplot(df_filtered2, aes(x = age_session_baseline, y = braintest_worst_ks, color = factor(Code))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE) +  # loess to explore shape
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), linetype = "dashed", se = FALSE) + # quadratic fit
  theme_minimal() +
  labs(
    x = "Baseline Age",
    y = "Worst KS Score",
    color = "Code"
  )

KS_model2 <- lmer(
  bt_worst_ks_rev ~ elapsed_time_z * Code + age_baseline_z + age_squared_z + gender + (1 | PID2),
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
ggplot(df_filtered2, aes(x = age_session_baseline, y = braintest_worst_at_adj, color = factor(Code))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE) +  # loess to explore shape
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), linetype = "dashed", se = FALSE) + # quadratic fit
  theme_minimal() +
  labs(
    x = "Baseline Age",
    y = "Worst AT Score",
    color = "Code"
  )


AT_model2 <- lmer(
  braintest_worst_at_adj ~ elapsed_time_z * Code + age_baseline_z + age_squared_z + gender + (1 | PID2),
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








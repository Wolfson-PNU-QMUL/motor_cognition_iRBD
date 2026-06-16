######## TMT sRBD vs controls longitudinal analysis ########
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

sRBD_control_TMT_data <- read_csv("sRBD_control_TMT_data.csv")
any(is.na(sRBD_control_TMT_data$PID2))
any(is.na(sRBD_control_TMT_data$TokenFinal)) 


## remove those with essential tremor and under 50
sRBD_control_TMT_data <- sRBD_control_TMT_data %>%
  filter(hx_essentialtremor != "Yes") %>% 
  filter(age_session_baseline >= 50)

sRBD_control_TMT_data %>%
  filter(Code == "Control", rx_anti_psychotics_use_final_ever != "Unknown") %>%
  count(rx_anti_psychotics_use_final_ever)

view(sRBD_control_TMT_data)

table(sRBD_control_TMT_data$PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline)
#rename the 'yes' 
sRBD_control_TMT_data <- sRBD_control_TMT_data %>%
  mutate(PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline = recode(
    PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline,
    "Yes, Survey + SurveyBRAIN + SurveySmell + SurveyBRAINSmell + SurveyAll" = "Yes",
    "No" = "No"
  ))

table(sRBD_control_TMT_data$PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline, sRBD_control_TMT_data$Code)


#remove transgender
sRBD_control_TMT_data <- sRBD_control_TMT_data %>%
  filter(gender != "Transgender/Nonbinary/Other")


#need to reformat df so can run lme- year in separate rows rather than column
sRBD_control_TMT_data2 <- sRBD_control_TMT_data %>%
  pivot_longer(
    cols = c(
      tmt_trailA_accuracy_4, tmt_trailA_accuracy_5, tmt_trailA_accuracy_7, tmt_trailA_accuracy_8,
      tmt_trailA_reactiontime_4, tmt_trailA_reactiontime_5, tmt_trailA_reactiontime_7, tmt_trailA_reactiontime_8,
      tmt_trailAminustrailB_exefcn_final_4, tmt_trailAminustrailB_exefcn_final_5, tmt_trailAminustrailB_exefcn_final_7, 
      tmt_trailAminustrailB_exefcn_final_8, tmt_trailB_accuracy_4, tmt_trailB_accuracy_5, tmt_trailB_accuracy_7,
      tmt_trailB_accuracy_8, tmt_trailB_reactiontime_4, tmt_trailB_reactiontime_5, tmt_trailB_reactiontime_7,
      tmt_trailB_reactiontime_8
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


#get rid of '0'
sRBD_control_TMT_data2 <- sRBD_control_TMT_data2 %>%
  filter(
    abs(tmt_trailA_accuracy ) > 1e-6,
    abs(tmt_trailA_reactiontime) > 1e-6,
    abs(tmt_trailAminustrailB_exefcn_final ) > 1e-6,
    abs(tmt_trailB_accuracy) > 1e-6,
    abs(tmt_trailB_reactiontime) > 1e-6,
  )

summary(sRBD_control_TMT_data2)


## only take out baseline scores and most recent followup for each participant (2 timepoints)
# and remove participants with no follow up data

# Map Year 1–4 to Session4–Session8
sRBD_control_TMT_data2 <- sRBD_control_TMT_data2 %>%
  mutate(
    session_label = recode(as.character(Year),
                           `1` = "Session4",
                           `2` = "Session5",
                           `3` = "Session7",
                           `4` = "Session8")
  )


# Find IDs with both a baseline and follow-up session
valid_ids <- sRBD_control_TMT_data2 %>%
  group_by(PID2) %>%
  summarise(
    has_baseline = any(session_label == tmt_longitudinal_baseline, na.rm = TRUE),
    has_followup = any(session_label == tmt_longitudinal_lastfollowup, na.rm = TRUE)
  ) %>%
  filter(has_baseline & has_followup) %>%
  pull(PID2)


# Now filter rows to just baseline and follow-up sessions *per participant*
sRBD_control_TMT_data2_clean <- sRBD_control_TMT_data2 %>%
  filter(PID2 %in% valid_ids) %>%
  filter(session_label == tmt_longitudinal_baseline | session_label == tmt_longitudinal_lastfollowup)

#check for NAs
sRBD_control_TMT_data2_clean %>%
  dplyr::select(tmt_trailA_reactiontime, tmt_trailAminustrailB_exefcn_final, tmt_trailB_reactiontime) %>% 
  summarise_all(~ sum(is.na(.)))
# no NAs


## remove PDs
table(sRBD_control_TMT_data2_clean$Code, sRBD_control_TMT_data2_clean$hx_neuro_parkinsons) # 40 RBD participants with PD
### make subgroup for those with PD
sRBD_and_PD_longitudinal_data <- sRBD_control_TMT_data2_clean %>%
  filter(hx_neuro_parkinsons == "Yes")
#remove those with PD
sRBD_control_TMT_data2_clean <- sRBD_control_TMT_data2_clean %>%
  anti_join(sRBD_and_PD_longitudinal_data, by = c("PID2", "TokenFinal"))


write.csv(sRBD_control_TMT_data2_clean, "Pilot/Longitudinal/sRBD vs controls/sRBD_con_long_TMT.csv", row.names = FALSE)


# inspect data for outliers
ggplot(sRBD_control_TMT_data2_clean, aes(x = factor(Year), y = tmt_trailA_reactiontime)) +
  geom_boxplot(fill = "#69b3a2", color = "#1f3552", outlier.color = "red", outlier.shape = 16, outlier.size = 2) +
  labs(
    title = "Distribution of TMT Scores by Year",
    x = "Year of Follow-Up",
    y = "Trail A Reaction Time"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )


ggplot(sRBD_control_TMT_data2_clean, aes(x = factor(Year), y = tmt_trailAminustrailB_exefcn_final)) +
  geom_boxplot(fill = "#69b3a2", color = "#1f3552", outlier.color = "red", outlier.shape = 16, outlier.size = 2) +
  labs(
    title = "Distribution of TMT Scores by Year",
    x = "Year of Follow-Up",
    y = "Executive Function"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )


ggplot(sRBD_control_TMT_data2_clean, aes(x = factor(Year), y = tmt_trailB_reactiontime)) +
  geom_boxplot(fill = "#69b3a2", color = "#1f3552", outlier.color = "red", outlier.shape = 16, outlier.size = 2) +
  labs(
    title = "Distribution of TMT Scores by Year",
    x = "Year of Follow-Up",
    y = "Trail B Reaction Time"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )



# have outliers check what cohort they belong to
Q1 <- quantile(sRBD_control_TMT_data2_clean$tmt_trailA_reactiontime, 0.25, na.rm = TRUE)
Q3 <- quantile(sRBD_control_TMT_data2_clean$tmt_trailA_reactiontime, 0.75, na.rm = TRUE)
IQR <- Q3 - Q1
lower <- Q1 - 1.5 * IQR
upper <- Q3 + 1.5 * IQR

outliers_overall <- sRBD_control_TMT_data2_clean %>%
  filter(tmt_trailA_reactiontime < lower | tmt_trailA_reactiontime > upper)

table(outliers_overall$Code)


Q1_2 <- quantile(sRBD_control_TMT_data2_clean$tmt_trailAminustrailB_exefcn_final, 0.25, na.rm = TRUE)
Q3_2 <- quantile(sRBD_control_TMT_data2_clean$tmt_trailAminustrailB_exefcn_final, 0.75, na.rm = TRUE)
IQR_2 <- Q3_2 - Q1_2
lower_2 <- Q1_2 - 1.5 * IQR_2
upper_2 <- Q3_2 + 1.5 * IQR_2

outliers_overall_2 <- sRBD_control_TMT_data2_clean %>%
  filter(tmt_trailAminustrailB_exefcn_final < lower_2 | tmt_trailAminustrailB_exefcn_final > upper_2)

table(outliers_overall_2$Code)


Q1_3 <- quantile(sRBD_control_TMT_data2_clean$tmt_trailB_reactiontime, 0.25, na.rm = TRUE)
Q3_3 <- quantile(sRBD_control_TMT_data2_clean$tmt_trailB_reactiontime, 0.75, na.rm = TRUE)
IQR_3 <- Q3_3 - Q1_3
lower_3 <- Q1_3 - 1.5 * IQR_3
upper_3 <- Q3_3 + 1.5 * IQR_3

outliers_overall_3 <- sRBD_control_TMT_data2_clean %>%
  filter(tmt_trailB_reactiontime < lower_3 | tmt_trailB_reactiontime > upper_3)

table(outliers_overall_3$Code)


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
sRBD_control_TMT_data2_clean <- sRBD_control_TMT_data2_clean %>%
  mutate(across(
    c(tmt_trailA_reactiontime, tmt_trailAminustrailB_exefcn_final, tmt_trailB_reactiontime),
    ~ if_else(Code == "Control", 
              winsorize_with_max_min(.), 
              .)
  ))


# check categories of data
str(sRBD_control_TMT_data2_clean)

# Check for missing values in specific columns
sum(is.na(sRBD_control_TMT_data2_clean$gender))
sum(is.na(sRBD_control_TMT_data2_clean$age_session_baseline))
sum(is.na(sRBD_control_TMT_data2_clean$tmt_date_baseline))
sum(is.na(sRBD_control_TMT_data2_clean$tmt_date_lastfollowup))




### check distribution of data- density plots
ggplot(sRBD_control_TMT_data2_clean, aes(x=age_session_baseline, group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(45,90)


TMT_baseline_male <- sRBD_control_TMT_data2_clean[sRBD_control_TMT_data2_clean$gender == "Male", ]

TMT_baseline_male %>%
  ggplot(aes(x = age_session_baseline, group = Code, fill = Code)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(45, 90)


TMT_baseline_female <- sRBD_control_TMT_data2_clean[sRBD_control_TMT_data2_clean$gender == "Female", ]

TMT_baseline_female %>%
  ggplot(aes(x = age_session_baseline, group = Code, fill = Code)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(45, 90)


df_filtered_sRBD <- sRBD_control_TMT_data2_clean %>%
  filter(session_label %in% c(tmt_longitudinal_baseline, tmt_longitudinal_lastfollowup)) %>%
  mutate(Timepoint = case_when(
    session_label == tmt_longitudinal_baseline ~ "Baseline",
    session_label == tmt_longitudinal_lastfollowup ~ "Follow-up"
  ))


ggplot(df_filtered, aes(x = tmt_trailA_accuracy, group = interaction(Code, Timepoint), fill = Code)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Timepoint) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(0,50)+
  labs(x = "TMT Trail A Accuracy", y = "Density", title = "Baseline vs Follow-up by Group")


ggplot(df_filtered, aes(x = tmt_trailA_reactiontime, group = interaction(Code, Timepoint), fill = Code)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Timepoint) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(200,3500) +
  labs(x = "TMT Trail A Reaction Time", y = "Density", title = "Baseline vs Follow-up by Group")


ggplot(df_filtered, aes(x = tmt_trailAminustrailB_exefcn_final, group = interaction(Code, Timepoint), fill = Code)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Timepoint) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(-1500,3700) +
  labs(x = "TMT Executive Function", y = "Density", title = "Baseline vs Follow-up by Group")


ggplot(df_filtered, aes(x = tmt_trailB_accuracy, group = interaction(Code, Timepoint), fill = Code)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Timepoint) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(0,50) +
  labs(x = "TMT Trail B Accuracy", y = "Density", title = "Baseline vs Follow-up by Group")


ggplot(df_filtered, aes(x = tmt_trailB_reactiontime, group = interaction(Code, Timepoint), fill = Code)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Timepoint) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_minimal() +
  xlim(200,5500) +
  labs(x = "TMT Trail B Reaction Time", y = "Density", title = "Baseline vs Follow-up by Group")




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


t.test(tmt_trailA_accuracy ~ Code, data = subset(df_filtered, Timepoint == "Baseline"))
t.test(tmt_trailA_accuracy ~ Code, data = subset(df_filtered, Timepoint == "Follow-up"))

t.test(tmt_trailA_reactiontime ~ Code, data = subset(df_filtered, Timepoint == "Baseline"))
t.test(tmt_trailA_reactiontime ~ Code, data = subset(df_filtered, Timepoint == "Follow-up"))

t.test(tmt_trailAminustrailB_exefcn_final ~ Code, data = subset(df_filtered, Timepoint == "Baseline"))
t.test(tmt_trailAminustrailB_exefcn_final ~ Code, data = subset(df_filtered, Timepoint == "Follow-up"))

t.test(tmt_trailB_accuracy ~ Code, data = subset(df_filtered, Timepoint == "Baseline"))
t.test(tmt_trailB_accuracy ~ Code, data = subset(df_filtered, Timepoint == "Follow-up"))

t.test(tmt_trailB_reactiontime ~ Code, data = subset(df_filtered, Timepoint == "Baseline"))
t.test(tmt_trailB_reactiontime ~ Code, data = subset(df_filtered, Timepoint == "Follow-up"))


### change education levels to secondary yes or no 
df_filtered <- df_filtered %>%
  mutate(education = ifelse(education_level %in% c("Secondary", "Further", "Higher"), "Yes", "No"))

education <- df_filtered %>%
  count(Code, education) %>%
  group_by(Code) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()
print(education)
fisher.test(table(df_filtered$Code, df_filtered$education))$p.value



##############   Linear Mixed Effect Regression    ##############
#assign cases and controls to data sheet
df_filtered2 <- df_filtered %>%
  mutate(Code = ifelse(Code == "sRBD", 1, 0))
df_filtered2$Code <- as.factor(df_filtered2$Code)
table(df_filtered2$Code, df_filtered2$Timepoint)


class(df_filtered2$tmt_date_baseline)
class(df_filtered2$tmt_date_lastfollowup)
df_filtered2 <- df_filtered2 %>%
  mutate(
    tmt_date_baseline = as.Date(tmt_date_baseline, format = "%d/%m/%Y"),
    tmt_date_lastfollowup = as.Date(tmt_date_lastfollowup, format = "%d/%m/%Y")
  )


#to avoid collinearity (similar data) produce elapsed_time- true time difference for each participant:
df_filtered2 <- df_filtered2 %>%
  group_by(PID2) %>%
  mutate(
    elapsed_time = as.numeric(difftime(tmt_date_lastfollowup, min(tmt_date_baseline), units = "days")) / 365.25
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



# Standardize reaction times
df_filtered2$TrailA_z <- scale(df_filtered2$tmt_trailA_reactiontime)
df_filtered2$TrailB_z <- scale(df_filtered2$tmt_trailB_reactiontime)

# Flip the direction so higher = better
df_filtered2$TrailA_negZ <- -1 * df_filtered2$TrailA_z
df_filtered2$TrailB_negZ <- -1 * df_filtered2$TrailB_z


# Trail A reaction time
ART_model <- lmer(
  TrailA_negZ ~ elapsed_time * Code + age_session_baseline + gender + education_level + (1 | PID2),
  data = df_filtered2
)
summary(ART_model)
# create a summary table with CIs
ART <- tidy(ART_model, effects = "fixed", conf.int = TRUE) %>%
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
print(ART)


# Plot ART across years for cases vs controls
ggplot(df_filtered2, aes(x = Timepoint, y = tmt_trailA_reactiontime, color = Code, group = Code)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3, shape = 21, fill = "white", stroke = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.8) +
  labs(
    title = "Trail A Reaction Time at Baseline and Follow-up",
    subtitle = "sRBD vs. Controls",
    x = "Timepoint",
    y = "Mean Reaction Time",
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
    plot.title = element_text(face = "bold", size = 14)
  )



# Trail B reaction time
BRT_model <- lmer(
  TrailB_negZ ~ elapsed_time * Code + age_session_baseline + gender + education_level + (1 | PID2),
  data = df_filtered2
)
summary(BRT_model)
# create a summary table with CIs
BRT <- tidy(BRT_model, effects = "fixed", conf.int = TRUE) %>%
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
print(BRT)


# Plot BRT across years for cases vs controls
ggplot(df_filtered2, aes(x = Timepoint, y = tmt_trailB_reactiontime, color = Code, group = Code)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3, shape = 21, fill = "white", stroke = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.8) +
  labs(
    title = "Trail B Reaction Time at Baseline and Follow-up",
    subtitle = "sRBD vs. Controls",
    x = "Timepoint",
    y = "Mean Reaction Time",
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
    plot.title = element_text(face = "bold", size = 14)
  )



# Executive function
EF_model <- lmer(
  tmt_trailAminustrailB_exefcn_final ~ elapsed_time * Code + age_session_baseline + gender + education_level + (1 | PID2),
  data = df_filtered2
)
summary(EF_model)
# create a summary table with CIs
EF <- tidy(EF_model, effects = "fixed", conf.int = TRUE) %>%
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
print(EF)


# Plot EF across years for cases vs controls
ggplot(df_filtered2, aes(x = Timepoint, y = tmt_trailAminustrailB_exefcn_final, color = Code, group = Code)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3, shape = 21, fill = "white", stroke = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.8) +
  labs(
    title = "Executive Function at Baseline and Follow-up",
    subtitle = "sRBD vs. Controls",
    x = "Timepoint",
    y = "Mean Executive Function",
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
    plot.title = element_text(face = "bold", size = 14)
  )



#### check linearity of age
#ART
ggplot(df_filtered2, aes(x = age_session_baseline, y = tmt_trailA_reactiontime, color = factor(Code))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE) +  # loess to explore shape
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), linetype = "dashed", se = FALSE) + # quadratic fit
  theme_minimal() +
  labs(
    x = "Baseline Age",
    y = "Trail A Reaction Time",
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


RTA_model2 <- lmer(
  TrailA_negZ ~ elapsed_time_z * Code + age_baseline_z + age_squared_z + gender + education_level + (1 | PID2),
  data = df_filtered2
)
summary(RTA_model2)
# create a summary table with CIs
RTA2 <- tidy(RTA_model2, effects = "fixed", conf.int = TRUE) %>%
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
print(RTA2)
anova(ART_model, RTA_model2)


#BRT
ggplot(df_filtered2, aes(x = age_session_baseline, y = tmt_trailB_reactiontime, color = factor(Code))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE) +  # loess to explore shape
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), linetype = "dashed", se = FALSE) + # quadratic fit
  theme_minimal() +
  labs(
    x = "Baseline Age",
    y = "Trail B Reaction Time",
    color = "Code"
  ) 


BRT_model2 <- lmer(
  TrailB_negZ ~ elapsed_time_z * Code + age_baseline_z + age_squared_z + gender + education_level + (1 | PID2),
  data = df_filtered2
)
summary(BRT_model2)
# create a summary table with CIs
BRT2 <- tidy(BRT_model2, effects = "fixed", conf.int = TRUE) %>%
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
print(BRT2)
anova(BRT_model, BRT_model2)



#EF
ggplot(df_filtered2, aes(x = age_session_baseline, y = tmt_trailAminustrailB_exefcn_final, color = factor(Code))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE) +  # loess to explore shape
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), linetype = "dashed", se = FALSE) + # quadratic fit
  theme_minimal() +
  labs(
    x = "Baseline Age",
    y = "Executive Function",
    color = "Code"
  ) 


EF_model2 <- lmer(
  tmt_trailAminustrailB_exefcn_final ~ elapsed_time_z * Code + age_baseline_z + age_squared_z + gender + education_level + (1 | PID2),
  data = df_filtered2
)
summary(EF_model2)
# create a summary table with CIs
EF2 <- tidy(EF_model2, effects = "fixed", conf.int = TRUE) %>%
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
print(EF2)
anova(EF_model, EF_model2)



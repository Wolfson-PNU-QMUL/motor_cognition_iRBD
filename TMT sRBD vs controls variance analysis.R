######## TMT sRBD vs controls longitudinal variance analysis ########
library(markdown)
library(data.table)
library(tidyverse)
library(readxl)
library(writexl)
library(e1071)
library(gridExtra)
library(broom.mixed)
library(boot)
library(readr)
library(ggdist) 
library(nlme)
library(lme4)
library(lmerTest)
library(nlme)
library(tidybayes)
library(moments)

getwd()

sink()

pdf()


#read in dataset
View(matched_sRBD_control_TMT_data)


table(matched_sRBD_control_TMT_data$PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline)
#rename the 'yes' 
matched_sRBD_control_TMT_data <- matched_sRBD_control_TMT_data %>%
  mutate(PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline = recode(
    PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline,
    "Yes, Survey + SurveyBRAIN + SurveySmell + SurveyBRAINSmell + SurveyAll" = "Yes",
    "No" = "No"
  ))


matched_sRBD_control_TMT_data %>%
  filter(Code == "Control", rx_anti_psychotics_use_final_ever != "Unknown") %>%
  count(rx_anti_psychotics_use_final_ever)

matched_sRBD_control_TMT_data %>%
  filter(Code == "Control", PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline != "Unknown") %>%
  count(PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline)

# Define the columns of interest
tmt_cols <- c(
  "tmt_trailA_reactiontime_4", "tmt_trailA_reactiontime_5", "tmt_trailA_reactiontime_7", "tmt_trailA_reactiontime_8",
  "tmt_trailB_reactiontime_4", "tmt_trailB_reactiontime_5", "tmt_trailB_reactiontime_7", "tmt_trailB_reactiontime_8",
  "tmt_trailBminustrailA_exefcn_final_4", "tmt_trailBminustrailA_exefcn_final_5", "tmt_trailBminustrailA_exefcn_final_7",
  "tmt_trailBminustrailA_exefcn_final_8"
)

# Remove participants with all test data missing
matched_sRBD_control_TMT_data2 <- matched_sRBD_control_TMT_data %>%
  filter(
    rowSums(is.na(select(., all_of(tmt_cols)))) <
      ncol(select(., all_of(tmt_cols)))
  )


#save file


######## Exploratory data analysis #########
#check columns
names(matched_sRBD_control_TMT_data2)[!sapply(matched_sRBD_control_TMT_data2, is.numeric)]

#make columns numeric
cols_to_convert <- c(
  "tmt_trailBminustrailA_exefcn_final_4", "tmt_trailBminustrailA_exefcn_final_5",
  "tmt_trailBminustrailA_exefcn_final_7", "tmt_trailBminustrailA_exefcn_final_8",
  "tmt_trailB_reactiontime_4", "tmt_trailB_reactiontime_5",
  "tmt_trailB_reactiontime_7", "tmt_trailB_reactiontime_8",
  "tmt_trailA_reactiontime_4", "tmt_trailA_reactiontime_5",
  "tmt_trailA_reactiontime_7",
  "tmt_trailA_reactiontime_8"
)

matched_sRBD_control_TMT_data2[cols_to_convert] <- lapply(matched_sRBD_control_TMT_data2[cols_to_convert], as.numeric)


#remove participants with ALL brain test data missing 
matched_sRBD_control_TMT_data2 <- matched_sRBD_control_TMT_data2 %>%
  filter(
    rowSums(is.na(select(., all_of(cols_to_convert)))) <
      ncol(select(., all_of(cols_to_convert)))
  )


# Pivot to long format but keep all non-brain test columns
sRBD_TMT_long <- matched_sRBD_control_TMT_data2 %>%
  pivot_longer(
    cols = all_of(cols_to_convert),
    names_to = c("Measure", "Session"),
    names_pattern = "(.*)_(\\d+)$",
    values_to = "Value"
  ) %>%
  mutate(Session = as.numeric(Session)) %>%  # ensure numeric, not factor
  group_by(PID2, Measure) %>%
  arrange(Session, .by_group = TRUE) %>%
  mutate(RelativeVisit = row_number()) %>%  # assign 1,2,3,... explicitly
  ungroup()


# check how many timepoints
sRBD_TMT_long %>%
  group_by(Measure, PID2) %>%
  summarise(n = n(), .groups = "drop") %>%
  count(Measure, n)

sRBD_TMT_long %>%
  filter(!is.na(Value)) %>%                       # only observed values
  group_by(Measure, PID2) %>%
  summarise(n_obs = n(), .groups = "drop") %>%    # number of observed timepoints
  filter(n_obs == 4) %>%                          # complete 4
  count(Measure, name = "n_IDs_complete_4tp")


#Variability-eligible dataset (maximising sample)
sRBD_TMT_long_all <- sRBD_TMT_long %>%
  filter(!is.na(Value)) %>%
  group_by(PID2, Code, Measure) %>%
  filter(n() >= 2) %>%                 
  arrange(Session, .by_group = TRUE) %>%
  mutate(RelativeVisit = row_number()) %>%   # recompute visit order
  ungroup()

sRBD_TMT_long_all %>%
  group_by(Measure, PID2) %>%
  summarise(n = n(), .groups = "drop") %>%
  count(Measure, n)


#check whether each participant has the same number of observations for TMT-A and TMT-B
sRBD_TMT_long_all %>%
  filter(Measure %in% c(
    "tmt_trailA_reactiontime",
    "tmt_trailB_reactiontime"
  )) %>%
  group_by(PID2, Measure) %>%
  summarise(n_visits = n_distinct(RelativeVisit), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = Measure,
    values_from = n_visits
  ) %>%
  count(
    tmt_trailA_reactiontime,
    tmt_trailB_reactiontime
  )

#check EF against A and B reacton times
sRBD_TMT_long_all %>%
  filter(Measure == "tmt_trailBminustrailA_exefcn_final") %>%
  group_by(PID2) %>%
  summarise(
    n_EF = n_distinct(RelativeVisit),
    .groups = "drop"
  ) %>%
  count(n_EF)

#check RT data
sRBD_TMT_long_all %>%
  filter(Measure %in% c(
    "tmt_trailA_reactiontime",
    "tmt_trailB_reactiontime"
  )) %>%
  group_by(Measure, Code) %>%
  summarise(
    n = sum(!is.na(Value)),
    min_value = min(Value, na.rm = TRUE),
    n_negative = sum(Value < 0, na.rm = TRUE),
    n_zero = sum(Value == 0, na.rm = TRUE),
    n_missing = sum(is.na(Value)),
    .groups = "drop"
  )
# Zero values in control group

#Identify participants and visits with zero values
sRBD_TMT_long_all %>%
  filter(
    Measure %in% c(
      "tmt_trailA_reactiontime",
      "tmt_trailB_reactiontime"
    ),
    Value == 0
  ) %>%
  select(PID2, TokenFinal, Code, RelativeVisit, Measure, Value)

sRBD_TMT_long %>%
  filter(
    Measure %in% c(
      "tmt_trailA_reactiontime",
      "tmt_trailB_reactiontime"
    ),
    Value == 0
  ) %>%
  select(PID2, TokenFinal, Code, RelativeVisit, Measure, Value)

matched_sRBD_control_TMT_data2 %>%
  filter(PID2 %in% c(4390, 15625, 21311)) %>%
  select(
    PID2,
    starts_with("tmt_trailA_reactiontime"),
    starts_with("tmt_trailB_reactiontime")
  ) %>%
  print(width = Inf)

# check which visits contain zeros
sRBD_TMT_long %>%
  filter(
    PID2 %in% c(4390, 15625, 21311),
    Measure %in% c(
      "tmt_trailA_reactiontime",
      "tmt_trailB_reactiontime"
    )
  ) %>%
  filter(Value == 0) %>%
  select(PID2, RelativeVisit, Measure, Value)


# QC for TMT-A and TMT-B
apply_tmt_qc <- function(sRBD_TMT_long_all,
                         value_col = "Value") {
  
  sRBD_TMT_long_all %>%
    mutate(
      Value_num = suppressWarnings(as.numeric(.data[[value_col]])),
      
      qc_flag = case_when(
        is.na(Value_num) ~ "missing",
        Measure %in% c(
          "tmt_trailA_reactiontime",
          "tmt_trailB_reactiontime"
        ) & Value_num <= 0 ~ "invalid_nonpositive",
        TRUE ~ "ok"
      ),
      
      Value_clean = case_when(
        qc_flag == "invalid_nonpositive" ~ NA_real_,
        TRUE ~ Value_num
      )
    )
}

sRBD_TMT_long_qc <- apply_tmt_qc(sRBD_TMT_long_all)


sRBD_TMT_long_qc %>%
  count(Measure, qc_flag)

sRBD_TMT_long_qc %>%
  filter(qc_flag == "invalid_nonpositive") %>%
  select(PID2, Code, RelativeVisit, Measure, Value, Value_clean)


# Raincloud plot
# Horizontal raincloud plots per visit
visits <- sort(unique(sRBD_TMT_long_qc$RelativeVisit))

for (v in visits) {
  p <- sRBD_TMT_long_qc %>%
    filter(RelativeVisit == v) %>%
    ggplot(aes(y = Code, x = Value_clean, fill = Code, color = Code)) +  # <- changed
    stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, orientation = "horizontal") +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, position = position_nudge(y = 0.3)) +
    geom_jitter(height = 0.1, size = 1.5, alpha = 0.8) +
    facet_wrap(~Measure, scales = "free_x") +
    scale_fill_manual(values = c("Control" = "#1f3552", "sRBD" = "#69b3a2")) +
    scale_color_manual(values = c("Control" = "#1f3552", "sRBD" = "#69b3a2")) +
    labs(
      y = "Group",
      x = "Score",
      title = paste("Raincloud Plot of TMT Measures – Visit", v)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "top"
    )
  
  print(p)
}



#check for outliers
sRBD_TMT_flagged <- sRBD_TMT_long_qc %>%
  group_by(Code, Measure, RelativeVisit) %>%
  mutate(
    Q1 = quantile(Value_clean, 0.25, na.rm = TRUE),
    Q3 = quantile(Value_clean, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    Lower = Q1 - 1.5 * IQR,
    Upper = Q3 + 1.5 * IQR,
    is_outlier = Value_clean < Lower | Value_clean > Upper
  ) %>%
  ungroup()


# Flag outliers and create summary table in one go
outlier_counts <- sRBD_TMT_long_qc %>%
  group_by(Code, Measure, RelativeVisit) %>%
  mutate(
    Q1 = quantile(Value_clean, 0.25, na.rm = TRUE),
    Q3 = quantile(Value_clean, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    Lower = Q1 - 1.5 * IQR,
    Upper = Q3 + 1.5 * IQR,
    is_outlier = Value_clean < Lower | Value_clean > Upper
  ) %>%
  ungroup() %>%
  group_by(Code) %>%
  summarise(
    total_outliers = sum(is_outlier, na.rm = TRUE),
    total_points = n(),
    perc_outliers = round(100 * total_outliers / total_points, 2)
  )

outlier_counts

# Identify outliers by Code, measure, and year
outlier_data <- sRBD_TMT_long_qc %>%
  group_by(Code, Measure, RelativeVisit) %>%
  mutate(
    Q1 = quantile(Value_clean, 0.25, na.rm = TRUE),
    Q3 = quantile(Value_clean, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    Lower = Q1 - 1.5 * IQR,
    Upper = Q3 + 1.5 * IQR,
    is_outlier = Value_clean < Lower | Value_clean > Upper
  ) %>%
  ungroup()

# Boxplot highlighting outliers
ggplot(outlier_data, aes(x = factor(RelativeVisit), y = Value_clean, fill = Code)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, outlier.color = "red", alpha = 0.7) +
  facet_wrap(~Measure, scales = "free_y") +
  labs(
    title = "Distribution of TMT Scores by Year and Group",
    x = "Year",
    y = "Score Value",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


sRBD_TMT_long_qc <- sRBD_TMT_long_qc %>%
  mutate(RelativeVisit = as.numeric(RelativeVisit))


#check if data is normally distributed 
ggplot(sRBD_TMT_long_qc %>% filter(Measure == "tmt_trailA_reactiontime"), 
       aes(x = Value_clean, fill = Code)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Distribution of Trail A Reaction Time by Group",
    x = "Trail A Reaction Time",
    y = "Density",
    fill = "Group"
  ) +
  xlim(0,3500) +
  theme_minimal(base_size = 14)


ggplot(sRBD_TMT_long_qc %>% filter(Measure == "tmt_trailB_reactiontime"), 
       aes(x = Value_clean, fill = Code)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Distribution of Trail B Reaction Time by Group",
    x = "Trail B Reaction Time",
    y = "Density",
    fill = "Group"
  ) +
  xlim(0,6500) +
  theme_minimal(base_size = 14)


ggplot(sRBD_TMT_long_qc %>% filter(Measure == "tmt_trailBminustrailA_exefcn_final"), 
       aes(x = Value_clean, fill = Code)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Distribution of Executive Function by Group",
    x = "Executive Function",
    y = "Density",
    fill = "Group"
  ) +
  xlim(-2000,4500) +
  theme_minimal(base_size = 14)


#test for normality
normality_results <- sRBD_TMT_long_qc %>%
  group_by(Measure, PID2) %>%       # assume ID is your participant identifier
  summarise(Value_clean = mean(Value_clean, na.rm = TRUE), .groups = "drop") %>%  # one value per person
  group_by(Measure) %>%
  summarise(
    n = sum(!is.na(Value_clean)),
    p_value = tryCatch(
      if (n >= 3) shapiro.test(Value_clean)$p.value else NA_real_,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(p_value) ~ "Not tested (too few values)",
      p_value > 0.05 ~ "Likely normal",
      TRUE ~ "Non-normal"
    )
  )

normality_results



#inspect visually
# Loop through each measure
for(m in unique(sRBD_TMT_long_qc$Measure)) {
  
  data_m <- sRBD_TMT_long_qc %>%
    filter(Measure == m, !is.na(Value_clean))
  
  cat(m, "skewness:", round(skewness(data_m$Value_clean), 2), "\n")
  
  gridExtra::grid.arrange(
    ggplot(data_m, aes(x = Value_clean)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "#69b3a2", alpha = 0.6) +
      geom_density(color = "#1f3552", linewidth = 1) +
      labs(title = paste(m, "- Histogram & Density"), x = "Value", y = "Density") +
      theme_minimal(base_size = 12),
    
    ggplot(data_m, aes(sample = Value_clean)) +
      stat_qq() + stat_qq_line(color = "red") +
      labs(title = paste(m, "- QQ Plot"), x = "Theoretical Quantiles", y = "Sample Quantiles") +
      theme_minimal(base_size = 12),
    
    ncol = 2
  )
}


# Log transform TMT-A and TMT-B
sRBD_TMT_long_qc <- sRBD_TMT_long_qc %>%
  mutate(
    Value_log = if_else(
      Measure %in% c(
        "tmt_trailA_reactiontime",
        "tmt_trailB_reactiontime"
      ),
      log(if_else(
        Measure %in% c(
          "tmt_trailA_reactiontime",
          "tmt_trailB_reactiontime"
        ),
        Value_clean,
        1
      )),
      NA_real_
    )
  )


#check
sRBD_TMT_long_qc %>%
  filter(Measure %in% c(
    "tmt_trailA_reactiontime",
    "tmt_trailB_reactiontime"
  )) %>%
  group_by(Measure) %>%
  summarise(
    n = sum(!is.na(Value_log)),
    min = min(Value_log, na.rm = TRUE),
    max = max(Value_log, na.rm = TRUE),
    mean = mean(Value_log, na.rm = TRUE),
    .groups = "drop"
  )


#check each participant has 1 TMT-A and TMT-B
sRBD_TMT_long_qc %>%
  filter(Measure %in% c(
    "tmt_trailA_reactiontime",
    "tmt_trailB_reactiontime"
  )) %>%
  count(PID2, RelativeVisit, Measure) %>%
  summarise(
    duplicate_combinations = sum(n > 1)
  )


# Calculate new log-transformed EF
TMT_EF <- sRBD_TMT_long_qc %>%
  filter(Measure %in% c(
    "tmt_trailA_reactiontime",
    "tmt_trailB_reactiontime"
  )) %>%
  select(
    PID2, TokenFinal, RelativeVisit, Code,
    age_session_baseline, gender,
    Measure, Value_log
  ) %>%
  tidyr::pivot_wider(
    names_from = Measure,
    values_from = Value_log
  ) %>%
  mutate(
    Value_log_EF = tmt_trailB_reactiontime -
      tmt_trailA_reactiontime
  ) %>%
  filter(
    !is.na(tmt_trailA_reactiontime),
    !is.na(tmt_trailB_reactiontime)
  )

# Add new log-transformed EF to long dataset
TMT_log <- TMT_EF %>%
  transmute(
    PID2,
    TokenFinal,
    RelativeVisit,
    Code,
    age_session_baseline,
    gender,
    Measure = "TMT_EF_log",
    Value = Value_log_EF
  )

# Rebuild final long dataset
sRBD_TMT_long_final <- bind_rows(
  sRBD_TMT_long_qc,
  TMT_log
)

# Check measure counts
sRBD_TMT_long_final %>%
  count(Measure)

# Check logged values
sRBD_TMT_long_final %>%
  filter(Measure %in% c(
    "tmt_trailA_reactiontime",
    "tmt_trailB_reactiontime",
    "TMT_EF_log"
  )) %>%
  group_by(Measure) %>%
  summarise(
    n_logged = sum(
      if_else(Measure == "TMT_EF_log", !is.na(Value), !is.na(Value_log))
    ),
    .groups = "drop"
  )


#check EF corresponds to participants with both valid RTs
TMT_EF %>%
  summarise(
    n_EF = n(),
    n_missing_EF = sum(is.na(Value_log_EF)),
    min_EF = min(Value_log_EF, na.rm = TRUE),
    max_EF = max(Value_log_EF, na.rm = TRUE)
  )

# Check age and gender for EF
sRBD_TMT_long_final %>%
  filter(Measure == "TMT_EF_log") %>%
  summarise(
    n = n(),
    missing_age = sum(is.na(age_session_baseline)),
    missing_gender = sum(is.na(gender))
  )


######### Variance Analysis ######### 
#create column for logged values
sRBD_TMT_long_final <- sRBD_TMT_long_final %>%
  mutate(
    Value_for_analysis = case_when(
      Measure %in% c(
        "tmt_trailA_reactiontime",
        "tmt_trailB_reactiontime"
      ) ~ Value_log,
      Measure == "TMT_EF_log" ~ Value,
      TRUE ~ NA_real_
    )
  )

#Step 1a — Compute mean, SD, and sample size per group × year × measure
variance_summary <- sRBD_TMT_long_final %>%
  filter(Measure != "tmt_trailBminustrailA_exefcn_final") %>%
  group_by(Measure, RelativeVisit, Code) %>%
  summarise(
    n = sum(!is.na(Value_for_analysis)),
    mean_value = mean(Value_for_analysis, na.rm = TRUE),
    sd_value = sd(Value_for_analysis, na.rm = TRUE),
    .groups = "drop"
  )

print(variance_summary)


#Step 1b — Bootstrap SD for confidence intervals (small sample size)
# Function to compute SD for bootstrapping
sd_fun <- function(data, indices) {
  return(sd(data[indices], na.rm = TRUE))
}

# Bootstrap SD function per group
bootstrap_sd_ci <- function(x, R = 1000, conf = 0.95) {
  boot_obj <- boot(data = x, statistic = sd_fun, R = R)
  ci <- boot.ci(boot_obj, type = "perc")$percent[4:5]  # lower and upper CI
  return(data.frame(sd = sd(x, na.rm = TRUE),
                    sd_lower = ci[1],
                    sd_upper = ci[2]))
}

# Apply to all groups
variance_boot <- sRBD_TMT_long_final %>%
  filter(Measure != "tmt_trailBminustrailA_exefcn_final") %>%
  group_by(Measure, RelativeVisit, Code) %>%
  summarise(
    n = sum(!is.na(Value_for_analysis)),
    boot_results = list(
      bootstrap_sd_ci(Value_for_analysis[!is.na(Value_for_analysis)])
    ),
    .groups = "drop"
  ) %>%
  unnest(cols = c(boot_results))

print(variance_boot, n = Inf, width = Inf)

#save file

#Step 1c — Plot mean ± SD with bootstrap CI
# Combine mean and bootstrap SD
plot_data <- variance_summary %>%
  left_join(variance_boot, by = c("Measure", "RelativeVisit", "Code"))

# Get unique measures
measures <- unique(plot_data$Measure)


# Named vector for measure labels
measure_labels <- c(
  "tmt_trailA_reactiontime" = "TMT-A",
  "tmt_trailB_reactiontime" = "TMT-B",
  "TMT_EF_log"              = "TMT Executive Function"
)

# Directionality info for subtitle
measure_direction <- c(
  "tmt_trailA_reactiontime" = "Higher = worse",
  "tmt_trailB_reactiontime" = "Higher = worse",
  "TMT_EF_log"              = "Higher = worse"
)

plot_data <- plot_data %>%
  mutate(Code = ifelse(Code == "sRBD", "Normosmic pRBD", Code))


# Loop through measures and plot with subtitle
for (m in measures) {
  
  p <- plot_data %>%
    filter(Measure == m) %>%
    ggplot(aes(x = RelativeVisit, y = mean_value, color = Code, group = Code, fill = Code)) +
    geom_ribbon(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(
      title = paste0("TMT test – ", measure_labels[m], " Over Time (Mean ± SD)"),
      subtitle = measure_direction[m],
      x = "Visit",
      y = "Mean score (log-transformed)",
      color = "Group",
      fill = "Group",
      caption = "Lines = group mean; shaded area = ±1 SD"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, face = "italic"),
      legend.position = "top"
    )
  
  print(p)
}

# Save the combined mean + SD + bootstrap CI

#Step 2: Model time-varying variance (heteroscedastic mixed model)
#make visit variable factor and control group the reference
sRBD_TMT_long_final <- sRBD_TMT_long_final %>%
  mutate(
    Code = relevel(factor(Code), ref = "Control"),
    RelativeVisit_f = factor(RelativeVisit)
  )



# Filter TMT-A scores
TMTA_data <- sRBD_TMT_long_final %>% 
  filter(Measure == "tmt_trailA_reactiontime")


# Visit-specific variance model
m_visit_TMTA <- nlme::lme(
  Value_for_analysis ~ Code * RelativeVisit + age_session_baseline + gender,
  random = ~ 1 | PID2,
  data = TMTA_data,
  na.action = na.exclude,
  weights = nlme::varIdent(form = ~ 1 | RelativeVisit_f)
)

#variance by visit and group (add one variance term for group)
m_group_TMTA <- update(
  m_visit_TMTA,
  weights = nlme::varComb(
    nlme::varIdent(form = ~ 1 | RelativeVisit_f),
    nlme::varIdent(form = ~ 1 | Code)
  )
)

anova(m_visit_TMTA, m_group_TMTA) #does variance differ between groups?

#direction of size and difference (95% CI & variance ratio)
ci_TMTA <- intervals(m_group_TMTA, which = "var-cov")$varStruct
ci_TMTA     # standard deviation ratio, cases vs controls, with 95% CI
ci_TMTA^2   # variance ratio, cases vs controls, with 95% CI


# Step 3: calculate residual-based MSSD
# Calculate residuals for TMTA
TMTA_data$resid_fixed <- TMTA_data$Value_for_analysis -
  fitted(m_visit_TMTA, level = 0)

# calculate residual-based MSSD for each participant
mssd_TMTA <- TMTA_data %>%
  arrange(PID2, RelativeVisit) %>%
  group_by(PID2, Code) %>%
  summarise(
    n_timepoints = sum(!is.na(resid_fixed)),
    mssd = if (sum(!is.na(resid_fixed)) >= 2) {
      mean(diff(resid_fixed[!is.na(resid_fixed)])^2)
    } else {
      NA_real_
    },
    .groups = "drop"
  )
#this calculates differences between consecutive observed visits

# Compare MSSD between groups
wilcox.test(mssd ~ Code, data = mssd_TMTA, exact = FALSE)


mssd_TMTA %>% count(Code, is.na(mssd))


# participants per group at each visit
TMTA_data %>%
  filter(!is.na(Value_for_analysis)) %>%
  with(table(Code, RelativeVisit))



#TMT-B reaction time
unique(sRBD_TMT_long_final$Measure)
TMTB_data <- sRBD_TMT_long_final %>% 
  filter(Measure == "tmt_trailB_reactiontime")


# Visit-specific variance model
m_visit_TMTB <- nlme::lme(
  Value_for_analysis ~ Code * RelativeVisit + age_session_baseline + gender,
  random = ~ 1 | PID2,
  data = TMTB_data,
  na.action = na.exclude,
  weights = nlme::varIdent(form = ~ 1 | RelativeVisit_f)
)

#variance by visit and group (add one variance term for group)
m_group_TMTB <- update(
  m_visit_TMTB,
  weights = nlme::varComb(
    nlme::varIdent(form = ~ 1 | RelativeVisit_f),
    nlme::varIdent(form = ~ 1 | Code)
  )
)

anova(m_visit_TMTB, m_group_TMTB) #does variance differ between groups?

#direction of size and difference (95% CI & variance ratio)
ci_TMTB <- intervals(m_group_TMTB, which = "var-cov")$varStruct
ci_TMTB     # standard deviation ratio, cases vs controls, with 95% CI
ci_TMTB^2   # variance ratio, cases vs controls, with 95% CI


# Step 3: calculate residual-based MSSD
# Calculate residuals for TMTB
TMTB_data$resid_fixed <- TMTB_data$Value_for_analysis -
  fitted(m_visit_TMTB, level = 0)

# calculate residual-based MSSD for each participant
mssd_TMTB <- TMTB_data %>%
  arrange(PID2, RelativeVisit) %>%
  group_by(PID2, Code) %>%
  summarise(
    n_timepoints = sum(!is.na(resid_fixed)),
    mssd = if (sum(!is.na(resid_fixed)) >= 2) {
      mean(diff(resid_fixed[!is.na(resid_fixed)])^2)
    } else {
      NA_real_
    },
    .groups = "drop"
  )
#this calculates differences between consecutive observed visits

# Compare MSSD between groups
wilcox.test(mssd ~ Code, data = mssd_TMTB, exact = FALSE)


mssd_TMTB %>% count(Code, is.na(mssd))


# participants per group at each visit
TMTB_data %>%
  filter(!is.na(Value_for_analysis)) %>%
  with(table(Code, RelativeVisit))



# Executive Function
unique(sRBD_TMT_long_final$Measure)
TMTEF_data <- sRBD_TMT_long_final %>% 
  filter(Measure == "TMT_EF_log")


# Visit-specific variance model
m_visit_TMTEF <- nlme::lme(
  Value_for_analysis ~ Code * RelativeVisit + age_session_baseline + gender,
  random = ~ 1 | PID2,
  data = TMTEF_data,
  na.action = na.exclude,
  weights = nlme::varIdent(form = ~ 1 | RelativeVisit_f)
)

#variance by visit and group (add one variance term for group)
m_group_TMTEF <- update(
  m_visit_TMTEF,
  weights = nlme::varComb(
    nlme::varIdent(form = ~ 1 | RelativeVisit_f),
    nlme::varIdent(form = ~ 1 | Code)
  )
)

anova(m_visit_TMTEF, m_group_TMTEF) #does variance differ between groups?

#direction of size and difference (95% CI & variance ratio)
ci_EF <- intervals(m_group_TMTEF, which = "var-cov")$varStruct
ci_EF     # standard deviation ratio, cases vs controls, with 95% CI
ci_EF^2   # variance ratio, cases vs controls, with 95% CI


# Step 3: calculate residual-based MSSD
# Calculate residuals for EF
TMTEF_data$resid_fixed <- TMTEF_data$Value_for_analysis -
  fitted(m_visit_TMTEF, level = 0)

# calculate residual-based MSSD for each participant
mssd_TMTEF <- TMTEF_data %>%
  arrange(PID2, RelativeVisit) %>%
  group_by(PID2, Code) %>%
  summarise(
    n_timepoints = sum(!is.na(resid_fixed)),
    mssd = if (sum(!is.na(resid_fixed)) >= 2) {
      mean(diff(resid_fixed[!is.na(resid_fixed)])^2)
    } else {
      NA_real_
    },
    .groups = "drop"
  )
#this calculates differences between consecutive observed visits

# Compare MSSD between groups
wilcox.test(mssd ~ Code, data = mssd_TMTEF, exact = FALSE)


mssd_TMTEF %>% count(Code, is.na(mssd))


# participants per group at each visit
table(TMTEF_data$Code[!is.na(TMTEF_data$Value_for_analysis)],
      TMTEF_data$RelativeVisit)


###### demographics ###### 
sRBD_TMT_long_all %>%
  distinct(PID2, Code) %>%
  count(Code)

#keep one row per participant for demographics
demo_df <- sRBD_TMT_long_all %>%
  arrange(PID2, RelativeVisit) %>% 
  distinct(PID2, .keep_all = TRUE)

t_test_age_result <- t.test(age_session_baseline ~ Code, data = demo_df)
print(t_test_age_result)

std_dev_ages <- demo_df %>% group_by(Code) %>% 
  summarize(std_dev_ages = sd(age_session_baseline, na.rm = TRUE))
print(std_dev_ages)

sex_percentage <- demo_df %>%
  count(Code, gender) %>%
  group_by(Code) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()
print(sex_percentage)

chi_test_result_sex <- demo_df %>%
  with(table(gender, Code)) %>%
  chisq.test()
print(chi_test_result_sex)


tmt_skewness <- sRBD_TMT_long_qc %>%
  filter(Measure %in% c(
    "tmt_trailA_reactiontime",
    "tmt_trailB_reactiontime",
    "tmt_trailBminustrailA_exefcn_final"
  )) %>%
  filter(!is.na(Value_clean)) %>%
  group_by(Measure, Code) %>%
  summarise(
    skewness = round(skewness(Value_clean), 2),
    .groups = "drop"
  )

print(tmt_skewness)

summary_means <- sRBD_TMT_long_qc %>%
  filter(Measure %in% c(
    "tmt_trailA_reactiontime",
    "tmt_trailB_reactiontime",
    "tmt_trailBminustrailA_exefcn_final"
  )) %>%
  filter(!is.na(Value_clean)) %>%
  group_by(Measure, RelativeVisit) %>%
  summarise(
    n_control = sum(Code == "Control"),
    median_control = median(Value_clean[Code == "Control"]),
    IQR_control = IQR(Value_clean[Code == "Control"]),
    
    n_sRBD = sum(Code == "sRBD"),
    median_sRBD = median(Value_clean[Code == "sRBD"]),
    IQR_sRBD = IQR(Value_clean[Code == "sRBD"]),
    
    p_value = wilcox.test(
      Value_clean[Code == "Control"],
      Value_clean[Code == "sRBD"],
      exact = FALSE
    )$p.value,
    
    .groups = "drop"
  ) %>%
  mutate(
    control_median_IQR = sprintf("%.2f (%.2f)", median_control, IQR_control),
    sRBD_median_IQR = sprintf("%.2f (%.2f)", median_sRBD, IQR_sRBD)
  ) %>%
  select(
    Measure, RelativeVisit,
    n_control, control_median_IQR,
    n_sRBD, sRBD_median_IQR,
    p_value
  )

print(summary_means, n = Inf, width = Inf)



# Close PDF
dev.off()

# Close console output
sink()

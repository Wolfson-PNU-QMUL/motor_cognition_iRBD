######## DFT RBD vs controls longitudinal variance analysis ########
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
View(matched_RBD_control_DFT_data)

table(matched_RBD_control_DFT_data$PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline)
#rename the 'yes' 
matched_RBD_control_DFT_data <- matched_RBD_control_DFT_data %>%
  mutate(PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline = recode(
    PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline,
    "Yes, Survey + SurveyBRAIN + SurveySmell + SurveyBRAINSmell + SurveyAll" = "Yes",
    "No" = "No"
  ))


matched_RBD_control_DFT_data %>%
  filter(Code == "Control", rx_anti_psychotics_use_final_ever != "Unknown") %>%
  count(rx_anti_psychotics_use_final_ever)

matched_RBD_control_DFT_data %>%
  filter(Code == "Control", PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline != "Unknown") %>%
  count(PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline)


# Define the columns of interest
dft_cols <- c(
  "onefingertaptest_worst_is_5", "onefingertaptest_worst_is_7",
  "onefingertaptest_worst_at_5", "onefingertaptest_worst_at_7",
  "onefingertaptest_worst_ks_5", "onefingertaptest_worst_ks_7"
)

# Remove participants with all test data missing
matched_RBD_control_DFT_data2 <- matched_RBD_control_DFT_data %>%
  filter(
    rowSums(is.na(select(., all_of(dft_cols)))) <
      ncol(select(., all_of(dft_cols)))
  )


#save file

######## Exploratory data analysis #########
# List of your BRAIN test columns
#check columns
names(matched_RBD_control_DFT_data2)[!sapply(matched_RBD_control_DFT_data2, is.numeric)]

#make columns numeric
cols_to_convert <- c( "onefingertaptest_worst_is_5", "onefingertaptest_worst_is_7",
                      "onefingertaptest_worst_ks_5", "onefingertaptest_worst_ks_7", 
                      "onefingertaptest_worst_at_5", "onefingertaptest_worst_at_7")

matched_RBD_control_DFT_data2[cols_to_convert] <- lapply(matched_RBD_control_DFT_data2[cols_to_convert], as.numeric)


#remove participants with ALL brain test data missing 
matched_RBD_control_DFT_data2 <- matched_RBD_control_DFT_data2 %>%
  filter(
    rowSums(is.na(select(., starts_with("onefingertaptest_worst")))) <
      ncol(select(., starts_with("onefingertaptest_worst")))
  )


# Pivot to long format but keep all non-brain test columns
RBD_DFT_long <- matched_RBD_control_DFT_data2 %>%
  pivot_longer(
    cols = all_of(dft_cols),                # only pivot DFT columns
    names_to = c("Measure", "Session"),
    names_pattern = "(.*)_(\\d+)$",
    values_to = "Value"
  ) %>%
  mutate(Session = as.numeric(Session)) %>%
  group_by(PID2, Measure) %>%             # define per-participant & per-measure
  arrange(PID2, Measure, Session) %>%
  mutate(RelativeVisit = row_number()) %>% # baseline = 1, next = 2, etc.
  ungroup()


# check how many timepoints
RBD_DFT_long %>%
  group_by(Measure, PID2) %>%
  summarise(n = n(), .groups = "drop") %>%
  count(Measure, n)

RBD_DFT_long %>%
  filter(!is.na(Value)) %>%
  group_by(Measure, PID2) %>%
  summarise(n_obs = n(), .groups = "drop") %>%
  filter(n_obs == 2) %>%
  count(Measure, name = "n_IDs_with_2_tp")


#Variability-eligible dataset (maximising sample)
RBD_DFT_long_2tp <- RBD_DFT_long %>%
  filter(!is.na(Value)) %>%
  group_by(PID2, Code, Measure) %>%
  filter(n() >= 2) %>%                 
  arrange(Session, .by_group = TRUE) %>%
  mutate(RelativeVisit = row_number()) %>%   # recompute visit order
  ungroup()

RBD_DFT_long_2tp %>%
  group_by(Measure, PID2) %>%
  summarise(n = n(), .groups = "drop") %>%
  count(Measure, n)



# Raincloud plot
# Horizontal raincloud plots per visit
visits <- sort(unique(RBD_DFT_long_2tp $RelativeVisit))

for (v in visits) {
  p <- RBD_DFT_long_2tp  %>%
    filter(RelativeVisit == v) %>%
    ggplot(aes(y = Code, x = Value, fill = Code, color = Code)) +  # swapped x/y
    stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, orientation = "horizontal") +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, position = position_nudge(y = 0.3)) +
    geom_jitter(height = 0.1, size = 1.5, alpha = 0.8) +
    facet_wrap(~Measure, scales = "free_x") +  # flipped facet axis
    scale_fill_manual(values = c("Control" = "#1f3552", "iRBD" = "#69b3a2")) +
    scale_color_manual(values = c("Control" = "#1f3552", "iRBD" = "#69b3a2")) +
    labs(
      y = "Group",
      x = "Score",
      title = paste("Raincloud Plot of DFT Measures – Visit", v)
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
RBD_DFT_flagged <- RBD_DFT_long_2tp %>%
  group_by(Code, Measure, RelativeVisit) %>%
  mutate(
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    Lower = Q1 - 1.5 * IQR,
    Upper = Q3 + 1.5 * IQR,
    is_outlier = Value < Lower | Value > Upper
  ) %>%
  ungroup()


# Flag outliers and create summary table in one go
outlier_counts <- RBD_DFT_long_2tp %>%
  group_by(Code, Measure, RelativeVisit) %>%
  mutate(
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    Lower = Q1 - 1.5 * IQR,
    Upper = Q3 + 1.5 * IQR,
    is_outlier = Value < Lower | Value > Upper
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
outlier_data <- RBD_DFT_long_2tp %>%
  group_by(Code, Measure, RelativeVisit) %>%
  mutate(
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    Lower = Q1 - 1.5 * IQR,
    Upper = Q3 + 1.5 * IQR,
    is_outlier = Value < Lower | Value > Upper
  ) %>%
  ungroup()


# Boxplot highlighting outliers
ggplot(outlier_data, aes(x = factor(RelativeVisit), y = Value, fill = Code)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, outlier.color = "red", alpha = 0.7) +
  facet_wrap(~Measure, scales = "free_y") +
  labs(
    title = "Distribution of DFT Scores by Year and Group",
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


RBD_DFT_long_2tp <- RBD_DFT_long_2tp %>%
  mutate(RelativeVisit = as.numeric(RelativeVisit))


#check IS values
RBD_DFT_long_2tp %>%
  filter(Measure == "onefingertaptest_worst_is") %>%
  summarise(
    min_IS = min(Value, na.rm = TRUE),
    n_negative = sum(Value < 0, na.rm = TRUE),
    n_zero = sum(Value == 0, na.rm = TRUE)
  )

RBD_DFT_long_2tp %>%
  group_by(Measure) %>%
  summarise(
    min_value = min(Value, na.rm = TRUE),
    n_negative = sum(Value < 0, na.rm = TRUE),
    n_zero = sum(Value == 0, na.rm = TRUE),
    .groups = "drop"
  )


#log transform
apply_dft_qc_and_transform <- function(RBD_DFT_long_2tp,
                                       measure_col = "Measure",
                                       value_col = "Value") {
  RBD_DFT_long_2tp %>%
    mutate(
      Value_num = suppressWarnings(as.numeric(.data[[value_col]])),
      
      qc_flag = case_when(
        is.na(Value_num) ~ "missing",
        Value_num <= 0 ~ "invalid_nonpositive",
        TRUE ~ "ok"
      ),
      
      Value_clean = case_when(
        qc_flag == "invalid_nonpositive" ~ NA_real_,
        TRUE ~ Value_num
      ),
      
      Value_clean_log = if_else(
        !is.na(Value_clean),
        log(Value_clean),
        NA_real_
      )
    )
}

RBD_DFT_long_qc <- apply_dft_qc_and_transform(RBD_DFT_long_2tp)

#check transformation worked
RBD_DFT_long_qc %>%
  group_by(Measure) %>%
  summarise(
    n_valid = sum(!is.na(Value_clean_log)),
    min_value = min(Value_clean_log, na.rm = TRUE),
    max_value = max(Value_clean_log, na.rm = TRUE),
    .groups = "drop"
  )

RBD_DFT_long_qc %>%
  count(Measure, qc_flag)

RBD_DFT_long_qc <- RBD_DFT_long_qc %>%
  mutate(RelativeVisit = as.numeric(RelativeVisit))




#check if data is normally distributed 
ggplot(RBD_DFT_long_2tp %>% filter(Measure == "onefingertaptest_worst_is"), 
       aes(x = Value, fill = Code)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Distribution of IS by Group",
    x = "Incoordination Score",
    y = "Density",
    fill = "Group"
  ) +
  xlim(3,12) +
  theme_minimal(base_size = 14)


ggplot(RBD_DFT_long_2tp %>% filter(Measure == "onefingertaptest_worst_ks"), 
       aes(x = Value, fill = Code)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Distribution of KS by Group",
    x = "Kinesia Score",
    y = "Density",
    fill = "Group"
  ) +
  xlim(-25,150) +
  theme_minimal(base_size = 14)

ggplot(RBD_DFT_long_2tp %>% filter(Measure == "onefingertaptest_worst_at"), 
       aes(x = Value, fill = Code)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Distribution of AT by Group",
    x = "Akinesia Time",
    y = "Density",
    fill = "Group"
  ) +
  xlim(-10,300) +
  theme_minimal(base_size = 14)

RBD_DFT_long_2tp %>%
  group_by(Measure) %>%
  summarise(
    p_value = shapiro.test(Value)$p.value
  )
#not normally distributed



#inspect visually
# Loop through each measure
for(m in unique(RBD_DFT_long_2tp$Measure)) {
  data_m <- RBD_DFT_long %>% filter(Measure == m, !is.na(Value))
  
  # Print skewness
  cat(m, "skewness:", round(skewness(data_m$Value), 2), "\n")
  
  # Plot histogram + QQ plot side by side
  gridExtra::grid.arrange(
    ggplot(data_m, aes(x = Value)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "#69b3a2", alpha = 0.6) +
      geom_density(color = "#1f3552", size = 1) +
      labs(title = paste(m, "- Histogram & Density"), x = "Value", y = "Density") +
      theme_minimal(base_size = 12),
    
    ggplot(data_m, aes(sample = Value)) +
      stat_qq() + stat_qq_line(color = "red") +
      labs(title = paste(m, "- QQ Plot"), x = "Theoretical Quantiles", y = "Sample Quantiles") +
      theme_minimal(base_size = 12),
    
    ncol = 2
  )
}



######### Variance Analysis ######### 
#Step 1a — Compute mean, SD, and sample size per group × year × measure
variance_summary <- RBD_DFT_long_qc %>%
  group_by(Measure, RelativeVisit, Code) %>%
  summarise(
    n = sum(!is.na(Value_clean_log)),
    mean_value = mean(Value_clean_log, na.rm = TRUE),
    sd_value = sd(Value_clean_log, na.rm = TRUE),
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
variance_boot <- RBD_DFT_long_qc %>%
  group_by(Measure, RelativeVisit, Code) %>%
  summarise(
    n = sum(!is.na(Value_clean_log)),
    boot_results = list(bootstrap_sd_ci(Value_clean_log)),
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


# Named vector for nice measure labels
measure_labels <- c(
  "onefingertaptest_worst_is" = "Incoordination Score",
  "onefingertaptest_worst_at" = "Akinesia Time",
  "onefingertaptest_worst_ks"     = "Kinesia Score"
)

# Directionality info for subtitle
measure_direction <- c(
  "onefingertaptest_worst_is" = "Higher = worse",
  "onefingertaptest_worst_at" = "Higher = worse",
  "onefingertaptest_worst_ks"     = "Lower = worse"
)

plot_data <- plot_data %>%
  mutate(Code = ifelse(Code == "iRBD", "PSG-iRBD", Code))

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
      title = paste0("DFT – ", measure_labels[m], " Over Time (Mean ± SD)"),
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
RBD_DFT_long_qc <- RBD_DFT_long_qc %>%
  mutate(
    Code = relevel(factor(Code), ref = "Control"),
    RelativeVisit_f = factor(RelativeVisit)
  )

# Filter IS scores
IS_data <- RBD_DFT_long_qc %>% 
  filter(Measure == "onefingertaptest_worst_is")


# Visit-specific variance model
m_visit_IS <- nlme::lme(
  Value_clean_log ~ Code * RelativeVisit + age_session_baseline + gender,
  random = ~ 1 | PID2,
  data = IS_data,
  na.action = na.exclude,
  weights = nlme::varIdent(form = ~ 1 | RelativeVisit_f)
)


#variance by visit and group (add one variance term for group)
m_group_IS <- update(
  m_visit_IS,
  weights = nlme::varComb(
    nlme::varIdent(form = ~ 1 | RelativeVisit_f),
    nlme::varIdent(form = ~ 1 | Code)
  )
)

anova(m_visit_IS, m_group_IS) #does variance differ between groups?

#direction of size and difference (95% CI & variance ratio)
ci_IS <- intervals(m_group_IS, which = "var-cov")$varStruct
ci_IS   # A.2 is the visit effect, ignore it

# Group row is B.Control (Control relative to iRBD), so invert to get iRBD vs Control
sd_ratio_IS <- 1 / ci_IS["B.Control", c("upper", "est.", "lower")]
names(sd_ratio_IS) <- c("lower", "est.", "upper")
sd_ratio_IS      # SD ratio, iRBD vs Control
sd_ratio_IS^2    # variance ratio, iRBD vs Control


# Step 3: calculate residual-based MSSD
# Calculate residuals for IS
IS_data$resid_fixed <- IS_data$Value_clean_log -
  fitted(m_visit_IS, level = 0)

# calculate residual-based MSSD for each participant
mssd_IS <- IS_data %>%
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
wilcox.test(mssd ~ Code, data = mssd_IS, exact = FALSE)


mssd_IS %>% count(Code, is.na(mssd))


# participants per group at each visit
table(
  IS_data$Code[!is.na(IS_data$Value_clean_log)],
  IS_data$RelativeVisit[!is.na(IS_data$Value_clean_log)]
)



## AT
unique(RBD_DFT_long_qc$Measure)
# Filter AT scores
AT_data <- RBD_DFT_long_qc %>% 
  filter(Measure == "onefingertaptest_worst_at")


# Visit-specific variance model
m_visit_AT <- nlme::lme(
  Value_clean_log ~ Code * RelativeVisit + age_session_baseline + gender,
  random = ~ 1 | PID2,
  data = AT_data,
  na.action = na.exclude,
  weights = nlme::varIdent(form = ~ 1 | RelativeVisit_f)
)


#variance by visit and group (add one variance term for group)
m_group_AT <- update(
  m_visit_AT,
  weights = nlme::varComb(
    nlme::varIdent(form = ~ 1 | RelativeVisit_f),
    nlme::varIdent(form = ~ 1 | Code)
  )
)

anova(m_visit_AT, m_group_AT) #does variance differ between groups?

#direction of size and difference (95% CI & variance ratio)
ci_AT <- intervals(m_group_AT, which = "var-cov")$varStruct
ci_AT

sd_ratio_AT <- 1 / ci_AT["B.Control", c("upper", "est.", "lower")]
names(sd_ratio_AT) <- c("lower", "est.", "upper")
sd_ratio_AT      # SD ratio, iRBD vs Control
sd_ratio_AT^2    # variance ratio, iRBD vs Control



# Step 3: calculate residual-based MSSD
# Calculate residuals for IS
AT_data$resid_fixed <- AT_data$Value_clean_log -
  fitted(m_visit_AT, level = 0)

# calculate residual-based MSSD for each participant
mssd_AT <- AT_data %>%
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
wilcox.test(mssd ~ Code, data = mssd_AT, exact = FALSE)


mssd_AT %>% count(Code, is.na(mssd))


# participants per group at each visit
table(
  AT_data$Code[!is.na(AT_data$Value_clean_log)],
  AT_data$RelativeVisit[!is.na(AT_data$Value_clean_log)]
)



## KS
unique(RBD_DFT_long_qc$Measure)
# Filter AT scores
KS_data <- RBD_DFT_long_qc %>% 
  filter(Measure == "onefingertaptest_worst_ks")


# Visit-specific variance model
m_visit_KS <- nlme::lme(
  Value_clean_log ~ Code * RelativeVisit + age_session_baseline + gender,
  random = ~ 1 | PID2,
  data = KS_data,
  na.action = na.exclude,
  weights = nlme::varIdent(form = ~ 1 | RelativeVisit_f)
)


#variance by visit and group (add one variance term for group)
m_group_KS <- update(
  m_visit_KS,
  weights = nlme::varComb(
    nlme::varIdent(form = ~ 1 | RelativeVisit_f),
    nlme::varIdent(form = ~ 1 | Code)
  )
)

anova(m_visit_KS, m_group_KS) #does variance differ between groups?

#direction of size and difference (95% CI & variance ratio)
ci_KS <- intervals(m_group_KS, which = "var-cov")$varStruct
ci_KS

sd_ratio_KS <- 1 / ci_KS["B.Control", c("upper", "est.", "lower")]
names(sd_ratio_KS) <- c("lower", "est.", "upper")
sd_ratio_KS      # SD ratio, iRBD vs Control
sd_ratio_KS^2    # variance ratio, iRBD vs Control


# Step 3: calculate residual-based MSSD
# Calculate residuals for IS
KS_data$resid_fixed <- KS_data$Value_clean_log -
  fitted(m_visit_KS, level = 0)

# calculate residual-based MSSD for each participant
mssd_KS <- KS_data %>%
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
wilcox.test(mssd ~ Code, data = mssd_KS, exact = FALSE)


mssd_KS %>% count(Code, is.na(mssd))


# participants per group at each visit
table(
  KS_data$Code[!is.na(KS_data$Value_clean_log)],
  KS_data$RelativeVisit[!is.na(KS_data$Value_clean_log)]
)




###### demographics ###### 
RBD_DFT_long_2tp %>%
  distinct(PID2, Code) %>%
  count(Code)

#keep one row per participant for demographics
demo_df <- RBD_DFT_long_2tp %>%
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

fisher.test(table(demo_df$gender, demo_df$Code))

RBD_DFT_long_qc %>%
  filter(!is.na(Value_clean)) %>%
  group_by(Measure, Code) %>%
  summarise(
    skewness = round(skewness(Value_clean), 2),
    .groups = "drop"
  )

summary_means <- RBD_DFT_long_qc %>%
  group_by(Measure, RelativeVisit) %>%
  summarise(
    n_control = sum(Code == "Control" & !is.na(Value_clean)),
    n_RBD = sum(Code == "iRBD" & !is.na(Value_clean)),
    
    median_control = median(
      Value_clean[Code == "Control"], na.rm = TRUE
    ),
    q1_control = quantile(
      Value_clean[Code == "Control"], 0.25, na.rm = TRUE
    ),
    q3_control = quantile(
      Value_clean[Code == "Control"], 0.75, na.rm = TRUE
    ),
    
    median_RBD = median(
      Value_clean[Code == "iRBD"], na.rm = TRUE
    ),
    q1_RBD = quantile(
      Value_clean[Code == "iRBD"], 0.25, na.rm = TRUE
    ),
    q3_RBD = quantile(
      Value_clean[Code == "iRBD"], 0.75, na.rm = TRUE
    ),
    
    p_value = wilcox.test(
      Value_clean ~ Code,
      data = pick(Value_clean, Code),
      exact = FALSE
    )$p.value,
    
    .groups = "drop"
  )

print(summary_means, n = nrow(summary_means))


# Close PDF
dev.off()

# Close console output
sink()


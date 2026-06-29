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


matched_RBD_control_DFT_data <- read_csv("DFT/RBD vs controls/matched_RBD_control_DFT_data.csv")
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


# log transform IS
cols_is <- c("onefingertaptest_worst_is_5", "onefingertaptest_worst_is_7")

matched_RBD_control_DFT_data2 <- matched_RBD_control_DFT_data2 %>%
  mutate(
    across(
      all_of(cols_is),
      ~ {
        x <- readr::parse_number(as.character(.x))
        if (any(x < 0, na.rm = TRUE)) {
          warning("Negative values found in IS column; log1p invalid for negatives.")
        }
        log1p(x)
      }
    )
  )


######## Exploratory data analysis #########
# List of test columns
cols_dft <- c(
  "onefingertaptest_worst_is_5", "onefingertaptest_worst_is_7",
  "onefingertaptest_worst_ks_5", "onefingertaptest_worst_ks_7", 
  "onefingertaptest_worst_at_5", "onefingertaptest_worst_at_7")

#check columns
names(matched_RBD_control_DFT_data2)[!sapply(matched_RBD_control_DFT_data2, is.numeric)]

#make columns numeric
cols_to_convert <- c( "onefingertaptest_worst_is_5", "onefingertaptest_worst_is_7",
                      "onefingertaptest_worst_ks_5", "onefingertaptest_worst_ks_7", 
                      "onefingertaptest_worst_at_5", "onefingertaptest_worst_at_7")

matched_RBD_control_DFT_data2[cols_to_convert] <- lapply(matched_RBD_control_DFT_data2[cols_to_convert], as.numeric)


#remove participants with ALL test data missing 
matched_RBD_control_DFT_data2 <- matched_RBD_control_DFT_data2 %>%
  filter(
    rowSums(is.na(select(., starts_with("onefingertaptest_worst")))) <
      ncol(select(., starts_with("onefingertaptest_worst")))
  )


# Pivot to long format but keep all non test columns
RBD_DFT_long <- matched_RBD_control_DFT_data2 %>%
  pivot_longer(
    cols = all_of(cols_dft),                # only pivot DFT columns
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
  xlim(0,150) +
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
  xlim(-0,300) +
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

#AT is skewed
# log transform AT
RBD_DFT_long_2tp <- RBD_DFT_long_2tp %>%
  mutate(
    Value_num = suppressWarnings(as.numeric(Value)),
    Value_clean_log = case_when(
      Measure == "onefingertaptest_worst_at" & !is.na(Value_num) & Value_num > 0 ~ log1p(Value_num),
      Measure == "onefingertaptest_worst_at" ~ NA_real_,   # AT <= 0 becomes NA
      TRUE ~ Value_num
    )
  )



######### Variance Analysis ######### 
#Step 1a — Compute mean, SD, and sample size per group × year × measure
variance_summary <- RBD_DFT_long_2tp %>%
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
variance_boot <- RBD_DFT_long_2tp %>%
  group_by(Measure, RelativeVisit, Code) %>%
  summarise(
    n = sum(!is.na(Value_clean_log)),
    boot_results = list(bootstrap_sd_ci(Value_clean_log)),
    .groups = "drop"
  ) %>%
  unnest(cols = c(boot_results))

print(variance_boot, n = Inf, width = Inf)


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
  "onefingertaptest_worst_ks" = "Lower = worse"
)


plot_data <- plot_data %>%
  mutate(Code = ifelse(Code == "iRBD", "PSG-iRBD", Code))

# Loop through measures and plot with subtitle
for (m in measures) {
  
  p <- plot_data %>%
    filter(Measure == m) %>%
    ggplot(aes(x = RelativeVisit, y = mean_value, color = Code, group = Code, fill = Code)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = mean_value - sd_lower, ymax = mean_value + sd_upper),
                alpha = 0.2, color = NA) +
    labs(
      title = paste0("DFT – ", measure_labels[m], " Over Time (Mean ± SD)"),
      subtitle = measure_direction[m],
      x = "Year",
      y = "Score",
      color = "Group",
      fill = "Group"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, face = "italic"),
      legend.position = "top"
    )
  
  print(p)
}


#Step 2: subject-level variance (“wobblers”) + wilcoxon rank
calc_subject_variability <- function(long_df,
                                     measures,
                                     id_col   = "PID2",
                                     group_col = "Code",
                                     time_col = "RelativeVisit",
                                     value_col = "Value_clean_log",
                                     exact_two_sided = FALSE) {
  
  stopifnot(all(c("Measure", id_col, group_col, time_col, value_col) %in% names(long_df)))
  stopifnot(all(measures %in% unique(long_df$Measure)))
  
  # compute SD, Range, MSSD for one subject (already ordered by time)
  subject_metrics <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(c(n_timepoints = 0, sd_subject = NA_real_, range_subject = NA_real_, mssd = NA_real_))
    }
    sdv <- if (length(x) >= 2) sd(x) else NA_real_
    rgv <- if (length(x) >= 1) max(x) - min(x) else NA_real_
    msv <- if (length(x) >= 2) mean(diff(x)^2) else NA_real_
    c(n_timepoints = length(x), sd_subject = sdv, range_subject = rgv, mssd = msv)
  }
  
  # One row per participant × measure with SD, Range, MSSD
  subj_tbl <- long_df %>%
    dplyr::filter(.data$Measure %in% measures) %>%
    dplyr::arrange(.data[[id_col]], .data$Measure, .data[[time_col]]) %>%
    dplyr::group_by(.data[[id_col]], .data[[group_col]], .data$Measure) %>%
    dplyr::summarise(metrics = list(subject_metrics(.data[[value_col]])),
                     .groups = "drop") %>%
    tidyr::unnest_wider(metrics)
  
  #Wilcoxon tests for SD, Range, MSSD per measure
  metrics_to_test <- c("sd_subject", "range_subject", "mssd")
  wilcox_rows <- list()
  
  for (m in unique(subj_tbl$Measure)) {
    dat_m <- subj_tbl %>% dplyr::filter(Measure == m)
    if (dplyr::n_distinct(dat_m[[group_col]]) < 2) next
    
    for (met in metrics_to_test) {
      d2 <- dat_m %>% dplyr::filter(!is.na(.data[[met]]))
      
      if (dplyr::n_distinct(d2[[group_col]]) < 2) {
        p <- NA_real_; W <- NA_real_
      } else {
        test <- wilcox.test(
          stats::as.formula(paste(met, "~", group_col)),
          data = d2,
          exact = exact_two_sided
        )
        p <- unname(test$p.value)
        W <- unname(test$statistic)
      }
      
      wilcox_rows[[length(wilcox_rows) + 1]] <- tibble::tibble(
        Measure = m,
        Metric = met,
        W = W,
        p_value = p,
        n_subjects_total = nrow(d2),
        n_groups = dplyr::n_distinct(d2[[group_col]])
      )
    }
  }
  
  wilcox_tbl <- dplyr::bind_rows(wilcox_rows) %>%
    dplyr::mutate(Metric = dplyr::recode(Metric,
                                         sd_subject = "SD",
                                         range_subject = "Range",
                                         mssd = "MSSD"))
  
  list(subject_metrics = subj_tbl, wilcox_tests = wilcox_tbl)
}


#Step 2b: run
DFT_measures <- c("onefingertaptest_worst_is", "onefingertaptest_worst_at", "onefingertaptest_worst_ks")

res_var <- calc_subject_variability(
  long_df   = RBD_DFT_long_2tp,
  measures  = DFT_measures,
  value_col = "Value_clean_log",
  exact_two_sided = FALSE
)

subject_variance_all <- res_var$subject_metrics
stat_results_all     <- res_var$wilcox_tests


#create boxplots
make_variability_boxplot <- function(subject_tbl, measure_name, plot_title) {
  subject_tbl %>%
    dplyr::filter(Measure == measure_name) %>%
    tidyr::pivot_longer(
      cols = c(sd_subject, range_subject, mssd),
      names_to = "Metric",
      values_to = "Value_clean_log"
    ) %>%
    dplyr::mutate(
      Metric = dplyr::recode(Metric,
                             sd_subject = "SD",
                             range_subject = "Range",
                             mssd = "MSSD")
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = Code, y = Value_clean_log, fill = Code)) +
    ggplot2::geom_boxplot(alpha = 0.5) +
    ggplot2::facet_wrap(~ Metric, scales = "free_y") +
    ggplot2::labs(title = plot_title, x = "Group", y = "Value") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(legend.position = "none",
                   strip.text = ggplot2::element_text(face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
}

p_IS <- make_variability_boxplot(subject_variance_all,
                                 "onefingertaptest_worst_is",
                                 "Within-Subject Variability of IS")

p_AT <- make_variability_boxplot(subject_variance_all,
                                 "onefingertaptest_worst_at",
                                 "Within-Subject Variability of AT")

p_KS <- make_variability_boxplot(subject_variance_all,
                                 "onefingertaptest_worst_ks",
                                 "Within-Subject Variability of KS")

print(p_IS); print(p_AT); print(p_KS)
subject_variance_all %>% count(Measure, n_timepoints)

wilcox_table <- stat_results_all %>%
  select(Measure, Metric, p_value) %>%
  tidyr::pivot_wider(names_from = Metric, values_from = p_value)

wilcox_table


#Step 3:	Model time-varying variance (heteroscedastic mixed model)
# Filter IS scores
IS_data <- RBD_DFT_long_2tp %>% 
  filter(Measure == "onefingertaptest_worst_is")
#make relative visit a factor 
IS_data <- IS_data %>%
  mutate(RelativeVisit_f = factor(RelativeVisit))

#Homoscedastic LME (residual variance constant)
lme_hom <- nlme::lme(
  Value_clean_log ~ Code * RelativeVisit + age_session_baseline + gender,
  random = ~ 1 | PID2,
  data = IS_data,
  na.action = na.exclude
)

# Visit-specific residual variance
lme_varident <- update(
  lme_hom,
  weights = nlme::varIdent(form = ~ 1 | RelativeVisit_f)
)

anova(lme_hom, lme_varident)

summary(lme_hom)


# Filter AT scores
AT_data <- RBD_DFT_long_2tp %>% 
  filter(Measure == "onefingertaptest_worst_at")
#make relative visit a factor 
AT_data <- AT_data %>%
  mutate(RelativeVisit_f = factor(RelativeVisit))


#Homoscedastic LME (residual variance constant)
lme_hom_AT <- nlme::lme(
  Value_clean_log ~ Code * RelativeVisit + age_session_baseline + gender,
  random = ~ 1 | PID2,
  data = AT_data,
  na.action = na.exclude
)

# Visit-specific residual variance
lme_varident_AT <- update(
  lme_hom_AT,
  weights = nlme::varIdent(form = ~ 1 | RelativeVisit_f)
)

anova(lme_hom_AT, lme_varident_AT)

summary(lme_hom_AT)


# Filter KS scores
KS_data <- RBD_DFT_long_2tp %>% 
  filter(Measure == "onefingertaptest_worst_ks")
#make relative visit a factor 
KS_data <- KS_data %>%
  mutate(RelativeVisit_f = factor(RelativeVisit))


#Homoscedastic LME (residual variance constant)
lme_hom_KS <- nlme::lme(
  Value_clean_log ~ Code * RelativeVisit + age_session_baseline + gender,
  random = ~ 1 | PID2,
  data = KS_data,
  na.action = na.exclude
)

# Visit-specific residual variance
lme_varident_KS <- update(
  lme_hom_KS,
  weights = nlme::varIdent(form = ~ 1 | RelativeVisit_f)
)

anova(lme_hom_KS, lme_varident_KS)

summary(lme_hom_KS)


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

chi_test_result_sex <- demo_df %>%
  with(table(gender, Code)) %>%
  chisq.test()
print(chi_test_result_sex)

summary_means <- RBD_DFT_long_2tp %>%
  group_by(Measure, RelativeVisit) %>%
  summarise(
    n_control = sum(Code == "Control" & !is.na(Value_clean_log)),
    n_RBD     = sum(Code == "iRBD" & !is.na(Value_clean_log)),
    
    mean_control = if (first(Measure) == "onefingertaptest_worst_is") {
      median(Value_clean_log[Code == "Control"], na.rm = TRUE)
    } else {
      mean(Value_clean_log[Code == "Control"], na.rm = TRUE)
    },
    
    sd_control = sd(Value_clean_log[Code == "Control"], na.rm = TRUE),
    
    mean_RBD = if (first(Measure) == "onefingertaptest_worst_is") {
      median(Value_clean_log[Code == "iRBD"], na.rm = TRUE)
    } else {
      mean(Value_clean_log[Code == "iRBD"], na.rm = TRUE)
    },
    
    sd_RBD = sd(Value_clean_log[Code == "iRBD"], na.rm = TRUE),
    
    p_value = wilcox.test(Value_clean_log ~ Code,
                          data = pick(Value_clean_log, Code),
                          exact = FALSE)$p.value,
    .groups = "drop"
  )


print(summary_means, n = nrow(summary_means))

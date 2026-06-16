##### TMT sRBD + hyposmia vs controls data analysis baseline #####
library(dplyr)
library(ggplot2)
library(readr)
library(readxl)
library(moments)
library(broom) # to tidy
library(ResourceSelection) # for Hosmer-Lemeshow test 

sRBD_hyp_control_TMT_data <- read_csv("sRBD_hyp_control_TMT_data.csv")

## remove those with essential tremor and under 50
sRBD_hyp_control_TMT_data <- sRBD_hyp_control_TMT_data %>%
  filter(hx_essentialtremor != "Yes") %>% 
  filter(age_session_baseline >= 50)

sRBD_hyp_control_TMT_data %>%
  filter(Code == "Control", rx_anti_psychotics_use_final_ever != "Unknown") %>%
  count(rx_anti_psychotics_use_final_ever)

table(sRBD_hyp_control_TMT_data$PREDICTPDGroupingsDetail_ParkinsonsDiseaseHighRisk_baseline,
      sRBD_hyp_control_TMT_data$Code)

### remove sessions other than baseline
sRBD_hyp_control_baseline_data <- sRBD_hyp_control_TMT_data %>%
  select(-tmt_trailA_accuracy_5, -tmt_trailA_accuracy_7, -tmt_trailA_reactiontime_5, -tmt_trailA_reactiontime_7,
         -tmt_trailAminustrailB_exefcn_final_5, -tmt_trailAminustrailB_exefcn_final_7, -tmt_trailB_accuracy_5,
         -tmt_trailB_accuracy_7, -tmt_trailB_reactiontime_5, -tmt_trailB_reactiontime_7)


# combine baseline data for session 4 and 8
## only use session 8 if session 4 is missing 
sRBD_hyp_control_baseline_data2 <- sRBD_hyp_control_baseline_data %>%
  mutate(tmt_trailA_accuracy_baseline = coalesce(tmt_trailA_accuracy_4, tmt_trailA_accuracy_8)) %>%
  select(-tmt_trailA_accuracy_4, -tmt_trailA_accuracy_8)  # Remove old columns if not needed

sRBD_hyp_control_baseline_data2 <- sRBD_hyp_control_baseline_data2 %>%
  mutate(tmt_trailA_reactiontime_baseline = coalesce(tmt_trailA_reactiontime_4, tmt_trailA_reactiontime_8)) %>%
  select(-tmt_trailA_reactiontime_4, -tmt_trailA_reactiontime_8)

sRBD_hyp_control_baseline_data2 <- sRBD_hyp_control_baseline_data2 %>%
  mutate(tmt_trailAminustrailB_exefcn_final_baseline = coalesce(tmt_trailAminustrailB_exefcn_final_4, tmt_trailAminustrailB_exefcn_final_8)) %>%
  select(-tmt_trailAminustrailB_exefcn_final_4, -tmt_trailAminustrailB_exefcn_final_8)

sRBD_hyp_control_baseline_data2 <- sRBD_hyp_control_baseline_data2 %>%
  mutate(tmt_trailB_accuracy_baseline = coalesce(tmt_trailB_accuracy_4, tmt_trailB_accuracy_8)) %>%
  select(-tmt_trailB_accuracy_4, -tmt_trailB_accuracy_8)

sRBD_hyp_control_baseline_data2 <- sRBD_hyp_control_baseline_data2 %>%
  mutate(tmt_trailB_reactiontime_baseline = coalesce(tmt_trailB_reactiontime_4, tmt_trailB_reactiontime_8)) %>%
  select(-tmt_trailB_reactiontime_4, -tmt_trailB_reactiontime_8)

table(sRBD_hyp_control_baseline_data2$Code)


#### baseline analysis sRBD hyposmia vs controls ####
summary(sRBD_hyp_control_baseline_data2)
#get rid of IDs with missing data
sRBD_hyp_control_baseline_data2 <- sRBD_hyp_control_baseline_data2 %>%
  filter(!is.na(tmt_trailA_reactiontime_baseline) & !is.na(tmt_trailB_accuracy_baseline) & !is.na(tmt_trailA_accuracy_baseline) &
           !is.na(tmt_trailAminustrailB_exefcn_final_baseline) & !is.na(tmt_trailB_reactiontime_baseline))

table(sRBD_hyp_control_baseline_data2$Code)


## make subgroup for those with PD
sRBD_hyp_and_PD_bl_data <- sRBD_hyp_control_baseline_data2 %>%
  filter(hx_neuro_parkinsons == "Yes")
#remove those with PD
sRBD_hyp_control_baseline_data2 <- sRBD_hyp_control_baseline_data2 %>%
  anti_join(sRBD_hyp_and_PD_bl_data, by = c("PID2", "TokenFinal"))

table(sRBD_hyp_control_baseline_data2$hx_neuro_parkinsons, sRBD_hyp_control_baseline_data2$Code)

write.csv(sRBD_hyp_control_baseline_data2, "Pilot/Baseline/sRBD hyposmia vs controls/sRBD_hyp_control_baseline_TMT_data.csv", row.names = FALSE)

any(sRBD_hyp_control_baseline_data2$tmt_trailAminustrailB_exefcn_final_baseline == 0)
any(sRBD_hyp_control_baseline_data2$tmt_trailB_reactiontime_baseline == 0)
any(sRBD_hyp_control_baseline_data2$tmt_trailA_reactiontime_baseline == 0)

# have 0s remove 
sRBD_hyp_control_baseline_data2 <- sRBD_hyp_control_baseline_data2 %>%
  filter(
    abs(tmt_trailAminustrailB_exefcn_final_baseline ) > 1e-6,
    abs(tmt_trailB_reactiontime_baseline) > 1e-6,
    abs(tmt_trailA_reactiontime_baseline ) > 1e-6,
    abs(tmt_trailA_accuracy_baseline ) > 1e-6,
    abs(tmt_trailB_accuracy_baseline ) > 1e-6
  )


table(sRBD_hyp_control_baseline_data2$Code)


# inspect data for outliers
ggplot(sRBD_hyp_control_baseline_data2, aes(y = tmt_trailAminustrailB_exefcn_final_baseline)) +
  geom_boxplot(fill = "#69b3a2", color = "#1f3552", outlier.color = "red", outlier.shape = 16, outlier.size = 2) +
  labs(
    title = "Distribution of TMT Scores at Baseline",
    y = "Executive Function"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )


ggplot(sRBD_hyp_control_baseline_data2, aes(y = tmt_trailB_reactiontime_baseline)) +
  geom_boxplot(fill = "#69b3a2", color = "#1f3552", outlier.color = "red", outlier.shape = 16, outlier.size = 2) +
  labs(
    title = "Distribution of Trail B Reaction Time at Baseline",
    y = "Reaction Time"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

ggplot(sRBD_hyp_control_baseline_data2, aes(y = tmt_trailA_reactiontime_baseline)) +
  geom_boxplot(fill = "#69b3a2", color = "#1f3552", outlier.color = "red", outlier.shape = 16, outlier.size = 2) +
  labs(
    title = "Distribution of Trail A Reaction Time at Baseline",
    y = "Reaction Time"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )


ggplot(sRBD_hyp_control_baseline_data2, aes(y = tmt_trailB_accuracy_baseline)) +
  geom_boxplot(fill = "#69b3a2", color = "#1f3552", outlier.color = "red", outlier.shape = 16, outlier.size = 2) +
  labs(
    title = "Distribution of Trail B Accuracy at Baseline",
    y = "Accurary"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )


# have outliers check what cohort they belong to
Q1 <- quantile(sRBD_hyp_control_baseline_data2$tmt_trailAminustrailB_exefcn_final_baseline, 0.25, na.rm = TRUE)
Q3 <- quantile(sRBD_hyp_control_baseline_data2$tmt_trailAminustrailB_exefcn_final_baseline, 0.75, na.rm = TRUE)
IQR <- Q3 - Q1
lower <- Q1 - 1.5 * IQR
upper <- Q3 + 1.5 * IQR

outliers_overall <- sRBD_hyp_control_baseline_data2 %>%
  filter(tmt_trailAminustrailB_exefcn_final_baseline < lower | tmt_trailAminustrailB_exefcn_final_baseline > upper)

table(outliers_overall$Code)

Q1_2 <- quantile(sRBD_hyp_control_baseline_data2$tmt_trailB_reactiontime_baseline, 0.25, na.rm = TRUE)
Q3_2 <- quantile(sRBD_hyp_control_baseline_data2$tmt_trailB_reactiontime_baseline, 0.75, na.rm = TRUE)
IQR_2 <- Q3_2 - Q1_2
lower_2 <- Q1_2 - 1.5 * IQR_2
upper_2 <- Q3_2 + 1.5 * IQR_2

outliers_overall_2 <- sRBD_hyp_control_baseline_data2 %>%
  filter(tmt_trailB_reactiontime_baseline < lower | tmt_trailB_reactiontime_baseline > upper)

table(outliers_overall_2$Code)

Q1_3 <- quantile(sRBD_hyp_control_baseline_data2$tmt_trailA_reactiontime_baseline, 0.25, na.rm = TRUE)
Q3_3 <- quantile(sRBD_hyp_control_baseline_data2$tmt_trailA_reactiontime_baseline, 0.75, na.rm = TRUE)
IQR_3 <- Q3_3 - Q1_3
lower_3 <- Q1_3 - 1.5 * IQR_3
upper_3 <- Q3_3 + 1.5 * IQR_3

outliers_overall_3 <- sRBD_hyp_control_baseline_data2 %>%
  filter(tmt_trailA_reactiontime_baseline < lower | tmt_trailA_reactiontime_baseline > upper)

table(outliers_overall_3$Code)


#remove outliers in controls 
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

sRBD_hyp_control_baseline_data2_clean <- sRBD_hyp_control_baseline_data2 %>%
  mutate(across(
    c(tmt_trailA_reactiontime_baseline, tmt_trailAminustrailB_exefcn_final_baseline, tmt_trailB_reactiontime_baseline),
    ~ if_else(Code == "Control", 
              winsorize_with_max_min(.), 
              .)
  ))

str(sRBD_hyp_control_baseline_data2_clean)


# Check for missing values in specific columns
sum(is.na(sRBD_hyp_control_baseline_data2_clean$gender))
sum(is.na(sRBD_hyp_control_baseline_data2_clean$age_session_baseline))




### check distribution of baseline data - density plots
ggplot(sRBD_hyp_control_baseline_data2_clean, aes(x=age_session_baseline, group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(45,90)

TMT_baseline_male <- sRBD_hyp_control_baseline_data2_clean[sRBD_hyp_control_baseline_data2_clean$gender == "Male", ]

ggplot(TMT_baseline_male, aes(x=age_session_baseline, group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(45,90)

TMT_baseline_female <- sRBD_hyp_control_baseline_data2_clean[sRBD_hyp_control_baseline_data2_clean$gender == "Female", ]

ggplot(TMT_baseline_female, aes(x=age_session_baseline, group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(45, 90)

ggplot(sRBD_hyp_control_baseline_data2_clean, aes(x=tmt_trailA_reactiontime_baseline , group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(0, 3750)
#check skewness
TrailA_rt <- sRBD_hyp_control_baseline_data2_clean$tmt_trailA_reactiontime_baseline
TrailA_rt <- rnorm(5802)
skew_value <- skewness(TrailA_rt)
print(paste("Skewness:", skew_value))
if (abs(skew_value) < 0.5) {
  print("Data is fairly symmetrical")
} else if (abs(skew_value) < 1) {
  print("Data is moderately skewed")
} else {
  print("Data is highly skewed")
}


ggplot(sRBD_hyp_control_baseline_data2_clean, aes(x=tmt_trailAminustrailB_exefcn_final_baseline , group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(-1000 , 3000)
#check skewness
TrailA_rt <- sRBD_hyp_control_baseline_data2_clean$tmt_trailAminustrailB_exefcn_final_baseline
TrailA_rt <- rnorm(5802)
skew_value <- skewness(TrailA_rt)
print(paste("Skewness:", skew_value))
if (abs(skew_value) < 0.5) {
  print("Data is fairly symmetrical")
} else if (abs(skew_value) < 1) {
  print("Data is moderately skewed")
} else {
  print("Data is highly skewed")
}


ggplot(sRBD_hyp_control_baseline_data2_clean, aes(x=tmt_trailB_accuracy_baseline , group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(20,30)
#check skewness
TrailA_rt <- sRBD_hyp_control_baseline_data2_clean$tmt_trailB_accuracy_baseline
TrailA_rt <- rnorm(5802)
skew_value <- skewness(TrailA_rt)
print(paste("Skewness:", skew_value))
if (abs(skew_value) < 0.5) {
  print("Data is fairly symmetrical")
} else if (abs(skew_value) < 1) {
  print("Data is moderately skewed")
} else {
  print("Data is highly skewed")
}


ggplot(sRBD_hyp_control_baseline_data2_clean, aes(x=tmt_trailB_reactiontime_baseline , group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(0,6000)
#check skewness
TrailA_rt <- sRBD_hyp_control_baseline_data2_clean$tmt_trailB_reactiontime_baseline
TrailA_rt <- rnorm(5802)
skew_value <- skewness(TrailA_rt)
print(paste("Skewness:", skew_value))
if (abs(skew_value) < 0.5) {
  print("Data is fairly symmetrical")
} else if (abs(skew_value) < 1) {
  print("Data is moderately skewed")
} else {
  print("Data is highly skewed")
}


ggplot(sRBD_hyp_control_baseline_data2_clean, aes(x=tmt_trailA_accuracy_baseline, group=Code, fill=Code)) +
  geom_density(alpha=.5)+ scale_fill_manual(values=c("#009E73","#E69F00")) + theme_minimal() +
  xlim(20,30)
#check skewness
TrailA_rt <- sRBD_hyp_control_baseline_data2_clean$tmt_trailA_accuracy_baseline
TrailA_rt <- rnorm(5802)
skew_value <- skewness(TrailA_rt)
print(paste("Skewness:", skew_value))
if (abs(skew_value) < 0.5) {
  print("Data is fairly symmetrical")
} else if (abs(skew_value) < 1) {
  print("Data is moderately skewed")
} else {
  print("Data is highly skewed")
}




## Demographics ##
table(sRBD_hyp_control_baseline_data2_clean$Code)
t_test_age_result <- t.test(age_session_baseline ~ Code, data = sRBD_hyp_control_baseline_data2_clean)
print(t_test_age_result)

std_dev_ages <- sRBD_hyp_control_baseline_data2_clean %>%
  group_by(Code) %>%
  summarize(std_dev_age = sd(age_session_baseline, na.rm = TRUE)) %>%
  ungroup()
print(std_dev_ages)

sex_percentage <- sRBD_hyp_control_baseline_data2_clean %>%
  count(Code, gender) %>%
  group_by(Code) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()
print(sex_percentage)

contingency_table <- table(sRBD_hyp_control_baseline_data2_clean$gender, sRBD_hyp_control_baseline_data2_clean$Code)
chi_test_result_sex <- chisq.test(contingency_table)
print(chi_test_result_sex)

t.test(tmt_trailA_accuracy_baseline ~ Code, data = sRBD_hyp_control_baseline_data2_clean)

t.test(tmt_trailA_reactiontime_baseline ~ Code, data = sRBD_hyp_control_baseline_data2_clean)

t.test(tmt_trailAminustrailB_exefcn_final_baseline ~ Code, data = sRBD_hyp_control_baseline_data2_clean)

t.test(tmt_trailB_accuracy_baseline ~ Code, data = sRBD_hyp_control_baseline_data2_clean)

t.test(tmt_trailB_reactiontime_baseline ~ Code, data = sRBD_hyp_control_baseline_data2_clean)


#change education levels to secondary yes or no 
sRBD_hyp_control_baseline_data2_clean <- sRBD_hyp_control_baseline_data2_clean %>%
  mutate(education = ifelse(education_level %in% c("Secondary", "Further", "Higher"), "Yes", "No"))

education <- sRBD_hyp_control_baseline_data2_clean %>%
  count(Code, education) %>%
  group_by(Code) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()
print(education)
fisher.test(table(sRBD_hyp_control_baseline_data2_clean$Code, sRBD_hyp_control_baseline_data2_clean$education))$p.value



#### Regression model ####
# convert code column to 0 (control) and 1 (iRBD) for glm
sRBD_hyp_control_baseline_data2_clean <- sRBD_hyp_control_baseline_data2_clean %>%
  mutate(Code = ifelse(Code == "sRBD_hyposmia", 1, 0))
sRBD_hyp_control_baseline_data2_clean$Code <- as.factor(sRBD_hyp_control_baseline_data2_clean$Code)
table(sRBD_hyp_control_baseline_data2_clean$Code)


#Reaction time A
ggplot(sRBD_hyp_control_baseline_data2_clean, aes(x=Code, y=tmt_trailA_reactiontime_baseline)) + 
  geom_boxplot()
# standardise and use z score
sRBD_hyp_control_baseline_data2_clean$tmt_trailA_reactiontime_baseline_z <- scale(sRBD_hyp_control_baseline_data2_clean$tmt_trailA_reactiontime_baseline)
model_reactiontime_trailA2 <- glm( Code ~ tmt_trailA_reactiontime_baseline_z  + gender + age_session_baseline + education, 
                                   data = sRBD_hyp_control_baseline_data2_clean,family = 'binomial')
summary(model_reactiontime_trailA2)
exp(cbind(coef(model_reactiontime_trailA2), confint(model_reactiontime_trailA2,level = 0.95)))


#Executive function
ggplot(sRBD_hyp_control_baseline_data2_clean, aes(x=Code, y=tmt_trailAminustrailB_exefcn_final_baseline)) + 
  geom_boxplot()
# standardise and use z score
sRBD_hyp_control_baseline_data2_clean$tmt_trailAminustrailB_exefcn_final_baseline_z <- scale(sRBD_hyp_control_baseline_data2_clean$tmt_trailAminustrailB_exefcn_final_baseline)
model_execfunc2 <- glm( Code ~ tmt_trailAminustrailB_exefcn_final_baseline_z + gender + age_session_baseline + education, 
                        data = sRBD_hyp_control_baseline_data2_clean,family = 'binomial')
summary(model_execfunc2)
exp(cbind(coef(model_execfunc2), confint(model_execfunc2,level = 0.95))) 


#Reaction time B
ggplot(sRBD_hyp_control_baseline_data2_clean, aes(x=Code, y=tmt_trailB_reactiontime_baseline)) + 
  geom_boxplot()
# standardise and use z score
sRBD_hyp_control_baseline_data2_clean$tmt_trailB_reactiontime_baseline_z <- scale(sRBD_hyp_control_baseline_data2_clean$tmt_trailB_reactiontime_baseline)
model_reactiontimeB2 <- glm( Code ~ tmt_trailB_reactiontime_baseline_z + gender + age_session_baseline + education, 
                             data = sRBD_hyp_control_baseline_data2_clean,family = 'binomial') 
summary(model_reactiontimeB2)
exp(cbind(coef(model_reactiontimeB2), confint(model_reactiontimeB2,level = 0.95))) 



### check age linearity
#Reaction time A
ggplot(sRBD_hyp_control_baseline_data2_clean, aes(x = age_session_baseline, y = fitted(model_reactiontime_trailA2))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue") +
  theme_minimal() +
  labs(title = "Check for Linearity in Log-Odds, RT A", x = "Age", y = "Predicted Probability of Code = 1")

sRBD_hyp_control_baseline_data2_clean$age_baseline_sq <- sRBD_hyp_control_baseline_data2_clean$age_session_baseline^2

model_age_quad <- glm(Code ~ tmt_trailA_reactiontime_baseline_z + gender + age_session_baseline + age_baseline_sq + education,
                      data = sRBD_hyp_control_baseline_data2_clean, family = "binomial")
summary(model_age_quad)
exp(cbind(coef(model_age_quad), confint(model_age_quad,level = 0.90)))
anova(model_reactiontime_trailA2, model_age_quad)
final_model_RTA <- model_reactiontime_trailA2


#Executive function
ggplot(sRBD_hyp_control_baseline_data2_clean, aes(x = age_session_baseline, y = fitted(model_execfunc2))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue") +
  theme_minimal() +
  labs(title = "Check for Linearity in Log-Odds, EF", x = "Age", y = "Predicted Probability of Code = 1")

model_age_quad_ef <- glm(Code ~ tmt_trailAminustrailB_exefcn_final_baseline_z + gender + age_session_baseline + age_baseline_sq + education,
                         data = sRBD_hyp_control_baseline_data2_clean, family = "binomial")
summary(model_age_quad_ef)
exp(cbind(coef(model_age_quad_ef), confint(model_age_quad_ef,level = 0.90)))
anova(model_execfunc2, model_age_quad_ef)
final_model_EF <- model_execfunc2


#Trail B reaction time
ggplot(sRBD_hyp_control_baseline_data2_clean, aes(x = age_session_baseline, y = fitted(model_reactiontimeB2))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue") +
  theme_minimal() +
  labs(title = "Check for Linearity in Log-Odds, RT B", x = "Age", y = "Predicted Probability of Code = 1")

model_age_quad_rtb <- glm(Code ~ tmt_trailB_reactiontime_baseline_z + gender + age_session_baseline + age_baseline_sq + education,
                          data = sRBD_hyp_control_baseline_data2_clean, family = "binomial")
summary(model_age_quad_rtb)
exp(cbind(coef(model_age_quad_rtb), confint(model_age_quad_rtb,level = 0.90)))
anova(model_reactiontimeB2, model_age_quad_rtb)
final_model_RTB <- model_reactiontimeB2



#### QQ plots ####
# List of models and labels
models <- list(
  RTA = final_model_RTA,
  EF = final_model_EF,
  RTB = final_model_RTB
)

df_all_residuals <- bind_rows(
  lapply(names(models), function(name) {
    data.frame(
      residuals = residuals(models[[name]], type = "deviance"),
      Outcome = name,
      Code = sRBD_hyp_control_baseline_data2_clean$Code
    )
  })
)

ggplot(df_all_residuals, aes(sample = residuals, color = factor(Code))) +
  stat_qq() +
  stat_qq_line(color = "black") +
  facet_wrap(~ Outcome, scales = "free") +
  labs(
    title = "QQ Plots of Deviance Residuals by Outcome",
    color = "Group",
    x = "Theoretical Quantiles",
    y = "Deviance Residuals"    
  ) +
  theme_minimal()



### tidy up final models
tidy(final_model_RTA, conf.int = TRUE, conf.level = 0.90, exponentiate = TRUE)

tidy(final_model_EF, conf.int = TRUE, conf.level = 0.90, exponentiate = TRUE)

tidy(final_model_RTB, conf.int = TRUE, conf.level = 0.90, exponentiate = TRUE)


#### Hosmer-Lemeshow test to check if the model fits ####
#RTA
sRBD_hyp_control_baseline_data2_clean$Code_numeric <- as.numeric(as.character(sRBD_hyp_control_baseline_data2_clean$Code))

hoslem.test(x = sRBD_hyp_control_baseline_data2_clean$Code_numeric,
            y = fitted(final_model_RTA),
            g = 10)


#EF
hoslem.test(x = sRBD_hyp_control_baseline_data2_clean$Code_numeric,
            y = fitted(final_model_EF),
            g = 10)


#RTB
hoslem.test(x = sRBD_hyp_control_baseline_data2_clean$Code_numeric,
            y = fitted(final_model_RTB),
            g = 10)



### Create Forest plots ###
#RTA
# Extract ORs and 95% CIs
OR_data <- as.data.frame(exp(cbind(OR = coef(final_model_RTA), confint(final_model_RTA))))
OR_data$Variable <- rownames(OR_data)

# Replace raw variable names with  labels
OR_data$Variable <- recode(OR_data$Variable,
                           "(Intercept)" = "Intercept",
                           "tmt_trailA_reactiontime_baseline_z" = "Trail A RT (z-score)",
                           "genderMale" = "Male (vs Female)",
                           "age_session_baseline" = "Age",
                           "age_baseline_sq" = "Age²",
                           "educationYes" = "Education (Yes)"
)

# Exclude Intercept for clarity
ggplot(OR_data[-1, ], aes(x = Variable, y = OR)) +
  geom_point(size = 3, color = "#2c7fb8") +   
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                width = 0.2, size = 1.2, color = "#2c7fb8") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", size = 1) +
  coord_flip() +
  labs(title = "Odds Ratios from Logistic Regression",
       y = "Odds Ratio (log scale)", x = "") +
  scale_y_log10() +
  theme_minimal(base_size = 14)


#EF
# Extract ORs and 95% CIs
OR_data_EF <- as.data.frame(exp(cbind(OR = coef(final_model_EF), confint(final_model_EF))))
OR_data_EF$Variable <- rownames(OR_data_EF)

# Replace raw variable names with  labels
OR_data_EF$Variable <- recode(OR_data_EF$Variable,
                              "(Intercept)" = "Intercept",
                              "tmt_trailAminustrailB_exefcn_final_baseline_z" = "Executive Function (z-score)",
                              "genderMale" = "Male (vs Female)",
                              "age_session_baseline" = "Age",
                              "age_baseline_sq" = "Age²",
                              "educationYes" = "Education (Yes)"
)


# Exclude Intercept for clarity
ggplot(OR_data_EF[-1, ], aes(x = Variable, y = OR)) +
  geom_point(size = 3, color = "#2c7fb8") +   # <- size, not linewidth
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                width = 0.2, size = 1.2, color = "#2c7fb8") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", size = 1) +
  coord_flip() +
  labs(title = "Odds Ratios from Logistic Regression",
       y = "Odds Ratio (log scale)", x = "") +
  scale_y_log10() +
  theme_minimal(base_size = 14)

#RTB
# Extract ORs and 95% CIs
OR_data_RTB <- as.data.frame(exp(cbind(OR = coef(final_model_RTB), confint(final_model_RTB))))
OR_data_RTB$Variable <- rownames(OR_data_RTB)

# Replace raw variable names with  labels
OR_data_RTB$Variable <- recode(OR_data_RTB$Variable,
                               "(Intercept)" = "Intercept",
                               "tmt_trailB_reactiontime_baseline_z" = "Trail B RT (z-score)",
                               "genderMale" = "Male (vs Female)",
                               "age_session_baseline" = "Age",
                               "age_baseline_sq" = "Age²",
                               "educationYes" = "Education (Yes)"
)


# Exclude Intercept for clarity
ggplot(OR_data_RTB[-1, ], aes(x = Variable, y = OR)) +
  geom_point(size = 3, color = "#2c7fb8") +   # <- size, not linewidth
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                width = 0.2, size = 1.2, color = "#2c7fb8") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", size = 1) +
  coord_flip() +
  labs(title = "Odds Ratios from Logistic Regression",
       y = "Odds Ratio (log scale)", x = "") +
  scale_y_log10() +
  theme_minimal(base_size = 14)


# collate into one forest plot
RTA_data <- tidy(final_model_RTA, conf.int = TRUE, conf.level = 0.90, exponentiate = TRUE) %>%
  filter(term == "tmt_trailA_reactiontime_baseline_z") %>%
  mutate(Variable = "Reaction Time (Trail A)")

EF_data <- tidy(final_model_EF, conf.int = TRUE, conf.level = 0.90, exponentiate = TRUE) %>%
  filter(term == "tmt_trailAminustrailB_exefcn_final_baseline_z") %>%
  mutate(Variable = "Executive Function")

RTB_data <- tidy(final_model_RTB, conf.int = TRUE, conf.level = 0.90, exponentiate = TRUE) %>%
  filter(term == "tmt_trailB_reactiontime_baseline_z") %>%
  mutate(Variable = "Reaction Time (Trail B)")

all_TMT_variables <- bind_rows(RTA_data, EF_data, RTB_data) %>%
  mutate(p_label = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)))

ggplot(all_TMT_variables, aes(x = estimate, y = Variable)) +
  geom_point(size = 3, color = "#2c7fb8") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, size = 1.2, color = "#2c7fb8") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", size = 1) +
  geom_text(aes(label = paste("p =", p_label), x = conf.high * 1.15), 
            hjust = 0, size = 4.2, color = "black") +
  scale_x_log10(limits = c(0.5, 3.5)) +  # ← sets fixed axis range
  labs(title = "Baseline Associations of TMT Scores with pRBD with hyposmia vs. Controls",
       subtitle = "Odds Ratios from Logistic Regression Adjusted for Age, Gender and Education",
       x = "Odds Ratio (log scale)", y = "") +
  theme_minimal(base_size = 14)



#extract data for overall forest plot
extract_row <- function(model, term, variable, test, comparison,
                        conf_level = 0.95) {
  tidy(model, conf.int = TRUE, conf.level = conf_level, exponentiate = TRUE) %>%
    filter(term == !!term) %>%
    mutate(
      Variable   = variable,
      Test       = test,
      Comparison = comparison
    ) %>%
    select(Test, Comparison, Variable, estimate, conf.low, conf.high, p.value)
}


TMT_df_sRBD_hyp <- bind_rows(
  extract_row(final_model_RTA,
              "tmt_trailA_reactiontime_baseline_z",
              "Reaction Time (Trail A)",
              test = "TMT", comparison = "pRBD + hyposmia vs control"),
  extract_row(final_model_EF,
              "tmt_trailAminustrailB_exefcn_final_baseline_z",
              "Executive Function",
              test = "TMT", comparison = "pRBD + hyposmia vs control"),
  extract_row(final_model_RTB,
              "tmt_trailB_reactiontime_baseline_z",
              "Reaction Time (Trail B)",
              test = "TMT", comparison = "pRBD + hyposmia vs control")
)

saveRDS(TMT_df_sRBD_hyp, "Pilot/Baseline/sRBD hyposmia vs controls/forest_rows_TMT_hyp.rds")

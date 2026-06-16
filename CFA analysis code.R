########### CFA analysis- one composite outcome ########### 
library(readr)
library(dplyr)
library(tidyr)
library(lme4)
library(lavaan)


#read in datasets
RBD_control_BT_data <- read_csv("C:/Users/hwx036/OneDrive - Queen Mary, University of London/PhD/2025- Year 2/Motor analysis/Finger tapping - harmonised data/Motor testing/BRAIN test/RBD_control_BT_data.csv")
View(RBD_control_BT_data)

sRBD_control_BT_data <- read_csv("C:/Users/hwx036/OneDrive - Queen Mary, University of London/PhD/2025- Year 2/Motor analysis/Finger tapping - harmonised data/Motor testing/BRAIN test/sRBD_control_BT_data.csv")
View(sRBD_control_BT_data)

sRBD_hyp_control_BT_data <- read_csv("C:/Users/hwx036/OneDrive - Queen Mary, University of London/PhD/2025- Year 2/Motor analysis/Finger tapping - harmonised data/Motor testing/BRAIN test/sRBD_hyp_control_BT_data.csv")
View(sRBD_hyp_control_BT_data)

RBD_control_TMT_data <- read_csv("C:/Users/hwx036/OneDrive - Queen Mary, University of London/PhD/2024- Year 1/Preliminary analysis/Cognitive analysis/TMT/TMT data extraction/RBD_control_TMT_data.csv")
View(RBD_control_TMT_data)

sRBD_control_TMT_data <- read_csv("C:/Users/hwx036/OneDrive - Queen Mary, University of London/PhD/2024- Year 1/Preliminary analysis/Cognitive analysis/TMT/TMT data extraction/sRBD_control_TMT_data.csv")
View(sRBD_control_TMT_data)

sRBD_hyp_control_TMT_data <- read_csv("C:/Users/hwx036/OneDrive - Queen Mary, University of London/PhD/2024- Year 1/Preliminary analysis/Cognitive analysis/TMT/TMT data extraction/sRBD_hyp_control_TMT_data.csv")
View(sRBD_hyp_control_TMT_data)

#create BRAIN test master data set
RBD_control_BT_data <- RBD_control_BT_data %>%
  mutate(PID2 = as.character(PID2))

sRBD_control_BT_data <- sRBD_control_BT_data %>%
  mutate(PID2 = as.character(PID2))

sRBD_hyp_control_BT_data <- sRBD_hyp_control_BT_data %>%
  mutate(PID2 = as.character(PID2))

BT_all <- bind_rows(
  RBD_control_BT_data,
  sRBD_control_BT_data,
  sRBD_hyp_control_BT_data
)

#check for duplicates
BT_all %>%
  count(PID2) %>%
  filter(n > 1)

#convert to long format
#harmonise KS columns
BT_all <- BT_all %>%
  mutate(
    ks_4 = braintest_worst_ks_4,
    ks_5 = braintest_worst_ks_5,
    ks_7 = braintest_worst_ks_7,
    ks_8 = braintest_worst_ks_adj_8
  )

names(BT_all)[grepl("ks", names(BT_all))]

BT_all <- BT_all %>%
  select(-starts_with("braintest_worst_ks_"))

#remove duplication
BT_all <- BT_all %>%
  distinct(PID2, braintest_age_4, .keep_all = TRUE)

BT_all %>%
  count(PID2, braintest_age_4) %>%
  filter(n > 1)

BT_long <- BT_all %>%
  pivot_longer(
    cols = c(
      braintest_worst_at_adj_4, braintest_worst_at_adj_5,
      braintest_worst_at_adj_7, braintest_worst_at_adj_8,
      
      braintest_worst_is_adj_4, braintest_worst_is_adj_5,
      braintest_worst_is_adj_7, braintest_worst_is_adj_8,
      
      ks_4, ks_5, ks_7, ks_8
    ),
    names_to = c("Measure", "Year"),
    names_pattern = "(.*)_(\\d+)$",
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

#wide format for latent variable analysis.
#designed for specifying, estimating, and evaluating latent variable models such as confirmatory factor analysis
BT_wide <- BT_long %>%
  pivot_wider(
    id_cols = c(PID2, Code),
    names_from = Year,
    values_from = c(
      braintest_worst_at_adj,
      braintest_worst_is_adj,
      ks
    ),
    names_glue = "{.value}_t{Year}"
  )

#rename variables for clarity
BT_wide <- BT_wide %>%
  rename(
    at_t1 = braintest_worst_at_adj_t1,
    at_t2 = braintest_worst_at_adj_t2,
    at_t3 = braintest_worst_at_adj_t3,
    at_t4 = braintest_worst_at_adj_t4,
    
    is_t1 = braintest_worst_is_adj_t1,
    is_t2 = braintest_worst_is_adj_t2,
    is_t3 = braintest_worst_is_adj_t3,
    is_t4 = braintest_worst_is_adj_t4
  )

###### CFA analysis ######
#now building a longitudinal CFA
#each time point has a latent factor
#the same construct is measured across time

#change ks direction
BT_wide <- BT_wide %>%
  mutate(across(
    c(ks_t1, ks_t2, ks_t3, ks_t4),
    ~ -1 * .
  ))

#standardise variables before running CFA
BT_wide_scaled <- BT_wide %>%
  mutate(across(
    c(at_t1, at_t2, at_t3, at_t4,
      is_t1, is_t2, is_t3, is_t4,
      ks_t1, ks_t2, ks_t3, ks_t4),
    scale
  ))

model <- '
F1 =~ l1*at_t1 + l2*is_t1 + l3*ks_t1
F2 =~ l1*at_t2 + l2*is_t2 + l3*ks_t2
F3 =~ l1*at_t3 + l2*is_t3 + l3*ks_t3
F4 =~ l1*at_t4 + l2*is_t4 + l3*ks_t4
'

fit <- cfa(model, data = BT_wide_scaled, missing = "fiml")
summary(fit, fit.measures = TRUE, standardized = TRUE)

## configure model
#to confirm the same factor structure fits across all timepoints before testing any equality constraints
model_configural <- '
F1 =~ at_t1 + is_t1 + ks_t1
F2 =~ at_t2 + is_t2 + ks_t2
F3 =~ at_t3 + is_t3 + ks_t3
F4 =~ at_t4 + is_t4 + ks_t4
'

fit_configural <- cfa(model_configural,
                      data = BT_wide_scaled,
                      missing = "fiml")
#some estimated variances are negative
#so the model is not statistically well-identified / not stable

#Run a 3-wave configural model
model_configural_3wave <- '
F1 =~ at_t1 + is_t1 + ks_t1
F2 =~ at_t2 + is_t2 + ks_t2
F3 =~ at_t3 + is_t3 + ks_t3
'

fit_configural_3wave <- cfa(model_configural_3wave,
                            data = BT_wide_scaled,
                            missing = "fiml")
#more secure, only warning left-> empty cases ignored

#metric invariance (3-wave)
model_metric_3wave <- '
F1 =~ l1*at_t1 + l2*is_t1 + l3*ks_t1
F2 =~ l1*at_t2 + l2*is_t2 + l3*ks_t2
F3 =~ l1*at_t3 + l2*is_t3 + l3*ks_t3
'

fit_metric_3wave <- cfa(model_metric_3wave,
                        data = BT_wide_scaled,
                        missing = "fiml")

#compare configural vs metric fit
summary(fit_configural_3wave, fit.measures = TRUE)
summary(fit_metric_3wave, fit.measures = TRUE)
#if CFA drops by <0.01 metric invariance holds, if >0.01 go to partial invariance

#The loadings are not equal across time, the construct is changing over time

#Move to partial invariance
model_partial <- '
F1 =~ l1*at_t1 + l2*is_t1 + l3*ks_t1
F2 =~ l1*at_t2 + l2*is_t2 + l3*ks_t2
F3 =~ l1*at_t3 + is_t3 + l3*ks_t3
'

fit_partial <- cfa(model_partial,
                   data = BT_wide_scaled,
                   missing = "fiml")

summary(fit_partial, fit.measures = TRUE)

#Configural CFI = 0.738
#Metric CFI = 0.719 (drop ≈ 0.019 → too big
#Partial CFI = 0.737 (back to ~configural)
# - partial invariance is supported, because freeing is_t3 fixed the drop in fit

# at and ks behave consistently across time
# is does NOT behave the same at time 3, so you correctly freed it
# composite is now valid enough to compare over time

#extract variables from partical invariance model
BT_scores <- lavPredict(fit_partial)
print(BT_scores)


# Same for TMT 
#create TMT master data set
RBD_control_TMT_data <- RBD_control_TMT_data %>%
  mutate(PID2 = as.character(PID2))

sRBD_control_TMT_data <- sRBD_control_TMT_data %>%
  mutate(PID2 = as.character(PID2))

sRBD_hyp_control_TMT_data <- sRBD_hyp_control_TMT_data %>%
  mutate(PID2 = as.character(PID2))

TMT_all <- bind_rows(
  RBD_control_TMT_data,
  sRBD_control_TMT_data,
  sRBD_hyp_control_TMT_data
)

#check for duplicates
TMT_all %>%
  count(PID2) %>%
  filter(n > 1)

names(TMT_all)


#convert to long format
TMT_long <- TMT_all %>%
  pivot_longer(
    cols = c(
      tmt_trailA_reactiontime_4, tmt_trailA_reactiontime_5,
      tmt_trailA_reactiontime_7, tmt_trailA_reactiontime_8,
      
      tmt_trailB_reactiontime_4, tmt_trailB_reactiontime_5,
      tmt_trailB_reactiontime_7, tmt_trailB_reactiontime_8,
      
      tmt_trailAminustrailB_exefcn_final_4,
      tmt_trailAminustrailB_exefcn_final_5,
      tmt_trailAminustrailB_exefcn_final_7,
      tmt_trailAminustrailB_exefcn_final_8
    ),
    names_to = c("Measure", "Year"),
    names_pattern = "(.*)_(\\d+)$",
    values_to = "Value"
  ) %>%
  mutate(
    Year = recode(Year, `4` = 1, `5` = 2, `7` = 3, `8` = 4),
    Year = as.numeric(Year)
  )

TMT_long <- TMT_long %>%
  mutate(
    Measure = case_when(
      grepl("trailA_reactiontime", Measure) ~ "tmtA",
      grepl("trailB_reactiontime", Measure) ~ "tmtB",
      grepl("AminustrailB", Measure) ~ "tmtAB"
    )
  )


#collapse duplicates
TMT_long <- TMT_long %>%
  group_by(PID2, Measure, Year) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

TMT_wide <- TMT_long %>%
  pivot_wider(
    id_cols = PID2,
    names_from = c(Measure, Year),
    values_from = Value,
    names_glue = "{Measure}_t{Year}"
  )

write.csv(TMT_wide, "TMT_yearly_output.csv", row.names = FALSE)
write.csv(BT_wide, "BT_yearly_output.csv", row.names = FALSE)

model_tmt <- '
F1 =~ tmtA_t1 + tmtB_t1 + tmtAB_t1
F2 =~ tmtA_t2 + tmtB_t2 + tmtAB_t2
F3 =~ tmtA_t3 + tmtB_t3 + tmtAB_t3
F4 =~ tmtA_t4 + tmtB_t4 + tmtAB_t4
'

fit_tmt <- cfa(model_tmt, data = TMT_wide, missing = "fiml")
#change to 3-wave


model_tmt_3wave <- '
F1 =~ tmtA_t1 + tmtB_t1 + tmtAB_t1
F2 =~ tmtA_t2 + tmtB_t2 + tmtAB_t2
F3 =~ tmtA_t3 + tmtB_t3 + tmtAB_t3
'

fit_tmt_3wave <- cfa(model_tmt_3wave,
                     data = TMT_wide,
                     missing = "fiml")
#TMT does not support a stable longitudinal latent factor model in your dataset

#check within-wave correlations
cor(TMT_wide[, c("tmtA_t1","tmtB_t1","tmtAB_t1")], use="pairwise.complete.obs")
#standardise tmtB and tmtAB (correlated whereas tmtA not correlated- cannot be same construct)
TMT_wide_scaled <- TMT_wide %>%
  mutate(across(
    c(tmtB_t1, tmtB_t2, tmtB_t3, tmtB_t4,
      tmtAB_t1, tmtAB_t2, tmtAB_t3, tmtAB_t4),
    ~ as.numeric(scale(.))
  ))

TMT_wide_scaled <- TMT_wide_scaled %>%
  mutate(
    tmt_comp_t1 = rowMeans(cbind(tmtB_t1, tmtAB_t1), na.rm = TRUE),
    tmt_comp_t2 = rowMeans(cbind(tmtB_t2, tmtAB_t2), na.rm = TRUE),
    tmt_comp_t3 = rowMeans(cbind(tmtB_t3, tmtAB_t3), na.rm = TRUE),
    tmt_comp_t4 = rowMeans(cbind(tmtB_t4, tmtAB_t4), na.rm = TRUE)
  )
summary(TMT_wide_scaled$tmt_comp_t1)

colSums(is.na(TMT_wide_scaled[, c("tmtB_t1", "tmtAB_t1")]))

colSums(is.na(TMT_wide_scaled[, c("tmtB_t2", "tmtAB_t2")]))
colSums(is.na(TMT_wide_scaled[, c("tmtB_t3", "tmtAB_t3")]))
colSums(is.na(TMT_wide_scaled[, c("tmtB_t4", "tmtAB_t4")]))


#combine TMT and BT data
#add PID2 to BT latent scores
BT_scores_df <- as.data.frame(BT_scores)

BT_scores_df$PID2 <- BT_wide_scaled$PID2
BT_scores_df$Code <- BT_wide_scaled$Code

BT_scores_df <- BT_scores_df %>%
  rename(
    bt_comp_t1 = F1,
    bt_comp_t2 = F2,
    bt_comp_t3 = F3
  )

#merge
final_data <- BT_scores_df %>%
  left_join(TMT_wide_scaled, by = "PID2")
colnames(final_data)

#check overlap at each timepoint
sum(complete.cases(final_data[, c("bt_comp_t1", "tmt_comp_t1")]))
sum(complete.cases(final_data[, c("bt_comp_t2", "tmt_comp_t2")]))
sum(complete.cases(final_data[, c("bt_comp_t3", "tmt_comp_t3")]))
sum(complete.cases(final_data[, c("bt_comp_t4", "tmt_comp_t4")]))


final_data <- final_data %>%
  mutate(
    overall_comp_t1 = rowMeans(cbind(bt_comp_t1, tmt_comp_t1), na.rm = TRUE),
    overall_comp_t2 = rowMeans(cbind(bt_comp_t2, tmt_comp_t2), na.rm = TRUE),
    overall_comp_t3 = rowMeans(cbind(bt_comp_t3, tmt_comp_t3), na.rm = TRUE)
  )


summary(final_data$overall_comp_t1)
summary(final_data$overall_comp_t2)
summary(final_data$overall_comp_t3)

colSums(is.na(final_data[, c("overall_comp_t1",
                             "overall_comp_t2",
                             "overall_comp_t3")]))
#final outcome is a hybrid: BT+TMT when possible, but mostly BT as missingness increases

#model how composite outcome changes over time (longitudinal analysis)
#reshape data to long format
final_long <- final_data %>%
  pivot_longer(
    cols = starts_with("overall_comp_t"),
    names_to = "time",
    values_to = "outcome"
  ) %>%
  mutate(
    time = as.numeric(gsub("overall_comp_t", "", time))
  ) %>%
  left_join(final_data %>% select(PID2, Code) %>% distinct(), by = "PID2")

#rename code column currently 'Code.y'
final_long <- final_long %>% 
  rename(Code = Code.y)

#run basic model of lme
model <- lmer(outcome ~ time + (1 | PID2), data = final_long)
summary(model)

#Does the composite outcome change over time differently across groups
model <- lmer(outcome ~ time * Code + (1 | PID2), data = final_long)
summary(model)

#include random slope for wobblers 
#drop random slope
model <- lmer(outcome ~ time * Code + (1 | PID2), data = final_long)
summary(model)

summary(model)
anova(model)
confint(model)


#get p values
library(lmerTest)

model <- lmer(outcome ~ time * Code + (1 | PID2), data = final_long)

summary(model)


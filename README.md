# Digital motor and cognitive variability as early markers of the prodromal Parkinson’s disease spectrum

# Project Overview
This repository contains the code used to analyse longitudinal digital cognitive and motor data from participants across the prodromal Parkinson's disease specturm (PSG-iRBD, pRBD with hyposmia and normosmic pRBD). THe primary aim of the project was to explore whether motor and cognitive deficits extend from iRBD to probable (pRBD) (with and without hyposmia), to determine whether longitudinal variability in performance can be detected in addition to changes in mean performance over time and to further test remote markers of motor and cognitive assessment. 

## Study Groups
Participants were classified into the following groups:
- Polysomnography-confirmed isolated REM sleep behaviour disorder (PSG-iRBD)
- Normosmic probable REM sleep behaviour disorder (pRBD)
- Probable REM sleep behaviour disorder with hyposmia (pRBD with hyposmia)
- Healthy controls

## Measures
## Digital Motor Assessments
- BRadykinesia Akinesia Incoordination (BRAIN) test
-  Distal Finger Tapping Test (DFT)

## Congitive Assessments
- Trail Making Test (TMT)

## Date cleaning
All data was cleaned and wrangled before performing the code in the scripts in this repository. 

## Statistical Analyses
The following analyses were performed:
- Baseline group comparisons
- Longitudinal linear mixed effects modelling
- Confirmatory factor analyis
- Intra-individual variability analysis
- Heteroscedastic mixed-effects modelling of variability

Intra-individual variability was quantified using Mean Squared Successive Differences (MSSD) across repeated assessments. 
Longitudinal trajectories were examined using linear mixed-effects models. Variability differences between groups and over time were assessed using heteroscedastic mixed-effects models, allowing residual variance to vary acorss groups and/or time.

## Analylis Scripts
The repository contains separate scripts for each primary analysis conducted in the manuscript.

| Script | Description |
|----------|-------------|
| 01_data_preparation.R | Data cleaning and derivation of analysis variables |
| 02_demographics.R | Demographic and clinical group comparisons |
| 03_BRAIN_baseline.R | Baseline BRAIN test analyses |
| 04_BRAIN_longitudinal.R | Longitudinal BRAIN test analyses |
| 05_DFT_baseline.R | Baseline DFT analyses |
| 06_DFT_longitudinal.R | Longitudinal DFT analyses |
| 07_TMT_baseline.R | Baseline TMT analyses |
| 08_TMT_longitudinal.R | Longitudinal TMT analyses |
| 09_figures_tables.R | Generation of manuscript figures and tables |



## Software Requirements
Analyses were performed using:
- R (verision 4.5.1)
- RStudio

Key packages include:

tidyverse

glm

lme4

lmerTest

nlme

emmeans

ggplot2

patchwork

lavaan

## Data Availability
Participant-level data is not publically available. This repository contains analysis code only.



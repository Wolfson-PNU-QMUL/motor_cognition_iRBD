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
| BRAIN test RBD vs controls baseline analysis.R | Baseline BRAIN test analysis in PSG-iRBD vs controls |
| DFT RBD vs controls baseline analysis.R | Baseline DFT analysis in PSG-iRBD vs controls |
| TMT RBD vs controls baseline analysis.R | Baseline TMT analysis in PSG-iRBD vs controls |
| BRAIN test sRBD vs controls baseline analysis.R | Baseline BRAIN test analysis in sRBD vs controls |
| DFT sRBD vs controls baseline analysis.R | Baseline DFT analysis in sRBD vs controls |
| TMT sRBD vs controls baseline analysis.R | Baseline TMT analysis in sRBD vs controls |
| BRAIN test sRBD hyposmia vs controls baseline analysis.R | Baseline BRAIN test analysis in sRBD + hyposmia vs controls |
| DFT sRBD hyposmia vs controls baseline analysis.R | Baseline DFT analysis in sRBD + hyposmia vs controls |
| TMT sRBD hyposmia vs controls baseline analysis.R | Baseline TMT analysis in sRBD + hyposmia vs controls |
| BRAIN test RBD vs controls LME model analysis.R | Longitudinal BRAIN test analyses in PSG-iRBD vs controls |
| DFT RBD vs controls LME model analysis.R | Longitudinal DFT analyses in PSG-iRBD vs controls |
| TMT RBD vs controls LME model analysis.R | Longitudinal TMT analyses in PSG-iRBD vs controls |
| BRAIN test sRBD vs controls LME model analysis.R | Longitudinal BRAIN test analyses in sRBD vs controls |
| DFT sRBD vs controls LME model analysis.R | Longitudinal DFT analyses in sRBD vs controls |
| TMT sRBD vs controls LME model analysis.R | Longitudinal TMT analyses in sRBD vs controls |
| BRAIN test sRBD hyposmia vs controls LME model analysis.R | Longitudinal BRAIN test analyses in sRBD + hypsomia vs controls |
| DFT sRBD hyposmia vs controls LME model analysis.R | Longitudinal DFT analyses in sRBD + hypsomia vs controls |
| TMT sRBD hyposmia vs controls LME model analysis.R | Longitudinal TMT analyses in sRBD + hypsomia vs controls |
| CFA analysis code.R | Confirmatory factor analysis |
| BRAIN test RBD vs controls variance analysis.R | Variance analysis BRAIN test in PSG-iRBD vs controls |
| DFT RBD vs controls variance analysis.R | Variance analysis DFT in PSG-iRBD vs controls |
| TMT RBD vs controls variance analysis.R | Variance analysis TMT in PSG-iRBD vs controls |
| BRAIN test sRBD vs controls variance analysis.R | Variance analysis BRAIN test in sRBD vs controls |
| DFT sRBD vs controls variance analysis.R | Variance analysis DFT in sRBD vs controls |
| TMT sRBD vs controls variance analysis.R | Variance analysis TMT in sRBD vs controls |
| BRAIN test sRBD hyposmia vs controls variance analysis.R | Variance analysis BRAIN test in sRBD + hyposmia vs controls |
| DFT sRBD hyposmia vs controls variance analysis.R | Variance analysis DFT in sRBD + hyposmia vs controls |
| TMTsRBD hyposmia vs controls variance analysis.R | Variance analysis TMT in sRBD + hyposmia vs controls |

## Software Requirements
Analyses were performed using:
- R (verision 4.5.1)
- RStudio

Key packages include:

- tidyverse
- glm
- lme4
- lmerTest
- nlme
- emmeans
- ggplot2
- patchwork
- lavaan

## Data Availability
Participant-level data is not publically available. This repository contains analysis code only.

## Author
Harneek Chohan h.chohan@qmul.ac.uk

# AMI-EDRS Web Application

## Overview

AMI-EDRS is a web-deployable dynamic risk stratification tool for early prognostic assessment after acute myocardial infarction (AMI). The application reconstructs early organ-domain dynamics from routinely available ICU data, assigns early trajectory phenotypes, estimates post-landmark mortality risk, and visualizes individualized risk profiles.

The AMI-EDRS framework integrates:

- Host vulnerability
- Clinical severity
- Early organ-domain dysfunction signals
- C1-C5 trajectory phenotype
- Elastic-net dynamic risk prediction
- External recalibration when applied to external cohorts

## Repository structure

```text
AMI-EDRS-WebApp/
├── app.R
├── R/
│   ├── preprocess.R
│   ├── score_calculation.R
│   ├── trajectory_assignment.R
│   ├── prediction.R
│   └── plot_functions.R
├── model/
│   ├── AMI_EDRS_model.rds
│   ├── AMI_EDRS_preprocessing_params.rds
│   ├── AMI_EDRS_trajectory_centers.rds
│   └── AMI_EDRS_recalibration_params.rds
├── data/
│   ├── input_template.csv
│   ├── demo_patient_lowrisk.csv
│   └── demo_patient_highrisk.csv
├── www/
├── README.md
├── LICENSE
└── .gitignore
```

## Installation

Install the required R packages:

```r
install.packages(c(
  'shiny', 'shinydashboard', 'tidyverse', 'data.table',
  'glmnet', 'xgboost', 'ranger', 'pROC', 'PRROC',
  'ggplot2', 'plotly', 'DT', 'openxlsx'
))
```

## Run the app

```r
shiny::runApp()
```

or from the repository root:

```r
shiny::runApp('AMI-EDRS-WebApp')
```

## Input data

The app supports single-patient input and batch prediction. A template file is provided:

```text
data/input_template.csv
```

Two simulated demonstration cases are provided:

```text
data/demo_patient_lowrisk.csv
data/demo_patient_highrisk.csv
```

These demo files contain simulated values and do not represent real patients.

## Output

The app returns:

- Organ-domain dysfunction scores
- C1-C5 trajectory phenotype assignment
- Individualized mortality risk estimate
- Risk stratum
- Risk gauge visualization
- Early organ-domain dynamic plots
- Batch-level prediction table when applicable

## Data privacy

This repository does not include MIMIC-IV, eICU, or any patient-level clinical data. Users should not upload protected health information or identifiable patient data to any public server.

## Research-use disclaimer

This application is intended for research demonstration only. It is not a medical device and should not be used as the sole basis for clinical decision-making. Clinical use requires prospective validation, local calibration, workflow integration, governance approval, and continuous performance monitoring.

## Citation

If you use AMI-EDRS, please cite the corresponding manuscript:

> AMI-EDRS: A Web-Deployable Trajectory-Enhanced Dynamic Risk Model for Early Prognostic Stratification After Acute Myocardial Infarction.

## License

This repository is released under the MIT License unless otherwise specified.

# Code and data for "Vector-borne and zoonotic disease risk twice as high for marginalized groups in high-income countries"

## Authors
Kyle M Shanebeck, Terrell J Sipin, Ashley Hazel, Skylar Hopkins, Andrew J MacDonald

Contact author: Kyle M Shanebeck (kyle.m.shanebeck@gmail.com)

## Abstract
As global change accelerates, the front lines of infectious disease risk are shifting into wealthy northern nations where entrenched social inequities compromise public health response. Identifying and addressing these vulnerabilities is a critical priority in this era of rapid disease emergence. We provide a global meta-analysis of differences between marginalized groups and the general population for vector-borne and zoonotic disease risks in wealthy, populous, high latitude countries. Marginalized groups have 2.2 times higher risk of exposure, 1.6 times higher risk of infection, and 1.8 times higher risk of worse health outcomes due to substandard healthcare. We quantify where and why disparities in risk occur, identify potential epicenters of current and future risk, and highlight targeted interventions to protect community health, economic stability, and national security.

## Workflow

### 1_rmaMV.R
Reads in study data built and cleaned by Kyle M Shanebeck, `Data/1_DataProcessing/cleaned_data_6_11_25_data_cleaned_collapse_risk_groups.csv`, and creates single and double factor multivariate linear mexed-effects models. Single models measure the effect sizes as a function of a single study-level variable, or moderator, as the predictor. On the other hand, double factor models measure the effect sizes as a function of the interaction of two moderators as the predictor. Both model types include the following random effects: a unique article ID (`RecordNo`), a unique ID associated with both the article and observation group (`groupID`), and odds ratio method (`OR_method`). Finally, each model is fitted via restricted maximum likelihood estimation (REML). Refer to `Data_variables_explanation.docx` for more details on the input data.

### 2_model_stats.R
For each model, computes fit statistics: 
* AIC
* $I^2$ (total, and for each random effect)
* $R^2$ (marginal and conditional)
* QE
* QM
* Mean $SE^2$
* $\tau^2$ (total, and for each random effect)
* Total variance

### 3_model_plots.R
Creates funnel plots and orchard plots for all models. Funnel plot types:
* Standard error
* Sampling variance
* Inverse standard error
* Inverse sampling variance

### 4_MS_orchards.R
Creates specific orchard plots with more customization that would be added to the main text or that were requested by Kyle M Shanebeck. Models include: 
* Marginalized group and group proxy
* Risk group and risk proxy
* Marginalized group and risk group
* Transmission and broad region (not in main text)

### 5_beta_CI_pval.R
Outputs individual and combined (2x3 and 3x2) log odds ratio plots for:

* Overall effect (intercept only) + marginalized group  
* Risk group
* Transmission
* Country
* Data source
* Scale

Also constructs a table of the point estimates, confidence intervals, and p-values for each group.

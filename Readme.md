Simulated example dataset and sample code for calculating the AUC in the presence of dependent censoring caused by treatment, using auxiliary information from multiple longitudinal risk factors. 

See analysis.R for application of the codes to a simulated dataset.

Zhang C, Ning J, Belle S H, Squires R H, Cai J, Li R. (2022+) Assessing Predictive Discrimination Performance of Biomarkers in the Presence of Treatment-induced Dependent Censoring. Journal of the Royal Statistical Society. Series C. Accepted.








AUCdep: calculate the net AUC for biomarkers after removing the disturbance of treatment-induced dependent censoring.


Description
------------

This is a simulated data example used to illustrate the idea and calculation of net AUC, sensitivity and specificity.

The user needs to download the 'Data.RData', the 'function_v1.0' and the 'jm_3Y2T.stan' in 'function' folder, and then run the 'analysis_v1.0' argument to have an intuitive understanding of the calculation. 

The user could modify the 'function_v1.0' and 'jm_3Y2T.stan' to accommodate different scenarios in joint modeling.

Variable list in Example.RData
------------
- id: Patient id

- time: Visit time

- x1: Baseline covariate, continuous

- x2: Baseline covariate, binary

- Y1: Continuous longitudinal risk facror

- Y2: Continous longitudinal risk factor

- Y3: Binary longitudinal risk factor

- obsT: Observed event time

- status: Outcome status

- delta: Type of first occuring event, 0 alive, 1 death and 2 dependent censoring



Analysis steps
------------

- Step 0: load requied packages. 

- Step 1: load required data and split as train and test set.

- Step 2: Prepare data for stan.

- Step 3: Jointly model three longitudinal risk factors, the survival outcome and the dependent censoring.

- Step 4: Estimate net AUC, sensitivity and specificity of Yvar in predicting death before T.end among subjects who are still alive at landmark time T.start.

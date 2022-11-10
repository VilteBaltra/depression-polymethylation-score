# Depression Polymethylation Score

This repository contains scripts for the derivation of depression polymethylation scores using BioMM. 

Data:
- ALSPAC children (phenotypic and genome-wide DNAm data)

Depression variable: 
- The Short Mood and Feelings Questionnaire (SMFQ) at 17.5y. SMFQ measures core depressive symptomology in children and adolescents (score range 0-26). 

Scripts:

1. '0.sample_preparation_forBioMM.R' selects individuals with SMFQ at 17.5y and DNAm at 17y. Saves selected sample with methylation in 'aries_out_17y_20221108.rds' and corresponding sample sheet as 'samplesheet_DNAm17y_SMFQ17.5y_20221108.csv'. 
2. '1.SVA_forBioMM.R' obtains SVs and saves them in 'pheno_SVs_20221108.rds' 
3. '2.mean_impute.R' removes CpGs with > 30% missingness and performs mean imputation on remaining CpGs. Saves imputed meth data as 'beta_mean_imputed_20221108.rds' 
4. '3.residualise_imp_data.R' regresses sex, age, cell type and 10 SVs from imputed DNAm data (beta_mean_imputed_20221108.rds) and saves regressed output in 'df.residuals_imp_data_20221108.rds'.
5. '4.format_data_forBioMM.R' formats data for BioMM and saves it as 'methylData_residualised.imputed.forBioMM_20221108.rds'.
6. '5.BioMM_DNAm_17y_SMFQ_17.5y.R' derives the polymethylation score (PMS) for depression and identifies implicated biological pathways.  the analysis on continuous SMFQ variable.


This project is a collaboration between COMMITMENT and the EarlyCause consortium.

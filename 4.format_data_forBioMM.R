# Date: 2022/11/08
# Author: VB
# Purpose: this script formats input data for BioMM 
# Output: 'methylData_residualised.imputed.forBioMM_20221108.rds' to be replaced with 'df.residuals_imp_data_nosexmismatch_20221108.rds'

# _______________________________________________#
#                    SET UP
# _______________________________________________#
# load libraries
library(tidyverse)
# set paths
setwd("/campaign/VB-FM5HPC-001/")
output = "/campaign/VB-FM5HPC-001/BioMM/nimbus_redo/output/"

# _______________________________________________#
#                   READ IN DATA
# _______________________________________________#
## read in imputed and residualised DNAm data
beta_imp <- readRDS(paste0(output, 'df.residuals_imp_data_nosexmismatch_20221108.rds')) # cpg in columns
# read in depression (SMFQ age 17.5y) data
pheno <- readRDS(paste0(output, 'pheno_SVs_20221108.rds')) # this will have 3 more individuals than beta_imp, 
# as in beta_imp sex mismatches were removed

# remove second-born siblings (mult == 1)
summary(as.factor(pheno$mult))
pheno <- pheno %>% filter(mult %in% c(0, NA))

# check dimentions of pheno data without siblings
dim(pheno)
summary(as.factor(pheno$mult))
# check dimentions of meth data
dim(beta_imp)
head(beta_imp[,1:4])

# _______________________________________________#
#                   FORMAT DATA
# _______________________________________________#
# merge methylation data with SMFQ variable
methylData<- merge(beta_imp, pheno[, c('Sample_Name', 'CCXD917')], by.x="row.names", by.y="Sample_Name")
dim(methylData) # two extra columns because Row.names was added as a column and CCXD917
methylData[1:6,1:6]

# rename CCXD917 to 'label' (expected by BioMM), add samples names as row names, and remove extra column with row.names "select(-1)"
row.names(methylData) <- methylData$Row.names
methylData <- methylData %>% rename(label = CCXD917)  %>%  select(-1)
methylData[1:6,1:6]
dim(methylData)

# change from factor to numeric
class(methylData$label)
methylData$label <- as.numeric(methylData$label)

summary(as.factor(methylData$label)) # ensure there are no NAs. If NAs present run:
# methylData <- methylData[methylData[,'label'] %in% c(0:26), ]
# dim(methylData)

# move label column to be fist column
methylData <- methylData %>% select(label, everything())
methylData[1:6,1:6]

# _______________________________________________#
#                    SAVE DATA
# _______________________________________________#

#Save formatted methyData 
saveRDS(methylData, file = paste0(output, 'methylData_residualised.imputed.forBioMM_20221108.rds'))

# ___________________end script__________________#


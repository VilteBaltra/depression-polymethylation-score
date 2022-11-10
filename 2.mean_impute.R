# Date: 2022/11/08
# Author: VB
# Purpose: This script removes CpGs with > 30% missingness and performs mean imputation on remaining CpGs. 
# The imputed methylation data is saved as 'beta_mean_imputed_20221108.rds'
# Output: mean imputed data is saved in 'beta_mean_imputed_20221108.rds'

# _______________________________________________#
#                    SET UP
# _______________________________________________#

# set paths
setwd("/campaign/VB-FM5HPC-001/")
output = "/campaign/VB-FM5HPC-001/BioMM/nimbus_redo/output/"

# _______________________________________________#
#                   READ IN DATA
# _______________________________________________#
# read in data
aries <- readRDS(paste0(output, 'aries_out_17y_20221108.rds'))
beta <- as.data.frame(aries$meth)


# _______________________________________________#
#                 PERCENT MISSING
# _______________________________________________#

### Calculate percent missing 
# per cpg
# added -1 to account for the NAs introduced as part of beta['miss_ind',] <- NA
percent_missing <- function(var) { (sum(is.na(var))-1) / (length(var)-1) * 100 } 

beta$miss_cpg <- NA

val = 1 # note: this loop takes way too long... there should be a quicker way to do it
for (i in row.names(beta)) {
  beta$miss_cpg[val] <- percent_missing(beta[i,])
  val=val+1
}

# per individual
IDs <- names(beta) 

beta['miss_ind',] <- NA

for (i in IDs) {
  beta['miss_ind',][i] <- percent_missing(beta[,i])
}

# remove individuals and cpgs with > 30% missingnesss

# number of CpG with missingness over 30%, 35%, 40%, 50%
sum(beta$miss_cpg > 30) # 221
sum(beta$miss_cpg > 35) # 207
sum(beta$miss_cpg > 40) # 184
sum(beta$miss_cpg > 50) # 82

# remove CpGs with >30% missingness
beta_nocpg <- beta[beta[,'miss_cpg'] <= 30,]
dim(beta_nocpg) # 450618   2342


# number of individuals with missingness over 30%
sum(beta_nocpg['miss_ind',] > 30) # 0 

# remove individuals with >30% missingness (does nothing, as none present)
beta_nocpg_noind <- beta_nocpg[,beta_nocpg['miss_ind',] <= 30]

# dimentions of pruned dataset for missingness 
dim(beta_nocpg_noind) # 450618   2342

nr_missing <- sum(is.na(beta_nocpg_noind))

nr_non_missing <- sum(!is.na(beta_nocpg_noind))

#percentage missingness
cat("total percentage of missing CpGs is", nr_missing / nr_non_missing * 100,"%")
# total percentage of missing CpGs is 0.3652105 % (i.e., < 1%) 
# _______________________________________________#
#                MEAN IMPUTATION
# _______________________________________________#

# mean imputation of multiple rows
sum(is.na(beta_nocpg_noind)) # 3840215
dim(beta_nocpg_noind) # 450618   2342

# remove individual missingness row 
beta_nocpg_noind <- beta_nocpg_noind[ !(row.names(beta_nocpg_noind) %in% c('miss_ind')), ]

# remove cpg missigness column
beta_nocpg_noind <- beta_nocpg_noind[, !(colnames(beta_nocpg_noind) %in% c('miss_cpg'))]

# should have one less row and one less column than before
dim(beta_nocpg_noind) # 450617   2341
sum(is.na(beta_nocpg_noind)) #  3840215

# impute data using row means
ind <- which(is.na(beta_nocpg_noind), arr.ind=TRUE)
beta_nocpg_noind[ind] <- rowMeans(beta_nocpg_noind, na.rm=TRUE)[ind[,1]]

# check nr of NAs
sum(is.na(beta_nocpg_noind)) # should have no NAs left 
# 0 
dim(beta_nocpg_noind) # should have same dimentions
# 450617   2341
# _______________________________________________#
#                   SAVE DATA
# _______________________________________________#

saveRDS(beta_nocpg_noind, file = paste0(output, 'beta_mean_imputed_20221108.rds'))

# ___________________end script__________________#


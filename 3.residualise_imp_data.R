# Date: 2022/11/08
# Author: VB
# Purpose: This script regresses sex, age, cell type and 10 SVs from imputed DNAm data ('beta_mean_imputed_20221108.rds') 
# Output: regressed output saved in 'df.residuals_imp_data_nosexmismatch_20221108.rds'

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
beta_imp <- readRDS(paste0(output, 'beta_mean_imputed_20221108.rds'))
pheno <- readRDS(paste0(output, 'pheno_SVs_20221108.rds'))

# check dimentions
dim(beta_imp)
dim(pheno)

# check no missingness
sum(is.na(beta_imp))
class(beta_imp)

# _______________________________________________#
#             REMOVE SEX MISMATCH
# _______________________________________________#

# remove phenotypic ('kz021') and DNAm sex ('Sex') mismatches 
# first recode Sex to numeric (1 = Males, 2 = Females; same as kz021)
pheno$Sex2 <- ifelse(pheno$Sex == 'M', 1, 2)
# identify sex mismatches
indices <- which(pheno$kz021 != pheno$Sex2) # there are three mismatches 
# get individuals IDs without the mismatches
ids <- pheno$Sample_Name[indices]
# remove them from pheno and beta_imp datasets
pheno <- pheno[!(pheno$Sample_Name %in% ids), ]
beta_imp <- beta_imp[, !(names(beta_imp) %in% ids)] # this should work if CpGs are in rows and IDs in columns

# check dimensions (should be 3 less rows in both dataset)
dim(pheno)
dim(beta_imp)

# _______________________________________________#
#                  RESIDUALISE
# _______________________________________________#

df.residuals <- apply(beta_imp,1,function(x){residuals(lm(unlist(x)~pheno$kz021+
                                                            pheno$CD8T+pheno$CD4T+pheno$NK+pheno$Bcell+pheno$Mono+pheno$Gran+
                                                            pheno$age+
                                                            pheno$V1 + pheno$V2 + pheno$V3 + pheno$V4 + pheno$V5 +
                                                            pheno$V6 + pheno$V7 + pheno$V8 + pheno$V9 + pheno$V10))})

# _______________________________________________#
#                   SAVE DATA
# _______________________________________________#

saveRDS(df.residuals, file = paste0(output, 'df.residuals_imp_data_nosexmismatch_20221108.rds'))

# ___________________end script__________________#




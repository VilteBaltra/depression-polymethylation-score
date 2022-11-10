# Date: 2022/11/08
# Author: VB
# Purpose: This script derives surrogate variables (SVs) to be used for BioMM 
# Output: SVs are saved in 'pheno_SVs_20221108.rds'

# _______________________________________________#
#                   SVA SET UP
# _______________________________________________#
# install.packages("remotes")
# remotes::install_github("perishky/meffil")

# source("http://bioconductor.org/biocLite.R")
# install.packages("devtools") # if the devtools package is not installed
# library(devtools)
# install_github("perishky/meffil")

library(meffil) #for SVA

# set paths
output = "/campaign/VB-FM5HPC-001/BioMM/nimbus_redo/output/"

# _______________________________________________#
#                   READ IN DATA
# _______________________________________________#

# read in samplesheet obtained in '0.sample_preparation_forBioMM.R' script
samplesheet <- read.csv(paste0(output, "samplesheet_DNAm17y_SMFQ17.5y_20221108.csv"))
#aries <- readRDS("/campaign/VB-FM5HPC-001/BioMM/methylData_residualised.imputed.sub_forBioMM.rds")
aries <- readRDS(paste0(output, 'aries_out_17y_20221108.rds'))


impute.matrix <- function(x, margin=1, fun=function(x) mean(x, na.rm=T)) {
  if (margin == 2) x <- t(x)
  
  idx <- which(is.na(x) | !is.finite(x), arr.ind=T)
  if (length(idx) > 0) {
    na.idx <- unique(idx[,"row"])
    v <- apply(x[na.idx,,drop=F],1,fun) ## v = summary for each row
    v[which(is.na(v))] <- fun(v)      ## if v[i] is NA, v[i] = fun(v)
    x[idx] <- v[match(idx[,"row"],na.idx)] ##
    stopifnot(all(!is.na(x)))
  }
  
  if (margin == 2) x <- t(x)
  x
}

random.seed <- set.seed(83)

# meth is your methylation data (CpGs in rows)
beta <- aries$meth
beta[1:10, 1:10]
dim(beta)

#smfq (CCXD917) is your outcome (depr)
variable <- samplesheet[,c("CCXD917")]

# add all your covariates here, which are in an object called samplesheet
# we probably want: age at DNAm, sex, cell type (any panel will do; whatever is alreayd available)
# covariates <- pheno[,c("salas.Bcell","salas.CD4T","salas.CD8T","salas.Gran","salas.Mono","salas.NK",
#                        "age.DNAm","child.sex")]
covariates <- samplesheet[,c("Bcell","CD4T","CD8T","Gran","Mono","NK",
                             "age","kz021")] # age = age.DNAm, Sex = sex.DNAm, kz021 = phenotypic sex. 
#I will use phenotypic sex but then check alignment with DNAm sex and remove mismatches

summary(as.factor(samplesheet$kz021))
#1    2 NA's 
#1084 1256    1
# we have 1084 males, 1256 females and 1 NA   1

summary(as.factor(samplesheet$Sex))
# F    M 
# 1255 1086
#according to DNAm Sex we have 1255 females and 1086 males 
# we can remove the 1 mismatch and 1 missing one

# I think you need to reduce down to probes common to both the 450k and epic array (age 17y only)
featureset="common"
features <- meffil.get.features(featureset)

stopifnot(length(rownames(beta)) > 0 && all(rownames(beta) %in% features$name))
stopifnot(ncol(beta) == length(variable))
stopifnot(is.null(covariates) || is.data.frame(covariates) && nrow(covariates) == ncol(beta))

original.variable <- variable
original.covariates <- covariates


sample.idx <- which(!is.na(variable))
if (!is.null(covariates))
  sample.idx <- intersect(sample.idx, which(apply(!is.na(covariates), 1, all)))

cat("Removing", ncol(beta) - length(sample.idx), "missing case(s).")


beta <- beta[,sample.idx]
variable <- variable[sample.idx]

covariates <- covariates[sample.idx,,drop=F]

# convert sex from character to numeric [already numeric]
# covariates$Sex <- ifelse(covariates$Sex == 'M', 1, 2) # male = 1, female = 2

pos.var.idx <- which(apply(covariates, 2, var, na.rm=T) > 0)
cat("Removing", ncol(covariates) - length(pos.var.idx), "covariates with no variance.")
covariates <- covariates[,pos.var.idx, drop=F]


covariate.sets <- list(none=NULL)
if (!is.null(covariates))
  covariate.sets$all <- covariates


surrogates.ret <- NULL 
beta.sva <- beta

autosomal.sites <- meffil.get.autosomal.sites(featureset)
autosomal.sites <- intersect(autosomal.sites, rownames(beta.sva))

most.variable <- length(autosomal.sites)

beta.sva <- beta.sva[autosomal.sites,]

var.idx <- order(rowVars(beta.sva, na.rm=T), decreasing=T)[1:most.variable]
beta.sva <- impute.matrix(beta.sva[var.idx,,drop=F])

cov.frame <- model.frame(~., data.frame(covariates, stringsAsFactors=F), na.action=na.pass)
mod0 <- model.matrix(~., cov.frame)
mod <- cbind(mod0, variable)


set.seed(random.seed)
sva.ret <- sva(beta.sva, mod=mod, mod0=mod0, n.sv=10)


#Check SVs aren't associated with Trait (you need to remove those SVs that associate with depr)
SVs<-sva.ret$sv
SV_check<-apply(SVs,2,function(x) summary(lm(x~variable))$coef[2,]) # none
SV_check #[updated output included below, no SV sig. associates with depr]
# [,1]         [,2]         [,3]         [,4]         [,5]
# Estimate   -3.312006e-05 2.019516e-05 3.760080e-05 6.396659e-05 9.771601e-06
# Std. Error  8.186794e-05 8.186974e-05 8.186711e-05 8.186011e-05 8.187055e-05
# t value    -4.045547e-01 2.466743e-01 4.592907e-01 7.814135e-01 1.193543e-01
# Pr(>|t|)    6.858419e-01 8.051819e-01 6.460682e-01 4.346385e-01 9.050049e-01
# [,6]          [,7]          [,8]          [,9]
# Estimate   -3.386447e-05 -4.297830e-05 -3.016877e-05 -1.092242e-04
# Std. Error  8.186781e-05  8.186598e-05  8.186842e-05  8.183962e-05
# t value    -4.136482e-01 -5.249836e-01 -3.685031e-01 -1.334613e+00
# Pr(>|t|)    6.791697e-01  5.996444e-01  7.125315e-01  1.821332e-01
# [,10]
# Estimate   -1.068957e-04
# Std. Error  8.184094e-05
# t value    -1.306139e+00
# Pr(>|t|)    1.916337e-01

row.names(SVs)=samplesheet[sample.idx,c("Sample_Name")]
pheno=merge(samplesheet,SVs, by.x="Sample_Name", by.y="row.names", all.x=T)
# should the be all.x = F to remove individuals with all NA for SVs?

# convert sex from character to numeric (this applies if we use DNAm 'Sex' instead of 'kz021')
#pheno$Sex <- ifelse(pheno$Sex == 'M', 1, 2) # male = 1, female = 2

# _______________________________________________#
#                   SAVE DATA
# _______________________________________________#

saveRDS(pheno, file = paste0(output, 'pheno_SVs_20221108.rds'))

# ___________________end script__________________#




# Date: 2022/11/08
# Author: VB
# Purpose: this script selects ALSPAC sample with DNAm at 17y and SMFQ at 17.5y
# Output: 'aries_out_17y_20221108.rds' and 'samplesheet_DNAm17y_SMFQ17.5y_20221108.csv'

# useful links:
# https://github.com/MRCIEU/aries
# http://htmlpreview.github.io/?https://github.com/MRCIEU/aries/blob/master/docs/tutorial/tutorial.html

# _______________________________________________#
#                       SET UP
# _______________________________________________#

# load packages 
# library(devtools)
# install_github("MRCIEU/aries")
library(aries)
library(foreign)
library(tidyverse)

# set paths
setwd("/campaign/VB-FM5HPC-001/")
aries.dir="woc_release"
scripts = "/campaign/VB-FM5HPC-001/BioMM/nimbus_redo/scripts/"
output = "/campaign/VB-FM5HPC-001/BioMM/nimbus_redo/output/"
pheno = "/campaign/VB-FM5HPC-001/ALSPAC_pheno_data/"

# source functions
source(paste0(scripts, 'aries.function.R'))


# _______________________________________________#
#                READ IN METH DATA
# _______________________________________________#

aries <- aries.select.mod(aries.dir, time.point="15up", featureset = 'common') # 'F7' for age 7 data and '15up' for age 17 
# (without meth for now)

# _______________________________________________#
#  SELECT SAMPLE WITH DNAm at 17y and SMFQ 17.5y
# _______________________________________________#

# get samples with SMQF data
samplesheet = aries$samples
samplesheet$cidB2957.qlet=paste0(samplesheet$cidB2957, samplesheet$QLET)
dim(samplesheet)
# 2857   31

# pheno is your dataframe, pruned down to those with depr info
pheno <- read.spss(paste0(pheno,'EarlyCause_AHupdated_CIDB2957_21OCT21.sav'), to.data.frame = TRUE, use.value.labels = FALSE)

#Select variables of Interest (i.e., ID, sibling status, twin status, depression outcome)
pheno_subset <- pheno[, c('cidB2957', 'qlet', 'mult', 'CCXD917', 'kz021', 'tc9991a')] 
rm(pheno)
# CCXD917 = SMFQ at age 17.5y (DV: Moods and Feelings total score)
# kz021 = sex, 1 = Male; 2 = Female.
# tc9991a = DV: Age of study teenager at completion (months)

pheno_subset$cidB2957.qlet <- paste0(pheno_subset$cidB2957, pheno_subset$qlet)
pheno_subset$cidB2957.qlet <- gsub(" ", "", pheno_subset$cidB2957.qlet) # remove space
dim(pheno_subset)
# 15645     7

# subset to those with SMFQ 17.5y info
pheno_subset_SMFQ <- pheno_subset %>% filter(!is.na(CCXD917)) 
dim(pheno_subset_SMFQ)
# 4498    7

# view SMFQ distribution at 17.5y
summary(as.factor((pheno_subset_SMFQ$CCXD917))) # no NA left 
# 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19 
# 308 344 410 475 415 367 334 284 237 180 171 155 116 287  70  71  47  40  26  27 
# 20  21  22  23  24  25  26 
# 29  25  19  17  16  11  17 


samplesheet <- merge(samplesheet, pheno_subset_SMFQ, by="cidB2957.qlet", all = F)
dim(samplesheet)
# 2369   37

samplesheet = samplesheet[samplesheet$qlet.x == "A",] # only keep first born child of twin pair
dim(samplesheet)
# 2341   37

# run once again to only read in sample with DNAm at age 17 and depr
aries <- aries.select.mod(aries.dir, time.point="15up", featureset = 'common', sample.names=samplesheet$Sample_Name)

# re-create samplesheet with pheno and covars based on the new DNAm data
rm(samplesheet)
samplesheet = aries$samples
samplesheet$cidB2957.qlet=paste0(samplesheet$cidB2957, samplesheet$qlet)
samplesheet = merge(samplesheet, pheno_subset_SMFQ, by="cidB2957.qlet", all = F)
samplesheet = merge(samplesheet, aries$cell.counts$`blood-gse35069`,
                    by.x="Sample_Name", by.y="row.names",
                    all = F)

# read in the methylation data for the selected 2341 individuals
aries$meth <- aries.methylation(aries)

# _______________________________________________#
#                   SAVE DATA
# _______________________________________________#
# save samplesheet for 2341 individuals with 17y DNAm and 17.5y SMFQ data
write.csv(samplesheet, file = paste0(output,'samplesheet_DNAm17y_SMFQ17.5y_20221108.csv'), row.names = F, quote = F)
# to read it in type: samplesheet <- read.csv("samplesheet_DNAm17y_SMFQ17.5y_20221108.csv")

# save 17y meth data together with the samplesheet
saveRDS(aries, file = paste0(output, 'aries_out_17y_20221108.rds'))

# ___________________end script__________________#



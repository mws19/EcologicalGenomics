# Morgan Southgate
# Ecological Genomics HW 1
# March 3, 2018

################################# Set-Up for Analysis#############################
setwd("~/EcologicalGenomics/Transcriptomics/HW1")

## load two necessary packages
library("DESeq2")
library("ggplot2")

## load data - noIT indicates Italy (native range) has been deleted from dataset
countsTable <- read.delim('allcountsdataRN_noIT.txt', head=T, stringsAsFactors=T, row.names=1)
head(countsTable)
tail(countsTable)

## convert data to matrix and check conversion
countData <- as.matrix(countsTable)
head(countData)
tail(countData)

## read in data about 48 samples
conds <- read.delim("cols_data_noIT.txt", header=T, stringsAsFactors = T, row.names=1)
head(conds)

## turn it into data frame
colData <- as.data.frame(conds)
head(colData)

########## Model 2: Sex-Specific Gene Expression Analysis ##################################### 

## add a new columns called group, that is the population + devstage + sex
## add columns to a data frame using the paste function
colData$group <- factor(paste0(colData$population, "-", colData$devstage, "-", colData$sex))
head(colData)

# create a new factor called group - so each one will be treated uniquely
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design= ~groupMA)
dds <- dds[rowSums(counts(dds))> 1, ]

dds <- DESeq(dds, parallel=T)

resultsNames(dds)
#[1] "Intercept"     "groupNC.AD4.F" "groupNC.AD4.M" "groupNC.L3L.F" "groupNC.L3L.M"
# [6] "groupNC.PD1.F" "groupNC.PD1.M" "groupNC.PP1.F" "groupNC.PP1.M" "groupWA.AD4.F"
# [11] "groupWA.AD4.M" "groupWA.L3L.F" "groupWA.L3L.M" "groupWA.PD1.F" "groupWA.PD1.M"
# [16] "groupWA.PP1.F" "groupWA.PP1.M"

#############Females############################################################
# Contrast females between two populations
res_pop.F <- results(dds,contrast=list(c("groupWA.L3L.F","groupWA.PP1.F","groupWA.PD1.F","groupWA.AD4.F"),c("groupNC.L3L.F","groupNC.PP1.F","groupNC.PD1.F","groupNC.AD4.F")),listValues = c(1/2,-1/2), alpha=0.05)

res_pop.F <- res_pop.F[order(res_pop.F$padj),]
head(res_pop.F)


#OTAU014798-RA  22.29407      3.6826741 0.6045166  6.091933 1.115557e-09
# OTAU014765-RA 257.12354      1.9643143 0.3374607  5.820868 5.854292e-09
# OTAU008015-RA 988.76557     -5.8136338 1.0347763 -5.618252 1.928994e-08
# OTAU013946-RA 946.74326      0.9005703 0.1609662  5.594778 2.209043e-08
# OTAU010901-RA  54.68902     -2.0663231 0.3727966 -5.542763 2.977350e-08
# OTAU017482-RA 126.29161     -6.8142306 1.2484546 -5.458132 4.811687e-08
# padj
# <numeric>
#   OTAU014798-RA 1.616219e-05
# OTAU014765-RA 4.240849e-05
# OTAU008015-RA 8.001152e-05
# OTAU013946-RA 8.001152e-05
# OTAU010901-RA 8.627168e-05
# OTAU017482-RA 1.161862e-04

#look at which ones are significant
# summary(res_pop.F)
# out of 16851 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 49, 0.29% 
# LFC < 0 (down)   : 48, 0.28% 
# outliers [1]     : 135, 0.8% 
# low counts [2]   : 2228, 13% 
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results



################ MA Plot ########################### 
plotMA(res_pop.F, main="WA vs NC Females", ylim=c(-2,2))
abline(h=c(-1,1), col="blue",lwd=2)

### ################ Write .csv outputs ##########
write.csv(res_pop.F,file="dge_females.csv", row.names=T, quote=F)
head(res_pop.F)

##### change pvalue to -log(pvalue) and export as .csv with rownames
neglogpval <- as.matrix(-log(res_pop.F$pvalue))
head(neglogpval)

res_pop_negpval <- cbind(row.names(res_pop.F), neglogpval)
head(res_pop_negpval)

colnames(res_pop_negpval)=c("gene","neglogpval")

write.csv(res_pop_negpval,file="dge_f_neglogpval.csv",row.names=F, quote=F)

############### Heat Maps #############
# pull out significance genes to be able to make a heat map
sig_popF <- res_pop.F[which(res_pop.F$padj <0.05), ] # great way to subset data
dim(sig_popF) # check

sig_popF_df<- as.data.frame(sig_popF)
sig_popF_df$Row.names <- rownames(sig_popF_df)
dim(sig_popF_df)

genesOfInterest_popF <- c(sig_popF_df$Row.names)
length(genesOfInterest_popF)

# variance stabilization
vsd <- vst(dds, blind=F)

#####
dds$combined = factor(paste0(dds$population, "-", dds$devstage, "-", dds$sex))
dds$combined <- factor(dds$combined, levels=c("WA-L3L-F","WA-L3L-M","WA-PP1-F","WA-PP1-M","WA-PD1-F","WA-PD1-M","WA-AD4-F","WA-AD4-M","NC-L3L-F","NC-L3L-M","NC-PP1-F","NC-PP1-M","NC-PD1-F","NC-PD1-M","NC-AD4-F","NC-AD4-M"), labels=c("WA-L3L-F","WA-L3L-M","WA-PP1-F","WA-PP1-M","WA-PD1-F","WA-PD1-M","WA-AD4-F","WA-AD4-M","NC-L3L-F","NC-L3L-M","NC-PP1-F","NC-PP1-M","NC-PD1-F","NC-PD1-M","NC-AD4-F","NC-AD4-M"))

# calculating means 
baseMeanPerGrp <- sapply( levels(dds$combined), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$combined == lvl] ) )

head(baseMeanPerGrp)
dim(baseMeanPerGrp)

# pulls out normalized counts (avg of 3 reps) for all of our significant genes 
m <- baseMeanPerGrp[genesOfInterest_popF, c("WA-L3L-F","WA-PP1-F","WA-PD1-F","WA-AD4-F", "NC-L3L-F","NC-PP1-F","NC-PD1-F","NC-AD4-F")] 
head(m)
dim(m)

# use the apply function to scale the matrix m: apply(object, rows, function)
mat_scaled = t(apply(m, 1, scale))
head(mat_scaled)

# make a heat map!
library(pheatmap)
pheatmap(mat_scaled, labels_col=c("WA-L3L-F","WA-PP1-F","WA-PD1-F","WA-AD4-F", "NC-L3L-F","NC-PP1-F","NC-PD1-F","NC-AD4-F", cluster_cols=T, cluster_rows=T))

################# GO Analysis ################
# input ="dge_f_neglogpval.csv"
# goAnnotations="gene_annotation_only.tab" 
# goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
# goDivision="MF" # either MF, or BP, or CC
# source("gomwu.functions.R")
# 
# gomwuStats(input, goDatabase, goAnnotations, goDivision,
#            perlPath="perl",largest=0.1,smallest=5,clusterCutHeight=0.25,Alternative="g")
# 
# # Plotting results
# library(ape)
# quartz()
# gomwuPlot(input,goAnnotations,goDivision,
#           absValue=-log(0.05,10),level1=0.1,level2=0.05,level3=0.01,
#           txtsize=1.2, treeHeight=0.5)
# 

####################### Males ##################################################################
res_pop.M <- results(dds, contrast=list(c("groupWA.PP1.M","groupWA.AD4.M","groupWA.L3L.M","groupWA.PD1.M"),c("groupNC.PP1.M","groupNC.AD4.M","groupNC.L3L.M","groupNC.PD1.M")),listValues = c(1/2,-1/2), alpha=0.05)

res_pop.M <- res_pop.M[order(res_pop.M$padj),]
head(res_pop.M)
# #                baseMean log2FoldChange     lfcSE      stat       pvalue
# <numeric>      <numeric> <numeric> <numeric>    <numeric>
#   OTAU017482-RA  126.2916      -8.256266 1.2419933 -6.647593 2.979252e-11
# OTAU012716-RA  188.8777       7.491769 1.2183487  6.149117 7.791530e-10
# OTAU005403-RA  340.2645      -1.158669 0.1935572 -5.986182 2.148235e-09
# OTAU014528-RA   29.9594      -7.290107 1.2715902 -5.733063 9.863294e-09
# OTAU009178-RA  113.8810      -3.366984 0.6303575 -5.341389 9.223731e-08
# OTAU002976-RA  998.4232      -1.656983 0.3242820 -5.109699 3.226720e-07
# padj
# <numeric>
#   OTAU017482-RA 4.124476e-07
# OTAU012716-RA 5.393297e-06
# OTAU005403-RA 9.913388e-06
# OTAU014528-RA 3.413686e-05
# OTAU009178-RA 2.553867e-04
# OTAU002976-RA 7.445118e-04

summary(res_pop.M)
# out of 16851 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 50, 0.3% 
# LFC < 0 (down)   : 46, 0.27% 
# outliers [1]     : 135, 0.8% 
# low counts [2]   : 2872, 17% 
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

############### MA PLot ##################### 
plotMA(res_pop.M, main="WA vs NC Males", ylim=c(-2,2))
abline(h=c(-1,1), col="blue",lwd=2)

############ Write .csv outputs ##################
write.csv(res_pop.M, file="dge_males.csv",row.names=T, quote=F)
head(res_pop.M)

# change pvalue to -log(pvalue) and export as .csv with rownames
neglogpval <- as.matrix(-log(res_pop.M$pvalue))
head(neglogpval)

res_pop_negpval <- cbind(row.names(res_pop.M), neglogpval)
head(res_pop_negpval)

colnames(res_pop_negpval)=c("gene","neglogpval")

write.csv(res_pop_negpval,file="dge_m_neglogpval.csv",row.names=F, quote=F)

########### Now males
# pull out significance genes to be able to make a heat map
sig_popM <- res_pop.M[which(res_pop.M$padj <0.05), ] # great way to subset data
dim(sig_popM) # check

sig_popM_df<- as.data.frame(sig_popM)
sig_popM_df$Row.names <- rownames(sig_popM_df)
dim(sig_popM_df)

genesOfInterest_popM <- c(sig_popM_df$Row.names)
length(genesOfInterest_popM)

# variance stabilization
vsd <- vst(dds, blind=F)

#####
# dds$combined = factor(paste0(dds$population, "-", dds$devstage, "-", dds$sex))
# dds$combined <- factor(dds$combined, levels=c("WA-L3L-M","WA-PP1-F","WA-PP1-M","WA-PD1-F","WA-PD1-M","WA-AD4-F","WA-AD4-M","NC-L3L-F","NC-L3L-M","NC-PP1-F","NC-PP1-M","NC-PD1-F","NC-PD1-M","NC-AD4-F","NC-AD4-M","WA-L3L-F","WA-L3L-M","WA-PP1-F","WA-PP1-M","WA-PD1-F","WA-PD1-M","WA-AD4-F","WA-AD4-M","NC-L3L-F","NC-L3L-M","NC-PP1-F","NC-PP1-M","NC-PD1-F","NC-PD1-M","NC-AD4-F","NC-AD4-M"))

# # calculating means 
# baseMeanPerGrp <- sapply( levels(dds$combined), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$combined == lvl] ) )
# 
# head(baseMeanPerGrp)
# dim(baseMeanPerGrp)

# pulls out normalized counts (avg of 3 reps) for all of our significant genes 
m <- baseMeanPerGrp[genesOfInterest_popM, c("WA-L3L-M","WA-PP1-M","WA-PD1-M","WA-AD4-M", "NC-L3L-M","NC-PP1-M","NC-PD1-M","NC-AD4-M")] 
head(m)
dim(m)

# use the apply function to scale the matrix m: apply(object, rows, function)
mat_scaled = t(apply(m, 1, scale))
head(mat_scaled)

# make a heat map!
library(pheatmap)

pheatmap(mat_scaled, labels_col=c("WA-L3L-M","WA-PP1-M","WA-PD1-M","WA-AD4-M", "NC-L3L-M","NC-PP1-M","NC-PD1-M","NC-AD4-M", cluster_cols=T, cluster_rows=T))

# ########## GO analysis
# input="dge_m_pop.csv" 
# goAnnotations="gene_annotation_only.tab"
# goDatabase="go.obo"
# goDivision="MF"
# source("gomwu.functions.R")
# 
# gomwuStats(input, goDatabase, goAnnotations, goDivision,
#            perlPath="perl",largest=0.1,smallest=5,clusterCutHeight=0.25,Alternative="g")
# 
# # Plotting results
# library(ape)
# quartz()
# gomwuPlot(input,goAnnotations,goDivision,
#           absValue=-log(0.05,10),level1=0.1,level2=0.05,level3=0.01,
#           txtsize=1.2, treeHeight=0.5)
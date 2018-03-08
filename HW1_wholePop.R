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

########## Model 1: Gene Expression Analysis between NC & WA populations############

## design call in DESeqDataSetFromMatrix function models "population effect" controlling for differences in devstage and sex
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design= ~devstage+sex+population)

## look at dimensions of dds
dim(dds)
## [1] 17483    48

## when you sum across all these columns, looking for at least one read
dds <- dds[rowSums(counts(dds)) >1,]

dim(dds)
## [1] 16851    48   okay, so only lost about 600 genes

## run DESeq function
dds <- DESeq(dds,modelMatrixType = "standard")

resultsNames(dds)
#[1] "Intercept"           "devstage_L3L_vs_AD4" "devstage_PD1_vs_AD4" "devstage_PP1_vs_AD4" "sex_M_vs_F"          "population_WA_vs_NC"
## Why not all comparisons between devstages?

res <- results(dds)

# order matrix by p adjusted values
res <- res[order(res$padj),]
head(res)

#look at which ones are significant
summary(res)

# Draw out WA vs NC specific comparison - res_pop
res_pop <- results(dds, name= "population_WA_vs_NC", alpha = 0.05)

res_pop <- res_pop[order(res_pop$padj),]
head(res_pop)
# baseMean log2FoldChange     lfcSE      stat       pvalue         padj
# <numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
#   OTAU008667-RA 231.87111     -0.7803753 0.1132661 -6.889755 5.588875e-12 4.129899e-08
# OTAU012562-RA 251.77436     -0.7828048 0.1135355 -6.894800 5.394059e-12 4.129899e-08
# OTAU011160-RA  10.24152     -1.0863185 0.1836927 -5.913781 3.343429e-09 1.438959e-05
# OTAU012716-RA 188.87768      1.0137238 0.1721500  5.888609 3.894604e-09 1.438959e-05
# OTAU002976-RA 998.42315     -0.6159655 0.1075972 -5.724735 1.035950e-08 3.062062e-05
# OTAU014686-RA 603.66091     -0.8168345 0.1468336 -5.562995 2.651835e-08 6.531911e-05
# > 

summary(res_pop)
# out of 16851 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 159, 0.94% 
# LFC < 0 (down)   : 147, 0.87% 
# outliers [1]     : 133, 0.79% 
# low counts [2]   : 1939, 12% 
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

##########################MA Plot ############## 
plotMA(res_pop, main="WA vs NC", ylim=c(-2,2))
abline(h=c(-1,1), col="blue",lwd=2)

####################### Save .csv outputs #####################################
#save outputs for population-wide comparison of differential gene expression
write.csv(res_pop,file="dgeWAvsNC_pop.csv", row.names=T, quote=F)
head(res_pop)

# change pvalue to -log(pvalue) and export as .csv with rownames
neglogpval <- as.matrix(-log(res_pop$pvalue))
head(neglogpval)

res_pop_negpval <- cbind(row.names(res_pop), neglogpval)
head(res_pop_negpval)

colnames(res_pop_negpval)=c("gene","neglogpval")

write.csv(res_pop_negpval,file="dge_WAvsNC_pop_negpval.csv",row.names=F, quote=F)

#########  Genes Of Signifcance with Heat Maps###########################################
####pull out significance genes to be able to make a heat map
sig_pop <- res_pop[which(res_pop$padj <0.05), ] # great way to subset data
dim(sig_pop) # check

sig_pop_df <- as.data.frame(sig_pop)
sig_pop_df$Row.names <- rownames(sig_pop_df)
dim(sig_pop_df)

genesOfInterest_pop <- c(sig_pop_df$Row.names)
length(genesOfInterest_pop)

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
m <- baseMeanPerGrp[genesOfInterest_pop, c("WA-L3L-F","WA-L3L-M","WA-PP1-F","WA-PP1-M","WA-PD1-F","WA-PD1-M","WA-AD4-F","WA-AD4-M","NC-L3L-F","NC-L3L-M","NC-PP1-F","NC-PP1-M","NC-PD1-F","NC-PD1-M","NC-AD4-F","NC-AD4-M", "WA-L3L-F","WA-L3L-M","WA-PP1-F","WA-PP1-M","WA-PD1-F","WA-PD1-M","WA-AD4-F","WA-AD4-M","NC-L3L-F","NC-L3L-M","NC-PP1-F","NC-PP1-M","NC-PD1-F","NC-PD1-M","NC-AD4-F","NC-AD4-M")]
head(m)
dim(m)

# use the apply function to scale the matrix m: apply(object, rows, function)
mat_scaled = t(apply(m, 1, scale))
head(mat_scaled)

# make a heat map!
library(pheatmap)
pheatmap(mat_scaled, labels_col=c("WA-L3L-F","WA-L3L-M","WA-PP1-F","WA-PP1-M","WA-PD1-F","WA-PD1-M","WA-AD4-F","WA-AD4-M","NC-L3L-F","NC-L3L-M","NC-PP1-F","NC-PP1-M","NC-PD1-F","NC-PD1-M","NC-AD4-F","NC-AD4-M", "WA-L3L-F","WA-L3L-M","WA-PP1-F","WA-PP1-M","WA-PD1-F","WA-PD1-M","WA-AD4-F","WA-AD4-M","NC-L3L-F","NC-L3L-M","NC-PP1-F","NC-PP1-M","NC-PD1-F","NC-PD1-M","NC-AD4-F","NC-AD4-M"), cluster_cols=T, cluster_rows=T)

################################## GO Analysis##################################################
input="dge_WAvsNC_pop.csv" 
goAnnotations="gene_annotation_only.tab"
goDatabase="go.obo"
goDivision="BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",largest=0.1,smallest=5,clusterCutHeight=0.25,Alternative="g")

# Plotting results
library(ape)
quartz()
gomwuPlot(input,goAnnotations,goDivision,
          absValue=-log(0.05,10),level1=0.1,level2=0.05,level3=0.01,
          txtsize=1.2, treeHeight=0.5)

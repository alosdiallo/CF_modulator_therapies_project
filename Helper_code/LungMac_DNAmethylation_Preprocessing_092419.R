# =========================================================================
# DNA methylation Pre-processing - Dave's Lung Macs
# Code by: SL, MM
# 09/24/2019
# =========================================================================
# Set working directory
setwd('/Users/sarahminkyunglee/OneDrive - Dartmouth College/Dave_LungMacs')

# Load required packages
library(minfi)
library(ENmix)
library(FlowSorted.Blood.EPIC)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(ggplot2)
library(ggpubr)

# Read in SNP annotation (Zhou et al, 2016)
ZhouAnnot <- read.delim("Zhou.EPIC.anno.GRCh38.tsv")

# Location of idat files
idat <- "/Users/sarahminkyunglee/OneDrive - Dartmouth College/Dave_LungMacs/idat_files"

# Read in methylation .idat files and covariate data using minfi package
targets <- read.metharray.sheet(idat)
RGset <- read.metharray.exp(targets=targets, extended= TRUE) #RedGreen channel set

SampleNames <- pData(RGset)$SampleName
Group <- pData(RGset)$Beadchip

# Load covariate data
covariates <- read.csv("/Users/sarahminkyunglee/OneDrive - Dartmouth College/Dave_LungMacs/idat_files/SampleSheet_10012019.csv", header= TRUE)

# Recast relevant variable as factors
covariates$Beadchip <- as.factor(covariates$Sentrix_ID)
covariates$CompleteBarcode <- as.factor(covariates$Filename)

# =========================================================================
# Quality Control
# =========================================================================
# Set working directory for QC plots
setwd('/Users/sarahminkyunglee/OneDrive - Dartmouth College/Dave_LungMacs/Methylation_QualityControl')

# Get QC
ENmix::plotCtrl(RGset)

# Increase sample threshold to 10%
qcinfo <- ENmix::QCinfo(RGset, samplethre= 0.1) #0 samples 

# Display lower quality samples
covariates[covariates$CompleteBarcode %in% qcinfo$badsample,c('SampleID', 'Beadchip', 'Sentrix_Position', 'Age', 'Gender')]

# No low quality samples removed as they are all above the sample threshold

# =========================================================================
# Estimate Cell Counts
# =========================================================================
estcellcounts <- estimateCellCounts2(RGset, processMethod = "preprocessFunnorm",
                                     referencePlatform = "IlluminaHumanMethylationEPIC", 
                                     IDOLOptimizedCpGs = IDOLOptimizedCpGs,)
estcellcounts$counts

library(reshape2)
celltypes <- melt(estcellcounts$counts)

celltypes_plot <- ggplot(celltypes, aes(x= Var2, y= value)) + geom_boxplot() +
  labs(x= "Cell types", y= "Cell type proportion") + 
  theme(axis.title = element_text(size= 23), axis.text = element_text(size= 20), 
        panel.background = element_rect(color= "#FFFFFF"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size= 0.6, fill= NA))
celltypes_plot

ggsave("DNAm_QC/EstimateCellCounts.png", plot = celltypes_plot, height= 8, width = 8)

# =========================================================================
# Quality Control
# =========================================================================
# Set correct working directory
setwd('/Users/sarahminkyunglee/OneDrive - Dartmouth College/Dave_LungMacs')

# Run functional normalization
funnormEPIC <- preprocessFunnorm(RGset)

# Filter probes below the specified detection p-value
funnormEPIC <- funnormEPIC[!row.names(funnormEPIC) %in% qcinfo$badCpG,]

# Indicate number of probes removed - 9432 CpGs
paste(formatC(length(qcinfo$badCpG), big.mark= ','), 'CpGs removed due to low detection p value')

# Filter probes on Y chromosome
cpgRemove <- ZhouAnnot$probeID[ZhouAnnot$chrm == "chrY"]
cpgRemove <- cpgRemove[cpgRemove %in% row.names(funnormEPIC)]
funnormEPIC <- funnormEPIC[!row.names(funnormEPIC) %in% cpgRemove]

# Indicate number of probes removed - 153
paste(formatC(length(cpgRemove), big.mark= ','), 'loci on Y chromosome removed')

# Filter probes on X chromosome
cpgRemove <- ZhouAnnot$probeID[ZhouAnnot$chrm == "chrX"]
cpgRemove <- cpgRemove[cpgRemove %in% row.names(funnormEPIC)]
funnormEPIC <- funnormEPIC[!row.names(funnormEPIC) %in% cpgRemove]

# Indicate number of probes removed - 18,911
paste(formatC(length(cpgRemove), big.mark= ','), 'loci on X chromosome removed')

# Filter SNP/cross-hybridizing probes
cpgRemove <- ZhouAnnot$probeID[ZhouAnnot$MASK.general == T]
cpgRemove <- cpgRemove[cpgRemove %in% row.names(funnormEPIC)]
funnormEPIC <- funnormEPIC[!row.names(funnormEPIC) %in% cpgRemove]

# Indicate number of probes removed - 78,379
paste(formatC(length(cpgRemove), big.mark= ','), 'SNP associated or cross-hybridizing CpGs removed')

# Save initial results from functional normalization
save(list= c("funnormEPIC", "covariates"), file= '/Users/sarahminkyunglee/OneDrive - Dartmouth College/Dave_LungMacs/DNAm_ProcessedData/LungMac_FunNorm_EPIC.RData')

# =========================================================================
# Beta values
# =========================================================================
# Extract beta-values for full dataset
betas <- getBeta(funnormEPIC)

# Save funnorm Beta values
save(list= c("betas", "covariates"), file= '/Users/sarahminkyunglee/OneDrive - Dartmouth College/Dave_LungMacs/DNAm_ProcessedData/LungMac_FunNorm_Betas.Rdata')

# =========================================================================
# PCA
# =========================================================================
# Look at PCA for overall betas
pca_overall <- princomp(betas)

# Merge first 5 PCs with covariate data
pca_overall <- cbind(pca_overall$loadings[,1:5],covariates[covariates$CompleteBarcode %in% colnames(betas),])
colnames(pca_overall)[1:5] <- paste('PC',seq(1:5),sep = '')

# Generate output plots
pdf('DNAm_QC/NormalizedBetas_PCA.pdf')

# Plot by array chip
ggscatter(data = pca_overall,x = 'PC1',y = 'PC2',color = 'Beadchip')
ggscatter(data = pca_overall,x = 'PC2',y = 'PC3',color = 'Beadchip')

#  Plot by subject
ggscatter(data = pca_overall,x = 'PC1',y = 'PC2',color = 'Subject')
ggscatter(data = pca_overall,x = 'PC2',y = 'PC3',color = 'Subject')

# Plot by gender
ggscatter(data = pca_overall,x = 'PC1',y = 'PC2',color = 'Gender')
ggscatter(data = pca_overall,x = 'PC2',y = 'PC3',color = 'Gender')

# Plot by treatment
ggscatter(data = pca_overall,x = 'PC1',y = 'PC2',color = 'Treatment')
ggscatter(data = pca_overall,x = 'PC2',y = 'PC3',color = 'Treatment')
dev.off()

# =========================================================================
# DNA methylation Heatmap - Dave's Lung Macs
# Code by: SL, MM, DC
# 09/24/2019
# =========================================================================
# Set working directory
setwd('/Users/sarahminkyunglee/OneDrive - Dartmouth College/Dave_LungMacs')

# Load required packages
library(pheatmap)
library(RColorBrewer)
require(grid)
fixInNamespace("draw_colnames", "pheatmap") #vjust = 1, hjust = 0.5, rot = 0

# Load beta values
load('/Users/sarahminkyunglee/OneDrive - Dartmouth College/Dave_LungMacs/DNAm_ProcessedData/LungMac_FunNorm_Betas.Rdata')

# Clean covariate data for annotation
row.names(covariates) <- covariates$CompleteBarcode

# Change colnames of betas 
betas2 <- betas
colnames(betas2) <- covariates$SampleName

covariates2 <- covariates
rownames(covariates2) <- covariates2$SampleName


#Make the annotation dataframe for the heatmap
heat_annot <- data.frame(
  row.names= covariates2$SampleName,
  Gender= covariates2$Gender,
  Treatment= covariates2$Treatment)


#Define colors for annotation colorbars
ann_colors <- list(Gender= c('F'= "grey", 
                             'M'= "black"), 
                   Treatment= c('Control'= "grey", 
                                'Treatment'= "black"))

# Identify most variable CpGs
# Ensure beta values and covariate data are in the same order
all(colnames(betas2) == row.names(covariates2))

# Calculate variable of each CpG across all samples(all subjects, all tissue types)
library(matrixStats)
CpG_Var <- matrixStats::rowVars(betas2)
rankVar <- data.frame('order' = rank(-CpG_Var), 'var'= CpG_Var)
rankVar$top10000 <- ifelse(rankVar$order <= 10000, 'Yes', 'No')
ggpubr::ggscatter(rankVar, x= 'order', y= 'var', color= 'top10000', palette= c('black', 'red') )

# Get rank of CpGs by
rankvar= rank(-CpG_Var)

# top 10,000 cpgs
data.topvar <- betas2[rankVar$top10000 == 'Yes',]
data.topvar <- as.matrix(data.topvar)

# =========================================================================
# Plot
# =========================================================================
#Generate heatmap
pheatmap(
  data.topvar,
  annotation_names_col= T,
  show_rownames= FALSE, #CpGs
  show_colnames= TRUE, #Samples
  annotation_col= heat_annot, 
  #annotation_row= row_annot,
  annotation_colors= ann_colors,
  cluster_cols= FALSE,
  color= colorRampPalette(c("yellow", "blue"))(1024),
  clustering_distance_rows= "manhattan",
  clustering_distance_colors= "manhattan",
  clustering_method= "average",
  border_color= NA,
  fontsize= 13
)

pheatmap(
  data.topvar,
  annotation_names_col= T,
  show_rownames= FALSE, #CpGs
  show_colnames= TRUE, #Samples
  annotation_col= heat_annot, 
  #annotation_row= row_annot,
  annotation_colors= ann_colors,
  cluster_cols= FALSE,
  color= colorRampPalette(c("yellow", "blue"))(1024),
  clustering_distance_rows= "manhattan",
  clustering_distance_colors= "manhattan",
  clustering_method= "average",
  file= "DNAm_Figures/Top10000CpGs_Heatmap.png",
  border_color= NA,
  fontsize= 13
)


# =========================================================================
# Supervised clustering
# =========================================================================
library(RPMM)
library(doParallel)

OUT_FILE_NAME <- "~/Dave_LungMacs/DNAm_Files/RPMM_cluster_membership_final.csv"; 

load_Yinv <- function(n.CpG=10000) {
  #'@description Load, subset, and transform data for RPMM
  load('/Users/sarahminkyunglee/OneDrive - Dartmouth College/Dave_LungMacs/DNAm_ProcessedData/LungMac_FunNorm_Betas.Rdata');
  covariates$CompleteBarcode <- paste(covariates$Sentrix_ID, covariates$Sentrix_Position,sep="_");
  if(all(covariates$CompleteBarcode == colnames(betas))){ #checkpoint
    print("Samples already matched, so proceed to name switching...");
    colnames(betas) <- covariates$SampleID; 
  } else {
    print("Matching sample names...");
    betas <- betas[ , match(covariates$CompleteBarcode, colnames(betas))]; 
    print("Now ready to name switching...");
    colnames(betas) <- covariates$SampleID; 
  }
  sele <- order(rowVars(betas), decreasing=TRUE)[1:n.CpG]; 
  Y_inv <- t(betas[sele, ]); 
  assign("Y_inv", Y_inv, envir=.GlobalEnv);
}

getRPMMClustLabels <- function(rpmmObject, Y_inv=NULL) {
  #'@description Extracts RPMM hard cluster labels
  #'@param rpmmObject RPMM object
  #'@param Y_inv Optional. Input matrix for RPMM computation. If provided, the sample name will be updated
  
  hardLabels <- blcTreeLeafClasses(rpmmObject);
  hardLabels <- as.data.frame(hardLabels);
  colnames(hardLabels) <- "RPMMClusters";
  if(! is.null(Y_inv)) rownames(hardLabels) <- hardLabels$Sample_Name <- rownames(Y_inv);
  return(hardLabels);
}

getRPMMSampOrder <- function(rpmmClusters, Y_inv) {
  #'@description Retrieves sample orders fo heat map visualization
  #'@describeIn Hinoue et al. online tutorial
  #'@param rpmmClusters data.frame with row.names = sample ID & 1 column named RPMM with cluster assignments
  #'@param Y_inv Input matrix for RPMM computation
  
  sampOrder <- c();
  for(r in names(table(rpmmClusters$RPMM))) {
    samps <- rownames(rpmmClusters)[rpmmClusters$RPMM == r];
    clu <- t(Y_inv[samps, ]);
    s_i <- seriation::seriate(clu, margin=2);
    so_i <- seriation::get_order(s_i);
    sampOrder <- c(sampOrder, samps[so_i]);
  }
  sampOrder <- data.frame(
    Sample_Name = sampOrder,
    RPMMSampleOrder = 1:length(sampOrder)
  );
  return(sampOrder);
}

main <- function() {
  ## RPMM with max_level of 2:
  load_Yinv();
  sup_rpmm <- blcTree(Y_inv, verbose=1, maxlevel=2-1); 
  print(sup_rpmm);
  
  ## Extract RPMM cluster labels:
  rpmmClusters <- getRPMMClustLabels(sup_rpmm, Y_inv);
  
  ## Retrieve RPMM sample order:
  sampOrders <- getRPMMSampOrder(rpmmClusters, Y_inv);  
  
  ## Export:
  dat_sup <- merge(rpmmClusters, sampOrders, by="Sample_Name");
  write.csv(dat_sup, file= "DNAm_Files/RPMM_clustering_results.csv")
}

main(); 

# Add in RPMM clustering information 
RPMM <- read.csv("DNAm_Files/RPMM_clustering_results.csv", header=TRUE)
RPMM <- RPMM[,-1]

covariates <- merge(covariates, RPMM, by.x= "SampleID", by.y= "Sample_Name")
covariates <- covariates[, -13]

# Change colnames of betas 
betas2 <- betas
colnames(betas2) <- covariates$SampleName
  
covariates2 <- covariates
rownames(covariates2) <- covariates2$SampleName

#Make the annotation dataframe for the heatmap
heat_annot <- data.frame(
  row.names= covariates2$SampleName,
  Gender= covariates2$Gender,
  Treatment= covariates2$Treatment)

#Define colors for annotation colorbars
ann_colors <- list(Gender= c('F'= "grey", 
                             'M'= "black"), 
                   Treatment= c('Control'= "grey", 
                                'Treatment'= "black"))

# Identify most variable CpGs
# Ensure beta values and covariate data are in the same order
all(colnames(betas2) == row.names(covariates2))

# Calculate variable of each CpG across all samples(all subjects, all tissue types)
library(matrixStats)
CpG_Var <- matrixStats::rowVars(betas2)
rankVar <- data.frame('order' = rank(-CpG_Var), 'var'= CpG_Var)
rankVar$top10000 <- ifelse(rankVar$order <= 10000, 'Yes', 'No')
ggpubr::ggscatter(rankVar, x= 'order', y= 'var', color= 'top10000', palette= c('black', 'red') )

# Get rank of CpGs by
rankvar= rank(-CpG_Var)

# top 10,000 cpgs
data.topvar <- betas2[rankVar$top10000 == 'Yes',]
data.topvar <- as.matrix(data.topvar)

sampOrders <- covariates2$CompleteBarcode[order(covariates2$Treatment, decreasing=TRUE)];

# =========================================================================
# Plot
# =========================================================================
#Generate heatmap
pheatmap(
  data.topvar[ , sampOrders],
  show_rownames = FALSE, #CpGs
  show_colnames = TRUE, #samples
  cluster_cols = FALSE, 
  annotation_col = heat_annot,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow", "blue"))(1024),
  border_color = NA, 
  file= "DNAm_Figures/Treatment_Clustering_Heatmap.png",
  fontsize = 13
)
dev.off()

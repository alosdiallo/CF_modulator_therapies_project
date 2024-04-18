# =========================================================================
# DNA methylation EWAS - Dave's Lung Macs
# Code by: SL, MM
# 09/25/2019
# =========================================================================
# Set working directory
setwd('/Users/sarahminkyunglee/OneDrive - Dartmouth College/Dave_LungMacs')

# Load packages
library(limma)
library(qvalue)
library(ggpubr)
library(minfi)

# Load beta values
load('DNAm_ProcessedData/Subset_betas_var.001.rdata')

# Add rownames to covariates
rownames(covariates) <- covariates$CompleteBarcode

# Make sample ID a factor
covariates$SampleID <- as.factor(covariates$SampleID)
covariates$CompleteBarcode <- as.character(covariates$CompleteBarcode)

# Convert to M values
M.Values <- logit2(betas_subset)

# =========================================================================
# Calculate correlation structure - unadjusted models
# =========================================================================
# Estimate correlation for control vs treatment
# Define model matrix
XX <- model.matrix(~Treatment, data= covariates)

# Calculate correlation structure - takes a little bit of time
corfit <- duplicateCorrelation(M.Values, XX, block= covariates$Subject)

# Look at consensus correlation
corfit$consensus.correlation #0.3891053

# =========================================================================
# Calculate correlation structure - adjusted models
# =========================================================================
# Estimate correlation for control vs treatment
#Define model matrix
XX <- model.matrix(~Treatment + Gender + Age, data= covariates)

#Calculate correlation structure - takes a little bit of time
adj_corfit <- duplicateCorrelation(M.Values, XX, block= covariates$Subject)

#Look at consensus correlation
adj_corfit$consensus.correlation #0.8569234

# =========================================================================
# Model specifications
# =========================================================================
# Run limma model
results_limma <- eBayes(lmFit(M.Values,XX,correlation = adj_corfit$consensus.correlation))

# Tumor Comparison
q.values <- qvalue(results_limma$p.value[, "TreatmentTreatment"])
summary(q.values)

# Generate results table
results_table <- data.frame('beta' = results_limma$coefficients[, "TreatmentTreatment"],
                                     'pVal' = q.values$pvalues,
                                     'log10.pVal' = -log10(q.values$pvalues),
                                     'qVal' = q.values$qvalues)

# Save 
save(results_limma, results_table, file= "DNAm_Files/Subset_EWAS_Results.Rdata")

# =========================================================================
# Volcano plots
# =========================================================================
# Volcano plots colored by control vs treatment
# Generate list of significant sites of interest
cpgList <- row.names(results_table)[results_table$qVal < 0.01]

# Save list of CpGs qVal < 0.01
write.csv(cpgList, "DNAm_Files/Subset_DifferentiallyMethylatedCpGs_qval.01.csv")

# Generate list of significant sites of interest - qval < 0.05
cpgList_5 <- row.names(results_table)[results_table$qVal < 0.05]

# Save list of CpGs qVal < 0.01
write.csv(cpgList_5, "DNAm_Files/Subset_DifferentiallyMethylatedCpGs_qval.05.csv")

# Create a new annotation column for sites significant (q < 0.01) in tumor vs cont normal comparison
results_table$label <- factor(ifelse(row.names(results_table) %in% cpgList_5, 1, 0))

# Reorder to plot sig points on top
results_table <- results_table[order(results_table$label,decreasing = F),]

# Plot
# Establish Q-value cutoffs
qcutHigh <- min(results_table$log10.pVal[results_table$qVal < 0.01])
qcutLow <- min(results_table$log10.pVal[results_table$qVal < 0.05])

# Calculate total num CpGs in each range for annotation
cpgHigh <- nrow(results_table[results_table$qVal < 0.01, ])
cpgLow <- nrow(results_table[results_table$qVal >= 0.01 & results_table$qVal < 0.05,])
cpgUnder <- nrow(results_table[results_table$qVal >= 0.05,])

# Check all CpGs included
cpgHigh + cpgLow + cpgUnder == nrow(results_table)

# Plot results
png(filename = 'DNAm_Figures/Subset_VolcanoPlots_DifferentialMethylation.png',width = 800,height = 800)
p <- ggscatter(results_table,x = 'beta',y = 'log10.pVal',color = 'label',
               alpha = 0.75,size = 0.75,title = 'Control vs Treatment\n',
               xlab = 'Beta Coefficient',ylab = expression('-log'[10]*'p value'),
               xlim = c(-4,4),ylim = c(0,12))
plot <- ggpar(p,legend = 'none',palette = c('black','red3'),
                      font.x = 30,font.y = 30,font.main = 34,font.tickslab = 30) +
  # Add lines to indicate Q value cutoffs
  geom_hline(yintercept = qcutHigh, color = "red", linetype = "dashed") + 
  geom_hline(yintercept = qcutLow, color = "blue", linetype = "dashed") + 
  # Add annotations of number of CpGs in each bin
  annotate(geom = 'text',x = -3.35,y = qcutHigh + 0.5,label = paste(format(cpgHigh,big.mark = ','),'CpGs'),size = 8) +
  annotate(geom = 'text',x = -3.35,y = qcutLow + 0.5,label = paste(format(cpgLow,big.mark = ','),'CpGs'),size = 8) +
  annotate(geom = 'text',x = -3.35,y = qcutLow - 0.5,label = paste(format(cpgUnder,big.mark = ','),'CpGs'),size =8)
plot
dev.off()

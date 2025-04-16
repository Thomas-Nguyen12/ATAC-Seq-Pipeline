library(DiffBind) 
library(tidyverse) 
args = commandArgs(trailingOnly=TRUE)
# I can specify which argument via indexing
# for example, sample_sheet <- args[1], output_dir <- args[2]
out_dir <- args[2]
sample_sheet <- read_csv(args[1])


# setting the working directory to avoid setting too many directories
setwd(out_dir)
print ("---------------")
print ("Files in this Directory...")
list.files()
print ("---------------")


sample_sheet <- sample_sheet %>% mutate(
  bamReads=paste(out_dir, "/", SampleID, "_without_blacklist.bam", sep=""),
  Peaks=paste(out_dir, "/", SampleID, "_peak_peaks.narrowPeak", sep=""),
  PeakCallers='narrow'
)

placenta <- dba(sampleSheet=sample_sheet)
print ("Counting reads...")
# counting reads
placenta_counts <- dba.count(placenta)
pdf(file = paste(out_dir, "unnormalised_counts.pdf", sep='/'), width = 8, height = 11) # defaults to 7 x 7 inches
plot(placenta_counts)
dev.off()
saveRDS(placenta_counts, paste(out_dir, 'raw_sample_counts.rds', sep='/'))

print ("Normalising...")
placenta_norm <- dba.normalize(placenta_counts, method=DBA_NORM_LIB)
pdf(file = paste0(out_dir, paste(out_dir, "normalised_counts.pdf", sep='/'), width = 8, height = 11) # defaults to 7 x 7 inches
plot(placenta_norm)
dev.off()
saveRDS(placenta_norm, paste(out_dir, 'placenta_norm.rds', sep='/'))


norm <- dba.normalize(placenta_norm, bRetrieve=TRUE)

print ("Finished normalising...")
## checking norm factors 

print ("normalisation method...")
print (norm)
saveRDS(norm, paste(out_dir, "norm.rds", sep='/'))
normlibs <- cbind(FullLibSize=norm$lib.sizes, NormFacs=norm$norm.factors, NormLibSize=round(norm$lib.sizes/norm$norm.factors)) 
saveRDS(normlibs, paste(out_dir, 'normlibs.rds', sep='/'))

print ("Establishing model Design and contrast...")
placenta_contrast <- dba.contrast(placenta_norm, categories=DBA_CONDITION, reorderMeta=list(Condition="Stem"))
saveRDS(placenta_contrast, paste(out_dir, "placenta_contrast.rds", sep='/')
        
print ("Initiating the Differential Analysis...")
# Performing the differential analysis
## This runs the default DESeq2 analysis	
placenta_analyse <- dba.analyze(placenta_contrast)

print ("Saving the analysis output...")
saveRDS(placenta_analyse, paste(out_dir, "placenta_analyse.rds", sep='/'))
dba.show(placenta_analyse, bContrasts=TRUE) 

pdf(file = paste(out_dir, "sample_analyse.pdf", sep='/'), width = 8, height = 11) # defaults to 7 x 7 inches
plot(placenta_analyse, contrast=1)
dev.off()


pdf(file = paste(out_dir, "correlation.pdf", sep='/'), width = 8, height = 11) # defaults to 7 x 7 inches
dba.plotHeatmap(placenta_analyse, contrast=1, correlations=TRUE, attributes=DBA_CONDITION)
dev.off()



placenta_report <- dba.report(placenta_analyse)
saveRDS(placenta_report, paste(out_dir, "samplereport.rds", sep="/"))

## CREATING THE PLOTS 

## PCA plots: plot_parameters --> (placenta_analyse, plot title, categories)
pdf(file = paste(out_dir, "sample_analyze_PCA_diagram.pdf", sep='/'), width = 8, height = 11)
dba.plotPCA(placenta_analyse, DBA_NAME, label=DBA_CONDITION, method=DBA_DESEQ2)
def.off()

pdf(file = paste(out_dir, "sample_analyze_MA_plot.pdf", sep='/'), width = 8, height = 11)
dba.plotMA(placenta_analyse, method=DBA_DESEQ2)
dev.off()

pdf(file = paste(out_dir, "sample_analyze_volcano_diagram.pdf", sep='/'), width = 8, height = 11)
dba.plotVolcano(placenta_analyse, method=DBA_DESEQ2)
dev.off()

pdf(file = paste(out_dir, "sample_analyze_boxplot.pdf", sep='/'), width = 8, height = 11)
dba.plotBox(placenta_analyse, method=DBA_DESEQ2)
dev.off()

# profile plot
profiles <- dba.plotProfile(placenta_analyse, merge=c(DBA_TISSUE, DBA_REPLICATE, DBA_CONDITION))

pdf(file=paste(out_dir, 'profile_plot.pdf', sep='/'), width=8, height=11)
dba.plotProfile(profiles, DBA_DESEQ2)
dev.off()



hmap <- colorRampPalette(c("red", 'black', 'green'))
pdf(file = paste(out_dir, "sample_analyze_heatmap.pdf", sep='/'), width = 8, height = 11)
dba.plotHeatmap(placenta_analyse, contrast=1, correlations=FALSE, scale="row", colScheme=hmap, method=DBA_DESEQ2) 
dev.off()



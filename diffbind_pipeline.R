library(DiffBind) 
library(tidyverse) 
args = commandArgs(trailingOnly=TRUE)
# I can specify which argument via indexing
# for example, sample_sheet <- args[1], output_dir <- args[2]
out_dir <- args[2]
sample_sheet <- read_csv(args[1])
print ("Output Directory: ")
print (out_dir)

# To avoid the whole directory issue, I can set the working directory to the out_dir so all plots are automatically saved there


# setting the working directory to avoid setting too many directories
setwd(out_dir)
print ("---------------")
print ("Files in this Directory...")
list.files()
print ("---------------")
print ("Loading the sample sheet...")


sample_sheet <- sample_sheet %>% mutate(
  bamReads=paste(out_dir, "/", SampleID, "_without_blacklist.bam", sep=""),
  Peaks=paste(out_dir, "/", SampleID, "_peak_peaks.narrowPeak", sep=""),
  PeakCaller='narrow'
)



print ("Sample Sheet...")
print (sample_sheet)
print ("------------------")
print (" ")
print ("Passing Sample Sheet Through Diffbind...")

placenta <- dba(sampleSheet=sample_sheet)
print ("Counting reads...")
# counting reads
placenta_counts <- dba.count(placenta)


pdf(file = "unnormalised_counts.pdf") # defaults to 7 x 7 inches
plot(placenta_counts)
dev.off()


print ("Normalising...")
placenta_norm <- dba.normalize(placenta_counts)
norm <- dba.normalize(placenta_norm, bRetrieve=TRUE)
print ("Normalisation methods...")
print (norm)


pdf(file = "normalised_counts.pdf") # defaults to 7 x 7 inches
plot(placenta_norm)
dev.off()

saveRDS(placenta_norm, 'placenta_norm.rds')


norm <- dba.normalize(placenta_norm, bRetrieve=TRUE)

print ("Finished normalising...")
## checking norm factors 

print ("normalisation method...")
print (norm)

normlibs <- cbind(FullLibSize=norm$lib.sizes, NormFacs=norm$norm.factors, NormLibSize=round(norm$lib.sizes/norm$norm.factors)) 
saveRDS(normlibs, 'normlibs.rds')

print ("Establishing model Design and contrast...")
placenta_contrast <- dba.contrast(placenta_norm, categories=DBA_CONDITION, reorderMeta=list(Condition="Stem"))
saveRDS(placenta_contrast, "placenta_contrast.rds")
        
print ("Initiating the Differential Analysis...")
# Performing the differential analysis
## This runs the default DESeq2 analysis	
placenta_analyse <- dba.analyze(placenta_contrast)

print ("Saving the analysis output...")
saveRDS(placenta_analyse, "placenta_analyse.rds")
dba.show(placenta_analyse, bContrasts=TRUE) 

pdf(file = "sample_analyse.pdf") # defaults to 7 x 7 inches
plot(placenta_analyse, contrast=1)
dev.off()


pdf(file = "correlation.pdf") # defaults to 7 x 7 inches
dba.plotHeatmap(placenta_analyse, contrast=1, correlations=TRUE, attributes=DBA_CONDITION)
dev.off()



placenta_report <- dba.report(placenta_analyse)
saveRDS(placenta_report, "diffbind_report.rds")

## CREATING THE PLOTS 

## PCA plots: plot_parameters --> (placenta_analyse, plot title, categories)
pdf(file = "sample_analyze_PCA_diagram.pdf")
dba.plotPCA(placenta_analyse, label=DBA_CONDITION, method=DBA_DESEQ2)
dev.off()

pdf(file = "sample_analyze_MA_plot.pdf")
dba.plotMA(placenta_analyse, method=DBA_DESEQ2)
dev.off()

pdf(file = "sample_analyze_volcano_diagram.pdf")
dba.plotVolcano(placenta_analyse, method=DBA_DESEQ2)
dev.off()

pdf(file = "sample_analyze_boxplot.pdf")
dba.plotBox(placenta_analyse, method=DBA_DESEQ2)
dev.off()

# profile plot
profiles <- dba.plotProfile(placenta_analyse, merge=c(DBA_TISSUE, DBA_REPLICATE, DBA_CONDITION))

pdf(file='profile_plot.pdf')
dba.plotProfile(profiles, doPlot=TRUE)
dev.off()



pdf(file = "sample_analyze_heatmap.pdf")
dba.plotHeatmap(placenta_analyse, correlations=FALSE, scale="row",
                contrast=1, attributes=DBA_CONDITION, method=DBA_DESEQ2,
                ColAttributes=DBA_CONDITION)

dev.off()



library(DiffBind) 
library(tidyverse) 





# reading in the samples 
placenta <- dba(sampleSheet='~/diff_sample_sheet.csv')
print ("Counting reads...")
# counting reads
placenta_counts <- dba.count(placenta)
saveRDS(placenta_counts, "~/scratch/diffbind/placenta_counts.rds")

pdf(file = "~/scratch/diffbind/unnormalised_counts.pdf", width = 8, height = 11) # defaults to 7 x 7 inches
plot(placenta_counts)
dev.off()





print ("Normalising...")
placenta_norm <- dba.normalize(placenta_counts, method=DBA_DESEQ2)

pdf(file = "~/scratch/diffbind/normalised_counts.pdf", width = 8, height = 11) # defaults to 7 x 7 inches
plot(placenta_norm)
dev.off()



saveRDS(placenta_norm, '~/scratch/diffbind/placenta_norm.rds')


norm <- dba.normalize(placenta_norm, bRetrieve=TRUE)

print ("Finished normalising...")
## checking norm factors 

print ("normalisation method...")
print (norm)
saveRDS(norm, "~/scratch/diffbind/norm.rds")
normlibs <- cbind(FullLibSize=norm$lib.sizes, NormFacs=norm$norm.factors, NormLibSize=round(norm$lib.sizes/norm$norm.factors)) 
saveRDS(normlibs, '~/scratch/diffbind/normlibs.rds')
# establishing a model design and contrast
## you should be comparing stem vs EVT

placenta_contrast <- dba.contrast(placenta_norm, categories=DBA_CONDITION, reorderMeta=list(Condition="Stem"))
saveRDS(placenta_contrast, "~/scratch/diffbind/placenta_contrast.rds")
print ("Analysing...")
# Performing the differential analysis
## This runs the default DESeq2 analysis	
placenta_analyse <- dba.analyze(placenta_contrast)
saveRDS(placenta_analyse, "~/scratch/diffbind/placenta_analyse.rds")
dba.show(placenta_analyse, bContrasts=TRUE) 

## Plotting the results

pdf(file = "~/scratch/diffbind/placenta_analyse.pdf", width = 8, height = 11) # defaults to 7 x 7 inches
plot(placenta_analyse, contrast=1)
dev.off()


# Retrieving the differentially bound sites
placenta_report <- dba.report(placenta_analyse)
saveRDS(placenta_report, "~/scratch/diffbind/placenta_report.rds")

# Plots can be generated afterwards 

## venn Diagram
pdf(file = "~/scratch/diffbind/placenta_analyze_venn_diagram.pdf", width = 8, height = 11)
dba.plotVenn(placenta_analyse, contrast=1, bDB=TRUE, )
def.off()
## PCA plots: plot_parameters --> (placenta_analyse, plot title, categories)
pdf(file = "~/scratch/diffbind/placenta_analyze_PCA_diagram_db.pdf", width = 8, height = 11)
dba.plotPCA(placenta_analyse, label=DBA_CONDITION, method=DBA_DESEQ2,
            contrast=1)
dev.off()

pdf(file = "~/scratch/diffbind/placenta_analyze_PCA_diagram.pdf", width = 8, height = 11)
dba.plotPCA(placenta_analyse, label=DBA_CONDITION)
dev.off()






## MA plot
pdf(file = "~/scratch/diffbind/placenta_analyze_MA_plot.pdf", width = 8, height = 11)
dba.plotMA(placenta_analyse, method=DBA_DESEQ2)
dev.off()
# Volcano plots
pdf(file = "~/scratch/diffbind/placenta_analyze_volcano_diagram.pdf", width = 8, height = 11)
dba.plotVolcano(placenta_analyse, method=DBA_DESEQ2, bLabels=TRUE
                )
dev.off()
# Boxplot
pdf(file = "~/scratch/diffbind/placenta_analyze_boxplot.pdf", width = 8, height = 11)
dba.plotBox(placenta_analyse, attribute=DBA_CONDITION)
dev.off()
# Heatmaps 

hmap <- colorRampPalette(c("red", 'black', 'green'))
pdf(file = "~/scratch/diffbind/placenta_analyze_heatmap.pdf", width = 8, height = 11)
dba.plotHeatmap(placenta_analyse, correlations=FALSE, scale="row",
                contrast=1, attributes=DBA_CONDITION, method=DBA_DESEQ2,
                ColAttributes=DBA_CONDITION)
dev.off()




pdf(file='/gpfs/home/wmp21vtu/scratch/diffbind/placenta_analyse_venn.pdf')
dba.plotVenn(placenta_analyse$masks$Stem, placenta_analyse$masks$EVT, contrast=1, method=DBA_DESEQ2, labelAttributes=DBA_CONDITION)
dev.off()





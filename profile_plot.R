library(DiffBind)


placenta_analyse <- readRDS('/gpfs/home/wmp21vtu/scratch/diffbind/placenta_analyse.rds')
profiles <- dba.plotProfile(placenta_analyze)


placenta_analyse$config$RunParallel <- TRUE
# EVT vs Stem
mask.MCF7 <- placenta_analyse$masks$MCF7
mask.stem <- placenta_analyse$masks$Stem
mask.evt <- placenta_analyse$masks$EVT

pdf(file='~/scratch/diffbind/profile_plot.pdf', width=8, height=10) 
profiles <- dba.plotProfile(profiles)
dev.off()

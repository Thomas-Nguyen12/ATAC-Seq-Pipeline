# ATAC-Seq-Pipeline Overview
This project aims to create a pipeline for ATAC-Seq Differential Chromatin Accessibility analysis as part of my MScR in Biomolecular Science (Bioinformatics). This project uses both UNIX and R scripting languages depending on the tools/packages required. Hopefully, this tool will be useful to those needing Bioinformatic analysis of their ATAC-Seq experiment.

# Introduction

ATAC-Seq is a method to identify regions of open chromatin (DNA + histone molecules) across the genome. ATAC-Seq has been involved in studies investigating transcriptional regulators within autism (Fazel Darbandi et al., 2024) and Enhancer elements in the Placenta (Abdulghani, Jain and Tuteja, 2019), and as such, has major applications towards understanding development and disease 

![image](https://github.com/user-attachments/assets/6b1cd157-c989-42ed-996e-b6d9a320b5a6)
_Figure 1: ATAC-Seq schematic. Tn5 Tranposases cleave open DNA and simultaneously attach barcode sequencing primers in a process known as tagmentation. DNA libraries are then sequenced using next-generation sequencing technologies such as Illumina and subject to Bioinformatic Analysis (Dillinger, 2021)_

Compared to alternatives such as DNAse-Seq (Song and Crawford, 2010) and Faire-Seq (Giresi et al., 2007), ATAC-Seq requires fewer input samples and is much more time efficient due to a hyperactive Tn5 Transposase, requiring only 500-50,000 cells (Buenrostro et al., 2013).

To date, however, there are few tools specifically designed for ATAC-Seq analysis (Angarica and Del Sol, 2017), with many CHiP-Seq tools being used in place, assuming similar data properties (Chang et al., 2018). 

Several pipelines incorporate various tools and libraries for end-to-end analysis, such as PEPATAC (Smith et al., 2021) and CATCH-UP (Riva et al., 2023) that have varying focuses on serial alignments or bulk samples. However, none of them account for cell-type differences in chromatin accessibility. Consequently, I have developed a pipeline for paired-end sequencing data that creates ATAC peaksets and assesses the chromatin accessibility between two cell-types. Using this pipeline, users can investigate the functional implications of where DNA is differentially open or closed.

This workflow has seven steps:
1. Importing data
2. quality control using htstream v1.4.1
3. alignment using bwa mem v0.7.17
4. alignment shifting using deeptools v3.5.6
5. Removal of the ENCODE blacklist regions using bedtools v2.31.1
6. Peak calling using MACS3 v3.0.2
7. Differential chromatin analysis using DiffBind v3.16

This can be split into two major stages (Figure 2). The first is a quality control preprocessing stage that generates ATAC-Peaks and bigwig files for visualisation in a genome browser. The second uses DiffBind (Rory Stark<Rory.Stark@Cruk.Cam.Ac.Uk>, 2017) to interrorgate the cell-type differences. 

<img width="577" alt="Screenshot 2025-04-18 at 15 15 49" src="https://github.com/user-attachments/assets/9492946f-924f-4fe2-bfb0-05ac6ae8f99f" />

_Figure 2: Schematic of ATAC-Seq pipeline. (A) Paired-end sequenced samples are inputted into a preprocessing script that cleans and produces peaksets (narrowPeak) for further analysis and bigwig files for visualisation in a genome browser. (B) Using a user-provided sample sheet,  a second script applies DiffBind towards each sample peakset and contrasts them based on their cell-type (made using Biorender)_



# How to use: 

## Prerequisites and Caveats

1. This script will only accept paired-end sequencing folders (in fastq.gz format)
   
2. Make sure you have a reference genome installed in .fa format
   
3. A sample-sheet (in .csv format) is required with the columns: SampleID, should be the same as your sample name; Tissue, the tissue used; Condition the different conditions of your experiment; Replicate, the number of replicates of your experiment; Factor, used to separate groups. This is for the diffbind_pipeline.R script to know which samples to look for and compare. For more information, visit this website: https://www.rdocumentation.org/packages/DiffBind/versions/2.0.2/topics/dba.peakset.

4. The program will use the ENCODE blacklist regions HG38v2 co-ordinates to filter problematic regions. For more information, look here: https://zenodo.org/records/1491733

5. In creating and testing the code, I used bwa v0.07.17. However, i placed bwa v0.7.19 into the requirements file due to compatability issues. Both versions produce identical alignments as confirmed by the creator, Heng Li (https://github.com/lh3/bwa/releases)

## Installation
All required packages and libraries can be installed via the **Requirements.txt** file. For example

1. Clone the repository

2. Create an anaconda environment
``` conda create -n <environment> ```

3. Activate the environment
``` conda activate <environment> ```

4. Install the libraries

``` conda install -r requirements.txt  ```


## Protocol 

There are five required parameters for the script to be ran. 
1. -s = sample sheet (in .csv format)
2. -f = sample filename
3. -r = reference genome (the full path in .fa format)
4. -i = input directory
5. -o = output directory


The script can be ran from terminal using:

``` bash <path to script>/atac_pipeline.sh -i <input directory> -o <output directory> -s <sample sheet> -f <sample filename> -r <reference genome> ``` 

If you have a SLURM manager in a HPC environment, you can execute the script using a submission script. For example: _run_atac_pipeline.sub_:
```sbatch run_atac_pipeline.sub```

Where _run_atac_pipeline.sub_ contains the above code in the previous example

# Acknowledgements
I would like to thank my Primary Supervisor, Professor David Monk and my Bioinformatics mentor, Leighton Folkes, for their support in this work


# References
Angarica, V.E. and Del Sol, A. (2017) ‘Bioinformatics Tools for Genome-Wide Epigenetic Research’, Advances in experimental medicine and biology, 978, pp. 489–512. Available at: https://doi.org/10.1007/978-3-319-53889-1_25.

Buenrostro, J.D. et al. (2013) ‘Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position’, Nature Methods, 10(12), pp. 1213–1218. Available at: https://doi.org/10.1038/nmeth.2688.

Chang, P. et al. (2018) ‘Computational Methods for Assessing Chromatin Hierarchy’, Computational and Structural Biotechnology Journal, 16, pp. 43–53. Available at: https://doi.org/10.1016/j.csbj.2018.02.003.

Dillinger, S. (2021) What is ATAC-Seq & How Does it Work?, www.activemotif.com. Available at: https://www.activemotif.com/blog-atac-seq.

Giresi, P.G. et al. (2007) ‘FAIRE (Formaldehyde-Assisted Isolation of Regulatory Elements) isolates active regulatory elements from human chromatin’, Genome Research, 17(6), pp. 877–885. Available at: https://doi.org/10.1101/gr.5533506.

Riva, S.G. et al. (2023) ‘CATCH-UP: A High-Throughput Upstream-Pipeline for Bulk ATAC-Seq and ChIP-Seq Data’, Journal of Visualized Experiments [Preprint], (199). Available at: https://doi.org/10.3791/65633.

Rory Stark<Rory.Stark@Cruk.Cam.Ac.Uk>, G.B.C. (2017) ‘DiffBind’. Bioconductor. Available at: https://doi.org/10.18129/B9.BIOC.DIFFBIND.

Smith, J.P. et al. (2021) ‘PEPATAC: an optimized pipeline for ATAC-seq data analysis with serial alignments’, NAR genomics and bioinformatics, 3(4). Available at: https://doi.org/10.1093/nargab/lqab101.

Song, L. and Crawford, G.E. (2010) ‘DNase-seq: A High-Resolution Technique for Mapping Active Gene Regulatory Elements across the Genome from Mammalian Cells’, Cold Spring Harbor Protocols, 2010(2), p. pdb.prot5384-pdb.prot5384. Available at: https://doi.org/10.1101/pdb.prot5384.


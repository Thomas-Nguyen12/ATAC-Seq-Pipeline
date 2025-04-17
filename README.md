# ATAC-Seq-Pipeline
This project aims to create a pipeline for Differential Chromatin Accessibility analysis using ATAC-Seq. 

ATAC-Seq is a method to identify regions of open chromatin (DNA + histone molecules) across the genome. 

![image](https://github.com/user-attachments/assets/6b1cd157-c989-42ed-996e-b6d9a320b5a6)
_Figure 1: ATAC-Seq schematic. Tn5 Tranposases cleave open DNA and simultaneously attach barcode sequencing primers in a process known as tagmentation. DNA libraries are then sequenced using next-generation sequencing technologies such as Illumina and subject to Bioinformatic Analysis (Dillinger, 2021)_

Compared to alternatives such as DNAse-Seq (Song and Crawford, 2010) and Faire-Seq (Giresi et al., 2007), ATAC-Seq requires fewer input samples and is much more time efficient due to a hyperactive Tn5 Transposase, requiring only 500-50,000 cells (Buenrostro et al., 2013).

To date, however, there are few tools specifically designed for ATAC-Seq analysis (Angarica and Del Sol, 2017), with many CHiP-Seq tools being used in place, assuming similar data properties (Chang et al., 2018). 

Several pipelines incorporate various tools and libraries for end-to-end analysis, such as PEPATAC (Smith et al., 2021) and CATCH-UP (Riva et al., 2023) that have varying focuses on serial alignments or bulk samples. However, none of them account for cell-type differences in chromatin accessibility. Consequently, I have developed a pipeline for paired-end sequencing data that creates ATAC peaksets and assesses the chromatin accessibility between two cell-types. Using this pipeline, users can investigate the functional implications of where DNA is differentially open or closed.

The workflow has two stages (Figure 2). The first is a quality control preprocessing stage that generates ATAC-Peaks and bigwig files for visualisation in a genome browser. The second uses DiffBind to interrorgate the cell-type differences. 

<img width="457" alt="Screenshot 2025-04-17 at 13 21 21" src="https://github.com/user-attachments/assets/2bf93d83-18a5-47ec-a49c-94f28986704c" />
_Figure 2: Schematic of ATAC-Seq pipeline. (A) Paired-end sequenced samples are inputted into a preprocessing script that cleans and produces peaksets and bigwig files for visualisation. (B) Using a user-provided sample sheet,  a second script applies DiffBind towards each sample peakset and contrasts them based on their cell-type (made using Biorender)._



# How to use: 

## Installation
All required packages and libraries can be installed via the **Requirements.txt** file. For example

1. Create an anaconda environment
``` conda create -n <environment> ```

2. Activate the environment
``` conda activate <environment> ```

3. Install the libraries

``` conda install -r requirements.txt  ```


## Protocol 

There are five required parameters for the script to be ran. 
1. -s = sample sheet (in .csv format)
2. -f = sample filename
3. -r = reference genome (in .fa format)
4. -i = input directory
5. -o = output directory


The script can be ran from terminal using:

``` bash <path to script>/atac_pipeline.sh -i <input directory> -o <output directory> -s <sample sheet> -f <sample filename> -r <reference genome> ``` 

If you have a SLURM manager in a HPC environment, you can execute the script using a submission script. For example: _run_atac_pipeline.sub_:
```sbatch run_atac_pipeline.sub```

Where _run_atac_pipeline.sub_ contains the above code in the previous example

# References



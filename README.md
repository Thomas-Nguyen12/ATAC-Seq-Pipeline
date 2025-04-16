# ATAC-Seq-Pipeline
This project aims to create a pipeline for Differential Chromatin Accessibility analysis using ATAC-Seq. 

ATAC-Seq is a method to identify regions of Euchromatin (open chromatin) across the genome. 

![image](https://github.com/user-attachments/assets/6b1cd157-c989-42ed-996e-b6d9a320b5a6)
_Figure 1: ATAC-Seq schematic. Tn5 Tranposases cleave open DNA and simultaneously attach barcode sequencing primers in a process known as tagmentation. DNA libraries are then sequenced using next-generation sequencing technologies such as Illumina and subject to Bioinformatic Analysis (Dillinger, 2021)_



Libraries and tools were inspired by studies investigating the epigenome of Human Cytotrophoblasts and Extravillous Trophoblasts (Varberg et al., 2023).

Processing takes place across various stages. 



# How to use: 

## Installation
All required packages and libraries can be installed via the **Requirements.txt** file found within this repository. 



There are five required parameters for the script to be ran. 
1. -s = sample sheet (in .csv format)
2. -f = sample filename
3. -r = reference genome (in .fa format)
4. -i = input directory
5. -o = output directory


The script can be ran from terminal using
``` bash <path to script>/atac_pipeline.sh -i <input directory> -o <output directory> -s <sample sheet> -f <sample filename> -r <reference genome> ``` 

If you have a SLURM manager in a HPC environment, you can execute the script using a submission script. For example: _run_script.sub_:
```sbatch run_script.sub```

Where _run_script.sub_ contains the above code in the previous example
# References



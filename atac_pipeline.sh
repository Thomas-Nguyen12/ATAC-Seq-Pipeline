#!/bin/bash 
# I need to include an option to parse the samplesheet for diffbind

while getopts ':i:o:r:f:s:' OPTION; do 

	case "$OPTION" in 
		i)
			
			input_directory="$OPTARG"
			echo "Input Directory: $input_directory"
			;;
		o)
			output_directory="$OPTARG"
			echo "Output Directory: $output_directory"
			;;
		r)
			reference="$OPTARG"
			echo "Reference: $reference"
			;;
		s)
			sample_sheet="$OPTARG"
			echo "Sample Sheet: $sample_sheet"
			;;
		f)
			sample_name="$OPTARG"
			echo "Sample Filename: $sample_name"
			;;
		?)
			echo "This tool is used for Comprehensive Atac-Seq analysis..."
			echo "This tool is designed for paired end fastq files..."
			echo "Usage: $(basename $0) [-i input_directory] [-o output_directory] [-r Genome Reference PATH (including file name)] [-f sample filename (in csv format)] [-s sample sheet PATH]"
			exit 1 
			;;
	esac
done

echo "Start time: "
date 



echo "----------------------------"



## 1.1: htstream preprocessing
echo "Running the Preprocessing Pipeline..."

cd $output_directory
echo "Seq Screener..."	
hts_SeqScreener -1 "${input_directory}/${sample_name}_1.fastq.gz" -2  "${input_directory}/${sample_name}_2.fastq.gz" -f "${output_directory}/1.${sample_name}_SeqScreen"
echo "SuperDeduper..."
hts_SuperDeduper -1 "${output_directory}/1.${sample_name}_SeqScreen_R1.fastq.gz" -2 "${output_directory}/1.${sample_name}_SeqScreen_R2.fastq.gz" -e 250000 -f "${output_directory}/2.${sample_name}_SuperDeduper"

echo "Adapter Trimmer..."
hts_AdapterTrimmer -1 "${output_directory}/2.${sample_name}_SuperDeduper_R1.fastq.gz" -2 "${output_directory}/2.${sample_name}_SuperDeduper_R2.fastq.gz" -p 4 -f "${output_directory}/3.${sample_name}_AdapterTrimmer"

echo "Removing unknown nucleotides (NTrimmer)"
hts_NTrimmer -1 "${output_directory}/3.${sample_name}_AdapterTrimmer_R1.fastq.gz" -2 "${output_directory}/3.${sample_name}_AdapterTrimmer_R2.fastq.gz" -f "${output_directory}/4.${sample_name}_NTrimmer"

echo "Base Quality Trimming..."
hts_QWindowTrim -1 "${output_directory}/4.${sample_name}_NTrimmer_R1.fastq.gz" -2 "${output_directory}/4.${sample_name}_NTrimmer_R2.fastq.gz" -q 20 -w 10 -f "${output_directory}/5.${sample_name}_QWindowTrim"

echo "Short Read Trimming..."
hts_LengthFilter -1 "${output_directory}/5.${sample_name}_QWindowTrim_R1.fastq.gz" -2 "${output_directory}/5.${sample_name}_QWindowTrim_R2.fastq.gz" -n -m 50 -f "${output_directory}/6.${sample_name}_LengthFilter" 

echo "Preprocessing Successful"	




# Stage 2: Alignment
module add bwa/0.7.17	
echo "Starting alignment process..."
echo "Checking if bwa indexes exist." 	

## 2.1: Building BWA indexes

# Checking if the bwa indexes already exist

if find "$output_directory" -name "bwa_index*" | grep -q .; then
    echo "Found existing bwa indexes..."
else
    echo "No bwa indexes found. Building bwa indexes"
    bwa index "$reference" -p "${output_directory}/bwa_index"
fi

echo "Aligning the samples..."
#What do the bwa indexes look like? 
module add samtools/1.21	
## 2.2: BWA alignment (remember the path of the bwa index) 
bwa mem ${output_directory}/bwa_index ${output_directory}/6.${sample_name}_LengthFilter_R1.fastq.gz ${output_directory}/6.${sample_name}_LengthFilter_R2.fastq.gz | samtools sort -o ${output_directory}/${sample_name}_sorted.bam



# Stage 3: Post Alignment Processing

## 3.1: indexing the bam files (this should produce .bai files)

echo "indexing the samples for alignment Sieve..."
samtools index ${output_directory}/${sample_name}_sorted.bam -o ${output_directory}/${sample_name}_sorted.bam.bai


## 3.2 Shifting the alignments using AlignmentSieve
echo "running alignment sieve..."
alignmentSieve -b ${output_directory}/${sample_name}_sorted.bam -o ${output_directory}/${sample_name}_alignmentSieve.bam

## 3.3 Removing the ENCODE blacklist regions
echo "Removing the ENCODE blacklist regions..."
encode_blacklist_regions=/gpfs/home/wmp21vtu/scratch/atac/encode_blacklist_regions/Boyle-Lab-Blacklist-f4a45ab/lists/hg38-blacklist.v2.bed.gz

module add bedtools/2.31.1
bedtools intersect -a ${output_directory}/${sample_name}_alignmentSieve.bam -b $encode_blacklist_regions -v > ${output_directory}/${sample_name}_without_blacklist.bam

# Stage 4: Peak Calling using MACS3
echo "Calling Peaks..."
macs3 callpeak -t ${output_directory}/${sample_name}_without_blacklist.bam -f BAMPE -g hs -q 0.01 --outdir $output_directory --name ${sample_name}_peak

## calling BIGWIGS
# include the bam files that are filtered (without the blacklist)
# bamCoverage is a tool found in deeptools
echo "creating bigwig files..."

# I need to index the bam file first

samtools index ${output_directory}/${sample_name}_without_blacklist.bam 
bamCoverage -b ${output_directory}/${sample_name}_without_blacklist.bam -o ${output_directory}/${sample_name}_bigwig.bw



## RUNNING DIFFBIND
# I need to make sure that all the files are available 
# Maybe I can read the sample sheet and find if the sampleIDs have a corresponding narrowPeak file

# Reading the sampleSheet 

# 1. Read the sample sheet 
## the sampleID column will be the first column 
echo "Checking the sample sheet..."

# Read the first column of the sample sheet, ignoring header
count=0
sample_sheet_length=$(wc -l < $sample_sheet)
while IFS= read -r sample; do
    if find "$output_directory" -name "*${sample}_peak_peaks.narrowPeak" | grep -q .; then
        echo "Peakset File found for $sample"
        ((count++))
    else
        echo "Peakset File not found for $sample"
    fi
done < <(awk -F',' 'NR>1 {print $1}' "$sample_sheet")

# Print total count of found files
echo "Total files found: $count"

# Check if all expected files exist
if [ "$count" -eq "$sample_sheet_length" ]; then
    echo "All files found. Starting DiffBind..."
    Rscript /gpfs/home/wmp21vtu/pipeline/atac_seq/diffbind_pipeline.R "$sample_sheet" "$output_directory"
else
    echo "Sample Sheet not ready for DiffBind..."
    echo "DiffBind Pending until all samples are processed"
fi

echo "Sample Preprocessing Completed at:"
date
echo "-------------------------------------"




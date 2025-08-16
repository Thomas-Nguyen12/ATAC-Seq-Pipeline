#!/bin/bash 
# I need to include an option to parse the samplesheet for diffbind
# I can include a parameter that specifies which contrasts are wanted for DiffBind (I can also do this in the diffbind_pipeline.R script) 
# If no contrast is specified, then contrast 1 is assumed 
d=0



while getopts ':i:o:r:f:s:' OPTION; do 

	case "$OPTION" in 
		i)
			
			input_directory="$OPTARG"
			;;
		o)
			output_directory="$OPTARG"
			;;
		r)
			reference="$OPTARG"
			;;
		s)
			sample_sheet="$OPTARG"
			;;
		f)
			sample_name="$OPTARG"
			;;
		d)
			diffbind_enabled="$OPTARG"
			if [[ $diffbind_enabled -ne 0 && $diffbind_enabled -ne 1 ]]; then 
				echo "Error. <-d 0 = DiffBind DISABLED>... <-d 1 = DiffBind ENABLED>"
				exit 1 
			fi
			;; 
		?)
			echo "This tool is used for Comprehensive Atac-Seq analysis..."
			echo "This tool is designed for paired end fastq files..."
			echo "Usage: $(basename $0) [-i input_directory] [-o output_directory] [-r Genome Reference PATH (including file name)] [-f sample filename (in csv format)] [-s sample sheet PATH] [-d enable_diffbind]"
			exit 1 
			;;
	esac
done

start_time=$(date)
echo "Start time: ${start_time}"

echo "----- SETTINGS -----"
echo "input directory --> ${input_directory}"
echo "output directory --> ${output_directory}"
echo "reference --> ${reference}"
echo "sample sheet --> ${sample_sheet}"
echo "sample name --> ${sample_name}"

echo "----------------------------"
# installing the ENCODE blacklist regionsi
echo "Checking if the ENCODE blacklist regions are installed in your input directory..." 
cd $input_directory
DIRECTORY=Boyle-Lab-Blacklist-f4a45ab
if [ ! -d "$DIRECTORY" ]; then
	echo "$DIRECTORY does not exist. Installing ENCODE blacklist regions..."
	echo "This program will use the HG38v2 co-ordinates..."
	curl https://zenodo.org/records/1491733/files/Boyle-Lab/Blacklist-v2.0.zip?download=1 --output ENCODE.zip
	unzip ENCODE.zip
	download_blacklist_exit_code=$?
	if [ $download_blacklist_exit_code -ne 0 ]; then
		echo "There was an error downloading the blacklist regions. Exiting..."
		exit 1
	fi

else
	echo "$DIRECTORY does exist. Continuing on..."

  
fi


encode_blacklist_regions=${input_directory}/Boyle-Lab-Blacklist-f4a45ab/lists/hg38-blacklist.v2.bed.gz

echo "The encode blacklist regions are installed in ${encode_blacklist_regions}"

## 1.1: htstream preprocessing
echo "Running the Preprocessing Pipeline..."

cd $output_directory

# checking if the input files are in _1.fastq.gz or _R1.fastq.gz


if find "$input_directory" -name "${sample_name}_R1.fastq.gz" | grep -q .; then
    echo "Found the input files in _R1.fastq.gz format"
    echo "SeqScreener..."
    hts_SeqScreener -1 "${input_directory}/${sample_name}_R1.fastq.gz" \
                    -2 "${input_directory}/${sample_name}_R2.fastq.gz" \
                    -f "${output_directory}/1.${sample_name}_SeqScreen"
else
    echo "Found input files in alternate format"
    echo "SeqScrenner..."
    hts_SeqScreener -1 "${input_directory}/${sample_name}_1.fastq.gz" \
                    -2 "${input_directory}/${sample_name}_2.fastq.gz" \
                    -f "${output_directory}/1.${sample_name}_SeqScreen"
fi

SeqScreener_exit_code=$?

if [ $SeqScreener_exit_code -ne 0 ]; then 
	echo "There was an error with the SeqScreener. Exiting... " 
	exit 1
fi




## The rest of the code should follow as normal

echo "SuperDeduper..."
hts_SuperDeduper -1 "${output_directory}/1.${sample_name}_SeqScreen_R1.fastq.gz" -2 "${output_directory}/1.${sample_name}_SeqScreen_R2.fastq.gz" -e 250000 -f "${output_directory}/2.${sample_name}_SuperDeduper"


SuperDeduper_exit_code=$?

if [ $SuperDeduper_exit_code -ne 0 ]; then
	echo "There was an error with the SuperDeduper. Exiting..." 
	exit 1
fi


echo "Adapter Trimmer..."
hts_AdapterTrimmer -1 "${output_directory}/2.${sample_name}_SuperDeduper_R1.fastq.gz" -2 "${output_directory}/2.${sample_name}_SuperDeduper_R2.fastq.gz" -p 4 -f "${output_directory}/3.${sample_name}_AdapterTrimmer"

AdapterTrimmer_exit_code=$?
if [ $AdapterTrimmer_exit_code -ne 0 ]; then
	echo "There was an error with the Adapter Trimmer. Exiting..."
	exit 1
fi


echo "Removing unknown nucleotides (NTrimmer)"
hts_NTrimmer -1 "${output_directory}/3.${sample_name}_AdapterTrimmer_R1.fastq.gz" -2 "${output_directory}/3.${sample_name}_AdapterTrimmer_R2.fastq.gz" -f "${output_directory}/4.${sample_name}_NTrimmer"

NTrimmer_exit_code=$?
if [ $NTrimmer_exit_code -ne 0 ]; then
	echo "There was an error with the NTrimmer. Exiting..." 
	exit 1
fi



echo "Base Quality Trimming..."
hts_QWindowTrim -1 "${output_directory}/4.${sample_name}_NTrimmer_R1.fastq.gz" -2 "${output_directory}/4.${sample_name}_NTrimmer_R2.fastq.gz" -q 20 -w 10 -f "${output_directory}/5.${sample_name}_QWindowTrim"

QWindow_trim_exit_code=$?
if [ $QWindow_trim_exit_code -ne 0 ]; then
	echo "There was an error with the QWindowTrim. Exiting..." 
	exit 1
fi




echo "Short Read Trimming..."
hts_LengthFilter -1 "${output_directory}/5.${sample_name}_QWindowTrim_R1.fastq.gz" -2 "${output_directory}/5.${sample_name}_QWindowTrim_R2.fastq.gz" -n -m 50 -f "${output_directory}/6.${sample_name}_LengthFilter" 


LengthFilter_exit_code=$?
if [ $LengthFilter_exit_code -ne 0 ]; then
	echo "There was an error with the LengthFilter. Exiting..." 
	exit 1
fi


echo "Preprocessing Successful"	




# Stage 2: Alignment
echo "Starting alignment process..."
echo "Checking if bwa indexes exist." 	

## 2.1: Building BWA indexes
# bwa should be a part of the anaconda environment 'pipeline' 

# Checking if the bwa indexes already exist

if find "$output_directory" -name "bwa_index*" | grep -q .; then
    echo "Found existing bwa indexes..."
else
    echo "No bwa indexes found. Building bwa indexes"
    bwa index $reference -p ${output_directory}/bwa_index
fi

bwa_index_exit_code=$?
if [ $bwa_index_exit_code -ne 0 ]; then
	echo "There was an error with the bwa index. Check if Bwa indexes already exist in your output directory"
	exit 1
fi


echo "Aligning the samples..."
#What do the bwa indexes look like? 
## 2.2: BWA alignment (remember the path of the bwa index) 
bwa mem ${output_directory}/bwa_index ${output_directory}/6.${sample_name}_LengthFilter_R1.fastq.gz ${output_directory}/6.${sample_name}_LengthFilter_R2.fastq.gz | samtools sort -o ${output_directory}/${sample_name}_sorted.bam
bwa_mem_exit_code=$?
if [ $bwa_mem_exit_code -ne 0 ]; then
	echo "There was an error with the alignment. Exiting..." 
	exit 1
fi



# Stage 3: Post Alignment Processing

## 3.1: indexing the bam files (this should produce .bai files)

echo "indexing the samples for alignment Sieve..."
samtools index ${output_directory}/${sample_name}_sorted.bam -o ${output_directory}/${sample_name}_sorted.bam.bai

samtools_index_exit_code=$?
if [ $samtools_index_exit_code -ne 0 ]; then 
	echo "There was an error indexing the alignment file. Exiting..." 
	exit 1
fi



## 3.2 Shifting the alignments using AlignmentSieve
echo "running alignment sieve..."
alignmentSieve -b ${output_directory}/${sample_name}_sorted.bam -o ${output_directory}/${sample_name}_alignmentSieve.bam

alignmentSieve_exit_code=$?
if [ $alignmentSieve_exit_code -ne 0 ]; then
	echo "There was an error with the alignmentSieve. Exiting..." 
	exit 1
fi



## 3.3 Removing the ENCODE blacklist regions
echo "Removing the ENCODE blacklist regions..."

bedtools intersect -a ${output_directory}/${sample_name}_alignmentSieve.bam -b $encode_blacklist_regions -v > ${output_directory}/${sample_name}_without_blacklist.bam

bedtools_intersect_exit_code=$?

if [ $bedtools_intersect_exit_code -ne 0 ]; then
	echo "There was an error with removing the ENCODE blacklist regions. Exiting..." 
	exit 1
fi



# Stage 4: Peak Calling using MACS3
echo "Calling Peaks..."
macs3 callpeak -t ${output_directory}/${sample_name}_without_blacklist.bam -f BAMPE -g hs -q 0.01 --outdir $output_directory --name ${sample_name}_peak

macs3_exit_code=$?
if [ $macs3_exit_code -ne 0 ]; then
	echo "There was an error calling the Peak-Sets from the alignment files. Exiting..." 
	exit 1
fi


## calling BIGWIGS
# include the bam files that are filtered (without the blacklist)
# bamCoverage is a tool found in deeptools
echo "creating bigwig files..."

# I need to index the bam file first

samtools index ${output_directory}/${sample_name}_without_blacklist.bam
samtools_index_bigwig_exit_code=$?

if [ $samtools_index_bigwig_exit_code -ne 0 ]; then
	echo "There was an error indexing the alignment file before creating the bigWig files. Exiting..." 
	exit 1
fi
 
bamCoverage -b ${output_directory}/${sample_name}_without_blacklist.bam -o ${output_directory}/${sample_name}_bigwig.bw

bamCoverage_exit_code=$?

if [ $bamCoverage_exit_code -ne 0 ]; then
	echo "There was an error creating the bigWig files. Exiting..." 
	exit 1
fi


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

echo "========================= CHECKING ALL PEAKSETS =========================" 

while IFS= read -r sample; do
    if find "$output_directory" -name "*${sample}_peak_peaks.narrowPeak" | grep -q .; then
	
        echo "[COMPLETE]    $sample"
        ((count++))
    else
        echo "[PENDING]    $sample"
    fi
done < <(awk -F',' 'NR>1 {print $1}' "$sample_sheet")

# Print total count of found files
# I can also include the total number of files as a pseudo-progress meter
total_files=$(wc -l < ${sample_sheet}) 
 
echo "Total samples processed: $count out of $total_files"
echo "${count} / ${total_files} complete" 
# Check if all expected files exist
if [ "$count" -eq "$sample_sheet_length" ]; then
    echo "All files processed. Proceeding..."
    if [ $diffbind_enabled -eq 1 ]; then
	
	echo "DiffBind enabled. Executing..."	
	Rscript /gpfs/home/wmp21vtu/pipeline/atac_seq/diffbind_pipeline.R "$sample_sheet" "$output_directory"
    else
        echo "DiffBind Disabled. Exiting..."
else
    echo "Sample Sheet not ready for DiffBind..."
    echo "DiffBind Pending until all samples are processed"
fi

echo "Sample Preprocessing Completed at:"
date
echo "-------------------------------------"




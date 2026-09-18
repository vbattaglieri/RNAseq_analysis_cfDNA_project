#!/usr/bin/bash

set -e

### RNA-seq pipeline ###

#QC and rsem quantification on GRCh38 with gencode v44
#setup variables and env

log_file="rna-seq_"$(date +"%d-%m-%Y__%H.%M")".log"

# ---- User configuration: edit these paths for your environment ----
BASE_DIR="/path/to/install"           # root holding miniconda3/ and references/
CONDA_ENV="rnaseq"                    # conda env with fastqc, multiqc, rsem, rseqc, samtools, star
# -------------------------------------------------------------------

source ${BASE_DIR}/miniconda3/bin/activate ${CONDA_ENV}

echo -e "\n\n\t=====> RNA-seq QC - rsem - GRCh38 - gencode v44  <=====\n\n" > $log_file

if [ $# -ne 1 ]; then
	echo "===============================================================" >> $log_file
	echo -e "\n$0\n\nNeeded 1 argument: folder file containing the fastq.gz files (path without the final slash)\n" >> $log_file
	echo "===============================================================" >>  $log_file
	exit 1
fi

#Argument dir
folder_fastq="$1"

#define variables
housekeeping_genes="${BASE_DIR}/references/hg38/hg38.HouseKeepingGenes.bed"

#executable path (resolved from the active conda env)
env_bin="${BASE_DIR}/miniconda3/envs/${CONDA_ENV}/bin"
fastqc_exe="${env_bin}/fastqc"
multiqc_exe="${env_bin}/multiqc"
rsem_calc_expr_exe="${env_bin}/rsem-calculate-expression"
rseqc_exe="${env_bin}"


#def functions
fastqc() {
	local folder=$1
	local fastq_array=("$folder"/*fastq.gz) # store all fastq.gz file in an array
	mkdir -p fastqc # folder for results
	$fastqc_exe ${fastq_array[@]} --threads 40 --outdir fastqc >> $log_file

}

multiqc() {
	local folder=$1
	mkdir -p quality_control/multiqc
	$multiqc_exe -f -o quality_control/multiqc/ $folder >> $log_file
}	


rsem_gencode_v44() {
	mkdir -p rsem
	local folder=$1
	local fastq_array_R1=("${folder}"/*"1.fastq.gz") # store all R1 fastq.gz file in an array
	local name_array=()
	
	for i in "${fastq_array_R1[@]}";
	do
		local name=$(basename "$i" _1.fastq.gz)
		name_array+=("$name")
	done

	for i in "${name_array[@]}"
	do
		echo "${folder}/${i}_1.fastq.gz"
		"$rsem_calc_expr_exe" \
			--star \
			--star-gzipped-read-file \
			--paired-end \
			--strandedness reverse \
			--star-output-genome-bam \
			--num-threads 40 \
			"${folder}/${i}_1.fastq.gz" \
			"${folder}/${i}_2.fastq.gz" \
			${BASE_DIR}/references/hg38/gencode_v44/GRCh38_gencode.v44 \
			rsem/${i}

	done
}

sort_and_index() {
	local folder=$1
	
	for i in "${folder}"/*"genome.bam";
	do
		echo $i
		samtools sort \
			-@ 10 \
			-m 3G \
			-o "${folder}/$(basename "$i" .bam).sorted.bam" \
			$i
		samtools index \
			-@10 \
			"${folder}/$(basename "$i" .bam).sorted.bam"
done

}

rseQC() {
	mkdir -p rseQC
	local folder=$1
	local bam_array=("${folder}"/*".STAR.genome.sorted.bam") # store all sorted bam file in an array
	local bai_array=("${folder}"/*".STAR.genome.sorted.bam.bai") # store all sorted bai file in an array
	
	for i in "${bam_array[@]}"
	do
		name_sample=$(basename "$i" .STAR.genome.sorted.bam)
		"${rseqc_exe}/bam_stat.py" \
	 		-i "$i" > rseQC/${name_sample}"_bam_stat.out"

		"${rseqc_exe}/infer_experiment.py" \
	 		-i "$i" \
	 		-r $housekeeping_genes > rseQC/${name_sample}"_infer_experiment.out"
	done
	
	mkdir -p bam_temp
	for i in "${bam_array[@]}"
	do
		ln -s "$(readlink -f "$i")" "bam_temp/$(basename "$i")"
	done
	
	for i in "${bai_array[@]}"
	do
		ln -s "$(readlink -f "$i")" "bam_temp/$(basename "$i")"
	done

	"${rseqc_exe}/geneBody_coverage.py" \
	 	-r $housekeeping_genes \
	 	-i bam_temp \
	 	-o rseQC/all_geneBody_coverage

		
}



### START SCRIPT ###

echo -e "--> FastQC...\n" >> $log_file
fastqc $folder_fastq
fastqc_folder="fastqc" #store output path in a variable
echo -e "--> DONE\n\n" >> $log_file


#echo -e "--> MultiQC...\n" >> $log_file
#mkdir -p quality_control
#multiqc $fastqc_folder
#echo -e "--> DONE\n" >> $log_file


echo -e "--> rsem...\n" >> $log_file
rsem_gencode_v44 $folder_fastq
echo -e "--> DONE\n" >> $log_file


echo -e "--> sort_and_index...\n" >> $log_file
sort_and_index rsem
echo -e "--> DONE\n" >> $log_file


echo -e "--> rseQC...\n" >> $log_file
rseQC rsem
echo -e "--> DONE\n" >> $log_file

echo "ANALYSIS COMPLETED"


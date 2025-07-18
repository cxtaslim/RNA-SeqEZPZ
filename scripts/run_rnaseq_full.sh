#!/bin/bash

# print linenumber when set -x is set
export PS4='+${LINENO}:${BASH_SOURCE}: '

# script connecting all the individual scripts to do full rnaseq analysis
# How to run
# cd <project_dir>
# bash scripts/run_rnaseq_full.sh &> run_rnaseq_full.out &
# DO NOT change the name or location of run_rnaseq_full.out
# Examples:
# bash scripts/run_rnaseq_full.sh &> run_rnaseq_full.out &
# to run with specific time limit:
# bash scripts/run_rnaseq_full.sh \
#        	time=DD-HH:MM:SS &> run_rnaseq_full.out &
#
# by default, alignment is done to human reference genome hg19 unless specified using 'genome=hg38':
# bash scripts/run_rnaseq_full.sh genome=hg38 &> run_rnaseq_full.out &
#
# or to do nothing but echo all commands:
# bash scripts/run_rnaseq_full.sh run=echo &> run_rnaseq_full.out &
#
# or to change a bunch of parameters at once:
# bash scripts/run_rnaseq_full.sh run=echo padj=1 time=DD-HH:MM:SS \
# &> run_rnaseq_full.out &
#
# or to run and printing all trace commands (i.e. set -x):
# bash scripts/run_rnaseq_full.sh run=debug &> run_rnaseq_full.out &


#set -x
set -e
# set 'run' to echo to simply echoing all commands
# set to empty to run all commands
# clear variable used for optional arguments
unset run time PYTHONPATH
# get command line arguments
while [[ "$#" -gt 0 ]]; do
	if [[ $1 == "run"* ]];then
		run=$(echo $1 | cut -d '=' -f 2)
		shift
	fi
	if [[ $1 == "padj"* ]];then
		padj=$(echo $1 | cut -d '=' -f 2)
		shift
	fi
	if [[ $1 == "time"* ]];then
		time=$(echo $1 | cut -d '=' -f 2)
		shift
	fi
	if [[ $1 == "genome"* ]];then
		ref_ver=$(echo $1 | cut -d '=' -f 2)
		shift
	fi
	if [[ $1 == "ref_fa"* ]];then
                ref_fa=$(echo $1 | cut -d '=' -f 2)
                shift
        fi
	if [[ $1 == "ref_gtf"* ]];then
                ref_gtf=$(echo $1 | cut -d '=' -f 2)
                shift
        fi
	if [[ $1 == "batch_adjust"* ]];then
                batch_adjust=$(echo $1 | cut -d '=' -f 2)
                shift
        fi
	if [[ $1 == "ncpus_trim"* ]];then
                ncpus_trim=$(echo $1 | cut -d '=' -f 2)
                shift
        fi
	if [[ $1 == "ncpus_star"* ]];then
                ncpus_star=$(echo $1 | cut -d '=' -f 2)
                shift
        fi
	if [[ $1 == "help" ]];then
		echo ""
		echo 'usage: bash scripts/run_rnaseq_full.sh [OPTION] &> run_rnaseq_full.out &'
		echo ''
		echo DESCRIPTION
		echo -e '\trun full RNA-seq analysis: Quality Control, alignment, and differential analysis'
		echo ''
		echo OPTIONS
		echo ''
		echo help 
		echo -e '\tdisplay this help and exit'
		echo run=echo
		echo -e "\tdo not run, echo all commands. Default is running all commands"
		echo -e "\tif set to "debug", it will run with "set -x""
		echo -e "genome=hg19"
		echo -e "\tset reference genome. Default is hg19. Other option: hg38"
		echo -e "\tif using genome other than hg19 or hg38, need to specify both ref_fa and ref_gtf."
                echo -e "ref_fa=/path/to/ref.fa"
                echo -e "\tif using genome other than hg19 or hg38, need to specify ref_fa with path to fasta file"
                echo -e "\tof the reference genome."
                echo -e "ref_gtf=/path/to/ref.gtf"
                echo -e "\tif using genome other than hg19 or hg38, need to specify ref_gtf with path to gtf file"
                echo -e "\tof the reference genome."
		echo padj=0.05
		echo -e "\tset FDR of differential genes (as calculated by DESeq2) < 0.05. Default=0.05"
		echo -e "time=1-00:00:00"
		echo -e "\tset SLURM time limit time=DD-HH:MM:SS, where ‘DD’ is days, ‘HH’ is hours, etc."
		echo -e "\tDefault is 1 day.\n"
		echo batch_adjust=yes
                echo -e "\tby default differential analysis was done with replicate batch adjustment."
                echo -e "\tto turn off batch adjustment, set to no.\n"
		echo ncpus_trim=4
                echo -e "\tset the number of cpus for trimming"
                echo -e "\tDefault is 4 cpus."
                echo -e "\tSet to 1 if pigz error occurs."
		echo ncpus_star=20
                echo -e "\tspecifies the number of cpus used for STAR, bamCompare and feature counts."
                echo -e "\tDefault is 20 cpus."
                echo -e "\tthe number of cpus will be automatically adjusted to max number of cpus if it is"
                echo -e "\tless than the max available cpus."
	exit
	fi
done
date
# converting samples.txt to unix format to remove any invisible extra characters
dos2unix -k samples.txt &> /dev/null

# set default parameters
# note run parameter is set differently here.
if [[ -z "$run" ]];then
	run=
fi
if [[ -z "$padj" ]];then
	padj=0.05
fi
if [[ -z "$time" ]];then
	time=1-00:00:00
fi
if [[ -z "$ref_ver" ]];then
	ref_ver=hg19
fi
if [[ -z "$batch_adjust" ]];then
        batch_adjust=yes
fi
if [[ -z "$ncpus_trim" ]];then
        ncpus_trim=4
fi
if [[ -z "$ncpus_star" ]];then
        ncpus_star=20
fi

# singularity image directory
# found based on location of this script
img_dir=$(dirname $(dirname $(readlink -f $0)))

genome_dir=$img_dir/ref/$ref_ver
# set ref_fa to fasta file in genome_dir if variable ref_fa is not defined
if [[ -z $ref_fa ]];then
	# adding -L to avoid error when symbolic link is used
    	fasta_file=${genome_dir}/$(find -L $genome_dir -name *.fasta -o -name *.fa | xargs basename)
else
	# else set to ref_fa
        fasta_file=$ref_fa
fi

# set genome_dir and star_index_dir accordingly
genome_dir=$(dirname $fasta_file)
star_index_dir=$genome_dir/STAR_index

# set gtf file to ref_gtf in genome_dir if variable ref_gtf is not defined
if [[ -z $ref_gtf ]];then
	# adding -L to avoid error when symbolic link is used
    	gtf_file=${genome_dir}/$(find -L $genome_dir -name *.gtf | xargs basename)
else
        gtf_file=$ref_gtf
fi

# throws an error if gtf_file or fasta_file doesn't exist
if [[ ! -f $gtf_file || ! -f $fasta_file ]];then
        echo "Please check your genome. Either fasta file or gtf file is not found\n"
        exit 1
fi

# run parameter needs to be set differently here
if [[ $run == "debug"* ]];then
        set -x
	run=
	# this is for passing to individual script
	run_debug=debug
fi

# project directory
proj_dir=$(pwd)
cd $proj_dir


echo -e "\nRunning full RNA-seq analysis\n"
echo -e "Options used to run:"
echo padj="$padj"
echo time="$time"
echo genome="$ref_ver"
echo batch_adjust="$batch_adjust"
echo ncpus_trim="$ncpus_trim"
echo ""

echo -e "\nUsing singularity image and scripts in:" ${img_dir} "\n"
# IMPORTANT: It is assumed that:
# scripts to run analysis are in $img_dir/scripts
# reference to run analysis are in $img_dir/ref

# getting SLURM configuration
source $img_dir/scripts/slurm_config_var.sh

work_dir=$proj_dir/outputs
echo -e "All outputs will be stored in $work_dir\n"
if [ ! -d $work_dir ]; then
        mkdir $work_dir
fi
log_dir=$work_dir/logs
if [ ! -d $log_dir ]; then
        mkdir $log_dir
fi
echo -e "Logs and scripts ran will be stored in $log_dir\n"

# copying this script for records
if [[ ! -d $log_dir/scripts ]]; then
	mkdir -p $log_dir/scripts
fi
$(cp $img_dir/scripts/run_rnaseq_full.sh $log_dir/scripts/)

# getting samples info from samples.txt
$(sed -e 's/[[:space:]]*$//' samples.txt | sed 's/"*$//g' | sed 's/^"*//g' > samples_tmp.txt)
$(mv samples_tmp.txt samples.txt)
groupname_array=($(awk '!/#/ {print $1}' samples.txt))
repname_array=($(awk '!/#/ {print $3}' samples.txt))
email=$(awk '!/#/ {print $5;exit}' samples.txt | tr -d '[:space:]')
path_to_r1_fastq=($(awk '!/#/ {print $6}' samples.txt))
path_to_r2_fastq=($(awk '!/#/ {print $7}' samples.txt))

### check and run star index if it doesn't exist
echo ""
echo Check and generate STAR index if genome has not been indexed yet....
echo "See run_star_index.out for more info"
echo ""
cp run_rnaseq_full.out $log_dir/

. $img_dir/scripts/run_star_index.sh run=$run_debug time=$time genome=$ref_ver ref_fa=$fasta_file \
	ref_gtf=$gtf_file ncpus_star=$ncpus_star \
	&> run_star_index.out
# also copy to log_dir
cp run_star_index.out $log_dir/
cp run_rnaseq_full.out $log_dir/

# skip checking job if not generating star index.
if [[ $skip_run_star_index == 0 ]];then
	tmp0=$($run sbatch --dependency=$jid0 \
		--partition=$general_partition \
                --time=5:00 \
                --output=$log_dir/dummy_run_star_index.txt \
                --job-name=run_star_index \
        	--wait \
	        --export message="$message",proj_dir=$proj_dir \
 		$addtl_opt \
                --wrap "echo -e \"$message\" >> $proj_dir/run_rnaseq_full.out; \
						cp $proj_dir/run_rnaseq_full.out $log_dir/"| cut -f 4 -d' ')
	# message if jobs never satisfied
	check_jid0=$(echo $jid0 | sed 's/:/,/g')
	# check to make sure jobs are completed. Print messages if not.
	msg_ok="Done checking and generating STAR index as needed.\n"
	msg_fail="Checking and/or generating STAR index failed. Please check run_star_index.out\n"
	jid_to_check=$check_jid0,$tmp0
	out_file=$proj_dir/run_rnaseq_full.out
	check_run_star_index_jid=$($run sbatch \
        --partition=$general_partition \
        --output=$log_dir/check_run_star_index.out \
        --mail-type=END \
        --mail-user=$email \
        --wait \
        --time=$time \
        --parsable \
	$addtl_opt \
        --job-name=check_run_star_index \
        --export=out_file="$out_file",jid_to_check="$jid_to_check",msg_ok="$msg_ok",msg_fail="$msg_fail" \
        --wrap "bash $img_dir/scripts/check_job.sh")
else
	# message to run_rnaseq_full skipping run_star_index
	echo -e "run_star_index.sh was not run since genome index already exists.\n"
fi

### trimming and QC
date
echo -e "\nChecking whether trimming already ran to completion"
# set default to not skip trim
skip_run_trim_qc=0
cp $proj_dir/run_rnaseq_full.out $log_dir/
# check if trim_fastqc file exist
# enable nullglob for check_trim
shopt -s nullglob
# need to do this instead of $(ls ${proj_dir}/outputs/logs/trim_fastqc_*.out | wc -l)
# to avoid error when files not there
check_trim=(${proj_dir}/outputs/logs/trim_fastqc_*.out)
# disable nullglob
shopt -u nullglob
#echo ${#check_trim[@]}
if [[ ${#check_trim[@]} -gt 0 ]]; then
        # check to make sure all trimming were run to completion
        ncomplete=$(grep "Analysis complete" ${proj_dir}/outputs/logs/trim_fastqc*.out | wc -l)
        # divide by 2 for paired reads
	ncomplete=$((ncomplete/2))
	if [[ $ncomplete -eq ${#groupname_array[@]} ]];then
        	echo "Skip trimming since it's already done."
        	cp $proj_dir/run_rnaseq_full.out $log_dir/
		skip_run_trim_qc=1
	fi
fi

if [[ $skip_run_trim_qc == 0 ]]; then
	echo ""
        echo Trimming and QC.....see progress in run_trim_qc.out
        echo ""
        cp $proj_dir/run_rnaseq_full.out $log_dir/
        cd $proj_dir
        . $img_dir/scripts/run_trim_qc.sh run=$run_debug time=$time ncpus_trim=$ncpus_trim &> run_trim_qc.out
        cp run_trim_qc.out $log_dir/

        tmp1=$($run sbatch --dependency=afterok:$jid2 \
        	$addtl_opt \
	        --partition=$general_partition \
		--time=5:00 \
                --output=$log_dir/dummy_run_trim_qc.txt \
                --job-name=run_trim_qc \
                --export message="$message",proj_dir=$proj_dir \
                --wrap "echo -e \"$message\" >> $proj_dir/run_rnaseq_full.out; \
                cp $proj_dir/run_rnaseq_full.out $log_dir/"| cut -f 4 -d' ')
        # message if jobs never satisfied
        check_jid2=$(echo $jid2 | sed 's/:/,/g')
# check to make sure jobs are completed. Print messages if not.
msg_ok="\n\nDone trimming and QC.\n"
msg_ok="${msg_ok}See run_trim_qc.out for more detailed info."
msg_fail="Trimming and/or QC failed. Please check run_trim_qc.out\n"
jid_to_check=$check_jid2,$tmp1
out_file=$proj_dir/run_rnaseq_full.out
check_run_trim_jid=$($run sbatch \
        --partition=$general_partition \
        $addtl_opt \
	--output=$log_dir/check_run_trim.out \
        --mail-type=END \
        --mail-user=$email \
        --wait \
        --time=$time \
        --parsable \
        --job-name=check_run_trim \
        --export=out_file="$out_file",jid_to_check="$jid_to_check",msg_ok="$msg_ok",msg_fail="$msg_fail" \
        --wrap "bash $img_dir/scripts/check_job.sh")
fi

### skip alignment if STAR pass2 is done
date
echo -e "\nChecking whether alignment already ran to completion"
# set default to not skip alignment
skip_run_align_create_tracks_rna=0
# check if star_pass2 files exist
# enable nullglob for check_trim
shopt -s nullglob
# need to do this instead of $(ls ${proj_dir}/outputs/logs/star_pass2_*.out | wc -l)
# to avoid error when files not there
check_align=(${proj_dir}/outputs/logs/star_pass2_*.out)
# disable nullglob
shopt -u nullglob
if [[ ${#check_align[@]} -gt 0 ]]; then
        # check to make sure all star_pass2 were run to completion
        ncomplete=$(grep "finished successfully" ${proj_dir}/outputs/logs/star_pass2_*.out | wc -l)
        if [[ $ncomplete -eq ${#groupname_array[@]} ]];then
  	      	echo "Skip alignment since it's already done."
		skip_run_align_create_tracks_rna=1
	fi
fi

if [[ $skip_run_align_create_tracks_rna == 0 ]]; then
	# run alignment and create tracks
	echo -e "Aligning reads and creating tracks for visualization....\n"
	echo -e "See progress in run_align_create_tracks_rna.out.\n"
        cp $proj_dir/run_rnaseq_full.out $log_dir/
	cd $proj_dir
	. $img_dir/scripts/run_align_create_tracks_rna.sh run=$run_debug time=$time genome=$ref_ver \
        	ref_fa=$fasta_file ref_gtf=$gtf_file ncpus_star=$ncpus_star &> \
		run_align_create_tracks_rna.out
	cp $proj_dir/run_rnaseq_full.out $log_dir/

	message="\nDone aligning and creating tracks for visualization.\n"
	message="${message}See run_align_create_tracks_rna.out for more info.\n"
	tmp=$($run sbatch --dependency=afterok:$jid4c \
               	--partition=$general_partition \
		--time=5:00 \
               	--output=$log_dir/dummy_run_align_create_tracks_rna.txt \
               	--job-name=check_run_align \
		--wait \
		$addtl_opt \
               	--export message="$message",proj_dir=$proj_dir \
               	--wrap "echo -e \"$message\" >> $proj_dir/run_rnaseq_full.out" | cut -f 4 -d' ')
fi

### Running differential genes analysis
cp $proj_dir/run_rnaseq_full.out $log_dir/
cd $proj_dir
date
echo -e "Performing differential genes analysis.....\n"
echo -e "See progress in run_differential_analysis_rna.out.\n"
message=""
. $img_dir/scripts/run_differential_analysis_rna.sh run=$run_debug padj=$padj time=$time \
	batch_adjust=$batch_adjust &> run_differential_analysis_rna.out
cp $proj_dir/run_differential_analysis_rna.out $log_dir/

# Print messages when entire pipeline ran to completion.
message="Differential RNA-Seq analysis completed successfully.\n"
message="${message}See log run_differential_analysis_rna.out\n\n"
message="${message}Done running RNA-seq full analysis\n\n"
message="${message}Output files are in $work_dir\n"
message="${message}Output files for differential analysis are in\n"
message="${message}$work_dir/diff_analysis_rslt\n"
message="${message}RNA-seq_differential_analysis_report.html in \n"
message="${message}$work_dir/diff_analysis_rslt/\n"
message="${message}for full documentation of differential analysis\n"
message="${message}bw_files folder contains the track files for visualization in IGV.\n"
message="${message}STAR_2pass/Pass2/ contains the aligned bam files, and feature counts for each sample\n"
message="${message}fastqc_rslt contains the Quality Control files.\n"
message="${message}See multiqc_report.html for summary of all QC metrics\n"
message="${message}trim folder contains the trimmed fastq files\n\n"

tmp=$($run sbatch --dependency=afterok:$jid8 \
		$addtl_opt \
		--partition=$general_partition \
		--time=5:00 \
		--output=$log_dir/dummy_run_differential_analysis_rna.txt \
		--mail-type=END \
		--mail-user=$email \
		--job-name=run_rnaseq_full \
		--export message="$message",proj_dir=$proj_dir \
		--wrap "echo -e \"$message\"$(date) >> $proj_dir/run_rnaseq_full.out; \
		cp $proj_dir/run_rnaseq_full.out $log_dir/"| cut -f 4 -d' ')



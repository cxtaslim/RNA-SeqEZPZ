#!/bin/bash

# script to trim fastq files and run quality control
# How to run:
# cd <project_dir>
# bash scripts/run_trim_qc.sh &> run_trim_qc.out &
# DO NOT change the name or location of run_trim_qc.out
# Examples:
# cd project1
# bash scripts/run_trim_qc.sh &> run_trim_qc.out &
#
# or to run with specific time limit:
# bash scripts/run_trim_qc.sh time=DD-HH:MM:SS  &> run_trim_qc.out &
#
# or to run and printing all trace commands (i.e. set -x):
# bash scripts/run_trim_qc.sh run=debug &> run_trim_qc.out &
#


#set -x
#set -e

# clear python path to prevent mix up of python packages
unset PYTHONPATH
# get command line arguments
while [[ "$#" -gt 0 ]]; do
	if [[ $1 == "run"* ]];then
		run=$(echo $1 | cut -d '=' -f 2)
		shift
	fi
	if [[ $1 == "time"* ]];then
		time=$(echo $1 | cut -d '=' -f 2)
		shift
	fi
	if [[ $1 == "ncpus_trim"* ]];then
                ncpus_trim=$(echo $1 | cut -d '=' -f 2)
                shift
        fi

	if [[ $1 == "help" ]];then
		echo ''
		echo 'usage: bash scripts/run_trim_qc.sh [OPTION] &> run_trim_qc.out &'
		echo ''
		echo DESCRIPTION
		echo -e '\trun trim and qc for RNA-seq samples'
		echo ''
		echo OPTIONS
		echo ''
		echo help 
		echo -e '\tdisplay this help and exit'
		echo run=echo
		echo -e "\tdo not run, echo all commands. Default is running all commands"
                echo -e "\tif set to "debug", it will run with "set -x""
		echo -e "time=1-00:00:00"
		echo -e "\tset SLURM time limit time=DD-HH:MM:SS, where ‘DD’ is days, ‘HH’ is hours, etc."
		echo -e "\tDefault is 1 day.\n"
		echo ncpus_trim=4
		echo -e "\tset the number of cpus for trimming"
		echo -e "\tDefault is 4 cpus."
		echo -e "\tSet to 1 if pigz error occurs."
		exit
	fi
done

date
echo -e "\nTimming reads for low quality/adapter sequences then run Quality Control\n"
# converting samples.txt to unix format to remove any invisible extra characters
dos2unix -k samples.txt &> /dev/null

# set default parameters
debug=0
if [[ -z "$run" ]];then
	run=
fi
if [[ -z "$time" ]];then
	time=1-00:00:00
fi
if [[ $run == "debug"* ]];then
        set -x
        run=
        debug=1
fi
if [[ -z "$ncpus_trim" ]];then
        ncpus_trim=4
fi


echo -e "Options used to run:"
echo time="$time"
# automatically change number of cpus for star, bamCompare and bigwigCompare
# if greater than max available cpus
max_cpu=$(lscpu | grep 'CPU(s):' | head -n 1 | awk '{print $2}')
if [[ $max_cpu -lt $ncpus_trim ]]; then
        ncpus_trim=$max_cpu
fi

echo ncpus_trim="$ncpus_trim"
echo ""

# project directory
proj_dir=$(pwd)
cd $proj_dir

# this is where quality control files generated in this script will be stored
work_dir=$proj_dir/outputs
echo "all outputs will be stored in $work_dir"
if [ ! -d $work_dir ]; then
	mkdir $work_dir
fi
log_dir=$work_dir/logs
if [ ! -d $log_dir ]; then
	mkdir $log_dir
fi
echo "See $proj_dir/run_quality_control.out to check the analysis progress"
echo "All other logs and scripts ran will be stored in $log_dir"
echo ""

# singularity image directory
# find based on location of this script
img_dir=$(dirname $(dirname $(readlink -f $0)))

# singularity image name
img_name=rnaseq-pipe-container.sif
# IMPORTANT: It is assumed that:
# scripts to run analysis are in $img_dir/scripts
# reference to run analysis are in $img_dir/ref

echo -e "\nUsing singularity image and scripts in:" ${img_dir} "\n"
REPO=${img_dir}  # set REPO
GIT_SHA=$(git -C "$REPO" rev-parse --short=12 HEAD 2>/dev/null || echo "unknown")
GIT_DESC=$(git -C "$REPO" describe --always --dirty --long --tags --abbrev=12 2>/dev/null || echo "unknown")
echo "git_sha=${GIT_SHA}"
echo "git_version=${GIT_DESC}"
echo ""

# getting SLURM configuration
source $img_dir/scripts/slurm_config_var.sh

# copying scripts ran for records
if [[ ! -d $log_dir/scripts ]];then
	mkdir -p $log_dir/scripts
fi
$(cp $img_dir/scripts/run_trim_qc.sh $log_dir/scripts)
$(cp $img_dir/scripts/trim_fastqc_simg.sbatch $log_dir/scripts/)
$(cp $proj_dir/run_trim_qc.out $log_dir/)
$(cp $img_dir/scripts/multiqc_simg.sbatch $log_dir/scripts)

# getting samples info from samples.txt
$(sed -e 's/[[:space:]]*$//' samples.txt | sed 's/"*$//g' | sed 's/^"*//g' > samples_tmp.txt)
$(mv samples_tmp.txt samples.txt)
groupname_array=($(awk '!/#/ {print $1}' samples.txt))
repname_array=($(awk '!/#/ {print $3}' samples.txt))
email=$(awk '!/#/ {print $5;exit}' samples.txt | tr -d '[:space:]')
path_to_r1_fastq=($(awk '!/#/ {print $6}' samples.txt))
path_to_r2_fastq=($(awk '!/#/ {print $7}' samples.txt))
# note:
# not sure why the usual echo ${filename_string_array[@]} doesn't work
# but this works
# printf "%s\n" "${filename_string_array[@]}"
# trim adapter and low quality bases

# initialize jobid strings
jid=
cd $proj_dir
# added this to prevent error when using symbolic link that is upstream of
# home directory
# NOTE: for this to work, complete path should be used when creating symbolic link
# example:
# this WORKS:
# cd <project-dir>
# ln -s /gpfs0/home1/gdlessnicklab/lab/data/tmp_data_ct/fastq .
# this does NOT work:
# # cd <project-dir>
# ln -s /home1/gdlessnicklab/lab/data/tmp_data_ct/fastq .
# the path should be the path that is returned by 'readlink -f'

target_link=$(readlink -f fastq)

date
# merge fastq files
echo "Merging fastq files"
echo ""

if [[ ! -d $work_dir/merged_fastq ]];then
	mkdir -p $work_dir/merged_fastq
fi
len_row=${#path_to_r1_fastq[@]}
for  (( row = 0; row <= len_row-1; row++ ));do
	echo Merging ${path_to_r1_fastq[$row]//,/ }
	echo and save them as ${work_dir}/merged_fastq/${groupname_array[$row]}_${repname_array[$row]}_r1_fastq.gz
	echo ""
	cat ${path_to_r1_fastq[$row]//,/ } > ${work_dir}/merged_fastq/${groupname_array[$row]}_${repname_array[$row]}_r1_fastq.gz
        echo Merging ${path_to_r2_fastq[$row]//,/ } 
        echo and save them as ${work_dir}/merged_fastq/${groupname_array[$row]}_${repname_array[$row]}_r2_fastq.gz
        echo ""
	cat ${path_to_r2_fastq[$row]//,/ } > ${work_dir}/merged_fastq/${groupname_array[$row]}_${repname_array[$row]}_r2_fastq.gz
done
date

# note: trim_galore is set to --cores 4 which actually translate to 15 cores
# see trim_galore help for more info
for (( idx =0; idx <= len_row-1; idx++ ));do
	echo $idx
	groupname=${groupname_array[$idx]}
	echo $groupname
	repname=${repname_array[$idx]}
	echo $repname
	prefix=${groupname}_${repname}
	echo $prefix
	path_r1=${work_dir}/merged_fastq/${prefix}_r1_fastq.gz
	echo $path_r1
	path_r2=${work_dir}/merged_fastq/${prefix}_r2_fastq.gz
	echo $path_r2
			# echo $prefix $run \
			# $log_dir $email $proj_dir \
			# $groupname $repname \
			# execute inside singularity
			tmp_jid=$(SINGULARITYENV_PYTHONPATH= \
			SINGULARITYENV_run=$run \
			SINGULARITYENV_prefix=${prefix} \
			SINGULARITYENV_path_r1=$path_r1 \
			SINGULARITYENV_path_r2=$path_r2 \
			SINGULARITYENV_ncpus_trim=$ncpus_trim \
				$run sbatch --output=$log_dir/trim_fastqc_${prefix}.out \
					--partition=$general_partition \
					--job-name=trim_fastqc \
					--time=$time \
					--mail-type=FAIL \
					--mail-user=$email \
					--cpus-per-task=$ncpus_trim \
					$addtl_opt \
					--wrap "singularity exec \
					--bind $img_dir/scripts:/scripts \
					--bind $proj_dir:/mnt \
					--bind $path_r1 \
					--bind $path_r2 \
					$img_dir/$img_name \
					/bin/sh /scripts/trim_fastqc_simg.sbatch" | cut -f 4 -d' ')
			echo "Trimming and QCing $path_r1 $path_r2 Job id: $tmp_jid"
			echo "" 
			cp $proj_dir/run_trim_qc.out $log_dir/
			if [ -z "$jid" ]; then
				jid=$tmp_jid
			else
				jid=${jid}:${tmp_jid}
			fi
			#set +x
done

##### run multiqc
jid2=$(SINGULARITYENV_PYTHONPATH= \
	SINGULARITYENV_run=$run \
	SINGULARITYENV_proj_dir=$proj_dir \
	SINGULARITYENV_input_dir=/mnt/outputs/fastqc_rslt \
	$run sbatch --output=$log_dir/multiqc.out \
		--partition=$general_partition \
		--job-name=multiqc \
		--mail-type=FAIL \
		--mail-user=$email \
		--time=$time \
		--dependency=afterok:$jid \
		$addtl_opt \
		--wrap "singularity exec \
			--bind $proj_dir:/mnt \
			--bind $img_dir/scripts:/scripts \
			$img_dir/$img_name \
				/bin/bash /scripts/multiqc_simg.sbatch"| cut -f 4 -d' ')
echo -e "Running multiqc to combine all the quality control files Job id: $jid2\n"

# message if jobs never satisfied
check_jid=$(echo $jid | sed 's/:/,/g')
# check to make sure jobs are completed. Print messages if not.
msg_ok="Done trimming reads and quality control.\n\n"
msg_ok="${msg_ok}Summary of quality control result is in $work_dir/trim/fastqc_rslt/multiqc_report.html\n"
msg_ok="${msg_ok}Trimmed fastq files are in $work_dir/trim\n"
msg_ok="$log_dir/multiqc.out contains the commands ran"
msg_ok="multiqc was completed successfully.\n"
msg_fail="Either trimming or fastqc or multiqc failed. Please check trim_fastqc_* files in $log_dir\n"
msg_fail="Please check multiqc.out in $log_dir\n"
jid_to_check=$check_jid,$jid2
out_file=$proj_dir/run_trim_qc.out
check_trim_jid=$($run sbatch \
        --partition=$general_partition \
        --output=$log_dir/check_trim.out \
        --mail-type=END \
        --mail-user=$email \
        --wait \
        --time=$time \
        --parsable \
        --job-name=check_trim \
        $addtl_opt \
        --wrap "out_file='${out_file}' jid_to_check='${jid_to_check}' msg_ok='${msg_ok}' msg_fail='${msg_fail}' \
	bash $img_dir/scripts/check_job.sh")

cp $proj_dir/run_trim_qc.out $log_dir/


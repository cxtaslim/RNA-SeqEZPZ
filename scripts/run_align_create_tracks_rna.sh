#!/bin/bash

# How to run:
# cd <project_dir>
# bash scripts/run_align_create_tracks_rna.sh &> run_align_create_tracks_rna.out &
# DO NOT change the name or location of run_align_create_tracks_rna.out
# Examples:
# cd project1
# bash scripts/run_align_create_tracks_rna.sh &> run_align_create_tracks_rna.out &
#
# by default, alignment is done to human reference genome hg19 unless specified genome=hg38:
# bash scripts/run_align_create_tracks_rna.sh genome=hg38 &> run_align_create_tracks_rna.out &
#
# or to run with specific time limit:
# bash scripts/run_align_create_tracks_rna.sh time=DD-HH:MM:SS &> run_align_create_tracks_rna.out &
#
# or to run and printing all trace commands (i.e. set -x):
# bash scripts/run_align_create_tracks_rna.sh run=debug &> run_align_create_tracks_rna.out &

#set -x
set -e

# clear python path to avoid reading in user's package sites
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
	if [[ $1 == "ncpus_star"* ]];then
                ncpus_star=$(echo $1 | cut -d '=' -f 2)
                shift
        fi
	if [[ $1 == "help" ]] ;then
		echo ""
		echo 'usage: bash scripts/run_align_create_tracks_rna.sh [OPTION] &> run_align_create_tracks_rna.out &'
		echo ''
		echo DESCRIPTION
		echo -e '\trun alignment and tracks creation'
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
		echo -e "time=1-00:00:00"
		echo -e "\tset SLURM time limit time=DD-HH:MM:SS, where ‘DD’ is days, ‘HH’ is hours, etc."
		echo -e "\tDefault is 1 day\n"
		echo -e "ncpus_star=20"
		echo -e "\tspecifies the number of cpus used for STAR, bamCompare and feature counts."
		echo -e "\tDefault is 20 cpus."
		echo -e "\tthe number of cpus will be automatically adjusted to max number of cpus if ncpus_star"
		echo -e "\tis less than the max available cpus."
		exit
	fi
done

date
# converting samples.txt to unix format to remove any invisible extra characters
dos2unix -k samples.txt &> /dev/null

# set default parameters
# parameter to turn set -x on inside scripts
debug=0
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
if [[ $run == "debug"* ]];then
        set -x
        run=
        debug=1
fi
if [[ -z "$ncpus_star" ]];then
        ncpus_star=20
fi

echo -e "\nAligning reads to $ref_ver and create tracks for visualization\n"
echo -e "Options used to run:"
echo padj="$padj"
echo time="$time"
echo genome="$ref_ver"
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
echo "See $proj_dir/run_align_create_tracks_rna.out to check the analysis progress"
echo "All other logs including scripts ran will be stored in $log_dir"
echo ""

# automatically change number of cpus for star, bamCompare and bigwigCompare
# if greater than max available cpus
max_cpu=$(lscpu | grep 'CPU(s):' | head -n 1 | awk '{print $2}')
if [[ $max_cpu -lt $ncpus_star ]]; then
        ncpus_star=$max_cpu
fi

# singularity image directory
# find based on location of this script
img_dir=$(dirname $(dirname $(readlink -f $0)))

# singularity image name
img_name=rnaseq-pipe-container.sif
# IMPORTANT: It is assumed that:
# scripts to run analysis are in $img_dir/scripts
# reference to run analysis are in $img_dir/ref

echo -e "\nUsing singularity image and scripts in:" ${img_dir} "\n"

# getting SLURM configuration
source $img_dir/scripts/slurm_config_var.sh

# copying scripts ran for records
if [[ ! -d $log_dir/scripts ]];then
	mkdir -p $log_dir/scripts
fi
$(cp $img_dir/scripts/run_align_create_tracks_rna.sh ${log_dir}/scripts/)
$(cp $img_dir/scripts/star_pass1_simg.sbatch ${log_dir}/scripts/)
$(cp $img_dir/scripts/star_pass2_simg.sbatch ${log_dir}/scripts/)
$(cp $img_dir/scripts/combinebw_simg.sbatch ${log_dir}/scripts/)

### specify reference genome
# if ref genome is not in $img_dir/ref, set genome_dir to  project dir
genome_dir=$img_dir/ref/$ref_ver
# set fasta file to ref_fa if exist
if [[ -f $ref_fa ]];then
        fasta_file=$ref_fa
else
    	fasta_file=${genome_dir}/$(find -L $genome_dir -name *.fasta -o -name *.fa | xargs basename)
fi
# set genome_dir and star_index_dir accordingly
genome_dir=$(dirname $fasta_file)
star_index_dir=$genome_dir/STAR_index

# set gtf file to ref_gtf if exist
if [[ -f $ref_gtf ]];then
        gtf_file=$ref_gtf
else
    	gtf_file=${genome_dir}/$(find -L $genome_dir -name *.gtf | xargs basename)
fi
# throws an error if gtf_file or fasta_file doesn't exist
if [[ ! -f $gtf_file || ! -f $fasta_file ]];then
        echo "Please check your genome. Either fasta file or gtf file is not found\n"
	exit 1
fi

# if genome_dir doesn't exist
# set genome_dir to  project dir
# and generate star index and chrom sizes in project_dir
if [[ ! -d $genome_dir ]]; then
                genome_dir=$proj_dir/ref/$ref_ver
                mkdir -p $genome_dir
fi

# if star_index_dir not exist yet
# check if directory writeable if it is create the STAR_index dir
# if not set star_index_dir to proj_dir
if [[ ! -d $star_index_dir ]] ; then
        if [[ -w $(dirname $star_index_dir) ]]; then
                mkdir -p $star_index_dir
        else
	        # star_index will be in proj_dir
                star_index_dir=$proj_dir/ref/$ref_ver/STAR_index
                mkdir -p $star_index_dir
        fi
fi

chr_info_path=$(find -L $genome_dir -name *.chrom.sizes)
if [[ ! -z "$chr_info_path" ]]; then
        chr_info=$(basename $chr_info_path)
else
	echo chrom.sizes file not found. Please run run_star_index.sh to generate it.
	exit 1
fi

### placeholder for calculating effective genome size as suggested in deeptools docs
readlength=$(awk '{sum += $15; n++} END {if (n>0) print int(sum/(n-1));}' \
	 $proj_dir/outputs/fastqc_rslt/multiqc_data/multiqc_fastqc.txt)


# getting samples info from samples.txt
$(sed -e 's/[[:space:]]*$//' samples.txt | sed 's/"*$//g' | sed 's/^"*//g' > samples_tmp.txt)
$(mv samples_tmp.txt samples.txt)
groupname_array=($(awk '!/#/ {print $1}' samples.txt))
repname_array=($(awk '!/#/ {print $3}' samples.txt))
email=$(awk '!/#/ {print $5;exit}' samples.txt | tr -d '[:space:]')
path_to_r1_fastq=($(awk '!/#/ {print $6}' samples.txt))
path_to_r2_fastq=($(awk '!/#/ {print $7}' samples.txt))

if [[ $email == "NA" ]];then
        email=
fi

#### alignment ####
# STAR first pass
# initialize job ids
jid3=
for i in "${!groupname_array[@]}"; do 
	groupname=${groupname_array[$i]}
	repname=${repname_array[$i]}
	prefix=${groupname}_${repname}
	# get the files for the same groupname and replicate
	# (i.e from multiple lanes)
	cd $work_dir/trim
	# Since we don't need to search for filenames anymore
	# I need to change this so I don't require groupname to NOT be a substring of other groups!!!
	#file=($(ls *val_1.fq.gz | \
	#	awk -v groupname=$groupname -v repname=$repname \
	#	'$1 ~ groupname && $1 ~ repname {print $1}'))
	file=${prefix}_val_1.fq.gz
	# combining files from different lanes
	read1=$(printf ",%s" "${file[@]}")
	read1=${read1:1}
	# Since we don't need to search for filenames anymore
        # I need to change this so I don't require groupname to NOT be a substring of other groups!!!
	#file=($(ls *val_2.fq.gz | \
	#	awk -v groupname=$groupname -v repname=$repname \
	#	'$1 ~ groupname && $1 ~ repname {print $1}'))
	file=${prefix}_val_2.fq.gz
	read2=$(printf ",%s" "${file[@]}")
	read2=${read2:1}
	cd $proj_dir
	#echo $read1
	#echo $read2
	cd $proj_dir
	# set -x
	tmp_jid=$(SINGULARITYENV_PYTHONPATH= \
		SINGULARITYENV_read1=$read1 \
		SINGULARITYENV_read2=$read2 \
		SINGULARITYENV_run=$run \
		SINGULARITYENV_ncpus=$ncpus_star \
		SINGULARITYENV_prefix=$prefix \
		SINGULARITYENV_ref_ver=$ref_ver \
			$run sbatch --output=$log_dir/star_pass1_${prefix}.out \
			--cpus-per-task $ncpus_star \
			--partition=$high_mem_partition \
			--mail-type=FAIL \
			--mail-user=$email \
			--job-name=star_pass1 \
			--time=$time \
			$addtl_opt \
			--wrap "singularity exec \
				--bind $proj_dir:/mnt \
				--bind $img_dir/scripts:/scripts \
				--bind $star_index_dir:/star_index_dir \
				$img_dir/$img_name \
				/bin/bash /scripts/star_pass1_simg.sbatch"| cut -f 4 -d' ')
	echo "Running STAR first pass for $prefix job id: $tmp_jid"
	echo ""
	if [ -z "$jid3" ]; then
		jid3=$tmp_jid
	else
		jid3=${jid3}:${tmp_jid}
	fi
#set +x
done

# dummy job to wait for STAR first pass to finish before
# continuing on
tmp=$($run sbatch \
	--dependency=afterok:$jid3 \
	--output=$log_dir/dummy.txt \
	--job-name=after_star_pass1 \
	--time=5:00 \
	$addtl_opt \
	--wrap "echo dummy sbatch to wait for STAR 1st pass to be finished"| cut -f 4 -d' ')

# message if jobs never satisfied or cancelled
check_jid3=$(echo $jid3 | sed 's/:/,/g')
# check to make sure jobs are completed. Print messages.
msg_ok="STAR 1st pass completed successfully.\n"
msg_fail="STAR 1st pass run failed. Please check star_pass1_* files in $log_dir\n"
jid_to_check=$check_jid3,$tmp
out_file=$proj_dir/run_align_create_tracks_rna.out
check_run_star_pass1_jid=$($run sbatch \
        --partition=$general_partition \
        --output=$log_dir/check_run_star_pass1.out \
        --mail-type=END \
        --mail-user=$email \
        --wait \
        --time=$time \
        --parsable \
        --job-name=check_star_pass1 \
        --export=out_file="$out_file",jid_to_check="$jid_to_check",msg_ok="$msg_ok",msg_fail="$msg_fail" \
        $addtl_opt \
	--wrap "bash $img_dir/scripts/check_job.sh")

# STAR second pass
genomeForPass2=$work_dir/STAR_2pass/GenomeForPass2
mkdir -p $genomeForPass2
# get all the junctions files for input to 2nd pass alignment
cd $genomeForPass2
sjFiles=($(find `pwd` -name '*_SJ.out.tab.Pass1.sjdb'))
sjFiles=$(printf " %s" "${sjFiles[@]}")
sjFiles=${sjFiles:1}
sjFiles=${sjFiles//$genomeForPass2/..\/STAR_2pass\/GenomeForPass2}

cd $proj_dir
# initialize job ids
jid4=
for i in "${!groupname_array[@]}"; do
	groupname=${groupname_array[$i]}
	repname=${repname_array[$i]}
	prefix=${groupname}_${repname}
	# get the files for the same groupname and replicate
	# (i.e from multiple lanes)
	cd $work_dir/trim
	# Since we don't need to search for filenames anymore
        # I need to change this so I don't require groupname to NOT be a substring of other groups!!!
	#file=($(ls *val_1.fq.gz | \
	#	awk -v groupname=$groupname -v repname=$repname \
	#	'$1 ~ groupname && $1 ~ repname {print $1}'))
	file=${prefix}_val_1.fq.gz
	read1=$(printf ",%s" "${file[@]}")
	read1=${read1:1}
	# Since we don't need to search for filenames anymore
        # I need to change this so I don't require groupname to NOT be a substring of other groups!!!
	#file=($(ls *val_2.fq.gz | \
	#	awk -v groupname=$groupname -v repname=$repname \
	#	'$1 ~ groupname && $1 ~ repname {print $1}'))
	file=${prefix}_val_2.fq.gz
	read2=$(printf ",%s" "${file[@]}")
	read2=${read2:1}
	cd $proj_dir
	#echo $read1
	#echo $read2
	# set -x
	tmp_jid=$(SINGULARITYENV_PYTHONPATH= \
		SINGULARITYENV_read1=$read1 \
		SINGULARITYENV_read2=$read2 \
		SINGULARITYENV_run=$run \
		SINGULARITYENV_ncpus=$ncpus_star \
		SINGULARITYENV_sjFiles=$sjFiles \
		SINGULARITYENV_prefix=$prefix \
		SINGULARITYENV_ref_ver=$ref_ver \
		SINGULARITYENV_fasta_file=$fasta_file \
		SINGULARITYENV_gtf_file=$gtf_file \
		SINGULARITYENV_readlength=$readlength \
			$run sbatch --output=$log_dir/star_pass2_${prefix}.out \
				--cpus-per-task $ncpus_star \
				--partition=$high_mem_partition \
				--mail-type=FAIL \
				--dependency=afterok:$tmp \
				--mail-user=$email \
				--job-name=star_pass2 \
				--mem=${high_mem}G \
				--time=$time \
				$addtl_opt \
				--wrap "singularity exec \
					--bind $proj_dir:/mnt \
					--bind $img_dir/scripts:/scripts \
					--bind $star_index_dir:/star_index_dir \
					--bind $fasta_file \
					$img_dir/$img_name \
					/bin/bash /scripts/star_pass2_simg.sbatch"| cut -f 4 -d' ')
	echo "Running STAR second pass for $prefix job id: $tmp_jid"
	echo ""
	if [ -z "$jid4" ]; then
		jid4=$tmp_jid
	else
		jid4=${jid4}:${tmp_jid}
	fi
	#set +x
done

# if there are replicates, combine the bw files for each replicate
n_rep=$(printf "%s\n" "${repname_array[@]}" | sort -u | wc -l)
echo number of replicates $n_rep
if [[ n_rep -gt 1 ]]; then
	#echo There is replicates
	#set -x
	cd $work_dir
	# dummy sbatch waiting for previous jobs to finish
	tmp=$($run sbatch \
	--partition=$general_partition \
	--dependency=afterok:$jid4 \
	--output=$log_dir/dummy.txt \
	--time=5:00 \
	$addtl_opt \
	--wrap "echo dummy sbatch waiting for STAR 2nd pass to be finished"| cut -f 4 -d' ')

	#message if jobs never satisfied or cancelled
	#cancel jobs if dependency never satisfied
	check_jid4=$(echo $jid4 | sed 's/:/,/g')
	# check to make sure jobs are completed. Print messages if not.
	msg_ok="STAR 2nd pass completed successfully.\n"
	msg_fail="STAR 2nd pass run failed. Please check star_pass2_* files in $log_dir\n"
	jid_to_check=$check_jid4,$tmp
	out_file=$proj_dir/run_align_create_tracks_rna.out
	check_star_pass2_jid=$($run sbatch \
        	--partition=$general_partition \
        	--output=$log_dir/check_star_pass2.out \
        	--mail-type=END \
        	--mail-user=$email \
        	--wait \
        	--time=$time \
        	--parsable \
        	--job-name=check_star_pass2 \
        	--export=out_file="$out_file",jid_to_check="$jid_to_check",msg_ok="$msg_ok",msg_fail="$msg_fail" \
        	$addtl_opt \
		--wrap "bash $img_dir/scripts/check_job.sh")

	# initialize job ids
	jid4b=
	cd $work_dir/bw_files
	for file in *_${repname_array[0]}_norm.bw; do
		# set -x
		# combined the replicates of bigwig into one bigwig
		# echo $file
		rep_files_array=($(ls ${file/_${repname_array[0]}_norm.bw/*_norm.bw}))
		all_files_rep_cmd=$(printf " %s" "${rep_files_array[@]}")
		all_files_rep_cmd=${all_files_rep_cmd:1}
		echo $all_files_rep_cmd
		# find groupname
		cd $proj_dir
		#echo $chr_info
		groupname=$(awk -v file=$file 'file ~ $1 {print $1;exit}' samples.txt)
		tmp_jid=$(SINGULARITYENV_PYTHONPATH= \
		SINGULARITYENV_run=$run \
				SINGULARITYENV_all_files_rep_cmd=$all_files_rep_cmd \
				SINGULARITYENV_ncpus=$ncpus_star \
				SINGULARITYENV_groupname=$groupname \
				SINGULARITYENV_ref_ver=$ref_ver \
				SINGULARITYENV_fasta_file=$fasta_file \
				SINGULARITYENV_chr_info=$chr_info \
				$run sbatch --output=$log_dir/combinebw_${groupname}.out \
					--dependency=afterok:$tmp \
					--cpus-per-task $ncpus_star \
					--partition=$high_mem_partition \
					--mail-type=FAIL \
					--mail-user=$email \
					--job-name=combinebw \
					--time=$time \
					$addtl_opt \
					--wrap "singularity exec \
						--bind $proj_dir:/mnt \
						--bind $genome_dir:/ref \
						--bind $img_dir/scripts:/scripts \
						$img_dir/$img_name \
							/bin/bash /scripts/combinebw_simg.sbatch"| cut -f 4 -d' ')
		echo "Combining $out_prefix bigwig replicates job id: $tmp_jid"
		echo ""
		cd $work_dir/bw_files
		if [ -z "$jid4b" ]; then
			jid4b=$tmp_jid
		else
			jid4b=${jid4b}:${tmp_jid}
		fi
		# set +x
	done
	jid4c=$($run sbatch \
		--partition=$general_partition \
		--dependency=afterok:$jid4b \
		--output=$log_dir/dummy.txt \
		--mail-type=END \
		--mail-user=$email \
		--time=5:00 \
		--job-name=run_align_create_tracks_rna \
		$addtl_opt \
		--wrap "echo dummy sbatch after combinebw is finished"| cut -f 4 -d' ')
	# message if jobs never satisfied or canceled
	check_jid4b=$(echo $jid4b | sed 's/:/,/g')
	# check to make sure jobs are completed. Print messages if not.
	msg_ok="bw files were combined successfully.\n"
	msg_fail="Combining bw files failed. Please check combinebw_* files in $log_dir\n"
	jid_to_check=$check_jid4b,$jid4c
	out_file=$proj_dir/run_align_create_tracks_rna.out
	check_combine_jid=$($run sbatch \
        --partition=$general_partition \
        --output=$log_dir/check_combine.out \
        --mail-type=END \
        --mail-user=$email \
        --wait \
        --time=$time \
        --parsable \
        --job-name=check_combine \
        --export=out_file="$out_file",jid_to_check="$jid_to_check",msg_ok="$msg_ok",msg_fail="$msg_fail" \
        $addtl_opt \
	--wrap "bash $img_dir/scripts/check_job.sh")
else
	echo There are no replicates
	jid4c=$($run sbatch \
		--dependency=afterok:$jid4 \
		--partition=$general_partition \
		--output=$log_dir/dummy.txt \
		--mail-type=END \
		--mail-user=$email \
		--time=5:00 \
		--job-name=run_align_create_tracks_rna \
		$addtl_opt \
		--wrap "echo dummy sbatch after STAR 2nd pass is finished. There are no replicates"| cut -f 4 -d' ')
	# message if jobs never satisfied or cancelled
	check_jid4=$(echo $jid4 | sed 's/:/,/g')
	# check to make sure jobs are completed. Print messages if not.
	msg_ok="STAR 2nd run completed successfully.\n"
	msg_fail="STAR 2nd pass run failed. Please check star_pass2_* files in $log_dir\n"
	jid_to_check=$check_jid4,$jid4c
	out_file=$proj_dir/run_align_create_tracks_rna.out
	check_star_pass2=$($run sbatch \
        --partition=$general_partition \
        --output=$log_dir/check_star_pass2.out \
        --mail-type=END \
        --mail-user=$email \
        --wait \
        --time=$time \
        --parsable \
        --job-name=check_star_pass2 \
        --export=out_file="$out_file",jid_to_check="$jid_to_check",msg_ok="$msg_ok",msg_fail="$msg_fail" \
        $addtl_opt \
	--wrap "bash $img_dir/scripts/check_job.sh")
fi

message="Done alignment using STAR 2-pass approach and created bw files for visualization\n\n\
Tracks for each replicate and combined replicates are \"*norm.bw\" and \"comb.bw\" respectively\n\
in $work_dir\n\
Aligned bam files are in $work_dir/STAR_2pass/Pass2\n\n"

# remove STAR pass1 directory
rm -r $work_dir/STAR_2pass/Pass1
rm -r $work_dir/STAR_2pass/GenomeForPass2

tmp=$($run sbatch --dependency=afterok:$jid4c \
		--partition=$general_partition \
		--output=$log_dir/dummy.txt \
		--mail-type=END \
		--mail-user=$email \
		--time=5:00 \
		--job-name=run_align_create_tracks_rna \
		--export message="$message",proj_dir=$proj_dir \
		$addtl_opt \
		--wrap "echo -e \"$message\"$(date) >> $proj_dir/run_align_create_tracks_rna.out; \
			cp $proj_dir/run_align_create_tracks_rna.out $log_dir/run_align_create_tracks_rna.out"| \
			cut -f 4 -d' ')


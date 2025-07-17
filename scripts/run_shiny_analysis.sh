#!/bin/bash
#set -x
#set -e
# script to run shiny app in the beginning
# How to run
# cd <project_dir>
# bash scripts/run_shiny_analysis.sh
# Examples:
# bash scripts/run_shiny_analysis.sh
# <project_dir> is where the fastq and outputs directory are
# project directory where the fastq and outputs directory are
#
# to run with all commands printed
# bash scripts/run_shiny_analysis.sh run=debug
#
# Default time limit is 1 day to change to 1 day, 2 hours, 20 minutes and 30 sec
# bash scripts/run_shiny_analysis.sh time=1-2:20:30

# clear python path to prevent mixed up of python packages
date
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
	if [[ $1 == "filepath"* ]];then
                filepath=$(echo $1 | cut -d '=' -f 2)
                shift
        fi
        if [[ $1 == "max_nsamples"* ]];then
                max_nsamples=$(echo $1 | cut -d '=' -f 2)
                shift
        fi

        if [[ $1 == "help" ]];then
		echo ""
                echo 'usage: bash scripts/run_shiny_analysis.sh [OPTION]'
                echo ''
                echo DESCRIPTION
                echo -e '\trun shiny app to run full RNA-seq analysis'
                echo ''
                echo OPTIONS
                echo ''
		echo filepath=path/to/fastq/fa/gtf/files
		echo -e '\tpath to navigate to fastq,reference fasta, GTF files and project folder'
		echo -e '\tthis path will be bound and used inside singularity and shiny app to navigate to these files.'
		echo -e '\tif not provided, host root will be used'
                echo max_nsamples
                echo -e '\tset the maximum number of samples. Default is 50'
                echo -e '\tset max_nsamples if you have more than 50 samples.'
                echo -e '\t 50 samples correspond to number of rows in samples.txt'
                echo -e '\t NOT number of fastq files!'
                echo -e '\t note: that max_nsamples only needed when running shiny app to run analysis'
                echo run=echo
                echo -e "\tdo not run, echo all commands. Default is running all commands"
                echo -e "\tif set to "debug", it will run with "set -x""
		echo -e ""
                echo -e "time=1-00:00:00"
                echo -e "\tset SLURM time limit time=DD-HH:MM:SS, where ‘DD’ is days, ‘HH’ is hours, etc."
                echo -e ""
                echo help
                echo -e '\tdisplay this help and exit'
                exit
        fi
done

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
if [[ -z "$max_nsamples" ]];then
        max_nsamples=50
fi

proj_dir=$(pwd)
cd $proj_dir

# create a named pipe
# supress message if pipe exist
mkfifo $proj_dir/mypipe 2> /dev/null

work_dir=$proj_dir/outputs

# singularity image directory
# find based on location of this script
img_dir=$(dirname $(dirname $(readlink -f $0)))

# singularity image name
img_name=rnaseq-pipe-container.sif
# IMPORTANT: It is assumed that:
# scripts to run analysis are in $img_dir/scripts
# reference to run analysis are in $img_dir/ref

echo -e "\nUsing singularity image and scripts in:" ${img_dir} "\n"

# getting Nextflow configuration
#skip load_nextflow line
source <(grep -v "load_nextflow" scripts/nextflow_config_var.config )

echo -e "Options used to run:"
echo time="$time"
echo ""

# find available port
FROM=1035
TO=65535
port_num=$(singularity exec $img_dir/$img_name comm -23 \
<(seq "$FROM" "$TO" | sort) \
<(ss -Htan | awk '{print $4}' | cut -d':' -f2 | sort -u) \
| shuf | head -n 1)

# if filepath not specified then bind host root
if [[ -z "$filepath" ]];then
	bind_filepath="--bind /:/root"
	filepath=/
else
        # get the long path
        filepath=$(readlink -f $filepath)
	bind_filepath="--bind $filepath:/filepath"
fi

# if filepath doesn't exist, exit with error
if [[ ! -e "$filepath" ]]; then
        echo -e "\n\n Error: check your filepath.\n\n"
        exit 1
fi

# Activate environment where shiny is installed and go to "app.R" directory
export port_num filepath max_nsamples img_dir proj_dir img_name bind_filepath
jid=$(sbatch --time=$time \
	--partition=$general_partition \
	$addtl_opt \
	--export=port_num,filepath,max_nsamples,img_dir,proj_dir,img_name,bind_filepath \
	--output=run_shiny_analysis.out \
	--wrap "/bin/sh $img_dir/scripts/run_listen_app.sh"| cut -f 4 -d' ')

echo -e "\n\nYou need to have x11 display server such as Xming running.\n"

#status=$(squeue -j $jid -o "%t" -h)
#node=$(squeue -j $jid -o "%R" -h)
status=$(sacct -j $jid -Xn -Po state)
node=$(sacct -j $jid -Xn -Po nodelist)
while [[ $status != "RUNNING"* ]];do
       #status=$(squeue -j $jid -o "%t" -h)
       #node=$(squeue -j $jid -o "%R" -h)
       status=$(sacct -j $jid -Xn -Po state)
       node=$(sacct -j $jid -Xn -Po nodelist)
	echo -e "Please wait...\n"
	sleep 5
done

# Extra check to make sure shiny app is already at listening point
# Put this in because sometimes it takes time to load libraries
while true; do
	# get the last line of run_shiny_analysis.out
       if [ -f run_shiny_analysis.out ]; then
	  last_line=$(tail -n 1 run_shiny_analysis.out)
	# exit if shiny fail to load
	if grep -q "Execution halted" "run_shiny_analysis.out"; then exit 1; fi
	# check if it contains listening
	if [[ $last_line == *"Listening"* ]]; then
		break
	fi
	sleep 5
	else
	last_line=""
	fi
done

# adding configuration to send signal every four minutes (240 secs)
echo -e "Host *\n ServerAliveInterval 240" >> ~/.ssh/config

echo Please ignore \"Failed to open connection..\" message.
echo ""
echo If you received \"connection to $node closed\",
echo "please re-run the command line (i.e. run_shiny_analysis.sh)".
echo ""
echo -e "After typing in your password, please wait until firefox appears ....\n"
echo -e "It may take some time for firefox to be responsive as it tries to load all the data\n"
sleep 5
ssh -tX "$node" 'export port_num='"'$port_num'"'; 
        export img_dir='"'$img_dir'"'; \
        export img_name='"'$img_name'"'; \
	export proj_dir='"'$proj_dir'"'; \
	export run=$run; \
        bash $img_dir/scripts/run_firefox.sh'

## this should not be run until firefox is closed
scancel $jid
# deleting the last lines which is server alive text
head -n -1 ~/.ssh/config > temp && cp temp ~/.ssh/config && rm temp
# make permission stricter
chmod 600 ~/.ssh/config

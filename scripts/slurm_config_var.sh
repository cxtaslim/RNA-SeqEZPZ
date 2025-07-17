#!/usr/bin/env bash
# name of partition with high memory
high_mem_partition=himem
# name of  partition for general use
general_partition=general
# memory in GB for jobs that require high memory like STAR
high_mem=128
# memory in GB for jobs that require medium memory like feature counts
med_mem=32
# memory in GB for jobs that require very high memory like DESeq2 and SARTools
very_high_mem=150
# any additional sbatch options needed to run i.e. in OSC this is your -A <project-num>
addtl_opt=""


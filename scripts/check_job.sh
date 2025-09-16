#!/bin/bash

set -x

# script to print message and cancel if jobs never satisfied or cancelled
# change : to ;
jids=${jid_to_check//:/,}
# states of jobs that's not completed in sacct
state_jids=($(sacct -j $jids --format=JobID,State --noheader | grep -v "COMPLETED"))

# while jobs is not completed yet keep monitoring
while [ ${#state_jids[@]} -ne 0 ];do
        # get pending or running job ids
	jids_pending_arr=($(sacct -j $jids -Xn -Po jobid,state | \
		grep -E "PENDING|RUNNING" | awk -F "|" '{print $1}'))
	jids_pending=$(IFS=','; echo "${jids_pending_arr[*]}") 
	# if nothing is pending exit early
	if [[ -z "$jids_pending" ]]; then
		# jobs completed exit early
		break
	fi
        # only squeue pending to avoid error
	reason=($(squeue -j $jids_pending -o "%R" -h))
	state=($(sacct -j $jids --format=state --noheader ))
	# store jobname and jobids
	jobname=($(sacct -j $jids --format=jobname%100 --noheader ))
	jobid=($(sacct -j $jids --format=jobid --noheader ))
	# unique states of jids
	state_uni=($(echo "${state[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
        if [[ "${reason[*]}" == *"DependencyNeverSatisfied"* || "${state[*]}" == *"CANCELLED"* || \
                "${state[*]}" == *"FAILED"* || "${reason[*]}" == *"FAILED"* ]]; then
		# if dependency that never satisfied or if they are cancelled or failed
        	# cancel all queued jobs.
		scancel "$jids_pending"
	elif [[ "${#state_uni[@]}" -eq 1 ]] && [[ "${state_uni[0]}" == *"COMPLETED"* ]] ;then
                # jobs completed exit early
		break
        else
                sleep 10
        fi
	# update states with pending jobs
	state_jids=($(sacct -j $jids_pending --format=jobid,jobname%100,State --noheader | \
		grep -v -E "COMPLETED|CANCELLED"))
	echo ${state_jids[*]}
done

# jobs out of queue check and print error message 
state=($(sacct -j $jids --format=state --noheader ))
# store jobname and jobids
jobname=($(sacct -j $jids --format=jobname%100 --noheader ))
jobid=($(sacct -j $jids --format=jobid --noheader ))
# unique states of jids
state_uni=($(echo "${state[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
# if jobs failed for the following reasons, cancel all jobs in queue
# and write message to run_rnaseq_full.out if file exists
if [[ "${reason[*]}" == *"DependencyNeverSatisfied"* || "${state[*]}" == *"CANCELLED"* || \
	"${state[*]}" == *"FAILED"* || "${reason[*]}" == *"FAILED"* ]]; then
        for i in "${!state[@]}"; do
        	if [[ "${state[$i]}" = *"FAILED"* ]]; then
                	fail_idx=$i;
                        break
                fi
       done
	scancel "$jids"
	echo -e "$msg_fail" >> "${out_file}"
        fail_job_info="${jobname[$fail_idx]} jobid:${jobid[$fail_idx]} failed.\n"
        fail_job_info="${fail_job_info}See ${jobname[$fail_idx]}.out for more details.\n\n"
        echo -e "$fail_job_info" >> "${out_file}"
        if [[ -f "$(dirname $out_file)/run_rnaseq_full.out" ]];then
        	echo -e "$fail_job_info" >> "$(dirname $out_file)/run_rnaseq_full.out"
        fi
        cp "${out_file}" "$(dirname $out_file)/outputs/logs/$(basename $out_file)"
        date >> "${out_file}"
        exit 1
fi

# jobs should be successfully completed 
# print ok message
echo -e "$msg_ok" >> "${out_file}"
date >> "${out_file}"
echo -e "\n" >> "${out_file}"
cp "${out_file}" "$(dirname $out_file)/outputs/logs/$(basename $out_file)"

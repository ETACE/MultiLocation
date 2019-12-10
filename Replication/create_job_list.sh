#!/bin/bash
#Script to generate a set of run scripts (job_list_*.sh)
#Each script runs a set of jobs in sequence
#Jobs are assigned to the lists by a round-robin shuffle algorithm
 
echo '  Starting creation of job scripts in:'$BASE

#Cycle through the number of procs
let count=0

#Remove old scripts
rm -f job_list_*.sh


    for folder2 in $F2; do
        for folder3 in $F3; do
		for folder4 in $F4; do
			

        cd ./$F1/$folder2/$folder3/$folder4

        for run in $RUNS; do
            cd 'run_'$run

		if [ ! -f iters.db ]; then 

		#Add a line for this run
       		echo 'cd '$PWD>>$BASE/job_list_$count.sh
       		echo 'bash ./run.sh'>>$BASE/job_list_$count.sh
		echo 'cd -'>>$BASE/job_list_$count.sh
		echo 'Assigned job to list '$count

		let "count = $count + 1"
		if [ $count -eq $NUM_PROCS ]; then let count=0; fi


		  fi

           	 cd ..
		done
	
	




		
       
       		 cd $BASE

	done
	 
        done
  done

echo '  Finished creation of job scripts.'


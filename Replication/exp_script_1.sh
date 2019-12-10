#!/bin/bash

######### CREATION OF EXPERIMENT FOLDERS ###############################################
# Script to create the folder hierarchy
# - Create folders ./duration_x/intensity_y/frequency_z
# - Create folders for multiple batch runs: its_"batch"
# - Copy output.xml and run.sh from top-level to these run folders
########################################################################################
BASE=$PWD

echo '  Creating folder hierarchy...'
echo 'In exp_script_1.sh these values are used:'
echo 'Skill upgrade: [' $F1 ']'
echo 'gamma: [' $F2 ']'
echo 'delta: [' $F3 ']'

rm -f STATUS
#mkdir -p 'its'

for folder1 in $F1; do
   mkdir -p $folder1
    echo '    Created folder:' $folder1
    cd $folder1
    for folder2 in $F2; do
       mkdir -p $folder2
        echo '      Created folder:' $folder1/$folder2
        cd $folder2
        for folder3 in $F3; do
            mkdir -p $folder3
            echo '        Created folder:' $folder1/$folder2/$folder3
	    cd $folder3
           	 #echo $PWD
		for folder4 in $F4; do
           	 mkdir -p $folder4
           	 echo '        Created folder:' $folder1/$folder2/$folder3/$folder4
         	 cd $folder4
			
		    		for run in $RUNS; do
		    		   mkdir -p 'run_'$run
		    		   # echo '          Created folder:' $folder1/$folder2/$folder3/'run_'$run
		     	   	echo $folder1'/'$folder2'/'$folder3'/'$folder4'/run_'$run':CREATED'>> $BASE/STATUS
		      	 	# cp $BASE/run.sh ./'run_'$run
				
		done
		cd ..
          
        done
        cd ..
	done
        cd ..
    done
    cd ..
done


cd $BASE

echo '  Finished creation of folder hierarchy.'

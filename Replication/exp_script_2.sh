#!/bin/bash

######### CREATION OF THE SPECIFIC SETTINGS #########################################################
# Script to create the file specific.xml that contains the specific parameter settings for the cases
# - Enter folder ./duration_x/intensity_y/frequency_z
# - Write xml tags and parameter settings to the file specific.xml
#####################################################################################################
BASE=$PWD

echo '  Creating specific.xml files in folder hierarchy...'
echo 'In exp_script_2.sh these values are used:'
echo 'Skill upgrade: [' $F1_values ']'
echo 'Gamma: [' $F2_values ']'
echo 'Delta: [' $F3_values ']'

for f1 in $F1; do
    folder1=$f1
    for f2 in $F2_values; do
        folder2=$F2Name'_'$f2
        for f3 in $F3_values; do
		folder3=$F3Name'_'$f3
		for f4 in $F4_values; do
			folder4=$F4Name'_'$f4


			
           
                echo $PWD

                cd ./$folder1/$folder2/$folder3/$folder4/
               

		 		for run in $RUNS; do

				   cd  'run_'$run
		    		   
		    		 rm -f run.sh

				 echo '#!/bin/bash' >>run.sh

				      chmod u+x run.sh

				    
				
				    echo 'java -jar -Xmx64m ' $BASE'/BatchParallelModel.jar' $ITS $F2Name $f2 $F3Name $f3 $F4Name $f4 $otherParameters $parameters $PWD '>>/dev/null'  >>run.sh

				   #echo 'java -jar' $BASE'/BatchParallelModel.jar' $ITS $F2Name $f2 $F3Name $f3 $F4Name $f4 $otherParameters $parameters
			#	echo 'java -jar -Xmx64m ' $BASE'/BatchParallelAlternativeModel.jar' $ITS $F2Name $f2 $F3Name $f3 $F4Name $f4 $otherParameters $parameters $PWD '>>/dev/null'  >>run.sh

				
		     	   	cd ..
		      	 	
				
		done
		cd ..

	
                
                cd $BASE

		done

		done
        done
done

echo '  Finished creation of specific.xml files.'

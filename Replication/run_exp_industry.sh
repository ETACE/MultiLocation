#!/bin/bash

################################################################################################################
# This script is for setting up and running the experiment:
# Steps:
# - 1 : experiment settings
# - 2 : creation of experiment folder hierarchy
# - 3 : creation of the specific settings xml file
# - 4 : running the initial transient phase; stores a snapshot of the last iteration in transient.xml
# - 5 : running the experiment
################################################################################################################

echo 'Starting top-level experiment script...'

######### STEP 1: EXPERIMENT SETTINGS ##################################################################

#Set the base folder
export BASE=$PWD

echo $BASE

cd $BASE

#Iterations
export ITS=1000


#Set number of job processes to use
export NUM_PROCS=1



#Set number of batch runs
export TOTAL_RUNS=200

RUNS=''
for ((j=1; j<=TOTAL_RUNS; j++)); do
    export RUNS=$RUNS' '$j
done

 
echo 'Batch runs $RUNS '

export F1="its_industry"


export F2Name="fractionEffectivityImitators"
export F2="fractionEffectivityImitators_0.05 "
export F2_values="0.05"


export F3Name="rateLocationDecision"
export F3="rateLocationDecision_0.02"
export F3_values="0.02" 

export F4Name="locationCosts"
export F4="locationCosts_0.01"
export F4_values="0.01"


# Here one can set other parameters
export otherParameters="marketEntryExit true  strategyExperiment false  enteringCosts 3.0  marketEntryNu 3.0  fractionEffectivityInnovators 0.05"





######### STEP 1: CREATION OF EXPERIMENT FOLDER HIERARCHY 
bash ./exp_script_1.sh

######### STEP 2: CREATION OF THE SPECIFIC SETTINGS XML FILE 
bash ./exp_script_2.sh

wait
######### STEP 3: CREATING  JOB SCRIPTS 
bash ./create_job_list.sh



######### STEP 4: LAUNCHING  JOB SCRIPTS 
bash ./launch_job_list.sh

#bash ./clean_up.sh
:
wait

echo 'Finished top-level experiment script.'

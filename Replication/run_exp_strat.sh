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


#export F1="Scenario_No_Gap__HC_Tech_Gap Scenario_No_Gap__Tech_Gap Scenario_No_Gap__HC_Gap"
#export F2="region_cost_0.0 region_cost_0.1 region_cost_0.25 region_cost_0.5 region_cost_1.0 region_cost_10000.0"

export F1="its_strategy"




export F4Name="enteringCosts"
export F4="enteringCosts_3.0"
export F4_values="3.0"



export F3Name="strategyParameter"
export F3="strategyParameter_-0.15 strategyParameter_-0.125 strategyParameter_-0.1 strategyParameter_-0.075 strategyParameter_-0.05 strategyParameter_-0.025 strategyParameter_0.0 strategyParameter_0.025 strategyParameter_0.05 strategyParameter_0.075 strategyParameter_0.1 strategyParameter_0.125 strategyParameter_0.15 strategyParameter_0.175 strategyParameter_0.2 strategyParameter_0.225 strategyParameter_0.25 strategyParameter_0.275 strategyParameter_0.3"
export F3_values="-0.15 -0.125 -0.1 -0.075 -0.05 -0.025 0.0 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3"  





export F2Name="industryScenario"
export F2="industryScenario_1 industryScenario_2 industryScenario_3 industryScenario_4"
export F2_values="1 2 3 4"



export otherParameters="marketEntryExit false  strategyExperiment true  marketEntryNu 3.0  fractionEffectivityInnovators 0.05 fractionEffectivityImitators 0.05 rateLocationDecision 0.02 locationCosts 0.01"





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

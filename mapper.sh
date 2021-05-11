#!/bin/bash

# Initialize Modules
source /etc/profile

# Load Anaconda Module
module load anaconda/2020a

# Call your script as you would from the command line, passing in $1 and $2 as arugments
# Note that $1 and $2 are the arguments passed into this script
#python top5overall_map.py $1 $2
python data-import-map.py $1 $2
# LLMapReduce --mapper mapper.sh --input ./data/batch_input/ --output ./data/output --apptype=mimo --np=[1,4,1]

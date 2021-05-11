#!/bin/bash

# Loading required module and activate environment
source /etc/profile
module load anaconda/2020a
conda activate bio

# run data preprocess script
python data-import.py

#!/bin/bash

# Loading required module and activate environment
source /etc/profile

while getopts f: flag
do
	case "${flag}" in 
		f) infile=${OPTARG};;
	esac
done

# run data preprocess script
# python data-import-fasta.py
#python data-import-muscle.py $infile
python data-import-mpi.py $infile

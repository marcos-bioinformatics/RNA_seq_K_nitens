#! /bin/bash

## This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International$
## Author: Marcos Elizalde Horcada
## Date: February 2020
## Contact: Francisco José Romero Campero - fran@us.es

# This script needs to be launched with knitens.params as a parameter

# Specify the file containing loop parameters.
PARAM_FILE=$1

# Parameter reading
WD=$(grep working_directory: ${PARAM_FILE} | awk '{ print $2 }')  #el archivo param files ya se especifica como parámetro
SD=$(grep script_directory: ${PARAM_FILE} | awk '{ print $2 }' )

echo This is my working directory $WD.

# Generation of sample's folder
NS=$(grep number_samples: ${PARAM_FILE} | awk '{ print $2 }')

# Access working directory
cd $WD/samples

# Create the appropiate number of folders
for i in $(seq 1 1 $NS)
do
	mkdir sample_$i
done

# Access the scripts folder
cd $SD

# Loop through all samples and execute sample processing script
for i in $(seq 1 1 $NS)
do 
	ACC_NUM=$(grep accession_number_sample_$i: ${PARAM_FILE} | awk '{ print $2 }')
	qsub -N sample_$i -o sample_$i sample_processing.sh $WD/samples/sample_$i ${ACC_NUM} $i
	sleep 10m
done



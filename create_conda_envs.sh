#!/bin/bash -e

#Create a conda_envs directory if not existing
if ! [[ -d ./conda_envs ]]; then
	mkdir ./conda_envs
fi


# Define a list of conda envs to create: --> conda_env[env_name]="Tools to install in the env"
declare -A conda_env
conda_env[env-preprocessing]="trim-galore fastp trimmomatic multiqc"
conda_env[env-alignment]="bwa bowtie2 samtools"


# Iterate through each of the envs, create it if not already present in the conda_envs directory

for env_name in ${!conda_env[@]}; do
#	echo $env_name
	if ! [[ -d conda_envs/${env_name} ]]; then
		echo  "Creating ${env_name}..."
		conda_flags="-c bioconda -c conda-forge -c defaults --mkdir --yes --quiet"
		conda create ${conda_flags} --prefix conda_envs/${env_name} ${conda_env[${env_name}]}
		echo Created ${env_name}.
	else
		echo "${env_name} already exists, skipping..."
	fi
#	echo ${conda_env[${env_name}]}
done



#conda create -n trimming -c bioconda -c conda-forge trim-galore fastp trimmomatic

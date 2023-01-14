#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="QCParent"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --time=2:00:00
#SBATCH --mem=25G
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/Quinn/Scripts/SlurmOut/output_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/Quinn/Errors/error_%j.e

#---- Description ----
# This code takes a type, and then performs FASTQC on all reads in the Data directory corresponding to that type
#---------------------

# Load modules
module load UHTS/Quality_control/fastqc/0.11.9

#----- for holoclonal reads -----

# go to the data directory 
cd /data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Data/holo

for k in `ls -1 *.fastq.gz`;
do fastqc -t 6 ${k} -o /data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Results/1-QC/holo;
done

#----- for parental reads -----
# go to the data directory
cd ../parent

for k in `ls -1 *.fastq.gz`;
do fastqc -t 6 ${k} -o /data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Results/1-QC/parent;
done

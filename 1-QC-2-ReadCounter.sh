#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=02:00:00
#SBATCH --job-name=RScr12
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Scripts/1-QC/Output/out_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Scripts/1-QC/Errors/error_%j.e

#---- Description ----
# This script takes .fastq.gz files and outputs a .txt file with the count of the number of reads of each
#---------------------
 
# Set type (holo, para, mero, parent)
type="holo" 

# Make output text file, that can then be appended. 
echo "${type} reads:" > /data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Results/1-QC/${type}_reads.txt

# Move into the data directory
cd /data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Data/${type}

  # Since the files are zipped, we use zcat to read them and | to pipe it into the next command which is "wc -l" 
  # echo will then print the output (i.e. number of lines) followed by "/4". |bc sends "<number of lines>/4" to the 
  # basic calculator, which will do the calculation.
  # we then append this to the .txt file we made in the Results directory. 
  # this is looped over 

for i in `ls -1 *.fastq.gz`;
do echo $(zcat $i|wc -l)/4|bc >> /data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Results/1-QC/${type}_reads.txt
done

# ---- 2 ----
# Repeat the process for parental reads

type="parental"

  # Make output text file, that can then be appended.
echo "${type} reads:" > /data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Results/1-QC/${type}_reads.txt

  # Move into the data directory
cd /data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Data/${type}

  # Since the files are zipped, we use zcat to read them and | to pipe it into the next command which is "wc -l"
  # echo will then print the output (i.e. number of lines) followed by "/4". |bc sends "<number of lines>/4" to the
  # basic calculator, which will do the calculation.
  # we then append this to the .txt file we made in the Results directory.
  # this is looped over

for i in `ls -1 *.fastq.gz`;
do echo $(zcat $i|wc -l)/4|bc >> /data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Results/1-QC/${type}_reads.txt
done

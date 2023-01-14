#!/bin/bash

#SBATCH --job-name="index"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --mem=80G
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Scripts/2-Mapping/Output/output_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Scripts/2-Mapping/Errors/error_%j.e
#SBATCH --partition=pall

# This script creates an index for hisat2 and then runs hisat 2 to produce .sam output

#---- 0. Preparation ----

# load the module
module add UHTS/Aligner/hisat/2.2.1

# Store the path to the course directories as variables to improve readability
RefDir=/data/courses/rnaseq_course/lncRNAs/Project1/references
CourDir=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon

#---- 1. Indexing ----

# set the working directory
cd ${CourDir}/Results/2-Mapping/Index

# link the reference genome 
ln -s ${RefDir}/GRCh38.genome.fa .

# Use Hisat2-build to generate an index of the reference genome using default parameters
hisat2-build GRCh38.genome.fa GRCh38_index

#---- 2. Running Hisat2 on the Parental replicates ----
# Change into the SAM directory
cd ../SAM_Files

# loop through parental replicates 1,2 and 3
for i in {1..3}
# run hisat2 with the following arguments
    # -p 8 to use 8 threads 
    # --dta specifies to perform downstream transcriptome assembly, creates results compatible with StringTie
    # -x indicates the path to the reference index
    # -1 gets the m1 mate pairs 
    # -2 gets the m2 mate pairs
do hisat2 -p 8 --dta -x $ref_indx -1 ${CourDir}/Data/parent/P${i}*R1*.fastq.gz -2 ${CourDir}/Data/parent/P${i}*R2*.fastq.gz -S Parent_P${i}_L3.sam
done

#---- 3. Running Hisat2 on the holoclonal repicates ----
# loop through 1,2,5, and run Hisat2 with the same arguments as before
for i in {1,2,5}
do hisat2 -p 8 --dta -x $ref_indx -1 ${CourDir}/Data/holo/1_${i}*R1*.fastq.gz -2 ${CourDir}/Data/holo/1_${i}*R2*.fastq.gz -S Holo_1_${i}_L3.sam
done

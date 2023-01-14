#!/bin/bash

#SBATCH --job-name="sam-to-bam"
#SBATCH --nodes=2
#SBATCH --cpus-per-task=6
#SBATCH --time=08:00:00
#SBATCH --mem=40G
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Scripts/2-Mapping/Output
$#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Scripts/2-Mapping/Errors
#SBATCH --partition=pall

# This script uses SAMTools to convert the .sam files to .bam files and sorts them

# Store the path to the course directories as variables to improve readability
RefDir=/data/courses/rnaseq_course/lncRNAs/Project1/references
CourDir=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon

# load the module
module add UHTS/Analysis/samtools/1.10

# set the working directory
cd ${CourDir}/Results/2-Mapping/BAM_Files

# save the names of the SAM files to an array
replicates=("Holo_1_1_L3" "Holo_1_2_L3" "Holo_1_5_L3" "Parent_P1_L3" "Parent_P2_L3" "Parent_P3_L3")

# Use Samtools faidx to create an index
samtools faidx ${RefDir}/GRCh38.genome.fa > GRCh38.genome.fai

# Now use Samtools view to convert from SAM to BAM
# -b to generate a BAM output
# -t for a file which lists reference names and lengths

for i in {0..5};
do samtools view -b -t GRCh38.genome.fai ${CourDir}/Results/2-Mapping/SAM_Files/${replicates[$i]}.sam > ${CourDir}/Results/2-Mapping/BAM_Files/${replicates[$i]}_unsorted.bam
done

# Sort the unsorted BAM files 
   # -o writes output to a file instead of standard output
for i in {0..5};
do samtools sort -o ${replicates[$i]}_sorted.bam ${replicates[$i]}_unsorted.bam
done

#Index BAM
samtools index *_sorted.bam 

#!/bin/bash
#SBATCH --job-name=KalPrep
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --time=6:00:00
#SBATCH --mem=80G
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Scripts/4-Quantification/Output/output_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Scripts/4-Quantification/Errors/error_%j.e
#SBATCH --partition=pall

# move to the results directory 
cd /data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Results/4-Quantification 

# save some paths to variables
ref_dir="/data/courses/rnaseq_course/lncRNAs/Project1/references"
merged_gtf="/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Results/3-Assembly/all_merged.gtf"

# load the required modules
    # kallisto for indexing, and cufflinks to generate the fasta file from our reference and  merged gtf file 
module load UHTS/Analysis/kallisto/0.46.0
module load UHTS/Assembler/cufflinks/2.2.1

# Use cufflinks gffread to generate a fasta file with the spliced exons of each gff transcript 
    # -w specifies to generate a .fa file for the spliced exons
    # -g specifies the full path to the reference genome 
# gffread -w full_exons.fa -g ${ref_dir}/GRCh38.genome.fa ${merged_gtf}

# Use kallisto to create an index
    # -i specifies the name of the index file 
kallisto index -i kallisto_index.idx full_exons.fa 

 


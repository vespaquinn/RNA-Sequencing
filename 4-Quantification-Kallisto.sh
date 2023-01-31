#!/bin/bash
#SBATCH --job-name=KalPrep
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
#SBATCH --mem=80G
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Scripts/4-Quantification/Output/output_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Scripts/4-Quantification/Errors/error_%j.e
#SBATCH --partition=pall

#---- Description ----
# This script generates an index for Kallisto and then runs Kallisto Quant 
# to generate expression tables for all of the reads.
#---------------------
# move to the results directory 
cd /data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Results/4-Quantification 

# save some paths to variables
ref_dir="/data/courses/rnaseq_course/lncRNAs/Project1/references"
merged_gtf="/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Results/3-Assembly/all_merged.gtf"
data_dir="/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Data" 
# load the required modules
    # kallisto for indexing, and cufflinks to generate the fasta file from our reference and  merged gtf file 
module load UHTS/Analysis/kallisto/0.46.0
module load UHTS/Assembler/cufflinks/2.2.1

# Use cufflinks gffread to generate a fasta file with the spliced exons of each gff transcript 
    # -w specifies to generate a .fa file for the spliced exons
    # -g specifies the full path to the reference genome 
gffread -w full_exons.fa -g ${ref_dir}/GRCh38.genome.fa ${merged_gtf}

# Use kallisto to create an index
    # -i specifies the name of the index file 
kallisto index -i kallisto_index.idx full_exons.fa 

# Use kallisto quant to quantify expression in PARENTAL
    # -i indicates the index file
    # -b indicates the number of bootstrap samples, this is needed for later steps with sleuth 
    # -o indicates output directory to write output to 
    # -t indicates thread number 
    # --rf-stranded runs kallisto in strand specific mode... 
      # only fragments where the first read in the pair pseudoaligns
      # to the reverse strand of a transcript are processed.

for i in {1..3}
do kallisto quant -i kallisto_index.idx -b 100 -o Parental_${i} -t 8 --rf-stranded ${data_dir}/parent/P${i}*R1*.fastq.gz ${data_dir}/parent/P${i}*R2*.fastq.gz
done

# Now use kallisto quant with the same arguments as before but this time for holoclonal 
for i in {1,2,5}
do kallisto quant -i kallisto_index.idx -b 100 -o Holo_${i} -t 8 --rf-stranded ${data_dir}/holo/1_${i}*R1*.fastq.gz ${data_dir}/holo/1_${i}*R2*.fastq.gz
done

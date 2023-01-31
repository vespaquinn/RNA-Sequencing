#!/bin/bash
#SBATCH --job-name=GetFasta
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=20G
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Scripts/6-IntegrativeAnalysis/Output/output_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Scripts/6-IntegrativeAnalysis/Errors/error_%j.e
#SBATCH --partition=pall
#--------------------------------------------------------------------------------------------------------------------------
# ---- Setting up Variables ----
### First set up array of transcript names
transcripts=("MSTRG.5294.5" "MSTRG.30802.2" "MSTRG.24494.2" "MSTRG.5294.10" "MSTRG.8213.10" "MSTRG.14959.1" "MSTRG.30868.5" "MSTRG.7963.1" "MSTRG.28914.7" "MSTRG.12257.3" "MSTRG.30826.9" "MSTRG.10745.2" "MSTRG.22624.2" "MSTRG.20318.4" "MSTRG.2427.3" "MSTRG.17346.1" "MSTRG.6583.1" "MSTRG.8213.6" "MSTRG.5588.2" "MSTRG.6631.1")

### assign ref genome and BED dir
refgen="/data/courses/rnaseq_course/lncRNAs/Project1/references/GRCh38.genome.fa"
beddir="/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Results/6-Integrative_Analysis/BED_Files"
fastdir="/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Results/6-Integrative_Analysis/FASTA_Files"
# ---- Loading Modules ----
module load UHTS/Analysis/BEDTools/2.29.2

# ---- Extracting Fasta ----
for i in $(seq 0 19)
do bedtools getfasta -fo ${fastdir}/${transcripts[i]}.fasta -fi ${refgen} -bed ${beddir}/${transcripts[i]}.bed
done

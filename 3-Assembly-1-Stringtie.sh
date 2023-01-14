#!/bin/bash
#SBATCH --job-name=Stringtie
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --time=6:00:00
#SBATCH --mem=120G
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Scripts/3-Assembly/Output/output_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Scripts/3-Assembly/Errors/error_%j.e
#SBATCH --partition=pall

#load stringtie
module load UHTS/Aligner/stringtie/1.3.3b

#create a shortcut to the reference annotation
reference="/data/courses/rnaseq_course/lncRNAs/Project1/references/gencode.v21.chr_patch_hapl_scaff.annotation.gtf"

#create a shortcut to .bam file directory 
bamfiles="/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Results/2-Mapping/BAM_Files"

#set working directory
cd /data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Results/Assembly

# save the names of the different replicates to an array
replicates=("Holo_1_1_L3" "Holo_1_2_L3" "Holo_1_5_L3" "Parent_P1_L3" "Parent_P2_L3" "Parent_P3_L3")

# create a loop which goes through all 6 replicates and performs the stringtie function on each, generating separate .gtf transcriptome assembly files for each
# the arguments for stringtie are:
               # -o defines the name of output file
               # --rf indicates first-strand library caused by sequencing type
               # -G defines the reference annotation file 
for i in {0..5};
do stringtie -o ${replicates[$i]}_transcriptome_assembly.gtf --rf -G ${reference} ${bamfiles}/${replicates[$i]}_sorted.bam
done

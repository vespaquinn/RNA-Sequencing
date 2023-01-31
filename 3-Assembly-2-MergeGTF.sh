#!/bin/bash
#SBATCH --job-name=MergeString
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --time=6:00:00
#SBATCH --mem=120G
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Scripts/3-Assembly/Output/output_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Scripts/3-Assembly/Errors/error_%j.e
#SBATCH --partition=pall

#---- Description ----
# This script takes the GTF asssemblies generated in the 
# previous step and merges them into one meta assembly
#---------------------

#load stringtie
module load UHTS/Aligner/stringtie/1.3.3b

#create a shortcut to the reference annotation
reference="/data/courses/rnaseq_course/lncRNAs/Project1/references/gencode.v21.chr_patch_hapl_scaff.annotation.gtf"

#set working directory
cd /data/courses/rnaseq_course/lncRNAs/Project1/users/qcoxon/Results/Assembly

# save the names of all .gtf files to a .txt file for input
    # -1 lists one output per line 
ls -1 *.gtf > gtf_list.txt

#Use the stringtie --merge funtion to merge all 6 .gtf files into one 
# the arguments for stringtie are:
    # -o defines the name of output file
    # --merge defines the merge function
    # -G defines the reference annotation file 
    # --rf indicates data is stranded and to use the strand-specific mode        
stringtie --rf --merge -o all_merged.gtf -G ${reference} gtf_list.txt

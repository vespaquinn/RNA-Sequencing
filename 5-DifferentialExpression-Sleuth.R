# ------------------------------------------------------------------------------
#Differential Expression With Sleuth 
#-------------------------------------------------------------------------------
## 1. Loading libraries and set up paths to the data we want to imput
library(sleuth)
library(biomaRt)
library(annotables)

### save names of the kalisto results directories to a vairable 'sample id'
sample_id <-dir(file.path("C:","Users","quinn","OneDrive","Desktop","Uni","RNASequencing","Results","kallisto_outputs"))

### then save all the directories to 'kal_dirs'
kal_dirs <- file.path("C:","Users","quinn","OneDrive","Desktop","Uni","RNASequencing","Results","kallisto_outputs",sample_id)

### create an auciliary table detailing the relationship between the kallisto directories and the samples 
s2c <- data.frame(sample =sample_id, condition=(c(rep("holo",3),rep("control_parental",3))),path=kal_dirs)

#-------------------------------------------------------------------------------
## 2. allow conversion ENSEMBL transcript gene IDs to gene names 
tx2gene <- function(){
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                       ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  return(t2g)
}

t2g <- tx2gene()
#-------------------------------------------------------------------------------
## 3. Constructing the 'Sleuth Object'

#### initializing sleuth object, loading kallisto data into it
      # NB. use transformation function to get log2 results as by default sleuth uses natural log 
so <- sleuth_prep(s2c, target_mapping = t2g, extra_bootstrap_summary = TRUE,
                  transformation_function = function(x) log2(x + 0.5))

#### estimate parameters for the sleuth response error measurement (full) model
so <- sleuth_fit(so, ~condition, 'full')
#### estimate parameters for the sleuth (reduced) model
so <- sleuth_fit(so, ~1, 'reduced')
#### perform differential analysis (testing) using a wald test
so <- sleuth_wt(so, 'conditionholo', 'full')

#-------------------------------------------------------------------------------
### 3. Examining results
#### 
sleuth_table <- sleuth_results(so, 'conditionholo', 'wt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.000005)
head(sleuth_significant, 20)

# save sleuth significant to csv 
write.csv(sleuth_significant,"C:\\Users\\quinn\\OneDrive\\Desktop\\Uni\\RNASequencing\\Results\\SleuthTable.csv")

#### make a volcano plot
library(ggplot2)

plot_volcano(so, test = 'conditionholo', test_type = "wt", which_model = "full",
             sig_level = 0.05, point_alpha = 0.2, sig_color = "coral") +  
             xlab("log2FC") +
             ggtitle( "Holoclonal vs Parental")

#-------------------------------------------------------------------------------
### 4. Identifying Novel Genes and transcripts 
#### Novel genes/transcripts  are those whose trancript_id has the prefix:
####  'MSTRG'
sleuth_novel = cbind(sleuth_table, Novel=NA)
for (i in 1:nrow(sleuth_novel)){
  sleuth_novel$Novel[i] <- startsWith(sleuth_novel$target_id[i],"MSTRG")
}
sleuth_novel <- sleuth_novel[which(sleuth_novel$Novel),]

### Save the sleuth_novel tabe to a csv 
write.csv(sleuth_novel,"C:\\Users\\quinn\\OneDrive\\Desktop\\Uni\\RNASequencing\\Results\\SleuthNovel.csv")

#-------------------------------------------------------------------------------
### 5. Create a merged table with meta data from the GTF

gtf <- rtracklayer::import("..//all_merged.gtf")
gtf <- data.frame(gtf)
gtf_transcripts <- gtf[which(gtf$type == "transcript"),]

#### --- 5.1 for the novel transcripts ---
# loop through the transcript IDs in the sleuth_novel df and add them to vector
# 'positions'
positions <- c()
for (target in sleuth_novel$target_id){
  positions <- c(positions,which(gtf_transcripts$transcript_id == target))
}
# merge the gtf data with the sleuth qvalues
merged_novel <- data.frame(transcript = gtf_transcripts$transcript_id[positions],
                           start = gtf_transcripts$start[positions],
                           end = gtf_transcripts$end[positions],
                           qvalue = sleuth_novel$qval)
# save to csv 
write.csv(merged_novel,"MergedNovel.csv")
# Create table with formattable
formattable(merged_novel[1:20,])

#### --- 5.2 for the known transcripts ---
# loop through the transcript IDs in the sleuth_novel df and add them to vector
# 'positions'
positions <- c()
for (target in sleuth_known$target_id){
  positions <- c(positions,which(gtf_transcripts$transcript_id == target))
}
# merge the gtf data with the sleuth qvalues
merged_known <- data.frame(transcript = gtf_transcripts$transcript_id[positions],
                           start = gtf_transcripts$start[positions],
                           end = gtf_transcripts$end[positions],
                           qvalue = sleuth_known$qval)
# save to csv 
write.csv(merged_novel,"MergedKnown.csv")

# create table
formattable(merged_known)
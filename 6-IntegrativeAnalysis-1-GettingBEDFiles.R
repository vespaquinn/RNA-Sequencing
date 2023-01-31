#--------- Extracting Novel Transcripts and Creating Bed Files -----
sleuth_table <- read.table("SleuthTable.csv", sep =",", header = TRUE)

gtf <- rtracklayer::import("all_merged.gtf")
gtf_df=as.data.frame(gtf)

sleuth_novel = cbind(sleuth_table, Novel=NA)
for (i in 1:nrow(sleuth_novel)){
  sleuth_novel$Novel[i] <- startsWith(sleuth_novel$target_id[i],"MSTRG")
}
sleuth_novel <- sleuth_novel[which(sleuth_novel$Novel),]
sleuth_novel_max <- sleuth_novel[1:20,]
# now sleuth_novel_max is a list of the top 20 most significantly dif expressed
# novel transcripts
#### this loop adds a new row to merged_df for each of the IDs in sleuth_novel_max
#### it adds the first instance of it which will always be the whole transcript, as opposed
#### to one of the exons. 
merged_df <- gtf_df[0,]
for (i in 1:20){
  new_transcript <- gtf_df[which(gtf_df$transcript_id == sleuth_novel_max$target_id[i]),]
  merged_df[nrow(merged_df) + 1, ] <- new_transcript[1,]
}
##### now we can clean up the df by keeping the columns needed for the bed file
###### NB. the numbers 1,2,3,11,8,5, corresponds with chr, start, end, name, score,strand
merged_df <- merged_df[,c(1,2,3,11,8,5)]
#---- Now we can export this to a bed file 
for (i in 1:20){
  file_name <- paste(merged_df$transcript_id[i],'.bed', sep ='')
  write.table(merged_df[i,], file=file_name, quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
}
print(merged_df$transcript_id)

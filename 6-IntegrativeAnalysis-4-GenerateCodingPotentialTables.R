# load packages
library(formattable)
# load in the meta data
merged_novel <- read.csv("MergedNovel.csv", sep =",",header =T)
merged_novel <- merged_novel[,2:5]
#---- Create a list of transcript Names
transcripts=c("MSTRG.5294.5","MSTRG.30802.2","MSTRG.24494.2","MSTRG.5294.10",
              "MSTRG.8213.10","MSTRG.14959.1","MSTRG.30868.5","MSTRG.7963.1",
              "MSTRG.28914.7","MSTRG.12257.3","MSTRG.30826.9","MSTRG.10745.2",
              "MSTRG.22624.2","MSTRG.20318.4","MSTRG.2427.3","MSTRG.17346.1",
              "MSTRG.6583.1","MSTRG.8213.6","MSTRG.5588.2","MSTRG.6631.1")
#---- extract the relavent transcripts from the meta data and sort it by transcript
positions <- c()
for (transcript in transcripts){
  positions <- c(positions,which(merged_novel == transcript))
}
merged_extracted <- merged_novel[positions,]
merged_sorted <- merged_extracted[order(merged_extracted$transcript),]
# --- Create a blank df with the column names from the imported TSVs to append
df <- read.table("MSTRG.5294.5", sep = "\t")
df <- cbind(Transcript = NA, df)
df <- df[0,]
rownames(df) <- NULL

# --- Loop through the trancripts and append each to the empty df 
for (i in 1:20){
  newdf <- read.table(transcripts[i], sep = "\t")
  df[i,2:6] <- newdf[,1:5]
  df[i,1] <- transcripts[i]
}
rm(newdf)
#--- sort the df so they can be combined 
df_sorted <- df[order(df$Transcript),]
# --- combine the dfs and then sort by coding probability 
combined <- data.frame(merged_sorted,
                       coding_probability = df_sorted$coding_prob)
combined_sorted <- combined[order(combined$coding_probability),]
#--- with the suggested cutoff by CPAT of 0.364 save the positions of 
# transcripts with significant coding potential 
coding_potential <- which(combined_sorted$coding_probability > 0.364)
# Make the table
formattable(combined_sorted, align = c("l",rep("r", NCOL(combined_sorted) - 1)), list(
  `Indicator Name` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  area(row = coding_potential, col = 5) ~ color_tile("olivedrab3","olivedrab3")))
#--------------------------------------------------------------------------------
# Generate another table of the same data, sorted by qvalue

combined_sorted_qval <- combined[order(combined$qvalue),]
coding_potential2 <- which(combined_sorted_qval$coding_probability > 0.364)

formattable(combined_sorted_qval, align = c("l",rep("r", NCOL(combined_sorted) - 1)), list(
  `Indicator Name` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  area(row = coding_potential2, col = 5) ~ color_tile("olivedrab3","olivedrab3")))

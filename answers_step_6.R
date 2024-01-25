### Parse data from step 6

# Load dplyer
install.packages("dplyr")
library(dplyr)


## Read bedfiles
# these were generated on the cluster --> step_6/output_files
intersect_data <- read.table("C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/intersect.bed", 
                             sep = "\t")
poly_a_data <- read.table("C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/three_prime_intersect.bed", 
                          sep = "\t")
tss_data <- read.table("C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/tss_intersect.bed", 
                       sep = "\t")
cpat_data <- read.table("C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/cpat_out.dat", 
                        sep = "\t")
# also need file of all novel transcripts
novel_light_bed <- read.table("C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/novel_light.bed", 
                        sep = "\t")
transcript_data <- read.table("C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/sig_transcripts.tsv", 
                              sep = "\t", header = T)


## Create summary table for novel transcripts
novel_light_bed$intergenic <- NA
novel_light_bed$tss <- NA
novel_light_bed$poly_a <- NA
novel_light_bed$tss_poly_a <- NA
novel_light_bed$protein_potential <- NA


# Add tss info
for (i in 1:length(tss_data$V4)){
  target_index <- which(trimws(novel_light_bed$V4, "l") ==  trimws(tss_data$V4[i], "l"))
  novel_light_bed$tss[target_index] <- T
}

# Add intergenic info
for (i in 1:length(intersect_data$V4)){
  target_index <- which(trimws(novel_light_bed$V4, "l") ==  trimws(intersect_data$V4[i], "l"))
  novel_light_bed$intergenic[target_index] <- T
}

# Add poly-a info
for (i in 1:length(poly_a_data$V4)){
  target_index <- which(trimws(novel_light_bed$V4, "l") ==  trimws(poly_a_data$V4[i], "l"))
  novel_light_bed$poly_a[target_index] <- T
}

# Add info if both poly A AND tss are True
novel_light_bed$tss_poly_a[novel_light_bed$tss == novel_light_bed$poly_a] <- T



# Summary for novel transcripts
correct_tss <- sum(novel_light_bed$tss, na.rm = T)
correct_poly_a <- sum(novel_light_bed$poly_a, na.rm = T)
correct_5_3 <- sum(novel_light_bed$tss_poly_a, na.rm = T)
intergenic_transcripts <- sum(novel_light_bed$intergenic, na.rm = T)
total_novel <- length(novel_light_bed$V4)

# Calculate Ratios when compared to novel
tss_ratio <- correct_tss/total_novel
poly_a_ratio <- correct_poly_a/total_novel
corr_5_3_ratio <- correct_5_3/total_novel
intergenic_ratio <- intergenic_transcripts/total_novel 


# Create a table with all results and ratios --> still need to add coding potential
result_table <- data.frame("Category" = c("Correct TSS", "Correct Poly-A", "Correct TSS & Poly-A",
                                          "Intergenic Transcripts", "Total Novel Transcripts"),
                           "Number" = c(correct_tss, correct_poly_a, correct_5_3,
                                        intergenic_transcripts, total_novel),
                           "Ratio of Novel" = c(tss_ratio, poly_a_ratio, corr_5_3_ratio,
                                                intergenic_ratio, 1))



## Create summary table for significant transcripts
# Create new cols
transcript_data$intergenic <- NA
transcript_data$tss <- NA
transcript_data$poly_a <- NA
transcript_data$tss_poly_a <- NA
transcript_data$protein_potential <- NA

# Add tss info
for (i in 1:length(tss_data$V4)){
  target_index <- which(trimws(transcript_data$target_id, "l") ==  trimws(tss_data$V4[i], "l"))
  transcript_data$tss[target_index] <- T
}

# Add intergenic info
for (i in 1:length(intersect_data$V4)){
  target_index <- which(trimws(transcript_data$target_id, "l") ==  trimws(intersect_data$V4[i], "l"))
  transcript_data$intergenic[target_index] <- T
}

# Add poly-a info
for (i in 1:length(poly_a_data$V4)){
  target_index <- which(trimws(transcript_data$target_id, "l") ==  trimws(poly_a_data$V4[i], "l"))
  transcript_data$poly_a[target_index] <- T
}

# Add info if both poly A AND tss are True
transcript_data$tss_poly_a[transcript_data$tss == transcript_data$poly_a] <- T




## Summary for significant transcripts
num_tss_sig <- sum(transcript_data$tss, na.rm = T)
num_polyA_sig <- sum(transcript_data$poly_a, na.rm = T)
num_tss_polyA_sig <- sum(transcript_data$tss_poly_a, na.rm = T)
num_intergenic_sig <- sum(transcript_data$intergenic, na.rm = T)
total_significant <- length(transcript_data$target_id)
novel_significant <- sum(grepl("^MSTRG", transcript_data$target_id))

# Calculate Ratios when compared to significant transcripts
tss_ratio_sig <- num_tss_sig/novel_significant
poly_a_ratio_sig <- num_polyA_sig/novel_significant
corr_5_3_ratio_sig <- num_tss_polyA_sig/novel_significant
intergenic_ratio_sig <- num_intergenic_sig/novel_significant

# Create a table with all results and ratios --> still need to add coding potential
sig_result_table <- data.frame("Number" = c(num_tss_sig, num_polyA_sig, num_tss_polyA_sig,
                                 num_intergenic_sig, novel_significant),
                                "Ratio of Novel Significant" = c(tss_ratio_sig, poly_a_ratio_sig, corr_5_3_ratio_sig,
                                                      intergenic_ratio_sig, 1))
row.names(sig_result_table) <- c("Correct TSS", "Correct Poly-A", "Correct TSS & Poly-A",
                             "Intergenic Transcripts", "Total Novel Significant Transcripts")




## Create a table for significant lncRNAs
sig_lnc_data <- filter(transcript_data, biotype == "lncRNA")



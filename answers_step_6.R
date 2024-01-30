### Parse data from step 6

# Load dplyer
install.packages("dplyr")
library(dplyr)


## Read bedfiles from evaluation on cluster
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
# Also load sig gene data to make plots
gene_data <- read.table("C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/sig_genes.tsv", 
                              sep = "\t", header = T)


# Filter cpat output for transcripts with coding prob lower than 0.364 --> recommended on cpat doc page
filtered_cpat <- filter(cpat_data, coding_prob < 0.364)








## Create summary table for novel transcripts
novel_light_bed$intergenic <- NA
novel_light_bed$tss <- NA
novel_light_bed$poly_a <- NA
novel_light_bed$tss_poly_a <- NA
novel_light_bed$protein_potential_bool <- NA
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


# Add Info about coding probability
novel_light_bed$protein_potential <- cpat_data$coding_prob

for (i in 1:length(novel_light_bed$V1)){
  if (novel_light_bed$protein_potential[i] < 0.364){
    novel_light_bed$protein_potential_bool[i] <- T
  }
}



# Summary for novel transcripts
correct_tss <- sum(novel_light_bed$tss, na.rm = T)
correct_poly_a <- sum(novel_light_bed$poly_a, na.rm = T)
correct_5_3 <- sum(novel_light_bed$tss_poly_a, na.rm = T)
intergenic_transcripts <- sum(novel_light_bed$intergenic, na.rm = T)
is_protein_coding <- sum(is.na(novel_light_bed$protein_potential_bool))
total_novel <- length(novel_light_bed$V4)

# Calculate Ratios when compared to novel
tss_ratio <- correct_tss/total_novel
poly_a_ratio <- correct_poly_a/total_novel
corr_5_3_ratio <- correct_5_3/total_novel
intergenic_ratio <- intergenic_transcripts/total_novel
is_protein_coding_ratio <- is_protein_coding/total_novel


# Create a table with all results and ratios for step 6
result_table <- data.frame("Category" = c("Correct TSS", "Correct Poly-A", "Correct TSS & Poly-A",
                                          "Intergenic Transcripts", "Protein Coding Transcripts", 
                                          "Total Novel Transcripts"),
                           "Number" = c(correct_tss, correct_poly_a, correct_5_3,
                                        intergenic_transcripts, is_protein_coding, total_novel),
                           "Ratio of Novel" = c(tss_ratio, poly_a_ratio, corr_5_3_ratio,
                                                intergenic_ratio, is_protein_coding_ratio, 1))








## Create summary table for significant transcripts
# Create new cols
transcript_data$intergenic <- F
transcript_data$tss <- F
transcript_data$poly_a <- F
transcript_data$tss_poly_a <- F
transcript_data$protein_potential_bool <- F
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


# Add Info about coding probability
for (i in 1:length(novel_light_bed$V1)){
  if (!is.na(novel_light_bed$protein_potential_bool[i])){
    target_index <- which(trimws(transcript_data$target_id, "l") ==  trimws(novel_light_bed$V4[i], "l"))
    transcript_data$protein_potential_bool[target_index] <- T
  }
}

for (i in 1:length(transcript_data$target_id)){
  if (is.na(transcript_data$ens_gene[i])){
    target_index <- which(trimws(novel_light_bed$V4, "l") == trimws(transcript_data$target_id[i], "l"))
    if (!length(target_index) == 0){
      transcript_data$protein_potential[i] <- novel_light_bed$protein_potential[target_index]
    }
    else {
      next
    }
  }
}

# Add info if both poly A AND tss are True
transcript_data$tss_poly_a[transcript_data$tss & transcript_data$poly_a] <- T




## Summary for significant transcripts
num_tss_sig <- sum(transcript_data$tss, na.rm = T)
num_polyA_sig <- sum(transcript_data$poly_a, na.rm = T)
num_tss_polyA_sig <- sum(transcript_data$tss_poly_a, na.rm = T)
num_intergenic_sig <- sum(transcript_data$intergenic, na.rm = T)
total_significant <- length(transcript_data$target_id)
novel_significant <- sum(grepl("^MSTRG", transcript_data$target_id))
num_protein_coding <- novel_significant - sum(transcript_data$protein_potential_bool, na.rm = T)

# Calculate Ratios when compared to significant transcripts
tss_ratio_sig <- num_tss_sig/novel_significant
poly_a_ratio_sig <- num_polyA_sig/novel_significant
corr_5_3_ratio_sig <- num_tss_polyA_sig/novel_significant
intergenic_ratio_sig <- num_intergenic_sig/novel_significant
protein_coding_ratio <- num_protein_coding/novel_significant

# Create a table with all results and ratios --> still need to add coding potential
sig_result_table <- data.frame("Number" = c(num_tss_sig, num_polyA_sig, num_tss_polyA_sig,
                                 num_intergenic_sig, num_protein_coding, novel_significant),
                                "Ratio of Novel Significant" = c(tss_ratio_sig, poly_a_ratio_sig, corr_5_3_ratio_sig,
                                                      intergenic_ratio_sig, protein_coding_ratio, 1))
row.names(sig_result_table) <- c("Correct TSS", "Correct Poly-A", "Correct TSS & Poly-A",
                             "Intergenic Transcripts", "Protein Coding Transcripts", 
                             "Total Novel Significant Transcripts")



## Save both result tables for novel and novel significant transcripts
## Save complete transcript and gene expression table
result_path_novel <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/step6_answers_novel.tsv"
result_path_sig_novel <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/step6_answers_sig_novel.tsv"

# save files for upload later
write.table(sig_result_table, file = result_path_sig_novel , sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)
write.table(result_table, file = result_path_novel, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)











## Create a ordered table for a scoring with novel significant genes
ordered_results <- filter(transcript_data, grepl("^MST", target_id))

# Calculate a score for transcripts based on correct 3', 5', intergenic and protein coding potential
ordered_results$score <- 0

for (i in 1:length(ordered_results$target_id)){
  if (ordered_results$intergenic[i]){
    ordered_results$score[i] = ordered_results$score[i] + 2
  }
  if (ordered_results$tss_poly_a[i]){
    ordered_results$score[i] = ordered_results$score[i] + 3
  }
  if (ordered_results$protein_potential_bool[i]){
    ordered_results$score[i] = ordered_results$score[i] + 1
  }
}


# sort the ordered results
sorted_results <- ordered_results[order(ordered_results$score, decreasing = T),]



## Create a table for novel lncRNAs --> >200nt
lnc_table <- novel_light_bed
# Create column for orf length
lnc_table$ORF_length <- lnc_table$V3 - lnc_table$V2

lnc_table <- filter(lnc_table, ORF_length > 200, protein_potential < 0.364)


## Create a table for significant lncRNAs --> >200nt
lnc_table_sig <- transcript_data
# look if a novel significant transcript can be found in lnc_table, if yes add biotype lncRNA
for (i in 1:length(lnc_table_sig$biotype)){
  if (trimws(lnc_table_sig$target_id[i], "l") %in% trimws(lnc_table$V4, "l"))
    lnc_table_sig$biotype[i] <- "lncRNA"
}

# Filter for only lncRNAs
lnc_table_sig <- filter(lnc_table_sig, biotype == "lncRNA")
# Filter for novel
lnc_table_sig_novel <- filter(lnc_table_sig, is.na(ens_gene))
# Filter for intergenic
lnc_table_intergenic <- filter(lnc_table_sig_novel, intergenic == T)



sorted_results_path <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/step6_sorted_novel.tsv"
sig_lnc_path <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/step6_sig_lnc.tsv"
sig_novel_lnc_path <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/step6_sig_novel_lnc.tsv"


# save files for upload later
write.table(sorted_results, file = sorted_results_path , sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)
write.table(lnc_table_sig, file = sig_lnc_path, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)
write.table(lnc_table_sig_novel, file = sig_novel_lnc_path, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)





### Create a volcano plot for lncRNAs
# Volcano plot with Enhanced Volcano
# Create a DF for enhanced volcano
volcano_lnc_df <- data.frame(log2FC = lnc_table_sig$log2FC,
                                    q_vals = lnc_table_sig$qval)
row.names(volcano_lnc_df) <- lnc_table_sig$target_id

volcano_transcript_df <- data.frame(log2FC = transcript_data$log2FC,
                                    q_vals = transcript_data$qval)
row.names(volcano_transcript_df) <- transcript_data$target_id

volcano_gene_df <- data.frame(log2FC = gene_data$log2FC,
                              q_vals = gene_data$qval)
row.names(volcano_gene_df) <- gene_data$gene_id

# Plot Volcano
transcript_volcano <-EnhancedVolcano(volcano_transcript_df,
                                     lab = rownames(volcano_transcript_df),
                                     x = 'log2FC',
                                     y = 'q_vals',
                                     ylim = c(0, -log10(10e-200)),
                                     xlim = c(-8, 8),
                                     xlab = "",
                                     title = 'Transcript Level',
                                     pointSize = 3.0,
                                     labSize = 5.0,
)
gene_volcano <- EnhancedVolcano(volcano_gene_df,
                                lab = rownames(volcano_gene_df),
                                x = 'log2FC',
                                y = 'q_vals',
                                ylim = c(0, -log10(10e-200)),
                                xlim = c(-8, 8),
                                ylab = "",
                                title = 'Gene Level',
                                pointSize = 3.0,
                                labSize = 5.0
)


lnc_volcano <- EnhancedVolcano(volcano_lnc_df,
                              lab = rownames(volcano_lnc_df),
                              x = 'log2FC',
                              y = 'q_vals',
                              ylim = c(0, -log10(10e-200)),
                              xlim = c(-8, 8),
                              xlab = "",
                              ylab = "",
                              title = 'lncRNA level',
                              pointSize = 3.0,
                              labSize = 5.0
)




plot(transcript_volcano)

library(ggpubr)


### Plot the 4 Plots for Fig. 1
ggarrange(transcript_volcano, gene_volcano, lnc_volcano,
          ncol=3, nrow=1, 
          common.legend = TRUE, legend="bottom")

## Save all important files --> finish this tomorrow!
# sig_genes_path <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/sig_genes.tsv"
# sig_transcripts_path <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/sig_transcripts.tsv"
# 
# 
# write.table(sig_genes_table, file = sig_genes_path, sep = "\t", quote = FALSE,
#             col.names = TRUE, row.names = FALSE)
# write.table(sig_transcripts_table, file = sig_transcripts_path, sep = "\t", quote = FALSE,
#             col.names = TRUE, row.names = FALSE)




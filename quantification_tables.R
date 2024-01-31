### get gene level expression tables for all replicates

# load dplyr
library(dplyr)


# Load gene mapping created in step 5
gene_mapping <- read.table("C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/gene_mapping.tsv", 
                             sep = "\t", header = T)



# set working directory
setwd("C:/Users/nbaho/OneDrive/Desktop/Bioinf/R-studio/Kallisto_Data")


# Read kallisto data
holo_1_1 <- read.table("C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/Kallisto_Data/S1_1_L3/abundance.tsv", 
                        sep = "\t", header = T)
holo_1_2 <- read.table("C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/Kallisto_Data/S1_2_L3/abundance.tsv", 
                       sep = "\t", header = T)
holo_1_5 <- read.table("C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/Kallisto_Data/S1_3_L3/abundance.tsv", 
                       sep = "\t", header = T)
parental_1 <- read.table("C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/Kallisto_Data/P1_L3/abundance.tsv", 
                       sep = "\t", header = T)
parental_2 <- read.table("C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/Kallisto_Data/P2_L3/abundance.tsv", 
                       sep = "\t", header = T)
parental_3 <- read.table("C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/Kallisto_Data/P3_L3/abundance.tsv", 
                       sep = "\t", header = T)








for (i in 1:length(holo_1_1$target_id)){
  target_index <- which(gene_mapping$target_id == holo_1_1$target_id[i])
  if (!length(target_index) == 0){
    holo_1_1$gene_id[i] <- gene_mapping$gene_id[target_index]
  } else {
    holo_1_1$gene_id[i] <- NA
  }
}

holo_1_1f <- filter(holo_1_1, !is.na(gene_id))




for (i in 1:length(holo_1_2$target_id)){
  target_index <- which(gene_mapping$target_id == holo_1_2$target_id[i])
  if (!length(target_index) == 0){
    holo_1_2$gene_id[i] <- gene_mapping$gene_id[target_index]
  } else {
    holo_1_2$gene_id[i] <- NA
  }
}

# Filter just in case
holo_1_2_f <- filter(holo_1_2, !is.na(gene_id))


for (i in 1:length(holo_1_5$target_id)){
  target_index <- which(gene_mapping$target_id == holo_1_5$target_id[i])
  if (!length(target_index) == 0){
    holo_1_5$gene_id[i] <- gene_mapping$gene_id[target_index]
  } else {
    holo_1_5$gene_id[i] <- NA
  }
}

holo_1_5f <- filter(holo_1_5, !is.na(gene_id))
                   

for (i in 1:length(parental_1$target_id)){
  target_index <- which(gene_mapping$target_id == parental_1$target_id[i])
  if (!length(target_index) == 0){
    parental_1$gene_id[i] <- gene_mapping$gene_id[target_index]
  } else {
    parental_1$gene_id[i] <- NA
  }
}

parental_1f <- filter(parental_1, !is.na(gene_id))


for (i in 1:length(parental_2$target_id)){
  target_index <- which(gene_mapping$target_id == parental_2$target_id[i])
  if (!length(target_index) == 0){
    parental_2$gene_id[i] <- gene_mapping$gene_id[target_index]
  } else {
    parental_2$gene_id[i] <- NA
  }
}

parental_2f <- filter(parental_2, !is.na(gene_id))


for (i in 1:length(parental_3$target_id)){
  target_index <- which(gene_mapping$target_id == parental_3$target_id[i])
  if (!length(target_index) == 0){
    parental_3$gene_id[i] <- gene_mapping$gene_id[target_index]
  } else {
    parental_3$gene_id[i] <- NA
  }
}

parental_3f <- filter(parental_3, !is.na(gene_id))


                   
                          



# parental_3_gene <- gene_mapping
# for (i in 1:length(gene_mapping$target_id)){
#   target_index <- which(parental_3$target_id == gene_mapping$target_id[i])
#   if (!length(target_index) == 0){
#     parental_3_gene$tpm[i] <- parental_3$tpm[target_index]
#   } else {
#     parental_3_gene$tpm[i] <- 0
#   }
# }






# Sum up tpm, length, eff length and est_counts

holo_1_1_gene_sum <- holo_1_1f %>%
  group_by(gene_id) %>%
  summarise(length = sum(length), eff_length = sum(eff_length), 
            est_counts = sum(est_counts), tpm = sum(tpm)) %>%
  filter(tpm > 0)


holo_1_2_gene_sum <- holo_1_2_f %>%
  group_by(gene_id) %>%
  summarise(length = sum(length), eff_length = sum(eff_length), 
            est_counts = sum(est_counts), tpm = sum(tpm)) %>%
  filter(tpm > 0)

holo_1_5_gene_sum <- holo_1_5f %>%
  group_by(gene_id) %>%
  summarise(length = sum(length), eff_length = sum(eff_length), 
            est_counts = sum(est_counts), tpm = sum(tpm)) %>%
  filter(tpm > 0)

parental_1_gene_sum <- parental_1f %>%
  group_by(gene_id) %>%
  summarise(length = sum(length), eff_length = sum(eff_length), 
            est_counts = sum(est_counts), tpm = sum(tpm)) %>%
  filter(tpm > 0)

parental_2_gene_sum <- parental_2f %>%
  group_by(gene_id) %>%
  summarise(length = sum(length), eff_length = sum(eff_length), 
            est_counts = sum(est_counts), tpm = sum(tpm)) %>%
  filter(tpm > 0)

parental_3_gene_sum <- parental_3f %>%
  group_by(gene_id) %>%
  summarise(length = sum(length), eff_length = sum(eff_length), 
            est_counts = sum(est_counts), tpm = sum(tpm)) %>%
  filter(tpm > 0)





differential_expression_transcript <- data.frame(target_id = holo_1_1$target_id,
                                                 h1_1_tpm = holo_1_1$tpm,
                                                 h1_2_tpm = holo_1_2$tpm,
                                                 h1_5_tpm = holo_1_5$tpm,
                                                 p1_tpm = parental_1$tpm,
                                                 p2_tpm = parental_2$tpm,
                                                 p3_tpm = parental_3$tpm)




result_path_transcripts <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/differential_expression_transcripts.tsv"
result_path_gene_h1 <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/gene_diff_h1.tsv"
result_path_gene_h2 <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/gene_diff_h2.tsv"
result_path_gene_h5 <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/gene_diff_h5.tsv"
result_path_gene_p1 <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/gene_diff_p1.tsv"
result_path_gene_p2 <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/gene_diff_p2.tsv"
result_path_gene_p3 <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/gene_diff_p3.tsv"


# save files for upload later
write.table(differential_expression_transcript, file = result_path_transcripts, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)
write.table(holo_1_1_gene_sum, file = result_path_gene_h1, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)
write.table(holo_1_2_gene_sum, file = result_path_gene_h2, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)
write.table(holo_1_5_gene_sum, file = result_path_gene_h5, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)
write.table(parental_1_gene_sum, file = result_path_gene_p1 , sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)
write.table(parental_2_gene_sum, file = result_path_gene_p2, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)
write.table(parental_3_gene_sum, file = result_path_gene_p3, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)






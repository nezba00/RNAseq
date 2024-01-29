### Get answers for step 3
#--> how many exons, transcripts, and genes in meta assembly?
#--> how many are novel
#--> how many transcripts are composed of one single exon?


library(tidyr)
library(dplyr)



# Load merged stringtie output
gtf_data <- read.table("C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/all_merged.gtf", sep = "\t")
gtf_df <- data.frame(gtf_data)






# Count exons and transcripts
exon_ts_count <- table(gtf_df$V3)

exon_number <- exon_ts_count[1]
transcript_number <- exon_ts_count[2]






# get df without exons
no_exon_gtf <- filter(gtf_df, V3 == "transcript")


## Split last column of gtf_df into multiple columns
gtf_df_split <- separate_wider_delim(no_exon_gtf, V9, ";", names = c("gene_id","transcript_id"), too_many = "merge")
gtf_df_split <- separate_wider_delim(gtf_df_split, transcript_id, ";", names = c("transcript_id","rest"), too_many = "merge", too_few = "align_end")
gtf_df_split <- separate_wider_delim(gtf_df_split, rest, ";", names = c("gene_name","rest"), too_many = "merge", too_few = "align_end")


# remove "gene_id" from gene_id column
gtf_df_split$gene_id <- gsub("gene_id ", "", gtf_df_split$gene_id)
# remove "transcript_id" from transcript_id column
gtf_df_split$transcript_id <- gsub(" transcript_id ", "", gtf_df_split$transcript_id)







# Get number of genes (novel and annotated)
gene_number <- length(unique(gtf_df_split$gene_id))
#>66851

# Novel gene number
novel_genes <- length(unique(gtf_df_split$gene_id[is.na(gtf_df_split$gene_name)]))
#>6672

# Novel transcript number
novel_transcripts <- length(unique(grep("^MSTRG", gtf_df_split$transcript_id)))
#>13000






## Create an exon df 

exon_gtf <- filter(gtf_df, V3 == "exon")


## Split last column of gtf_df into multiple columns
exon_df_split <- separate_wider_delim(exon_gtf, V9, ";", names = c("gene_id","transcript_id"), too_many = "merge")
exon_df_split <- separate_wider_delim(exon_df_split, transcript_id, ";", names = c("transcript_id","rest"), too_many = "merge", too_few = "align_end")
exon_df_split <- separate_wider_delim(exon_df_split, rest, ";", names = c("exon_nr","rest"), too_many = "merge", too_few = "align_end")

# remove "gene_id" from gene_id column
exon_df_split$gene_id <- gsub("gene_id ", "", exon_df_split$gene_id)
# remove "transcript_id" from transcript_id column
exon_df_split$transcript_id <- gsub(" transcript_id ", "", exon_df_split$transcript_id)


# Get number of novel exons
novel_exons <- length(unique(grep("^MSTRG", exon_df_split$transcript_id)))


# filter df for rows that have exon 2
# --> total amount of genes/transcripts - amount of genes/transcripts with 2 exons
# --> = amount of genes/transcripts with only 1 exon

multi_exon_df <- filter(exon_df_split, exon_nr == " exon_number 2")


# Get number of multi exon genes (novel and annotated)
multi_exon_gene_number <- length(unique(multi_exon_df$gene_id))
#>41816

# Get number of multi exon transcripts (novel and annotated)
multi_exon_ts_number <- length(unique(multi_exon_df$transcript_id))
#>260176

# Calculate number of single exon transcripts/genes
single_exon_gene_number <- gene_number - multi_exon_gene_number
#> 25035
single_exon_ts_number <- transcript_number - multi_exon_ts_number
#>29592






#--> how many exons, transcripts, and genes in meta assembly?
#--> how many are novel
#--> how many transcripts are composed of one single exon?


# Create a df with all numbers 
results_3 <- data.frame(total = c(exon_number, transcript_number, gene_number),
                        novel = c(novel_exons, novel_transcripts, novel_genes),
                        single_exon = c(NA, single_exon_ts_number, single_exon_gene_number))
rownames(results_3) <- c("exons", "transcripts", "genes")




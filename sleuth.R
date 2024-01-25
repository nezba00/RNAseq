# Install rhdf5
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rhdf5")
library("rhdf5")
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
BiocManager::install("biomaRt")
library(biomaRt)


# Install sleuth
install.packages("devtools")
devtools::install_github("pachterlab/sleuth")

# Load Sleuth
library(sleuth)

# Load cowplot
install.packages("cowplot")
library(cowplot)

# Load dplyer
install.packages("dplyr")
library(dplyr)

library(tidyr)

# set working directory
setwd("C:/Users/nbaho/OneDrive/Desktop/Bioinf/R-studio")

# Get name of all kallisto dirs
id <- dir(file.path("..", "RNA_Seq_Project", "Kallisto_Data" ))

# Create path to the results
kallisto_dirs <- file.path("..", "RNA_Seq_Project", "Kallisto_Data", id, "abundance.h5")


# Load merged GTF file
gtf_data <- read.table("C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/all_merged.gtf", sep = "\t")
gtf_df <- data.frame(gtf_data)


# create metadata df
ref_df <- data.frame(sample = id,condition = c("parental", "parental", "parental", "HoloClonal", "HoloClonal", "HoloClonal"))

# append the path to the kallisto directories
ref_df$path <- kallisto_dirs


### Get Biotypes and make target mapping
# biomaRt::useMart() to connect to specified biomart database
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")

# Show possible attributes for getBM
attributes = listAttributes(mart)

# get information for sleuth table in the end
mapping <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name", "gene_biotype",
                                     "transcript_biotype"), mart = mart)

# For target mapping in sleuth we need at least one column "target id"!
mapping <- dplyr::rename(mapping, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
###building sleuth object
##so object contains info about experiment and details of the model for DE testing
#loading kallisto data into the so object




# Initialize Sleuth Object
#--> changed transformation function to log2FC
sleuth_object <- sleuth_prep(ref_df, extra_bootstrap_summary = TRUE, 
                             target_mapping = mapping,
                             transformation_function = function(x) log2(x + 0.5))


# Fit full model
sleuth_object <- sleuth_fit(sleuth_object, ~condition, 'full')

#---- Error message
# 2 NA values were found during variance shrinkage estimation due to mean observation values outside of the range used for the LOESS fit.
# The LOESS fit will be repeated using exact computation of the fitted surface to extrapolate the missing values.
# These are the target ids with NA values: ENST00000582591.1, ENST00000643908.1
#---



# Fit reduce model
sleuth_object <- sleuth_fit(sleuth_object, ~1, 'reduced')
#----
# 1 NA values were found during variance shrinkage estimation due to mean observation values outside of the range used for the LOESS fit.
# The LOESS fit will be repeated using exact computation of the fitted surface to extrapolate the missing values.
# These are the target ids with NA values: ENST00000643908.1
#----



# Perform Sleuth test XXX --> we should do it with wt
# sleuth_object <- sleuth_lrt(sleuth_object, 'reduced', 'full')
sleuth_object <- sleuth_wt(sleuth_object, "conditionparental")

# Read sleuth data into a table
# Code from pachterlab
sleuth_table <- sleuth_results(sleuth_object, 'conditionparental', 'wt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

# Create a new dataframe for results
sig_transcripts_table <- data.frame(target_id = sleuth_significant$target_id,
                                    ens_gene = sleuth_significant$ens_gene,
                                    ext_gene = sleuth_significant$ext_gene,
                                    biotype = sleuth_significant$transcript_biotype,
                                    log2FC = sleuth_significant$b,
                                    qval = sleuth_significant$qval)








### Do the same for gene level expression
# Create df with transcript id and gene name
no_exon_gtf <- filter(gtf_df, V3 == "transcript")
# select columns: CHR, star, end, transcript id, strand
mapping_precursor <- select(no_exon_gtf, "V1","V3","V4","V5","V7","V9")
# Split last column with many attributes separated by ";"
mapping_precursor_split <- separate_wider_delim(mapping_precursor, V9, ";", names = c("gene_id","transcript_id"), too_many = "merge")
mapping_split_transcripts <- separate_wider_delim(mapping_precursor_split, transcript_id, ";", names = c("transcript_id", "rest"), too_many = "merge")
mapping_split_transcripts <- separate_wider_delim(mapping_split_transcripts, rest, ";", names = c("gene_name", "ref_gene_id"), too_many = "merge", too_few = "align_end")

# remove "gene_id" from gene_id column
mapping_split_transcripts$gene_id <- gsub("gene_id ", "", mapping_split_transcripts$gene_id)
# remove "transcript_id" from transcript_id column
mapping_split_transcripts$transcript_id <- gsub("transcript_id ", "", mapping_split_transcripts$transcript_id)
# remove "gene_name" from transcript_id column
mapping_split_transcripts$gene_name <- gsub("gene_name ", "", mapping_split_transcripts$gene_name)


# df should contain target id and gene id
# novel_mapping <- data.frame("target_id" = trimws(mapping_split_transcripts$transcript_id, "l"),
#                             "ens_gene" = mapping_split_transcripts$gene_id,
#                             "ext_gene" = NA,
#                             "gene_biotype" = "novel",
#                             "transcript_biotype" = "novel")


## create the data frame with the transcript to gene_name/gene_id mapping
mapping_genes <- data.frame(matrix(ncol = 2, nrow = length(mapping_split_transcripts$gene_name)))

# rename cols
colnames(mapping_genes)[1] <- "target_id"
colnames(mapping_genes)[2] <- "gene_id"


# if there is no gene name --> transcript_id / gene_id
# if there is gene name --> transcript_id / gene_name
for (i in 1:length(mapping_split_transcripts$V1)){
  if (is.na(mapping_split_transcripts$gene_name[i])){
    mapping_genes$target_id[i] <- mapping_split_transcripts$transcript_id[i]
    mapping_genes$gene_id[i] <- mapping_split_transcripts$gene_id[i]
  }
  else {
    mapping_genes$target_id[i] <- mapping_split_transcripts$transcript_id[i]
    mapping_genes$gene_id[i] <- mapping_split_transcripts$gene_name[i]  
  }
}

# strip whitespace in front of target id and gene_id
mapping_genes$target_id <- trimws(mapping_genes$target_id, "l")
mapping_genes$gene_id <- trimws(mapping_genes$gene_id, "l")


# Add biotype column
mapping_genes$gene_biotype <- NA

# add biotypes from transcript mapping for fields that start with ENSG
for (i in 1:length(mapping_genes$gene_id)){
  if (startsWith(mapping_genes$gene_id[i], "ENSG")){
    # search index of gene in transcript mapping
    target_index <- which(mapping$ens_gene == mapping_genes$gene_id[i])
    mapping_genes$gene_biotype[i] <- mapping$gene_biotype[target_index[1]]
  }
  else {
    next
  }
}

# add biotypes from transcript mapping for fields where gene id is a real gene name
for (i in 1:length(mapping_genes$gene_id)){
  if (!startsWith(mapping_genes$gene_id[i], "ENSG") && !startsWith(mapping_genes$gene_id[i], "MSTR")){
    # search index of gene in transcript mapping
    target_index <- which(mapping$ext_gene == mapping_genes$gene_id[i])
    mapping_genes$gene_biotype[i] <- mapping$gene_biotype[target_index[1]]
  }
  else {
    next
  }
}



# Merge the novel mapping with the mapping from before
# gene_mapping <- mapping
# gene_mapping <- rbind(mapping, novel_mapping)


# Initialize Sleuth Object
#--> changed transformation function to log2FC
#--> aggregation column: is necessary if gene mode set to TRUE
sleuth_object_gene <- sleuth_prep(ref_df, extra_bootstrap_summary = TRUE, 
                             target_mapping = mapping_genes,
                             aggregation_column = "gene_id",
                             gene_mode = TRUE,
                             transformation_function = function(x) log2(x + 0.5))


# Fit full model
sleuth_object_gene <- sleuth_fit(sleuth_object_gene, ~condition, 'full')

#---- Error message
# 3 NA values were found during variance shrinkage estimation due to mean observation values outside of the range used for the LOESS fit.
# The LOESS fit will be repeated using exact computation of the fitted surface to extrapolate the missing values.
# These are the target ids with NA values:  GNAT1,  RNF152,  ZDHHC22
# computing variance of betas
#---

# Fit reduce model
sleuth_object_gene <- sleuth_fit(sleuth_object_gene, ~1, 'reduced')
#----
# 4 NA values were found during variance shrinkage estimation due to mean observation values outside of the range used for the LOESS fit.
# The LOESS fit will be repeated using exact computation of the fitted surface to extrapolate the missing values.
# These are the target ids with NA values:  GNAT1,  RNF152,  ZDHHC22,  MT-CO1
# computing variance of betas
#----



# Perform Sleuth test XXX --> we should do it with wt
# sleuth_object <- sleuth_lrt(sleuth_object, 'reduced', 'full')
sleuth_object_gene <- sleuth_wt(sleuth_object_gene, "conditionparental")

# Read sleuth data into a table
# Code from pachterlab
sleuth_gene_table <- sleuth_results(sleuth_object_gene, 'conditionparental', 'wt', show_all = TRUE)
sleuth_gene_significant <- dplyr::filter(sleuth_gene_table, qval <= 0.05)

# Create a new dataframe for results
sig_genes_table <- data.frame(gene_id = sleuth_gene_significant$target_id,
                              gene_biotype = sleuth_gene_significant$gene_biotype,
                                    log2FC = sleuth_gene_significant$b,
                                    qval = sleuth_gene_significant$qval)

#----


# Plot the data
plot_bootstrap(sleuth_object, "ENST00000371817.8", units = "est_counts", color_by = "condition")

# Volcano plot with Enhanced Volcano
#EnhancedVolcano()


# Volcano plot with sleuth 
plot_volcano(sleuth_object, test = "conditionparental", test_type = "wt", which_model = "full",
             sig_level = 0.1, point_alpha = 0.2, sig_color = "red",
             highlight = NULL)












#----
### Create New dataframe with necessary information for bedfile
lnc_RNA_data <- dplyr::filter(sleuth_significant, gene_biotype == "lncRNA")
novel_transcripts <- dplyr::filter(sleuth_significant, is.na(ens_gene))


dif_exp_results <- sleuth_significant$ 
dplyr::filter(gtf_df,  )


# Make 2 bed files
# --> one with ENST
# --> one with novel transcripts

no_exon_gtf <- filter(gtf_df, V3 == "transcript")
# select columns: CHR, star, end, transcript id, strand
bed_precursor <- select(no_exon_gtf, "V1","V3","V4","V5","V7","V9")
# Split last column with many attributes separated by ";"
bed_precursor_split <- separate_wider_delim(bed_precursor, V9, ";", names = c("gene_id","transcript_id"), too_many = "merge")
pre_bed_split_transcripts <- separate_wider_delim(bed_precursor_split, transcript_id, ";", names = c("transcript_id", "rest"), too_many = "merge")
pre_bed_split_transcripts <- separate_wider_delim(pre_bed_split_transcripts, rest, ";", names = c("gene_name", "ref_gene_id"), too_many = "merge", too_few = "align_end")

# rename other columns
colnames(pre_bed_split_transcripts)[1] <- "Chr"
colnames(pre_bed_split_transcripts)[3] <- "start"
colnames(pre_bed_split_transcripts)[4] <- "end"
colnames(pre_bed_split_transcripts)[5] <- "strand"


# remove "gene_id" from gene_id column
pre_bed_split_transcripts$gene_id <- gsub("gene_id ", "", pre_bed_split_transcripts$gene_id)
# remove "transcript_id" from transcript_id column
pre_bed_split_transcripts$transcript_id <- gsub("transcript_id ", "", pre_bed_split_transcripts$transcript_id)
# remove "gene_name" from gene_name column
pre_bed_split_transcripts$gene_name <- gsub("gene_name ", "", pre_bed_split_transcripts$gene_name)
# remove "ref_gene_id" from ref_gene_id column
pre_bed_split_transcripts$ref_gene_id <- gsub("ref_gene_id ", "", pre_bed_split_transcripts$ref_gene_id)
pre_bed_split_transcripts$ref_gene_id <- gsub(";", "", pre_bed_split_transcripts$ref_gene_id)

# split the dataframe in one for + stranded and one for - stranded
plus_pre_bed <- pre_bed_split_transcripts[pre_bed_split_transcripts$strand == "+",]
minus_pre_bed_inversed <- pre_bed_split_transcripts[pre_bed_split_transcripts$strand == "-",]

# if we got transcripts on - strand we need to switch start and stop position
# reorder columns in general
minus_pre_bed <- minus_pre_bed_inversed[,c(1,4,3,7,5,6,8,9,2)]

# Reorder Columns for plus_pre_bed
plus_pre_bed <- plus_pre_bed[,c(1,3,4,7,5,6,8,9,2)]


# Create a dataframe with only novel/known transcripts for + and - stranded transcripts
plus_pre_bed_novel <- plus_pre_bed[is.na(plus_pre_bed$gene_name),]
plus_pre_bed_known <- plus_pre_bed[!is.na(plus_pre_bed$gene_name),]

minus_pre_bed_novel <- minus_pre_bed[is.na(minus_pre_bed$gene_name),]
minus_pre_bed_known <- minus_pre_bed[!is.na(minus_pre_bed$gene_name),]

# rename cols from minus pre bed for merging later
#colnames(minus_pre_bed_novel)[colnames(minus_pre_bed_novel) == "start"] <- "stop"
#colnames(minus_pre_bed_novel)[colnames(minus_pre_bed_novel) == "end"] <- "start"
#colnames(minus_pre_bed_novel)[colnames(minus_pre_bed_novel) == "stop"] <- "end"




# merge plus/minus dfs for each novel and known transcripts in a df
full_bed_novel <- merge(plus_pre_bed_novel, minus_pre_bed_novel, all = T)
full_bed_known <- merge(plus_pre_bed_known, minus_pre_bed_known, all = T)


bed_path_novel <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/novel.bed"
bed_path_known <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/known.bed"

# save bed files for upload later
write.table(full_bed_novel, file = bed_path_novel, sep = "\t", quote = FALSE,
            col.names = FALSE, row.names = FALSE)
write.table(full_bed_known, file = bed_path_known, sep = "\t", quote = FALSE,
            col.names = FALSE, row.names = FALSE)


### Create Bed files for the 5' end and 3' end
# we want a window of +- 50bp
# --> start for - stranded transcripts is end??

five_prime_bed <- full_bed_novel 

# rename relevant cols
colnames(five_prime_bed)[colnames(five_prime_bed) == "start"] <- "start -50"
colnames(five_prime_bed)[colnames(five_prime_bed) == "end"] <- "start +50"

# replace start and end values with start +- 50
for (index in  1:nrow(five_prime_bed)){
  if (five_prime_bed[index, 5] == "+"){
    start <- five_prime_bed[index, 2]
    five_prime_bed[index, 2] <- start - 50
    five_prime_bed[index, 3] <- start + 50
  }
  else {
    start <- five_prime_bed[index, 3]
    five_prime_bed[index, 2] <- start - 50
    five_prime_bed[index, 3] <- start + 50
        
  }
} 

## Do the same for 3 prime

three_prime_bed <- full_bed_novel 

# rename relevant cols
colnames(three_prime_bed)[colnames(three_prime_bed) == "start"] <- "end -50"
colnames(three_prime_bed)[colnames(three_prime_bed) == "end"] <- "end +50"

# replace start and end values with end +- 50
for (index in  1:nrow(three_prime_bed)){
  if (three_prime_bed[index, 5] == "+"){
    end <- three_prime_bed[index, 3]
    three_prime_bed[index, 2] <- end - 50
    three_prime_bed[index, 3] <- end + 50
  }
  else {
    end <- three_prime_bed[index, 2]
    three_prime_bed[index, 2] <- end - 50
    three_prime_bed[index, 3] <- end + 50
    
  }
} 

# Check for negative values where we substracted 50
print(sum(five_prime_bed$`start -50` < 0))
#>0
print(sum(three_prime_bed$`end -50` < 0))
#>3

# exchange negative numbers for 0
three_prime_bed$`end -50`[three_prime_bed$`end -50` < 0] <- 0

# Check again for negative values where we subtracted 50
print(sum(three_prime_bed$`end -50` < 0))
#>0

## save data to bed file
bed_path_five_prime <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/five_prime.bed"
bed_path_three_prime <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/three_prime.bed"

# save bed files for upload later
write.table(five_prime_bed, file = bed_path_five_prime, sep = "\t", quote = FALSE,
            col.names = FALSE, row.names = FALSE)
write.table(three_prime_bed, file = bed_path_three_prime, sep = "\t", quote = FALSE,
            col.names = FALSE, row.names = FALSE)



## Make a light version of the bed files that has strand in 6th column
five_prime_light <- data.frame(matrix(ncol = 6, nrow = nrow(five_prime_bed)))
five_prime_light$X1 <- five_prime_bed$Chr
five_prime_light$X2 <- five_prime_bed$`start -50`
five_prime_light$X3 <- five_prime_bed$`start +50`
five_prime_light$X4 <- five_prime_bed$transcript_id
five_prime_light$X5 <- rep.int(1, nrow(five_prime_bed))
five_prime_light$X6 <- five_prime_bed$strand

three_prime_light <- data.frame(matrix(ncol = 6, nrow = nrow(three_prime_bed)))
three_prime_light$X1 <- three_prime_bed$Chr
three_prime_light$X2 <- three_prime_bed$`end -50`
three_prime_light$X3 <- three_prime_bed$`end +50`
three_prime_light$X4 <- three_prime_bed$transcript_id
three_prime_light$X5 <- rep.int(1, nrow(three_prime_bed))
three_prime_light$X6 <- three_prime_bed$strand

## save data to bed file
bed_path_five_prime_light <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/five_prime_light.bed"
bed_path_three_prime_light <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/three_prime_light.bed"

# save bed files for upload later
write.table(five_prime_light, file = bed_path_five_prime_light, sep = "\t", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(three_prime_light, file = bed_path_three_prime_light, sep = "\t", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)


## make a new bed that follows convention for bed_novel
bed_novel_light <- data.frame(matrix(ncol = 6, nrow = nrow(full_bed_novel)))
bed_novel_light$X1 <- full_bed_novel$Chr
bed_novel_light$X2 <- full_bed_novel$start
bed_novel_light$X3 <- full_bed_novel$end
bed_novel_light$X4 <- paste0('"', trimws(full_bed_novel$transcript_id, "l"),'"')
bed_novel_light$X5 <- rep.int(1, nrow(full_bed_novel))
bed_novel_light$X6 <- full_bed_novel$strand

## save data to bed file
bed_path_novel_light <- "C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/novel_light.bed"

# save bed files for upload later
write.table(bed_novel_light, file = bed_path_novel_light, sep = "\t", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)



# Read kallisto output for parental lines
#p1_l3 = read.table('C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/P1_L3/abundance.tsv', sep="\t", header=TRUE, as.is=1)
#p2_l3 = read.table('C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/P2_L3/abundance.tsv', sep="\t", header=TRUE, as.is=1)
#p3_l3 = read.table('C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/P3_L3/abundance.tsv', sep="\t", header=TRUE, as.is=1)

# Read kallisto output for holo lines
#s1_l3 = read.table('C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/S1_1_L3/abundance.tsv', sep="\t", header=TRUE, as.is=1)
#s2_l3 = read.table('C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/S1_2_L3/abundance.tsv', sep="\t", header=TRUE, as.is=1)
#s3_l3 = read.table('C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/S1_3_L3/abundance.tsv', sep="\t", header=TRUE, as.is=1)


# t2g target mappimg


loaded_packages <- installed.packages()[, "Package"]
print(loaded_packages)

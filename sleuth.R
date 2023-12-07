# Install rhdf5
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rhdf5")
library("rhdf5")
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)


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


# Initialize Sleuth Object
sleuth_object <- sleuth_prep(ref_df, extra_bootstrap_summary = TRUE)


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



# Perform Sleuth test
sleuth_object <- sleuth_lrt(sleuth_object, 'reduced', 'full')


# Read sleuth data into a table
# Code from pachterlab
sleuth_table <- sleuth_results(sleuth_object, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)


# Plot the data
plot_bootstrap(sleuth_object, "ENST00000371817.8", units = "est_counts", color_by = "condition")

# Volcano plot with Enhanced Volcano
EnhancedVolcano()


# Volcano plot with sleuth 
plot_volcano(sleuth_object,"LRT" , test_type = "wt", which_model = "full",
             sig_level = 0.1, point_alpha = 0.2, sig_color = "red",
             highlight = NULL)






# Read kallisto output for parental lines
#p1_l3 = read.table('C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/P1_L3/abundance.tsv', sep="\t", header=TRUE, as.is=1)
#p2_l3 = read.table('C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/P2_L3/abundance.tsv', sep="\t", header=TRUE, as.is=1)
#p3_l3 = read.table('C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/P3_L3/abundance.tsv', sep="\t", header=TRUE, as.is=1)

# Read kallisto output for holo lines
#s1_l3 = read.table('C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/S1_1_L3/abundance.tsv', sep="\t", header=TRUE, as.is=1)
#s2_l3 = read.table('C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/S1_2_L3/abundance.tsv', sep="\t", header=TRUE, as.is=1)
#s3_l3 = read.table('C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/S1_3_L3/abundance.tsv', sep="\t", header=TRUE, as.is=1)


# t2g target mappimg


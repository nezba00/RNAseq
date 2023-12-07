# Install rhdf5
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rhdf5")


# Install sleuth
install.packages("devtools")
devtools::install_github("pachterlab/sleuth")

# Load Sleuth
library(sleuth)

# Read kallisto output for parental lines
p1_l3 = read.table('C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/P1_L3/abundance.tsv', sep="\t", header=FALSE, as.is=1)
p2_l3 = read.table('C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/P2_L3/abundance.tsv', sep="\t", header=FALSE, as.is=1)
p3_l3 = read.table('C:/Users/nbaho/OneDrive/Desktop/Bioinf/RNA_Seq_Project/P3_L3/abundance.tsv', sep="\t", header=FALSE, as.is=1)

# Read kallisto output for holo lines




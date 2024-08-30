################################################################################
### Commands for server 

# cd /storage09/nawan2/Collapse_CpGs_IDwise/
# ls
# head g11.XX_C.3_2_23_0_3M.CpG_report.merged_CpG_evidence.cov

# (base) [nawan2@hodor01] (4)$ head g11.XX_C.3_2_23_0_3M.CpG_report.merged_CpG_evidence.cov
# 1       3050095 3050096 75.000000       18      6
# 1       3050195 3050196 87.878788       29      4
# 1       3050223 3050224 85.714286       30      5
# 1       3050262 3050263 54.838710       17      14
# 1       3050331 3050332 24.390244       10      31
# 1       3050403 3050404 76.086957       35      11
# 1       3050556 3050557 77.192982       44      13
# 1       3050666 3050667 65.217391       30      16
# 1       3051031 3051032 54.794521       40      33
# 1       3051200 3051201 72.058824       49      19

# NOTE:
# The coverage output looks like this (tab-delimited; 1-based genomic coords): 
#   <chromosome>  <start position>  <end position>  <methylation percentage>  <count methylated>  <count unmethylated>

# To see the number of lines (rows) in the .cov file
# wc -l your_file.cov
# wc -l g11.XX_C.3_2_23_0_3M.CpG_report.merged_CpG_evidence.cov
# 21124929 (result)

# # Creating a new directory called DMR_analysis_DSS
# cd /storage09/nawan2/
# mkdir DMR_analysis_DSS
# cd /storage09/nawan2/DMR_analysis_DSS
# pwd

# # To see the list of all screen sessions
# screen -ls

# # Creating a screen session
screen -S dmr
screen -r 254344.dmr
TMOUT=$((96 * 3600))
export TMOUT
# 
# # Start R:
# # conda activate r_env
# # R
# 
# screen -X -S 254344.dmr quit
# (Done!)



################################################################################


################################################################################
### Converting .cov files to .txt files to make a BS object

library(DSS)
?makeBSseqData
# The colnames of dat MUST BE "chr", "pos", "N", "X".
# For example, see head(dat1.1) in the example codes below.

# # EXAMPLE CODES THAT ARE GOING TO BE USEFUL FOR US:

# # Creating BS object:
# library(DSS)
# require(bsseq)
# path = file.path(system.file(package="DSS"), "extdata")
# dat1.1 = read.table(file.path(path, "cond1_1.txt"), header=TRUE)
# dat1.2 = read.table(file.path(path, "cond1_2.txt"), header=TRUE)
# dat2.1 = read.table(file.path(path, "cond2_1.txt"), header=TRUE)
# dat2.2 = read.table(file.path(path, "cond2_2.txt"), header=TRUE)
# BSobj = makeBSseqData( list(dat1.1, dat1.2, dat2.1, dat2.2),
#                        c("C1","C2", "N1", "N2") )[1:1000,]
# BSobj
# ## An object of type 'BSseq' with
# ##   1000 methylation loci
# ##   4 samples
# ## has not been smoothed
# ## All assays are in-memory

# # Creating BS design:
# Strain = rep(c("A", "B", "C"), 4)
# Sex = rep(c("M", "F"), each=6)
# design = data.frame(Strain,Sex)
# design

# Setting the current working directory to /storage09/nawan2/DMR_analysis_DSS
setwd("/storage09/nawan2/DMR_analysis_DSS")

# List of files
file_paths <- list.files("/storage09/nawan2/Collapse_CpGs_IDwise", pattern = "\\.cov$", full.names = TRUE)

# Function to process each file
process_cov_file <- function(file_path) {
  # Read the .cov file into R
  cov_data <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)
  
  # Rename the columns
  colnames(cov_data) <- c("chr", "start_position", "end_position", "methylation_percentage", "count_methylated", "count_unmethylated")
  
  # Creating new columns to match the requirement of the makeBSseqData()
  cov_data$pos <- cov_data$start_position
  cov_data$N <- cov_data$count_methylated + cov_data$count_unmethylated
  cov_data$X <- cov_data$count_methylated
  # chr=(1) Chromosome number; 
  # pos=(2) Genomic coordinates; 
  # N=(3) Read coverage of the position from BS-seq data; 
  # X=(4) Number of reads showing methylation of the position
  
  # Select required columns
  cov_data <- cov_data[, c("chr", "pos", "N", "X")]
  
  # Write the processed data to a tab-separated .txt file
  output_file <- gsub("\\.cov$", ".txt", basename(file_path))
  # Test: gsub("\\.cov$", ".txt", basename("/storage09/nawan2/Collapse_CpGs_IDwise/g11.XX_C.3_2_23_0_3M.CpG_report.merged_CpG_evidence.cov"))
  write.table(cov_data, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# Apply the function to each file
lapply(file_paths, process_cov_file)
################################################################################


################################################################################
### Creating the BSseqData object

# Setting the current working directory to /storage09/nawan2/DMR_analysis_DSS
setwd("/storage09/nawan2/DMR_analysis_DSS")

library(DSS)
require(bsseq)


### Make loading the following datasets faster:

# library(data.table)
# g43 <- fread("g43.XY_T.4_20_23_1_1M.CpG_report.merged_CpG_evidence.txt", header = TRUE)

# library(future.apply)
# plan(multicore)
# loaded_data <- future_lapply(list_of_files, function(file) fread(file))

g11 = read.table("g11.XX_C.3_2_23_0_3M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
# dim(g11)
# dim(na.omit(g11))
# > dim(g11)
# [1] 21124929        4
# > dim(na.omit(g11))
# [1] 21124929        4

g12 = read.table("g12.XX_C.3_23_23_1_3M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
# dim(g12)
# dim(na.omit(g12))
# > dim(g12)
# [1] 21129473        4
# > dim(na.omit(g12))
# [1] 21129473        4

g14 = read.table("g14.XX_C.2_23_23_0_5M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
# dim(g14)
# dim(na.omit(g14))
# > dim(g14)
# [1] 21072671        4
# > dim(na.omit(g14))
# [1] 21072671        4

g15 = read.table("g15.XX_C.6_24_23_0_2M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
# dim(g15)
# dim(na.omit(g15))
# > dim(g15)
# [1] 21077636        4
# > dim(na.omit(g15))
# [1] 21077636        4

g21 = read.table("g21.XX_T.4_24_23_0_2M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
# dim(g21)
# dim(na.omit(g21))
# > dim(g21)
# [1] 21098115        4
# > dim(na.omit(g21))
# [1] 21098115        4

g22 = read.table("g22.XX_T.4_24_23_0_4M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
# dim(g22)
# dim(na.omit(g22))
# > dim(g22)
# [1] 21065440        4
# > dim(na.omit(g22))
# [1] 21065440        4

g23 = read.table("g23.XX_T.4_20_23_0_3M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
# dim(g23)
# dim(na.omit(g23))
# > dim(g23)
# [1] 21029834        4
# > dim(na.omit(g23))
# [1] 21029834        4

g24 = read.table("g24.XX_T.4_20_23_1_3M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
# dim(g24)
# dim(na.omit(g24))
# > dim(g24)
# [1] 21103756        4
# > dim(na.omit(g24))
# [1] 21103756        4

g25 = read.table("g25.XX_T.4_06_23_0_3M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
# dim(g25)
# dim(na.omit(g25))
# > dim(g25)
# [1] 21116991        4
# > dim(na.omit(g25))
# [1] 21116991        4

g31 = read.table("g31.XY_C.4_24_23_1_4M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
# dim(g31)
# dim(na.omit(g31))
# > dim(g31)
# [1] 21332874        4
# > dim(na.omit(g31))
# [1] 21332874        4

g32 = read.table("g32.XY_C.1_18_23_0_4M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
# dim(g32)
# dim(na.omit(g32))
# > dim(g32)
# [1] 21313517        4
# > dim(na.omit(g32))
# [1] 21313517        4

g33 = read.table("g33.XY_C.3_9_23_0_5M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
# dim(g33)
# dim(na.omit(g33))
# > dim(g33)
# [1] 21295370        4
# > dim(na.omit(g33))
# [1] 21295370        4

g34 = read.table("g34.XY_C.6_24_23_0_4M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
dim(g34)
dim(na.omit(g34))
# > dim(g34)
# [1] 21302283        4
# > dim(na.omit(g34))
# [1] 21302283        4

g35 = read.table("g35.XY_C.6_24_23_0_3M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
# dim(g35)
# dim(na.omit(g35))
# > dim(g35)
# [1] 21283068        4
# > dim(na.omit(g35))
# [1] 21283068        4

g41 = read.table("g41.XY_T.4_24_23_0_5M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
# dim(g41)
# dim(na.omit(g41))
# > dim(g41)
# [1] 21281403        4
# > dim(na.omit(g41))
# [1] 21281403        4

g43 = read.table("g43.XY_T.4_20_23_1_1M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
# dim(g43)
# dim(na.omit(g43))
# > dim(g43)
# [1] 21311941        4
# > dim(na.omit(g43))
# [1] 21311941        4

g44 = read.table("g44.XY_T.4_24_23_0_8M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
# dim(g44)
# dim(na.omit(g44))
# > dim(g44)
# [1] 21269593        4
# > dim(na.omit(g44))
# [1] 21269593        4

g45 = read.table("g45.XY_T.6_24_23_1_2M.CpG_report.merged_CpG_evidence.txt", header=TRUE)
# dim(g45)
# dim(na.omit(g45))
# > dim(g45)
# [1] 21350393        4
# > dim(na.omit(g45))
# [1] 21350393        4

################################################################################

### Datasets with only common CpGs across the samples:

### Important step:
### Create a new column "chr_pos" by combining "chr" and "pos" columns
g11$chr_pos <- with(g11, paste("chr", chr, "_pos", pos, sep=""))
head(g11)
g12$chr_pos <- with(g12, paste("chr", chr, "_pos", pos, sep=""))
g14$chr_pos <- with(g14, paste("chr", chr, "_pos", pos, sep=""))
g15$chr_pos <- with(g15, paste("chr", chr, "_pos", pos, sep=""))

g21$chr_pos <- with(g21, paste("chr", chr, "_pos", pos, sep=""))
g22$chr_pos <- with(g22, paste("chr", chr, "_pos", pos, sep=""))
g23$chr_pos <- with(g23, paste("chr", chr, "_pos", pos, sep=""))
g24$chr_pos <- with(g24, paste("chr", chr, "_pos", pos, sep=""))
g25$chr_pos <- with(g25, paste("chr", chr, "_pos", pos, sep=""))

g31$chr_pos <- with(g31, paste("chr", chr, "_pos", pos, sep=""))
g32$chr_pos <- with(g32, paste("chr", chr, "_pos", pos, sep=""))
g33$chr_pos <- with(g33, paste("chr", chr, "_pos", pos, sep=""))
g34$chr_pos <- with(g34, paste("chr", chr, "_pos", pos, sep=""))
g35$chr_pos <- with(g35, paste("chr", chr, "_pos", pos, sep=""))

g41$chr_pos <- with(g41, paste("chr", chr, "_pos", pos, sep=""))
g43$chr_pos <- with(g43, paste("chr", chr, "_pos", pos, sep=""))
g44$chr_pos <- with(g44, paste("chr", chr, "_pos", pos, sep=""))
g45$chr_pos <- with(g45, paste("chr", chr, "_pos", pos, sep=""))

### Finding common CpGs across all samples

common_CpGs <- Reduce(intersect,(list(g11$chr_pos, g12$chr_pos, g14$chr_pos, g15$chr_pos,
              g21$chr_pos, g22$chr_pos, g23$chr_pos, g24$chr_pos, g25$chr_pos,
              g31$chr_pos, g32$chr_pos, g33$chr_pos, g34$chr_pos, g35$chr_pos,
              g41$chr_pos, g43$chr_pos, g44$chr_pos, g45$chr_pos)))

length(common_CpGs)
# 20786908

autosome_chr <- as.character(1:19)
class(autosome_chr)

g11_common <- subset(g11, chr_pos %in% common_CpGs)
g11_common <- subset(g11_common, select = -c(chr_pos))
dim(g11_common)
# [1] 20786908        4
g11_common_autosomes <- subset(g11_common, chr %in% autosome_chr) 
dim(g11_common_autosomes)
# [1] 19895623        4
# Save data frame as CSV file in the current working directory
write.csv(g11_common_autosomes, file = "g11_common_autosomes.csv", row.names = FALSE)


g12_common <- subset(g12, chr_pos %in% common_CpGs)
g12_common <- subset(g12_common, select = -c(chr_pos))
g12_common_autosomes <- subset(g12_common, chr %in% autosome_chr) 
dim(g12_common_autosomes)
# [1] 19895623        4
write.csv(g12_common_autosomes, file = "g12_common_autosomes.csv", row.names = FALSE)

g14_common <- subset(g14, chr_pos %in% common_CpGs)
g14_common <- subset(g14_common, select = -c(chr_pos))
g14_common_autosomes <- subset(g14_common, chr %in% autosome_chr) 
dim(g14_common_autosomes)
# [1] 19895623        4
write.csv(g14_common_autosomes, file = "g14_common_autosomes.csv", row.names = FALSE)

g15_common <- subset(g15, chr_pos %in% common_CpGs)
g15_common <- subset(g15_common, select = -c(chr_pos))
g15_common_autosomes <- subset(g15_common, chr %in% autosome_chr) 
dim(g15_common_autosomes)
# [1] 19895623        4
write.csv(g15_common_autosomes, file = "g15_common_autosomes.csv", row.names = FALSE)

g21_common <- subset(g21, chr_pos %in% common_CpGs)
g21_common <- subset(g21_common, select = -c(chr_pos))
g21_common_autosomes <- subset(g21_common, chr %in% autosome_chr) 
dim(g21_common_autosomes)
# [1] 19895623        4
write.csv(g21_common_autosomes, file = "g21_common_autosomes.csv", row.names = FALSE)

g22_common <- subset(g22, chr_pos %in% common_CpGs)
g22_common <- subset(g22_common, select = -c(chr_pos))
g22_common_autosomes <- subset(g22_common, chr %in% autosome_chr) 
dim(g22_common_autosomes)
# [1] 19895623        4
write.csv(g22_common_autosomes, file = "g22_common_autosomes.csv", row.names = FALSE)

g23_common <- subset(g23, chr_pos %in% common_CpGs)
g23_common <- subset(g23_common, select = -c(chr_pos))
g23_common_autosomes <- subset(g23_common, chr %in% autosome_chr) 
dim(g23_common_autosomes)
# [1] 19895623        4
write.csv(g23_common_autosomes, file = "g23_common_autosomes.csv", row.names = FALSE)

g24_common <- subset(g24, chr_pos %in% common_CpGs)
g24_common <- subset(g24_common, select = -c(chr_pos))
g24_common_autosomes <- subset(g24_common, chr %in% autosome_chr) 
dim(g24_common_autosomes)
# [1] 19895623        4
write.csv(g24_common_autosomes, file = "g24_common_autosomes.csv", row.names = FALSE)

g25_common <- subset(g25, chr_pos %in% common_CpGs)
g25_common <- subset(g25_common, select = -c(chr_pos))
g25_common_autosomes <- subset(g25_common, chr %in% autosome_chr) 
dim(g25_common_autosomes)
# [1] 19895623        4
write.csv(g25_common_autosomes, file = "g25_common_autosomes.csv", row.names = FALSE)

g31_common <- subset(g31, chr_pos %in% common_CpGs)
g31_common <- subset(g31_common, select = -c(chr_pos))
g31_common_autosomes <- subset(g31_common, chr %in% autosome_chr) 
dim(g31_common_autosomes)
# [1] 19895623        4
write.csv(g31_common_autosomes, file = "g31_common_autosomes.csv", row.names = FALSE)

g32_common <- subset(g32, chr_pos %in% common_CpGs)
g32_common <- subset(g32_common, select = -c(chr_pos))
g32_common_autosomes <- subset(g32_common, chr %in% autosome_chr) 
dim(g32_common_autosomes)
# [1] 19895623        4
write.csv(g32_common_autosomes, file = "g32_common_autosomes.csv", row.names = FALSE)

g33_common <- subset(g33, chr_pos %in% common_CpGs)
g33_common <- subset(g33_common, select = -c(chr_pos))
g33_common_autosomes <- subset(g33_common, chr %in% autosome_chr) 
dim(g33_common_autosomes)
# [1] 19895623        4
write.csv(g33_common_autosomes, file = "g33_common_autosomes.csv", row.names = FALSE)

g34_common <- subset(g34, chr_pos %in% common_CpGs)
g34_common <- subset(g34_common, select = -c(chr_pos))
g34_common_autosomes <- subset(g34_common, chr %in% autosome_chr) 
dim(g34_common_autosomes)
# [1] 19895623        4
write.csv(g34_common_autosomes, file = "g34_common_autosomes.csv", row.names = FALSE)

g35_common <- subset(g35, chr_pos %in% common_CpGs)
g35_common <- subset(g35_common, select = -c(chr_pos))
g35_common_autosomes <- subset(g35_common, chr %in% autosome_chr) 
dim(g35_common_autosomes)
# [1] 19895623        4
write.csv(g35_common_autosomes, file = "g35_common_autosomes.csv", row.names = FALSE)

g41_common <- subset(g41, chr_pos %in% common_CpGs)
g41_common <- subset(g41_common, select = -c(chr_pos))
g41_common_autosomes <- subset(g41_common, chr %in% autosome_chr) 
dim(g41_common_autosomes)
# [1] 19895623        4
write.csv(g41_common_autosomes, file = "g41_common_autosomes.csv", row.names = FALSE)

g43_common <- subset(g43, chr_pos %in% common_CpGs)
g43_common <- subset(g43_common, select = -c(chr_pos))
g43_common_autosomes <- subset(g43_common, chr %in% autosome_chr) 
dim(g43_common_autosomes)
# [1] 19895623        4
write.csv(g43_common_autosomes, file = "g43_common_autosomes.csv", row.names = FALSE)

g44_common <- subset(g44, chr_pos %in% common_CpGs)
g44_common <- subset(g44_common, select = -c(chr_pos))
g44_common_autosomes <- subset(g44_common, chr %in% autosome_chr) 
dim(g44_common_autosomes)
# [1] 19895623        4
write.csv(g44_common_autosomes, file = "g44_common_autosomes.csv", row.names = FALSE)

g45_common <- subset(g45, chr_pos %in% common_CpGs)
g45_common <- subset(g45_common, select = -c(chr_pos))
g45_common_autosomes <- subset(g45_common, chr %in% autosome_chr) 
dim(g45_common_autosomes)
# [1] 19895623        4
write.csv(g45_common_autosomes, file = "g45_common_autosomes.csv", row.names = FALSE)

################################################################################

### Checking frequency of chromosones (and names) in each sample

# List of dataset names
datasets <- c("g11", "g12", "g14", "g15",
              "g21", "g22", "g23", "g24", "g25",
              "g31", "g32", "g33", "g34", "g35",
              "g41", "g43", "g44", "g45")

# Load necessary libraries
library(data.table)
library(openxlsx)

# Create a new Excel workbook
wb <- createWorkbook()

# Loop through each dataset
for (dataset in datasets) {
  # Read the dataset
  data <- get(dataset)
  
  # Calculate frequency table for "chr" variable
  freq <- data.table(table(data$chr))
  
  # Rename columns
  setnames(freq, c("chr", "Freq"))
  
  # Add a new worksheet to the workbook
  addWorksheet(wb, sheetName = dataset)
  
  # Write frequency table to the worksheet
  writeDataTable(wb, sheet = dataset, x = freq, startCol = 1, startRow = 1)
}

# Save the workbook to a single Excel file
saveWorkbook(wb, file = "frequency_tables.xlsx", overwrite = TRUE)

################################################################################

### Checking frequency of chromosones (and names) in each sample with only common CpGs across all samples

# List of dataset names
datasets <- c("g11_common", "g12_common", "g14_common", "g15_common",
              "g21_common", "g22_common", "g23_common", "g24_common", "g25_common",
              "g31_common", "g32_common", "g33_common", "g34_common", "g35_common",
              "g41_common", "g43_common", "g44_common", "g45_common")

# Load necessary libraries
library(data.table)
library(openxlsx)

# Create a new Excel workbook
wb <- createWorkbook()

# Loop through each dataset
for (dataset in datasets) {
  # Read the dataset
  data <- get(dataset)
  
  # Calculate frequency table for "chr" variable
  freq <- data.table(table(data$chr))
  
  # Rename columns
  setnames(freq, c("chr", "Freq"))
  
  # Add a new worksheet to the workbook
  addWorksheet(wb, sheetName = dataset)
  
  # Write frequency table to the worksheet
  writeDataTable(wb, sheet = dataset, x = freq, startCol = 1, startRow = 1)
}

# Save the workbook to a single Excel file
saveWorkbook(wb, file = "frequency_tables_common_CpGs.xlsx", overwrite = TRUE)

################################################################################

################################################################################

### Frequency tables for 
# N<=1
# N<=2
# N<=3
# N<=4
# N<=5
# N<=10
# N<=15
# N<=20
# N<=30
# N<=40
# N<=50
# for each sample

# Load necessary libraries
library(data.table)
library(openxlsx)

# List of dataset names
datasets <- c("g11", "g12", "g14", "g15",
              "g21", "g22", "g23", "g24", "g25",
              "g31", "g32", "g33", "g34", "g35",
              "g41", "g43", "g44", "g45")

# Function to calculate frequency table for "N" variable
calculate_freq <- function(data) {
  thresholds <- c(1, 2, 3, 4, 5, 10, 15, 20, 30, 40, 50)
  freq <- sapply(thresholds, function(threshold) {
    sum(data$N <= threshold)
  })
  thresholds <- c(paste("N<=", thresholds))
  percent <- round(freq / nrow(data) * 100, 2)  # Adjusted calculation for percentage
  
  return(data.frame(N = thresholds, Freq = freq, Percent = percent))
}

# Create a new Excel workbook
wb <- createWorkbook()

# Loop through each dataset
for (dataset in datasets) {
  # Read the dataset
  data <- get(dataset)
  
  # Calculate frequency table for "N" variable
  freq <- calculate_freq(data)
  
  # Add a new worksheet to the workbook
  addWorksheet(wb, sheetName = dataset)
  
  # Write frequency table to the worksheet
  writeData(wb, sheet = dataset, freq, startCol = 1, startRow = 1)
}

# Save the workbook to a single Excel file
saveWorkbook(wb, file = "frequency_tables_N_thresholds_adjusted.xlsx", overwrite = TRUE)

################################################################################

################################################################################

### Frequency tables for 
# N<=1
# N<=2
# N<=3
# N<=4
# N<=5
# N<=10
# N<=15
# N<=20
# N<=30
# N<=40
# N<=50
# for each sample with datasets with common CpGs

# Load necessary libraries
library(data.table)
library(openxlsx)

# List of dataset names
datasets <- c("g11_common", "g12_common", "g14_common", "g15_common",
              "g21_common", "g22_common", "g23_common", "g24_common", "g25_common",
              "g31_common", "g32_common", "g33_common", "g34_common", "g35_common",
              "g41_common", "g43_common", "g44_common", "g45_common")

# Function to calculate frequency table for "N" variable
calculate_freq <- function(data) {
  thresholds <- c(1, 2, 3, 4, 5, 10, 15, 20, 30, 40, 50)
  freq <- sapply(thresholds, function(threshold) {
    sum(data$N <= threshold)
  })
  thresholds <- c(paste("N<=", thresholds))
  percent <- round(freq / nrow(data) * 100, 2)  # Adjusted calculation for percentage
  
  return(data.frame(N = thresholds, Freq = freq, Percent = percent))
}

# Create a new Excel workbook
wb <- createWorkbook()

# Loop through each dataset
for (dataset in datasets) {
  # Read the dataset
  data <- get(dataset)
  
  # Calculate frequency table for "N" variable
  freq <- calculate_freq(data)
  
  # Add a new worksheet to the workbook
  addWorksheet(wb, sheetName = dataset)
  
  # Write frequency table to the worksheet
  writeData(wb, sheet = dataset, freq, startCol = 1, startRow = 1)
}

# Save the workbook to a single Excel file
saveWorkbook(wb, file = "frequency_tables_N_common_CpGs.xlsx", overwrite = TRUE)

################################################################################

### Creating 18 new datasets with at least/minimum 10 reads (m10r)

g11_m10r <- subset(g11, N >= 10)
g12_m10r <- subset(g12, N >= 10)
g14_m10r <- subset(g14, N >= 10)
g15_m10r <- subset(g15, N >= 10)

g21_m10r <- subset(g21, N >= 10)
g22_m10r <- subset(g22, N >= 10)
g23_m10r <- subset(g23, N >= 10)
g24_m10r <- subset(g24, N >= 10)
g25_m10r <- subset(g25, N >= 10)

g31_m10r <- subset(g31, N >= 10)
g32_m10r <- subset(g32, N >= 10)
g33_m10r <- subset(g33, N >= 10)
g34_m10r <- subset(g34, N >= 10)
g35_m10r <- subset(g35, N >= 10)

g41_m10r <- subset(g41, N >= 10)
g43_m10r <- subset(g43, N >= 10)
g44_m10r <- subset(g44, N >= 10)
g45_m10r <- subset(g45, N >= 10)

################################################################################

### Making BS object with all CpGs

BSobj = makeBSseqData( list(g11, g12, g14, g15,
                            g21, g22, g23, g24, g25,
                            g31, g32, g33, g34, g35,
                            g41, g43, g44, g45),
                       c("3_2_23_0_3M","3_23_23_1_3M", "2_23_23_0_5M", "6_24_23_0_2M",
                         "4_24_23_0_2M", "4_24_23_0_4M", "4_20_23_0_3M", "4_20_23_1_3M", "4_06_23_0_3M",
                         "4_24_23_1_4M", "1_18_23_0_4M", "3_9_23_0_5M", "6_24_23_0_4M", "6_24_23_0_3M",
                         "4_24_23_0_5M", "4_20_23_1_1M", "4_24_23_0_8M", "6_24_23_1_2M") )
BSobj

# > BSobj
# An object of type 'BSseq' with
# 21609770 methylation loci
# 18 samples
# has not been smoothed
# All assays are in-memory

################################################################################

### Making BS object with CpGs that has at least 10 reads

BSobj_m10r = makeBSseqData( list(g11_m10r, g12_m10r, g14_m10r, g15_m10r,
                            g21_m10r, g22_m10r, g23_m10r, g24_m10r, g25_m10r,
                            g31_m10r, g32_m10r, g33_m10r, g34_m10r, g35_m10r,
                            g41_m10r, g43_m10r, g44_m10r, g45_m10r),
                       c("3_2_23_0_3M","3_23_23_1_3M", "2_23_23_0_5M", "6_24_23_0_2M",
                         "4_24_23_0_2M", "4_24_23_0_4M", "4_20_23_0_3M", "4_20_23_1_3M", "4_06_23_0_3M",
                         "4_24_23_1_4M", "1_18_23_0_4M", "3_9_23_0_5M", "6_24_23_0_4M", "6_24_23_0_3M",
                         "4_24_23_0_5M", "4_20_23_1_1M", "4_24_23_0_8M", "6_24_23_1_2M") )
BSobj_m10r
# Warning message:
#   In min(dd) : no non-missing arguments to min; returning Inf
# > BSobj_m10r
# An object of type 'BSseq' with
# 20872010 methylation loci
# 18 samples
# has not been smoothed
# All assays are in-memory

# Are there missing values in input data?

################################################################################



################################################################################

### Making BS object with data with only common CpGs across samples

BSobj_common = makeBSseqData( list(g11_common, g12_common, g14_common, g15_common,
                                 g21_common, g22_common, g23_common, g24_common, g25_common,
                                 g31_common, g32_common, g33_common, g34_common, g35_common,
                                 g41_common, g43_common, g44_common, g45_common),
                            c("3_2_23_0_3M","3_23_23_1_3M", "2_23_23_0_5M", "6_24_23_0_2M",
                              "4_24_23_0_2M", "4_24_23_0_4M", "4_20_23_0_3M", "4_20_23_1_3M", "4_06_23_0_3M",
                              "4_24_23_1_4M", "1_18_23_0_4M", "3_9_23_0_5M", "6_24_23_0_4M", "6_24_23_0_3M",
                              "4_24_23_0_5M", "4_20_23_1_1M", "4_24_23_0_8M", "6_24_23_1_2M") )
BSobj_common
# > BSobj_common
# An object of type 'BSseq' with
# 20786908 methylation loci
# 18 samples
# has not been smoothed
# All assays are in-memory

################################################################################



################################################################################
### Creating the design

type = rep(c("Control", "Test", "Control", "Test"), times=c(4, 5, 5, 4))
sex = rep(c("Female", "Male"), each=9)
design = data.frame(type, sex)
design
################################################################################


################################################################################
### Fit a linear model using DMLfit.multiFactor function
### For data without removing lower read counts

DMLfit = DMLfit.multiFactor(BSobj, design=design, formula=~type+sex+type:sex)
# > colnames(DMLfit$X)
# [1] "(Intercept)"      "typeTest"         "sexMale"          "typeTest:sexMale"

### Use DMLtest.multiFactor function to test the "type" effect
DMLtest.type = DMLtest.multiFactor(DMLfit, coef="typeTest")
DMLtest.type_COMPLETE <- na.omit(DMLtest.type) 
# Saving the resulting data frame containing chromosome number, CpG site position,
# test statistics, p-values (from normal distribution), and FDR.
# Rows are sorted by chromosome/position of the CpG sites
write.table(DMLtest.type, file = "DMLtest.for.type_sorted.by.chr.pos.txt", sep = "\t", row.names = FALSE)
# Saving the COMPLETE dataset
write.table(na.omit(DMLtest.type), file = "DMLtest.for.type_sorted.by.chr.pos_COMPLETE.txt", sep = "\t", row.names = FALSE)
# Obtain top ranked CpG sites, one can sort the data frame using following codes:
ix=sort(DMLtest.type[,"pvals"], index.return=TRUE)$ix
# head(DMLtest.type[ix,])
# Saving the resulting data frame containing chromosome number, CpG site position,
# test statistics, p-values (from normal distribution), and FDR.
# Rows are sorted by p-value of the CpG sites
write.table(DMLtest.type[ix,], file = "DMLtest.for.type_sorted.by.pval.txt", sep = "\t", row.names = FALSE)
# Saving the COMPLETE dataset
write.table(na.omit(DMLtest.type[ix,]), file = "DMLtest.for.type_sorted.by.pval_COMPLETE.txt", sep = "\t", row.names = FALSE)

### Use DMLtest.multiFactor function to test the "sex" effect
DMLtest.sex = DMLtest.multiFactor(DMLfit, coef="sexMale")
DMLtest.sex_COMPLETE <- na.omit(DMLtest.sex) 
# Saving the resulting data frame containing chromosome number, CpG site position,
# test statistics, p-values (from normal distribution), and FDR.
# Rows are sorted by chromosome/position of the CpG sites
write.table(DMLtest.sex, file = "DMLtest.for.sex_sorted.by.chr.pos.txt", sep = "\t", row.names = FALSE)
# Saving the COMPLETE dataset
write.table(na.omit(DMLtest.sex), file = "DMLtest.for.sex_sorted.by.chr.pos_COMPLETE.txt", sep = "\t", row.names = FALSE)
# Obtain top ranked CpG sites, one can sort the data frame using following codes:
ix=sort(DMLtest.sex[,"pvals"], index.return=TRUE)$ix
# head(DMLtest.type[ix,])
# Saving the resulting data frame containing chromosome number, CpG site position,
# test statistics, p-values (from normal distribution), and FDR.
# Rows are sorted by p-value of the CpG sites
write.table(DMLtest.sex[ix,], file = "DMLtest.for.sex_sorted.by.pval.txt", sep = "\t", row.names = FALSE)
# Saving the COMPLETE dataset
write.table(na.omit(DMLtest.sex[ix,]), file = "DMLtest.for.sex_sorted.by.pval_COMPLETE.txt", sep = "\t", row.names = FALSE)

### Use DMLtest.multiFactor function to test the "type:sex" interaction effect
DMLtest.interaction = DMLtest.multiFactor(DMLfit, coef="typeTest:sexMale")
DMLtest.interaction_COMPLETE <- na.omit(DMLtest.interaction) 
# Saving the resulting data frame containing chromosome number, CpG site position,
# test statistics, p-values (from normal distribution), and FDR.
# Rows are sorted by chromosome/position of the CpG sites
write.table(DMLtest.interaction, file = "DMLtest.for.interaction_sorted.by.chr.pos.txt", sep = "\t", row.names = FALSE)
# Saving the COMPLETE dataset
write.table(na.omit(DMLtest.interaction), file = "DMLtest.for.interaction_sorted.by.chr.pos_COMPLETE.txt", sep = "\t", row.names = FALSE)
# Obtain top ranked CpG sites, one can sort the data frame using following codes:
ix=sort(DMLtest.interaction[,"pvals"], index.return=TRUE)$ix
# head(DMLtest.type[ix,])
# Saving the resulting data frame containing chromosome number, CpG site position,
# test statistics, p-values (from normal distribution), and FDR.
# Rows are sorted by p-value of the CpG sites
write.table(DMLtest.interaction[ix,], file = "DMLtest.for.interaction_sorted.by.pval.txt", sep = "\t", row.names = FALSE)
# Saving the COMPLETE dataset
write.table(na.omit(DMLtest.interaction[ix,]), file = "DMLtest.for.interaction_sorted.by.pval_COMPLETE.txt", sep = "\t", row.names = FALSE)

# > dim(DMLtest.type)
# [1] 21609770        5
# > dim(na.omit(DMLtest.type))
# [1] 21144793        5
# 21609770-21144793
# > 21609770-21144793
# [1] 464977

################################################################################

### IMPORTANT: EXPLORING THE REASON/SOURCE OF MISSINGNESS  

# Subset the data frame to include only rows with any missing values
rows_with_na <- apply(is.na(DMLtest.type), 1, any)

# Subset the original data frame based on rows_with_na
DMLtest.type_with_na <- DMLtest.type[rows_with_na, ]

# Display the first few rows of the subset with missing values
head(DMLtest.type_with_na)

# > head(DMLtest.type_with_na)
# chr     pos stat pvals fdrs
# 2651   1 3570829   NA    NA   NA
# 2652   1 3570967   NA    NA   NA
# 2655   1 3571381   NA    NA   NA
# 2656   1 3571472   NA    NA   NA
# 2657   1 3571662   NA    NA   NA
# 2658   1 3571709   NA    NA   NA


################################################################################

# > dim(DMLtest.sex)
# [1] 21609770        5
# > dim(na.omit(DMLtest.sex))
# [1] 21144793        5

# > dim(DMLtest.interaction)
# [1] 21609770        5
# > dim(na.omit(DMLtest.interaction))
# [1] 21144793        5
################################################################################

################################################################################
### Fit a linear model using DMLfit.multiFactor function
### For data with common CpGs across samples (without smoothing)

DMLfit_common = DMLfit.multiFactor(BSobj_common, design=design, formula=~type+sex+type:sex)
# > colnames(DMLfit$X)
# [1] "(Intercept)"      "typeTest"         "sexMale"          "typeTest:sexMale"

### Use DMLtest.multiFactor function to test the "type" effect
DMLtest.type_common = DMLtest.multiFactor(DMLfit_common, coef="typeTest")
sum(DMLtest.type_common$fdrs < 0.05)
# 18808688
# There were 21609770  methylation loci
library(fdrtool)
pp.out.common <- fdrtool(DMLtest.type_common$pvals, statistic = "pvalue", plot = FALSE)
# > length(pp.out.common$lfdr)
# [1] 20786908
sum(pp.out.common$lfdr<0.05)
# 18416788

# Stats themselves
ss <- DMLtest.type_common$stat
# Z normalized stats
zz <- (ss - mean(ss)) / sd(ss)

ss.out <- fdrtool(ss, statistic = "normal", plot = FALSE)
zz.out <- fdrtool(zz, statistic = "normal", plot = FALSE)

# lFDR
sum(ss.out$lfdr<0.05)
sum(zz.out$lfdr<0.05)

# > sum(ss.out$lfdr<0.05)
# [1] 17679508
# > sum(zz.out$lfdr<0.05)
# [1] 539


DMLtest.type_common_COMPLETE <- na.omit(DMLtest.type_common) 
# Saving the resulting data frame containing chromosome number, CpG site position,
# test statistics, p-values (from normal distribution), and FDR.
# Rows are sorted by chromosome/position of the CpG sites
write.table(DMLtest.type_common, file = "DMLtest.for.type_common_sorted.by.chr.pos.txt", sep = "\t", row.names = FALSE)
# Saving the COMPLETE dataset
write.table(na.omit(DMLtest.type_common), file = "DMLtest.for.type_common_sorted.by.chr.pos_COMPLETE.txt", sep = "\t", row.names = FALSE)
# Obtain top ranked CpG sites, one can sort the data frame using following codes:
# ix=sort(DMLtest.type[,"pvals"], index.return=TRUE)$ix
# # head(DMLtest.type[ix,])
# # Saving the resulting data frame containing chromosome number, CpG site position,
# # test statistics, p-values (from normal distribution), and FDR.
# # Rows are sorted by p-value of the CpG sites
# write.table(DMLtest.type[ix,], file = "DMLtest.for.type_sorted.by.pval.txt", sep = "\t", row.names = FALSE)
# # Saving the COMPLETE dataset
# write.table(na.omit(DMLtest.type[ix,]), file = "DMLtest.for.type_sorted.by.pval_COMPLETE.txt", sep = "\t", row.names = FALSE)



### Use DMLtest.multiFactor function to test the "sex" effect
DMLtest.sex_common = DMLtest.multiFactor(DMLfit_common, coef="sexMale")
DMLtest.sex_common_COMPLETE <- na.omit(DMLtest.sex_common) 
# Saving the resulting data frame containing chromosome number, CpG site position,
# test statistics, p-values (from normal distribution), and FDR.
# Rows are sorted by chromosome/position of the CpG sites
write.table(DMLtest.sex_common, file = "DMLtest.for.sex_common_sorted.by.chr.pos.txt", sep = "\t", row.names = FALSE)
# Saving the COMPLETE dataset
write.table(na.omit(DMLtest.sex_common), file = "DMLtest.for.sex_common_sorted.by.chr.pos_COMPLETE.txt", sep = "\t", row.names = FALSE)
# # Obtain top ranked CpG sites, one can sort the data frame using following codes:
# ix=sort(DMLtest.sex[,"pvals"], index.return=TRUE)$ix
# # head(DMLtest.type[ix,])
# # Saving the resulting data frame containing chromosome number, CpG site position,
# # test statistics, p-values (from normal distribution), and FDR.
# # Rows are sorted by p-value of the CpG sites
# write.table(DMLtest.sex[ix,], file = "DMLtest.for.sex_sorted.by.pval.txt", sep = "\t", row.names = FALSE)
# # Saving the COMPLETE dataset
# write.table(na.omit(DMLtest.sex[ix,]), file = "DMLtest.for.sex_sorted.by.pval_COMPLETE.txt", sep = "\t", row.names = FALSE)

### Use DMLtest.multiFactor function to test the "type:sex" interaction effect
DMLtest.interaction_common = DMLtest.multiFactor(DMLfit_common, coef="typeTest:sexMale")
DMLtest.interaction_common_COMPLETE <- na.omit(DMLtest.interaction_common) 
# Saving the resulting data frame containing chromosome number, CpG site position,
# test statistics, p-values (from normal distribution), and FDR.
# Rows are sorted by chromosome/position of the CpG sites
write.table(DMLtest.interaction_common, file = "DMLtest.for.interaction_common_sorted.by.chr.pos.txt", sep = "\t", row.names = FALSE)
# Saving the COMPLETE dataset
write.table(na.omit(DMLtest.interaction_common), file = "DMLtest.for.interaction_common_sorted.by.chr.pos_COMPLETE.txt", sep = "\t", row.names = FALSE)
# # Obtain top ranked CpG sites, one can sort the data frame using following codes:
# ix=sort(DMLtest.interaction[,"pvals"], index.return=TRUE)$ix
# # head(DMLtest.type[ix,])
# # Saving the resulting data frame containing chromosome number, CpG site position,
# # test statistics, p-values (from normal distribution), and FDR.
# # Rows are sorted by p-value of the CpG sites
# write.table(DMLtest.interaction[ix,], file = "DMLtest.for.interaction_sorted.by.pval.txt", sep = "\t", row.names = FALSE)
# # Saving the COMPLETE dataset
# write.table(na.omit(DMLtest.interaction[ix,]), file = "DMLtest.for.interaction_sorted.by.pval_COMPLETE.txt", sep = "\t", row.names = FALSE)

# > dim(DMLtest.type)
# [1] 21609770        5
# > dim(na.omit(DMLtest.type))
# [1] 21144793        5
# 21609770-21144793
# > 21609770-21144793
# [1] 464977

################################################################################


################################################################################
### Fit a linear model using DMLfit.multiFactor function
### For data with common CpGs across samples (with smoothing)

DMLfit_common_smth = DMLfit.multiFactor(BSobj_common, design=design, formula=~type+sex+type:sex, smoothing=TRUE)
# > colnames(DMLfit$X)
# [1] "(Intercept)"      "typeTest"         "sexMale"          "typeTest:sexMale"

### Use DMLtest.multiFactor function to test the "type" effect
DMLtest.type_common_smth = DMLtest.multiFactor(DMLfit_common_smth, coef="typeTest")
sum(DMLtest.type_common_smth$fdrs < 0.05)
# 19118826
# There were 21609770  methylation loci
sum(is.na(DMLtest.type_common_smth$pvals))
# 0
library(fdrtool)
pp.out <- fdrtool(DMLtest.type_common_smth$pvals, statistic = "pvalue", plot = FALSE)
# > length(pp.out$lfdr)
# [1] 20786908
sum(pp.out$lfdr<0.05)
# 18767470

DMLtest.type_common_smth_COMPLETE <- na.omit(DMLtest.type_common_smth) 
# Saving the resulting data frame containing chromosome number, CpG site position,
# test statistics, p-values (from normal distribution), and FDR.
# Rows are sorted by chromosome/position of the CpG sites
write.table(DMLtest.type_common_smth, file = "DMLtest.for.type_common_smth_sorted.by.chr.pos.txt", sep = "\t", row.names = FALSE)
# Saving the COMPLETE dataset
write.table(na.omit(DMLtest.type_common_smth), file = "DMLtest.for.type_common_smth_sorted.by.chr.pos_COMPLETE.txt", sep = "\t", row.names = FALSE)
# Obtain top ranked CpG sites, one can sort the data frame using following codes:
# ix=sort(DMLtest.type[,"pvals"], index.return=TRUE)$ix
# # head(DMLtest.type[ix,])
# # Saving the resulting data frame containing chromosome number, CpG site position,
# # test statistics, p-values (from normal distribution), and FDR.
# # Rows are sorted by p-value of the CpG sites
# write.table(DMLtest.type[ix,], file = "DMLtest.for.type_sorted.by.pval.txt", sep = "\t", row.names = FALSE)
# # Saving the COMPLETE dataset
# write.table(na.omit(DMLtest.type[ix,]), file = "DMLtest.for.type_sorted.by.pval_COMPLETE.txt", sep = "\t", row.names = FALSE)

### Use DMLtest.multiFactor function to test the "sex" effect
DMLtest.sex_common_smth = DMLtest.multiFactor(DMLfit_common_smth, coef="sexMale")
DMLtest.sex_common_smth_COMPLETE <- na.omit(DMLtest.sex_common_smth) 
# Saving the resulting data frame containing chromosome number, CpG site position,
# test statistics, p-values (from normal distribution), and FDR.
# Rows are sorted by chromosome/position of the CpG sites
write.table(DMLtest.sex_common_smth, file = "DMLtest.for.sex_common_smth_sorted.by.chr.pos.txt", sep = "\t", row.names = FALSE)
# Saving the COMPLETE dataset
write.table(na.omit(DMLtest.sex_common_smth), file = "DMLtest.for.sex_common_smth_sorted.by.chr.pos_COMPLETE.txt", sep = "\t", row.names = FALSE)
# # Obtain top ranked CpG sites, one can sort the data frame using following codes:
# ix=sort(DMLtest.sex[,"pvals"], index.return=TRUE)$ix
# # head(DMLtest.type[ix,])
# # Saving the resulting data frame containing chromosome number, CpG site position,
# # test statistics, p-values (from normal distribution), and FDR.
# # Rows are sorted by p-value of the CpG sites
# write.table(DMLtest.sex[ix,], file = "DMLtest.for.sex_sorted.by.pval.txt", sep = "\t", row.names = FALSE)
# # Saving the COMPLETE dataset
# write.table(na.omit(DMLtest.sex[ix,]), file = "DMLtest.for.sex_sorted.by.pval_COMPLETE.txt", sep = "\t", row.names = FALSE)

### Use DMLtest.multiFactor function to test the "type:sex" interaction effect
DMLtest.interaction_common_smth = DMLtest.multiFactor(DMLfit_common_smth, coef="typeTest:sexMale")
DMLtest.interaction_common_smth_COMPLETE <- na.omit(DMLtest.interaction_common_smth) 
# Saving the resulting data frame containing chromosome number, CpG site position,
# test statistics, p-values (from normal distribution), and FDR.
# Rows are sorted by chromosome/position of the CpG sites
write.table(DMLtest.interaction_common_smth, file = "DMLtest.for.interaction_common_smth_sorted.by.chr.pos.txt", sep = "\t", row.names = FALSE)
# Saving the COMPLETE dataset
write.table(na.omit(DMLtest.interaction_common_smth), file = "DMLtest.for.interaction_common_smth_sorted.by.chr.pos_COMPLETE.txt", sep = "\t", row.names = FALSE)
# # Obtain top ranked CpG sites, one can sort the data frame using following codes:
# ix=sort(DMLtest.interaction[,"pvals"], index.return=TRUE)$ix
# # head(DMLtest.type[ix,])
# # Saving the resulting data frame containing chromosome number, CpG site position,
# # test statistics, p-values (from normal distribution), and FDR.
# # Rows are sorted by p-value of the CpG sites
# write.table(DMLtest.interaction[ix,], file = "DMLtest.for.interaction_sorted.by.pval.txt", sep = "\t", row.names = FALSE)
# # Saving the COMPLETE dataset
# write.table(na.omit(DMLtest.interaction[ix,]), file = "DMLtest.for.interaction_sorted.by.pval_COMPLETE.txt", sep = "\t", row.names = FALSE)

# > dim(DMLtest.type)
# [1] 21609770        5
# > dim(na.omit(DMLtest.type))
# [1] 21144793        5
# 21609770-21144793
# > 21609770-21144793
# [1] 464977

################################################################################


################################################################################
### Venn diagram of the common and unique CpGs significant (at 5% level) between type and sex

# # Creating a screen session
# screen -S dmr2
# screen -r 1837459.dmr2
# TMOUT=$((96 * 3600))
# export TMOUT
# 
# # Start R:
# # conda activate r_env
# # R
# 
# screen -X -S 1837459.dmr2 quit
# (Done!)

# Setting the current working directory to /storage09/nawan2/DMR_analysis_DSS
setwd("/storage09/nawan2/DMR_analysis_DSS")

DMLtest.for.type_sorted.by.chr.pos <- read.table("DMLtest.for.type_sorted.by.chr.pos_COMPLETE.txt", header = T)
DMLtest.for.sex_sorted.by.chr.pos <- read.table("DMLtest.for.sex_sorted.by.chr.pos_COMPLETE.txt", header = T)
DMLtest.for.interaction_sorted.by.chr.pos <- read.table("DMLtest.for.interaction_sorted.by.chr.pos_COMPLETE.txt", header = T)

# For m10r samples
DMLtest.for.type_m10r_sorted.by.chr.pos <- read.table("DMLtest.for.type_m10r_sorted.by.chr.pos_COMPLETE.txt", header = T)
DMLtest.for.sex_m10r_sorted.by.chr.pos <- read.table("DMLtest.for.sex_m10r_sorted.by.chr.pos_COMPLETE.txt", header = T)
DMLtest.for.interaction_m10r_sorted.by.chr.pos <- read.table("DMLtest.for.interaction_m10r_sorted.by.chr.pos_COMPLETE.txt", header = T)

# For common CpG samples
DMLtest.for.type_common_sorted.by.chr.pos <- read.table("DMLtest.for.type_common_sorted.by.chr.pos_COMPLETE.txt", header = T)
DMLtest.for.sex_common_sorted.by.chr.pos <- read.table("DMLtest.for.sex_common_sorted.by.chr.pos_COMPLETE.txt", header = T)
DMLtest.for.interaction_common_sorted.by.chr.pos <- read.table("DMLtest.for.interaction_common_sorted.by.chr.pos_COMPLETE.txt", header = T)


### Important step:
### Create a new column "chr_pos" by combining "chr" and "pos" columns
DMLtest.for.type_sorted.by.chr.pos$chr_pos <- with(DMLtest.for.type_sorted.by.chr.pos, paste("chr", chr, "_pos", pos, sep=""))
head(DMLtest.for.type_sorted.by.chr.pos)
DMLtest.for.sex_sorted.by.chr.pos$chr_pos <- with(DMLtest.for.sex_sorted.by.chr.pos, paste("chr", chr, "_pos", pos, sep=""))
head(DMLtest.for.sex_sorted.by.chr.pos)
DMLtest.for.interaction_sorted.by.chr.pos$chr_pos <- with(DMLtest.for.interaction_sorted.by.chr.pos, paste("chr", chr, "_pos", pos, sep=""))
head(DMLtest.for.interaction_sorted.by.chr.pos)

# For m10r samples
DMLtest.for.type_m10r_sorted.by.chr.pos$chr_pos <- with(DMLtest.for.type_m10r_sorted.by.chr.pos, paste("chr", chr, "_pos", pos, sep=""))
DMLtest.for.sex_m10r_sorted.by.chr.pos$chr_pos <- with(DMLtest.for.sex_m10r_sorted.by.chr.pos, paste("chr", chr, "_pos", pos, sep=""))
DMLtest.for.interaction_m10r_sorted.by.chr.pos$chr_pos <- with(DMLtest.for.interaction_m10r_sorted.by.chr.pos, paste("chr", chr, "_pos", pos, sep=""))

dim(DMLtest.for.type_m10r_sorted.by.chr.pos)
dim(DMLtest.for.sex_m10r_sorted.by.chr.pos)
dim(DMLtest.for.interaction_m10r_sorted.by.chr.pos)

# For common CpG samples
DMLtest.for.type_common_sorted.by.chr.pos$chr_pos <- with(DMLtest.for.type_common_sorted.by.chr.pos, paste("chr", chr, "_pos", pos, sep=""))
DMLtest.for.sex_common_sorted.by.chr.pos$chr_pos <- with(DMLtest.for.sex_common_sorted.by.chr.pos, paste("chr", chr, "_pos", pos, sep=""))
DMLtest.for.interaction_common_sorted.by.chr.pos$chr_pos <- with(DMLtest.for.interaction_common_sorted.by.chr.pos, paste("chr", chr, "_pos", pos, sep=""))

dim(DMLtest.for.type_common_sorted.by.chr.pos)
dim(DMLtest.for.sex_common_sorted.by.chr.pos)
dim(DMLtest.for.interaction_common_sorted.by.chr.pos)


# > dim(DMLtest.for.type_m10r_sorted.by.chr.pos)
# [1] 20412347        6
# > dim(DMLtest.for.sex_m10r_sorted.by.chr.pos)
# [1] 20412347        6
# > dim(DMLtest.for.interaction_m10r_sorted.by.chr.pos)
# [1] 20412347        6


# Is there any interaction that is significant (FDR<0.05)?
sum(DMLtest.for.interaction_sorted.by.chr.pos$fdrs < 0.05)
# NA
# But why NA??

sum(is.na(DMLtest.for.interaction_sorted.by.chr.pos$fdrs))
# > sum(is.na(DMLtest.for.interaction_sorted.by.chr.pos$fdrs))
# [1] 464977
sum(is.na(DMLtest.for.type_sorted.by.chr.pos$fdrs))
# > sum(is.na(DMLtest.for.type_sorted.by.chr.pos$fdrs))
# [1] 464977
sum(is.na(DMLtest.for.sex_sorted.by.chr.pos$fdrs))
# > sum(is.na(DMLtest.for.sex_sorted.by.chr.pos$fdrs))
# [1] 464977

## After deleting missing (using complete data):
# Only one CpG (locus) has significant interaction
# > sum(DMLtest.for.interaction_sorted.by.chr.pos$fdrs < 0.05)
# [1] 1



# Sorted based on FDR<0.05
sig.CpGs.type <- DMLtest.for.type_sorted.by.chr.pos$chr_pos[DMLtest.for.type_sorted.by.chr.pos$fdrs < 0.05]
sig.CpGs.sex <- DMLtest.for.sex_sorted.by.chr.pos$chr_pos[DMLtest.for.sex_sorted.by.chr.pos$fdrs < 0.05]
sig.CpGs.interaction <- DMLtest.for.interaction_sorted.by.chr.pos$chr_pos[DMLtest.for.interaction_sorted.by.chr.pos$fdrs < 0.05]

length(sig.CpGs.type)
length(sig.CpGs.sex)
length(sig.CpGs.interaction)

# > length(sig.CpGs.type)
# [1] 18969580
# > length(sig.CpGs.sex)
# [1] 32803
# > length(sig.CpGs.interaction)
# [1] 1

length(intersect(sig.CpGs.type, sig.CpGs.sex))
# > length(intersect(sig.CpGs.type, sig.CpGs.sex))
# [1] 23760

length(unique(sig.CpGs.type))
length(unique(sig.CpGs.sex))
length(unique(sig.CpGs.interaction))

# > length(unique(sig.CpGs.type))
# [1] 18969580
# > length(unique(sig.CpGs.sex))
# [1] 32803
# > length(unique(sig.CpGs.interaction))
# [1] 1

# For m10r samples

sig.CpGs.type_m10r <- DMLtest.for.type_m10r_sorted.by.chr.pos$chr_pos[DMLtest.for.type_m10r_sorted.by.chr.pos$fdrs < 0.05]
sig.CpGs.sex_m10r <- DMLtest.for.sex_m10r_sorted.by.chr.pos$chr_pos[DMLtest.for.sex_m10r_sorted.by.chr.pos$fdrs < 0.05]
sig.CpGs.interaction_m10r <- DMLtest.for.interaction_m10r_sorted.by.chr.pos$chr_pos[DMLtest.for.interaction_m10r_sorted.by.chr.pos$fdrs < 0.05]

length(sig.CpGs.type_m10r)
length(sig.CpGs.sex_m10r)
length(sig.CpGs.interaction_m10r)

# > length(sig.CpGs.type_m10r)
# [1] 18500180
# > length(sig.CpGs.sex_m10r)
# [1] 25202
# > length(sig.CpGs.interaction_m10r)
# [1] 4


# For m10r samples

sig.CpGs.type_common <- DMLtest.for.type_common_sorted.by.chr.pos$chr_pos[DMLtest.for.type_common_sorted.by.chr.pos$fdrs < 0.05]
sig.CpGs.sex_common <- DMLtest.for.sex_common_sorted.by.chr.pos$chr_pos[DMLtest.for.sex_common_sorted.by.chr.pos$fdrs < 0.05]
sig.CpGs.interaction_common <- DMLtest.for.interaction_common_sorted.by.chr.pos$chr_pos[DMLtest.for.interaction_common_sorted.by.chr.pos$fdrs < 0.05]

length(sig.CpGs.type_common)
length(sig.CpGs.sex_common)
length(sig.CpGs.interaction_common)

# > length(sig.CpGs.type_m10r)
# [1] 18500180
# > length(sig.CpGs.sex_m10r)
# [1] 25202
# > length(sig.CpGs.interaction_m10r)
# [1] 4


### But why some of them are not unique? Maybe I should involve the chr? YES!
# Now unique because of creating the chr_pos variable and deleting missing.

# if (!require(devtools)) install.packages("devtools")
# install.packages("htmltools")
# library(htmltools)
# devtools::install_github("yanlinlin82/ggvenn")

# install.packages("ggvenn")
# 
# library(ggvenn)

# ### Using the VennDiagram package (I did not use this!)
# 
# # Install and load necessary package
# install.packages("VennDiagram")
# library(VennDiagram)
# 
# # Create a list of sets for the Venn diagram
# venn_list <- list(
#   Type = sig.CpGs.type,
#   Sex = sig.CpGs.sex,
#   Interaction = sig.CpGs.interaction
# )
# 
# # Create the Venn diagram
# venn.plot <- VennDiagram::venn.diagram(
#   x = venn_list, file=NULL,
#   category.names = c("Type", "Sex", "Interaction")
# )
# 
# # Display the Venn diagram
# grid.draw(venn.plot)
# 
# # Save the Venn diagram as a PDF file in the current working directory
# pdf("venn_diagram_sig_CpGs.pdf", width = 8, height = 6)  # Adjust width and height as needed
# grid.draw(venn.plot, main = "Venn diagram of significant CpGs")
# dev.off()  # Close the PDF device

################################################################################

### Using ggvenn package (I used this!)

# Setting the current working directory to /storage09/nawan2/DMR_analysis_DSS
setwd("/storage09/nawan2/DMR_analysis_DSS")

# Install and load necessary packages
install.packages("ggvenn")
install.packages("gridExtra")
library(ggvenn)
library(gridExtra)

# Create Venn diagram with labels
venn <- ggvenn(list(
  Type = sig.CpGs.type,
  Sex = sig.CpGs.sex,
  Interaction = sig.CpGs.interaction
))

# Add title to the plot
venn <- venn + ggtitle("Venn diagram of significant CpGs")

# Save the plot to the current working directory as PDF
ggsave("venn_diagram_sig.at.0.05_CpGs.pdf", plot = venn, width = 6, height = 6)



# Create Venn diagram with labels
venn_m10r <- ggvenn(list(
  Type = sig.CpGs.type_m10r,
  Sex = sig.CpGs.sex_m10r,
  Interaction = sig.CpGs.interaction_m10r
))

# Add title to the plot
venn_m10r <- venn_m10r + ggtitle("Venn diagram of significant CpGs")

# Save the plot to the current working directory as PDF
ggsave("venn_m10r_diagram_sig.at.0.05_CpGs.pdf", plot = venn_m10r, width = 6, height = 6)

################################################################################


################################################################################
### Volcano plots

# Load the necessary packages
library(ggplot2)

# Create a volcano plot for type
volcano_plot_type <- ggplot(DMLtest.for.type_sorted.by.chr.pos, aes(x = stat, y = -log10(fdrs))) +
  geom_point(aes(color = ifelse(stat > 0, "Positive", "Negative")), size = 3) +
  scale_color_manual(values = c("Positive" = "blue", "Negative" = "red")) +
  labs(x = "Stat", y = "-log10(FDR)", color = "Effects") +
  ggtitle("Volcano plot for Type (main effect)") +
  theme_minimal()

# Save the volcano plot as a PNG file in the current working directory
ggsave("volcano_plot_for_type.png", plot = volcano_plot_type, width = 6, height = 6, dpi = 300)


# Create a volcano plot for sex
volcano_plot_sex <- ggplot(DMLtest.for.sex_sorted.by.chr.pos, aes(x = stat, y = -log10(fdrs))) +
  geom_point(aes(color = ifelse(stat > 0, "Positive", "Negative")), size = 3) +
  scale_color_manual(values = c("Positive" = "blue", "Negative" = "red")) +
  labs(x = "Stat", y = "-log10(FDR)", color = "Effects") +
  ggtitle("Volcano plot for Sex (main effect)") +
  theme_minimal()

# Save the volcano plot as a PNG file in the current working directory
ggsave("volcano_plot_for_sex.png", plot = volcano_plot_sex, width = 6, height = 6, dpi = 300)


# Create a volcano plot for interaction
volcano_plot_interaction <- ggplot(DMLtest.for.interaction_sorted.by.chr.pos, aes(x = stat, y = -log10(fdrs))) +
  geom_point(aes(color = ifelse(stat > 0, "Positive", "Negative")), size = 3) +
  scale_color_manual(values = c("Positive" = "blue", "Negative" = "red")) +
  labs(x = "Stat", y = "-log10(FDR)", color = "Effects") +
  ggtitle("Volcano plot for interaction") +
  theme_minimal()

# Save the volcano plot as a PNG file in the current working directory
ggsave("volcano_plot_for_interaction.png", plot = volcano_plot_interaction, width = 6, height = 6, dpi = 300)



### Volcano plots for m10r

# Create a volcano plot for type
volcano_plot_type_m10r <- ggplot(DMLtest.for.type_m10r_sorted.by.chr.pos, aes(x = stat, y = -log10(fdrs))) +
  geom_point(aes(color = ifelse(stat > 0, "Positive", "Negative")), size = 3) +
  scale_color_manual(values = c("Positive" = "blue", "Negative" = "red")) +
  labs(x = "Stat", y = "-log10(FDR)", color = "Effects") +
  ggtitle("Volcano plot for Type (main effect)") +
  theme_minimal()

# Save the volcano plot as a PNG file in the current working directory
ggsave("volcano_plot_m10r_for_type.png", plot = volcano_plot_type_m10r, width = 6, height = 6, dpi = 300)


# Create a volcano plot for sex
volcano_plot_sex_m10r <- ggplot(DMLtest.for.sex_m10r_sorted.by.chr.pos, aes(x = stat, y = -log10(fdrs))) +
  geom_point(aes(color = ifelse(stat > 0, "Positive", "Negative")), size = 3) +
  scale_color_manual(values = c("Positive" = "blue", "Negative" = "red")) +
  labs(x = "Stat", y = "-log10(FDR)", color = "Effects") +
  ggtitle("Volcano plot for Sex (main effect)") +
  theme_minimal()

# Save the volcano plot as a PNG file in the current working directory
ggsave("volcano_plot_m10r_for_sex.png", plot = volcano_plot_sex_m10r, width = 6, height = 6, dpi = 300)


# Create a volcano plot for interaction
volcano_plot_interaction_m10r <- ggplot(DMLtest.for.interaction_m10r_sorted.by.chr.pos, aes(x = stat, y = -log10(fdrs))) +
  geom_point(aes(color = ifelse(stat > 0, "Positive", "Negative")), size = 3) +
  scale_color_manual(values = c("Positive" = "blue", "Negative" = "red")) +
  labs(x = "Stat", y = "-log10(FDR)", color = "Effects") +
  ggtitle("Volcano plot for interaction") +
  theme_minimal()

# Save the volcano plot as a PNG file in the current working directory
ggsave("volcano_plot_m10r_for_interaction.png", plot = volcano_plot_interaction_m10r, width = 6, height = 6, dpi = 300)

# ### To make non-significant dots grey:
# 
# # Define color based on FDR threshold
# color <- ifelse(DMLtest.for.type_m10r_sorted.by.chr.pos$fdrs > 0.05, "Non-significant", ifelse(DMLtest.for.type_m10r_sorted.by.chr.pos$stat > 0, "Positive", "Negative"))
# 
# # Create a volcano plot for type
# volcano_plot_type_m10r <- ggplot(DMLtest.for.type_m10r_sorted.by.chr.pos, aes(x = stat, y = -log10(fdrs))) +
#   geom_point(aes(color = color), size = 3) +
#   scale_color_manual(values = c("Positive" = "blue", "Negative" = "red", "Non-significant" = "grey")) +
#   labs(x = "Stat", y = "-log10(FDR)", color = "Effects") +
#   ggtitle("Volcano plot for Type (main effect)") +
#   theme_minimal()
# 
# # Save the volcano plot as a PNG file in the current working directory
# ggsave("volcano_plot_m10r_for_type.png", plot = volcano_plot_type_m10r, width = 6, height = 6, dpi = 300)
# 
# 
# # Similar modifications for other volcano plots


################################################################################



################################################################################
### Plot to see how many CpGs are common across samples

### An example before we begin

# Define the vectors
s1 <- c(1, 2, 3, 4, 5)
s2 <- c(1, 2, 3, 5, 6)
s3 <- c(1, 2, 7, 8)
s4 <- c(1, 9, 10)

# Create a list of vectors
vectors <- list(s1, s2, s3, s4)

# Initialize a vector to store whether each vector has at least one common element
vector_has_common <- logical(length(vectors))

# Loop through each vector and check if it has at least one common element with any other vector
for (i in 1:length(vectors)) {
  has_common <- any(sapply(vectors[-i], function(x) length(intersect(vectors[[i]], x))) > 0)
  vector_has_common[i] <- has_common
}

# Count the number of vectors that have at least one common element
num_vectors_with_common <- sum(vector_has_common)

# # Output the count
# cat("Number of vectors with at least one common element:", num_vectors_with_common, "\n")


# Initialize a vector to store whether each vector has at least two common elements
vector_has_at_least_two_common <- logical(length(vectors))

# Loop through each vector and check if it has at least two common elements with any other vector
for (i in 1:length(vectors)) {
  has_at_least_two_common <- any(sapply(vectors[-i], function(x) length(intersect(vectors[[i]], x))) >= 2)
  vector_has_at_least_two_common[i] <- has_at_least_two_common
}

# Count the number of vectors that have at least two common elements
num_vectors_with_at_least_two_common <- sum(vector_has_at_least_two_common)

# # Output the count
# cat("Number of vectors with at least two common elements:", num_vectors_with_at_least_two_common, "\n")


# Initialize a vector to store whether each vector has at least three common elements
vector_has_at_least_three_common <- logical(length(vectors))

# Loop through each vector and check if it has at least three common elements with any other vector
for (i in 1:length(vectors)) {
  has_at_least_three_common <- any(sapply(vectors[-i], function(x) length(intersect(vectors[[i]], x))) >= 3)
  vector_has_at_least_three_common[i] <- has_at_least_three_common
}

# Count the number of vectors that have at least three common elements
num_vectors_with_at_least_three_common <- sum(vector_has_at_least_three_common)

# # Output the count
# cat("Number of vectors with at least three or more common elements:", num_vectors_with_at_least_three_common, "\n")


# Initialize a vector to store whether each vector has at least three common elements
vector_has_at_least_four_common <- logical(length(vectors))

# Loop through each vector and check if it has at least three common elements with any other vector
for (i in 1:length(vectors)) {
  has_at_least_four_common <- any(sapply(vectors[-i], function(x) length(intersect(vectors[[i]], x))) >= 4)
  vector_has_at_least_four_common[i] <- has_at_least_four_common
}

# Count the number of vectors that have at least three common elements
num_vectors_with_at_least_four_common <- sum(vector_has_at_least_four_common)

# # Output the count
# cat("Number of vectors with at least four or more common elements:", num_vectors_with_at_least_four_common, "\n")


at_least_common_count <- c(num_vectors_with_common,
                           num_vectors_with_at_least_two_common,
                           num_vectors_with_at_least_three_common,
                           num_vectors_with_at_least_four_common)

# Define the x-axis labels
x_labels <- c("1+", "2+", "3+", "4+")

# Plot the barplot with custom x-axis labels
barplot(at_least_common_count, 
        xlab = "Number of Common CpGs", 
        ylab = "Number of Samples", 
        main = "Samples with Common CpGs", 
        names.arg = x_labels,
        col = "grey")

################################################################################

### Now, creating the barplot for samples with common CpGs

# Setting the current working directory to /storage09/nawan2/DMR_analysis_DSS
setwd("/storage09/nawan2/DMR_analysis_DSS")

### Important step:
### Create a new column "chr_pos" by combining "chr" and "pos" columns
g11$chr_pos <- with(g11, paste("chr", chr, "_pos", pos, sep=""))
g12$chr_pos <- with(g12, paste("chr", chr, "_pos", pos, sep=""))
g14$chr_pos <- with(g14, paste("chr", chr, "_pos", pos, sep=""))
g15$chr_pos <- with(g15, paste("chr", chr, "_pos", pos, sep=""))
g21$chr_pos <- with(g21, paste("chr", chr, "_pos", pos, sep=""))
g22$chr_pos <- with(g22, paste("chr", chr, "_pos", pos, sep=""))
g22$chr_pos <- with(g22, paste("chr", chr, "_pos", pos, sep=""))
g24$chr_pos <- with(g24, paste("chr", chr, "_pos", pos, sep=""))
g25$chr_pos <- with(g25, paste("chr", chr, "_pos", pos, sep=""))
g31$chr_pos <- with(g31, paste("chr", chr, "_pos", pos, sep=""))
g32$chr_pos <- with(g32, paste("chr", chr, "_pos", pos, sep=""))
g32$chr_pos <- with(g32, paste("chr", chr, "_pos", pos, sep=""))
g34$chr_pos <- with(g34, paste("chr", chr, "_pos", pos, sep=""))
g35$chr_pos <- with(g35, paste("chr", chr, "_pos", pos, sep=""))
g41$chr_pos <- with(g41, paste("chr", chr, "_pos", pos, sep=""))
g43$chr_pos <- with(g43, paste("chr", chr, "_pos", pos, sep=""))
g44$chr_pos <- with(g44, paste("chr", chr, "_pos", pos, sep=""))
g45$chr_pos <- with(g45, paste("chr", chr, "_pos", pos, sep=""))

# Create a list of vectors
vectors <- list(g11$chr_pos, g12$chr_pos, g14$chr_pos, g15$chr_pos,
                g21$chr_pos, g22$chr_pos, g23$chr_pos, g24$chr_pos, g25$chr_pos,
                g31$chr_pos, g32$chr_pos, g33$chr_pos, g34$chr_pos, g35$chr_pos,
                g41$chr_pos, g43$chr_pos, g44$chr_pos, g45$chr_pos)

# Initialize a vector to store whether each vector has at least one common element
vector_has_common <- logical(length(vectors))

# Loop through each vector and check if it has at least one common element with any other vector
for (i in 1:length(vectors)) {
  has_common <- any(sapply(vectors[-i], function(x) length(intersect(vectors[[i]], x))) > 0)
  vector_has_common[i] <- has_common
}

# Count the number of vectors that have at least one common element
num_vectors_with_common <- sum(vector_has_common)

# # Output the count
# cat("Number of vectors with at least one common element:", num_vectors_with_common, "\n")


# Initialize a vector to store whether each vector has at least two common elements
vector_has_at_least_two_common <- logical(length(vectors))

# Loop through each vector and check if it has at least two common elements with any other vector
for (i in 1:length(vectors)) {
  has_at_least_two_common <- any(sapply(vectors[-i], function(x) length(intersect(vectors[[i]], x))) >= 2)
  vector_has_at_least_two_common[i] <- has_at_least_two_common
}

# Count the number of vectors that have at least two common elements
num_vectors_with_at_least_two_common <- sum(vector_has_at_least_two_common)

# # Output the count
# cat("Number of vectors with at least two common elements:", num_vectors_with_at_least_two_common, "\n")


# Initialize a vector to store whether each vector has at least three common elements
vector_has_at_least_three_common <- logical(length(vectors))

# Loop through each vector and check if it has at least three common elements with any other vector
for (i in 1:length(vectors)) {
  has_at_least_three_common <- any(sapply(vectors[-i], function(x) length(intersect(vectors[[i]], x))) >= 3)
  vector_has_at_least_three_common[i] <- has_at_least_three_common
}

# Count the number of vectors that have at least three common elements
num_vectors_with_at_least_three_common <- sum(vector_has_at_least_three_common)

# # Output the count
# cat("Number of vectors with at least three or more common elements:", num_vectors_with_at_least_three_common, "\n")


# Initialize a vector to store whether each vector has at least three common elements
vector_has_at_least_four_common <- logical(length(vectors))

# Loop through each vector and check if it has at least three common elements with any other vector
for (i in 1:length(vectors)) {
  has_at_least_four_common <- any(sapply(vectors[-i], function(x) length(intersect(vectors[[i]], x))) >= 4)
  vector_has_at_least_four_common[i] <- has_at_least_four_common
}

# Count the number of vectors that have at least three common elements
num_vectors_with_at_least_four_common <- sum(vector_has_at_least_four_common)

# # Output the count
# cat("Number of vectors with at least four or more common elements:", num_vectors_with_at_least_four_common, "\n")


at_least_common_count <- c(num_vectors_with_common,
                           num_vectors_with_at_least_two_common,
                           num_vectors_with_at_least_three_common,
                           num_vectors_with_at_least_four_common)

# Define the x-axis labels
x_labels <- c("1+", "2+", "3+", "4+")

# Open a PDF device
pdf("barplot_samples_with_common_CpGs.pdf")

# Plot the barplot with custom x-axis labels
barplot(at_least_common_count, 
        xlab = "Number of Common CpGs", 
        ylab = "Number of Samples", 
        main = "Samples with Common CpGs", 
        names.arg = x_labels,
        col = "grey")

# Close the PDF device
dev.off()
################################################################################


################################################################################
### PCA

### Preparing the dataset for PCA

# First check the chr variable carefully

table(g11$chr)

# Creating the percent methylation
g11$pct_methylation_g11 <- g11$X/g11$N

# Variation across all samples for CpG
# 5-10% of the most variable CpGs (Coleman's paper), probably separated by chromosome

# Calculate the threshold for the top 5% pct_methylation values
threshold_g11 <- quantile(g11$pct_methylation_g11, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g11 <- g11[g11$pct_methylation_g11 > threshold_g11, ]
# Selecting only chr_pos (CpG) and pct_methylation
g11.pca <- top_5_percent_g11[,c("chr_pos", "pct_methylation_g11")]

pca.data <- g11.pca

# Creating the percent methylation
g12$pct_methylation_g12 <- g12$X/g12$N
# Calculate the threshold for the top 5% pct_methylation values
threshold_g12 <- quantile(g12$pct_methylation_g12, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g12 <- g12[g12$pct_methylation_g12 > threshold_g12, ]
# Selecting only chr_pos (CpG) and pct_methylation
g12.pca <- top_5_percent_g12[,c("chr_pos", "pct_methylation_g12")]

pca.data <- merge(pca.data, g12.pca, by="chr_pos", all=TRUE)

# Creating the percent methylation
g14$pct_methylation_g14 <- g14$X/g14$N
# Calculate the threshold for the top 5% pct_methylation values
threshold_g14 <- quantile(g14$pct_methylation_g14, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g14 <- g14[g14$pct_methylation_g14 > threshold_g14, ]
# Selecting only chr_pos (CpG) and pct_methylation
g14.pca <- top_5_percent_g14[,c("chr_pos", "pct_methylation_g14")]

pca.data <- merge(pca.data, g14.pca, by="chr_pos", all=TRUE)

# Creating the percent methylation
g15$pct_methylation_g15 <- g15$X/g15$N
# Calculate the threshold for the top 5% pct_methylation values
threshold_g15 <- quantile(g15$pct_methylation_g15, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g15 <- g15[g15$pct_methylation_g15 > threshold_g15, ]
# Selecting only chr_pos (CpG) and pct_methylation
g15.pca <- top_5_percent_g15[,c("chr_pos", "pct_methylation_g15")]

pca.data <- merge(pca.data, g15.pca, by="chr_pos", all=TRUE)

# Creating the percent methylation
g21$pct_methylation_g21 <- g21$X/g21$N
# Calculate the threshold for the top 5% pct_methylation values
threshold_g21 <- quantile(g21$pct_methylation_g21, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g21 <- g21[g21$pct_methylation_g21 > threshold_g21, ]
# Selecting only chr_pos (CpG) and pct_methylation
g21.pca <- top_5_percent_g21[,c("chr_pos", "pct_methylation_g21")]

pca.data <- merge(pca.data, g21.pca, by="chr_pos", all=TRUE)

# Creating the percent methylation
g22$pct_methylation_g22 <- g22$X/g22$N
# Calculate the threshold for the top 5% pct_methylation values
threshold_g22 <- quantile(g22$pct_methylation_g22, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g22 <- g22[g22$pct_methylation_g22 > threshold_g22, ]
# Selecting only chr_pos (CpG) and pct_methylation
g22.pca <- top_5_percent_g22[,c("chr_pos", "pct_methylation_g22")]

pca.data <- merge(pca.data, g22.pca, by="chr_pos", all=TRUE)

# Creating the percent methylation
g23$pct_methylation_g23 <- g23$X/g23$N
# Calculate the threshold for the top 5% pct_methylation values
threshold_g23 <- quantile(g23$pct_methylation_g23, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g23 <- g23[g23$pct_methylation_g23 > threshold_g23, ]
# Selecting only chr_pos (CpG) and pct_methylation
g23.pca <- top_5_percent_g23[,c("chr_pos", "pct_methylation_g23")]

pca.data <- merge(pca.data, g23.pca, by="chr_pos", all=TRUE)

# Creating the percent methylation
g24$pct_methylation_g24 <- g24$X/g24$N
# Calculate the threshold for the top 5% pct_methylation values
threshold_g24 <- quantile(g24$pct_methylation_g24, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g24 <- g24[g24$pct_methylation_g24 > threshold_g24, ]
# Selecting only chr_pos (CpG) and pct_methylation
g24.pca <- top_5_percent_g24[,c("chr_pos", "pct_methylation_g24")]

pca.data <- merge(pca.data, g24.pca, by="chr_pos", all=TRUE)

# Creating the percent methylation
g25$pct_methylation_g25 <- g25$X/g25$N
# Calculate the threshold for the top 5% pct_methylation values
threshold_g25 <- quantile(g25$pct_methylation_g25, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g25 <- g25[g25$pct_methylation_g25 > threshold_g25, ]
# Selecting only chr_pos (CpG) and pct_methylation
g25.pca <- top_5_percent_g25[,c("chr_pos", "pct_methylation_g25")]

pca.data <- merge(pca.data, g25.pca, by="chr_pos", all=TRUE)

# Creating the percent methylation
g31$pct_methylation_g31 <- g31$X/g31$N
# Calculate the threshold for the top 5% pct_methylation values
threshold_g31 <- quantile(g31$pct_methylation_g31, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g31 <- g31[g31$pct_methylation_g31 > threshold_g31, ]
# Selecting only chr_pos (CpG) and pct_methylation
g31.pca <- top_5_percent_g31[,c("chr_pos", "pct_methylation_g31")]

pca.data <- merge(pca.data, g31.pca, by="chr_pos", all=TRUE)

# Creating the percent methylation
g32$pct_methylation_g32 <- g32$X/g32$N
# Calculate the threshold for the top 5% pct_methylation values
threshold_g32 <- quantile(g32$pct_methylation_g32, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g32 <- g32[g32$pct_methylation_g32 > threshold_g32, ]
# Selecting only chr_pos (CpG) and pct_methylation
g32.pca <- top_5_percent_g32[,c("chr_pos", "pct_methylation_g32")]

pca.data <- merge(pca.data, g32.pca, by="chr_pos", all=TRUE)

# Creating the percent methylation
g33$pct_methylation_g33 <- g33$X/g33$N
# Calculate the threshold for the top 5% pct_methylation values
threshold_g33 <- quantile(g33$pct_methylation_g33, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g33 <- g33[g33$pct_methylation_g33 > threshold_g33, ]
# Selecting only chr_pos (CpG) and pct_methylation
g33.pca <- top_5_percent_g33[,c("chr_pos", "pct_methylation_g33")]

pca.data <- merge(pca.data, g33.pca, by="chr_pos", all=TRUE)

# Creating the percent methylation
g34$pct_methylation_g34 <- g34$X/g34$N
# Calculate the threshold for the top 5% pct_methylation values
threshold_g34 <- quantile(g34$pct_methylation_g34, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g34 <- g34[g34$pct_methylation_g34 > threshold_g34, ]
# Selecting only chr_pos (CpG) and pct_methylation
g34.pca <- top_5_percent_g34[,c("chr_pos", "pct_methylation_g34")]

pca.data <- merge(pca.data, g34.pca, by="chr_pos", all=TRUE)

# Creating the percent methylation
g35$pct_methylation_g35 <- g35$X/g35$N
# Calculate the threshold for the top 5% pct_methylation values
threshold_g35 <- quantile(g35$pct_methylation_g35, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g35 <- g35[g35$pct_methylation_g35 > threshold_g35, ]
# Selecting only chr_pos (CpG) and pct_methylation
g35.pca <- top_5_percent_g35[,c("chr_pos", "pct_methylation_g35")]

pca.data <- merge(pca.data, g35.pca, by="chr_pos", all=TRUE)

# Creating the percent methylation
g41$pct_methylation_g41 <- g41$X/g41$N
# Calculate the threshold for the top 5% pct_methylation values
threshold_g41 <- quantile(g41$pct_methylation_g41, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g41 <- g41[g41$pct_methylation_g41 > threshold_g41, ]
# Selecting only chr_pos (CpG) and pct_methylation
g41.pca <- top_5_percent_g41[,c("chr_pos", "pct_methylation_g41")]

pca.data <- merge(pca.data, g41.pca, by="chr_pos", all=TRUE)

# Creating the percent methylation
g43$pct_methylation_g43 <- g43$X/g43$N
# Calculate the threshold for the top 5% pct_methylation values
threshold_g43 <- quantile(g43$pct_methylation_g43, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g43 <- g43[g43$pct_methylation_g43 > threshold_g43, ]
# Selecting only chr_pos (CpG) and pct_methylation
g43.pca <- top_5_percent_g43[,c("chr_pos", "pct_methylation_g43")]

pca.data <- merge(pca.data, g43.pca, by="chr_pos", all=TRUE)

# Creating the percent methylation
g44$pct_methylation_g44 <- g44$X/g44$N
# Calculate the threshold for the top 5% pct_methylation values
threshold_g44 <- quantile(g44$pct_methylation_g44, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g44 <- g44[g44$pct_methylation_g44 > threshold_g44, ]
# Selecting only chr_pos (CpG) and pct_methylation
g44.pca <- top_5_percent_g44[,c("chr_pos", "pct_methylation_g44")]

pca.data <- merge(pca.data, g44.pca, by="chr_pos", all=TRUE)

# Creating the percent methylation
g45$pct_methylation_g45 <- g45$X/g45$N
# Calculate the threshold for the top 5% pct_methylation values
threshold_g45 <- quantile(g45$pct_methylation_g45, probs = 0.95)
# Subset the dataset to select rows with pct_methylation values greater than the threshold
top_5_percent_g45 <- g45[g45$pct_methylation_g45 > threshold_g45, ]
# Selecting only chr_pos (CpG) and pct_methylation
g45.pca <- top_5_percent_g45[,c("chr_pos", "pct_methylation_g45")]

pca.data <- merge(pca.data, g45.pca, by="chr_pos", all=TRUE)

pca.data.complete <- na.omit(pca.data)

dim(pca.data)
dim(pca.data.complete)


################################################################################


################################################################################
###

################################################################################


################################################################################
###

################################################################################
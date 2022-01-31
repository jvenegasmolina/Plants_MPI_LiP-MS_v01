# DDA Statistical analysis to identify differential enrichment at peptide level upon LiP-MS experiments using DEP. 
# by Jhon Venegas-Molina, 14/09/2021
# This script is based on the original DEP script (https://bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/DEP.R) with several modifications to allow LiP-MS analysis. 
# Package: DEP 1.14.0
# http://bioconductor.org/packages/release/bioc/html/DEP.html

## ----install DEP--------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("DEP")
#  
#  library("DEP")
#  

## ----required packages--------------------------------------------------------
suppressPackageStartupMessages({
  library("BiocStyle")
  library("DEP")
  library("tidyr")
  library("dplyr")
})

## ----Prepare MaxQuant output files for DEP analysis---------------------------
# Load the "peptides.txt" output file from the MaxQuant analysis
data_peptides <- read.delim("peptides.txt")

# Avoid spaces in column names, therefore, replace spaces with a dot using make.names function
names(data_peptides) <- make.names(names(data_peptides), unique=TRUE)

# Also spaces in column names could be replaced with underscore 
# names(data_peptides) <- gsub(" ", "_", names(data_peptides))

# Create a new column called "name". In this column, concatenate the protein name + sequence (e.g. Q8VYG9_LNVDQNPLVVPPEEVVK)
data_edited <- data_peptides %>%
  unite("name", c(Proteins, Sequence), remove = FALSE)

# Change "id" column to "ID", this is important for the generation of a Summarized Experiment object downstream
names(data_edited)[names(data_edited) == "id"] <- "ID"

# Write a data frame with the experimental design. It should contain 3 columns. One for "label" with the name of each sample, second "condition" indicating the treatment corresponding to each sample, and third "replicate" indicating the number of replicate 
LiP_MS_exp_design <- data.frame(label = c("MOCK_1",
                           "MOCK_2",
                           "MOCK_3",
                           "MOCK_4",
                           "Phe_C1_1",
                           "Phe_C1_2",
                           "Phe_C1_3",
                           "Phe_C1_4",
                           "Phe_C2_1",
                           "Phe_C2_2",
                           "Phe_C2_3",
                           "Phe_C2_4"),
                 condition = c("MOCK",
                               "MOCK",
                               "MOCK",
                               "MOCK",
                               "Phe_C1",
                               "Phe_C1",
                               "Phe_C1",
                               "Phe_C1",
                               "Phe_C2",
                               "Phe_C2",
                               "Phe_C2",
                               "Phe_C2"),
                 replicate = c(1,2,3,4)
)

# Check that the experimental design data frame is correct
print (LiP_MS_exp_design)

# Example:
# label       condition   replicate
#    MOCK_1      MOCK         1
#    MOCK_2      MOCK         2
#    MOCK_3      MOCK         3
#    MOCK_4      MOCK         4
#   Phe_C1_1    Phe_C1         1
#   Phe_C1_2    Phe_C1         2
#   Phe_C1_3    Phe_C1         3
#   Phe_C1_4    Phe_C1         4
#   Phe_C2_1    Phe_C2         1
#   Phe_C2_2    Phe_C2         2
#   Phe_C2_3    Phe_C2         3
#   Phe_C2_4    Phe_C2         4

## ----Data processing----------------------------------------------------------
data <- data_edited

# Check that data was loaded correctly
head(data)
dim(data)

# Check that data in the intensities columns corresponds to the class "numeric"
colnames(data) # To identify the number of columns corresponding to intensities

# Example to check the 75th column corresponding to a Intensity of a treatment
class(data[,75])

# Filter contaminant proteins, decoy database hits, and no unique peptides, which are indicated by "+" in the columns "Potential.contaminants", "Reverse", and "no" in "Unique..Proteins" respectively 
data <- filter(data, Reverse != "+", Potential.contaminant != "+", Unique..Proteins. != "no")

# Check new dimensions of the data
dim(data)

# Check if are there any duplicated names?
data$name %>% duplicated() %>% any()

## ----Data preparation---------------------------------------------------------
# Display experimental design
design <- LiP_MS_exp_design
knitr::kable(design)

# Generate a SummarizedExperiment object using the experimental design
Intensity_columns <- grep("Intensity.", colnames(data)) # To get Intensity column numbers
experimental_design <- design
data_se <- make_se(data, Intensity_columns, experimental_design)

# Check the SummarizedExperiment object
data_se

## ----plot_data_noFilt, fig.width = 4, fig.height = 4--------------------------
# Plot a barplot of the peptide identification overlap between samples
plot_frequency(data_se)

## ----filter_missval-----------------------------------------------------------
# Filter for proteins that are identified in all replicates of at least one condition
data_filt_strict <- filter_missval(data_se, thr = 0)

# Less stringent filtering:
# Filter for proteins that are identified in 3 out of 4 replicates of at least one condition (thr = number of allowed missing values). Use this option for downstream analysis 
data_filt <- filter_missval(data_se, thr = 1)

## ----plot_data, fig.width = 4, fig.height = 4---------------------------------
# Plot a barplot of the number of identified peptides per samples
plot_numbers(data_filt)

## ----plot_data2, fig.width = 3, fig.height = 4--------------------------------
# Plot a barplot of the peptide identification overlap between samples
plot_coverage(data_filt)

## ----normalize data-----------------------------------------------------------
# Normalize the data
data_norm <- normalize_vsn(data_filt)

## ----plot_norm, fig.width = 4, fig.height = 5---------------------------------
# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)

## ----plot_missval, fig.height = 4, fig.width = 3------------------------------
# Plot a heatmap of peptides with missing values
plot_missval(data_filt)

## ----plot_detect, fig.height = 4, fig.width = 4-------------------------------
# Plot intensity distributions and cumulative fraction of peptides with and without missing values
plot_detect(data_filt)

## ----impute, results = "hide", message = FALSE, warning = FALSE, error = TRUE----
# All possible imputation methods are printed in an error, if an invalid function name is given.
impute(data_norm, fun = "")

# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
# data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)

## ----plot_imp, fig.width = 4, fig.height = 4----------------------------------
# Plot intensity distributions before and after imputation. The effect of the imputation on the distributions can be visualized here 
plot_imputation(data_norm, data_imp)
plot_imputation(data_norm, data_imp_man)

## ----statistics---------------------------------------------------------------
# Differential enrichment analysis  based on linear models and empirical Bayes statistics

# Test every sample versus control (define value of control)
data_diff <- test_diff(data_imp, type = "control", control = "MOCK")

# Test all possible comparisons of samples
# data_diff_all_contrasts <- test_diff(data_imp, type = "all")

# Test manually defined comparisons
# data_diff_manual <- test_diff(data_imp, type = "manual", 
#                              test = c("Ubi4_vs_Ctrl", "Ubi6_vs_Ctrl"))

## ----add_reject---------------------------------------------------------------
# Denote significant peptides based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.01, lfc = log2(2))

## ----pca, fig.height = 3, fig.width = 4---------------------------------------
# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

## ----corr, fig.height = 3, fig.width = 5--------------------------------------
# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

## ----heatmap, fig.height = 5, fig.width = 3-----------------------------------
# Plot a heatmap of all significant peptides with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))

## ----heatmap2, fig.height = 5, fig.width = 3----------------------------------
# Plot a heatmap of all significant peptides (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)

## ----volcano, fig.height = 5, fig.width = 5-----------------------------------
# Plot a volcano plot for the contrast "Phe_C1 vs MOCK""
plot_volcano(dep, contrast = "Phe_C1_vs_MOCK", label_size = 2, add_names = FALSE)
#plot_volcano(dep, contrast = "Phe_C1_vs_MOCK", label_size = 2, add_names = TRUE)

# Plot a volcano plot for the contrast "Phe_C2 vs MOCK""
plot_volcano(dep, contrast = "Phe_C2_vs_MOCK", label_size = 2, add_names = FALSE)
#plot_volcano(dep, contrast = "Phe_C2_vs_MOCK", label_size = 2, add_names = TRUE)

##
# Check important peptides to show as examples
## ----bar, fig.height = 3.5, fig.width = 3.5-----------------------------------
# Plot a barplot for P08312_ADHDTFWFDTT and P0AB91_ELLPPVAL
plot_single(dep, proteins = c("P08312_ADHDTFWFDTT", "P0AB91_ELLPPVAL"))

# Plot a barplot for the peptide P08312_ADHDTFWFDTT with the data centered
plot_single(dep, proteins = "P08312_ADHDTFWFDTT", type = "centered")

##

## ----overlap, fig.height = 4, fig.width = 6-----------------------------------
# Plot a frequency plot of significant proteins for the different conditions
plot_cond(dep)

## ----results_table------------------------------------------------------------
# Generate a results table
data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()

## ----results_table2-----------------------------------------------------------
# Column names of the results table
colnames(data_results)

## ----get_df-------------------------------------------------------------------
# Generate a wide data.frame
df_wide <- get_df_wide(dep)
# Generate a long data.frame
df_long <- get_df_long(dep)

## ----save_load, eval = FALSE--------------------------------------------------
#  # Save analyzed data
save(data_se, data_norm, data_imp, data_diff, dep, file = "data.RData")
#  # These data can be loaded in future R sessions using this command
#  load("data.RData")

# Generate tables for edition in Excel------------------------------------------
#The df_wide file can be used for further presentation of results 
#write.table(df_long, file = "df_long.txt", sep = "\t",
#            row.names = FALSE)
write.table(df_wide, file = "Results_df_wide.txt", sep = "\t",
            row.names = FALSE)

# End of analysis

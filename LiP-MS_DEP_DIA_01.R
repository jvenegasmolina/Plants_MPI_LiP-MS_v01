# DIA Statistical analysis to identify differential enrichment at peptide level upon LiP-MS experiments using DEP. 
# by Jhon Venegas-Molina, 14/09/2021
# This script is based on the original DEP script (https://bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/DEP.R) with several mofidications to allow LiP-MS analysis. 
# Package: DEP 1.12.0
# http://bioconductor.org/packages/release/bioc/html/DEP.html

## ----required packages, echo = FALSE, warning=FALSE, results="hide"-----------
suppressPackageStartupMessages({
  library("BiocStyle")
  library("DEP")
  library("dplyr")
})

## ----install, eval = FALSE----------------------------------------------------
#  
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("DEP")
#  
#  library("DEP")
#  

## ----prepare OpenSwath output files for DEP analysis-----------------------------
# Download data from Galaxy (file "peptide_signal.tabular") 
# Take out contaminant proteins (CON_), iRT proteins (sp|Biognosys_pep-a|iRT-Kit_WR_fusion_), and no unique peptides (peptides that were identified in more than one protein) using R or Galaxy tools
# Change name of column "ProteinName_FullPeptide" for "name" and add column "ID" (even if this column is empty)
# Make a .txt file with the experimental design (e.g "LiP_MS_exp_design.txt"). It should contain 3 columns. One for "label" with the name of each sample, second "condition" indicating the treatment of each sample, and third "replicate" indicating the number of replicate 
# Save file as .txt. and give a new name (e.g. "PyProphet_filtered_peptide_signal.txt")

## ----load data----------------------------------------------------------------
# Loading a package required for data handling
library("dplyr")

data <- read.delim("PyProphet_filtered_peptide_signal.txt")
head(data)
dim(data)


# Check that the data in the columns of intensities corresponds to the class "numeric". Following is an example to check the 12th column corresponding to a Intensity of a treatment
class(data[,12])
colnames(data) # To identify the number of columns corresponding to intensities

## ----new dimension----------------------------------------------------------------
dim(data)

## ----unique-------------------------------------------------------------------
# Are there any duplicated gene names?
data$ProteinName_FullPeptideName %>% duplicated() %>% any()

## ----expdesign, echo = FALSE--------------------------------------------------
# Display experimental design
design <- read.delim("LiP_MS_exp_design.txt")
knitr::kable(design)

## ----to_exprset---------------------------------------------------------------
# Generate a SummarizedExperiment object using an experimental design
LFQ_columns <- grep("LFQ.", colnames(data)) # To get LFQ column numbers
experimental_design <- design
data_se <- make_se(data, LFQ_columns, experimental_design)

# Check the SummarizedExperiment object
data_se

## ----plot_data_noFilt, fig.width = 4, fig.height = 4--------------------------
# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

## ----filter_missval-----------------------------------------------------------
# Filter for proteins that are identified in all replicates of at least one condition
data_filt_strict <- filter_missval(data_se, thr = 0)

# Less stringent filtering:
# Filter for proteins that are identified in 3 out of 4 replicates of at least one condition (thr = number of allowed missing values). Use this option for downstream analysis
data_filt <- filter_missval(data_se, thr = 1)

## ----plot_data, fig.width = 4, fig.height = 4---------------------------------
# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

## ----plot_data2, fig.width = 3, fig.height = 4--------------------------------
# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

## ----normalize----------------------------------------------------------------
# Normalize the data
data_norm <- normalize_vsn(data_filt)

## ----plot_norm, fig.width = 4, fig.height = 5---------------------------------
# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)

## ----plot_missval, fig.height = 4, fig.width = 3------------------------------
# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

## ----plot_detect, fig.height = 4, fig.width = 4-------------------------------
# Plot intensity distributions and cumulative fraction of proteins with and without missing values
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
# Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)

## ----statistics---------------------------------------------------------------
# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "MOCK")

# Test all possible comparisons of samples
# data_diff_all_contrasts <- test_diff(data_imp, type = "all")

# Test manually defined comparisons
# data_diff_manual <- test_diff(data_imp, type = "manual", 
#                              test = c("Ubi4_vs_Ctrl", "Ubi6_vs_Ctrl"))

## ----add_reject---------------------------------------------------------------
# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

## ----pca, fig.height = 3, fig.width = 4---------------------------------------
# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

## ----corr, fig.height = 3, fig.width = 5--------------------------------------
# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

## ----heatmap, fig.height = 5, fig.width = 3-----------------------------------
# Plot a heatmap of all significant proteins with the data centered per protein
# Check what happened with sample "SA_1_mM_3". I think I got an email that this sample was re-run?
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))

## ----heatmap2, fig.height = 5, fig.width = 3----------------------------------
# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)

## ----volcano, fig.height = 5, fig.width = 5-----------------------------------
# Plot a volcano plot for the contrast "ATP_10_mM vs MOCK""
plot_volcano(dep, contrast = "ATP_10_mM_vs_MOCK", label_size = 2, add_names = FALSE)

# Plot a volcano plot for the contrast "ATP_25_mM vs Ctrl""
plot_volcano(dep, contrast = "ATP_25_mM_vs_MOCK", label_size = 2, add_names = TRUE)

##
# Check important peptides to show as examples
## ----bar, fig.height = 3.5, fig.width = 3.5-----------------------------------
# Plot a barplot for USP15 and IKBKG
# plot_single(dep, proteins = c("USP15", "IKBKG"))

# Plot a barplot for the protein USP15 with the data centered
# plot_single(dep, proteins = "USP15", type = "centered")

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

#Generate tables for using in Excel
#The df_wide file can be used for further presentation of result 
#write.table(df_long, file = "df_long.txt", sep = "\t",
#            row.names = FALSE)
write.table(df_wide, file = "df_wide.txt", sep = "\t",
            row.names = FALSE)

# End of analysis

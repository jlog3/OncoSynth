# 
# integrating tcga with synthea data for EMR + cancer pathway interactive visualizer
# to extract plausible biomarkers 

# Identify Relevant Data: For BRCA (Breast Cancer), filter the datasets to obtain relevant demographic and RNASeq data.
# Download Data: Download the RNASeq data along with the demographic data (age, gender, race, etc.) which are crucial for your study.


library(TCGAbiolinks)


get_gdc_clinical_data <- function(project) {
    # Query clinical data for TCGA-BRCA project
    query_clinical <- GDCquery_clinic(project = project, type = "clinical", save.csv = FALSE)
    # query_clinical object already contains the clinical data you need, so you might not have to download anything extra. 
    # The GDCquery_clinic function in TCGAbiolinks seems to have fetched the clinical data and stored it directly in the query_clinical object in a tabular format.
  
    # View the clinical data
    # head(clinical_data)
    return(query_clinical)
}


get_gdc_rna_data <- function(project, sample_types, download = FALSE, customdir = NULL) {
  # browser()
  # Query for RNA-Seq data
  query_rnaseq <- GDCquery(project = project, 
                           data.category = "Transcriptome Profiling",
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "STAR - Counts",
                           sample.type = sample_types)
  
  # Download data if requested
  if (download) {
    if (is.null(customdir)) {
      GDCdownload(query_rnaseq, method = "api", files.per.chunk = 5)
    } else {
      GDCdownload(query_rnaseq, method = "api", files.per.chunk = 5, directory = customdir)
    }
  }
  
  # barcodes <- getResults(query_rnaseq, cols = c("submitter_id"))
  barcodes <- getResults(query_rnaseq, cols = c("cases"))
  submitter_ids <- getResults(query_rnaseq, cols = c("cases.submitter_id"))
  # Assuming 'barcodes' and 'submitter_ids' are vectors of equal length
  barcode_submitter_mapping <- data.frame(
    barcode = barcodes,
    submitter_id = submitter_ids,
    stringsAsFactors = FALSE
  )
  
  # ^---needed to merge on clinical data
  # return(list(barcodes=barcodes, submitter_ids=submitter_ids))
  return(barcode_submitter_mapping)
}

process_chunk <- function(barcodes, project, sample_type, customdir) {
  query <- GDCquery(
    project = project, 
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = sample_type,
    barcode = barcodes
  )
  
  GDCdownload(query, directory = customdir)
  data <- GDCprepare(query, directory = customdir)
  return(data)
}

# get matrix from the gene data
process_and_write_csv <- function(data_subset, filename, writerawcounts=FALSE) {
  # browser()
  # Convert to a regular data frame or matrix
  assay_data_subset <- assay(data_subset)
  #---^ what the user should be giving us -- the rawcounts
  # then we need to do the below preproc and then can get the log folds ourselves
  
  # Write to CSV raw counts
  if (writerawcounts) {
    write.csv(assay_data_subset, file = filename, row.names = TRUE)
  }
  head(assay_data_subset)
  
  # Remove version numbers from Ensembl IDs
  rownames(assay_data_subset) <- gsub("\\..*", "", rownames(assay_data_subset))
  
  # Averaging duplicates: when there are multiple versions of the same gene
  assay_data_subset <- aggregate(assay_data_subset, by = list(rownames(assay_data_subset)), FUN = mean)
  
  # Rename the rownames
  rownames(assay_data_subset) <- assay_data_subset$Group.1
  assay_data_subset <- assay_data_subset[-1]
  
  
  # Assuming assay_data_subset is your data frame or matrix
  
  # Check for duplicate row names
  duplicate_row_names <- rownames(assay_data_subset)[duplicated(rownames(assay_data_subset))]
  
  # See if there are any duplicates
  if (length(duplicate_row_names) > 0) {
    print("There are duplicate row names:")
    print(duplicate_row_names)
  } else {
    print("There are no duplicate row names.")
  }
  
  return(assay_data_subset)
}



# Example usage
# data_all_primary_tumor <- get_gdc_data("TCGA-BRCA", "Primary Tumor", download = TRUE, customdir = 'GDCdata_BRCA')
library(TCGAbiolinks)
library(SummarizedExperiment)
# Load necessary library
library(dplyr)

# write boths 
# rawcount file
# clinical file
# use: run for both tumor and normal tissue settings
# once we have those 4 files we can later merge if we need but we will be able to analyze however

main <- function(  
  out_dir = "/home/john/Desktop/cancer_vis",
  write_rawcount_path = "tcga_synthea_test/tcga_tumor_rawcount.csv",
  write_clinical_path = 'tcga_synthea_test/tcga_tumor_clinical.csv',
  project = "TCGA-BRCA", # param:get_gdc_rna_data 
  tissue = "Primary Tumor", # param:get_gdc_rna_data
  custom_dir = 'GDCdata_BRCA', # param:get_gdc_rna_data; 
  chunk_size = 100,  #  Adjust based on your memory capacity
  sample_size = 2 # less than chunk size
  ) {

  # browser()
  
  setwd(out_dir)
  
  clinical <- get_gdc_clinical_data(project)  # can remove download parameter since it will find them if custom dir is same as last downloaded Of the 1183 files for download 1183 already exist.   All samples have been already downloaded
  # only identifier for case is submitter_id  TCGA-E9-A1NF  
  # head(clinical[[1]][[1]], 2)
  
  barcode_submitter_mapping <- get_gdc_rna_data(project, tissue, download = TRUE, customdir = custom_dir)
  # add column to clinical data for TCGA submitter_id so we can extract just the subset we choose from the RNAseq data
  clinical <- merge(clinical, barcode_submitter_mapping, by = "submitter_id")
  
  # barcode_chunks <- split(barcodes$submitter_id, ceiling(seq_along(barcodes$submitter_id) / chunk_size))
  # barcode_chunks <- split(barcodes, ceiling(seq_along(barcodes) / chunk_size))
  barcode_chunks <- split(barcode_submitter_mapping$barcode, ceiling(seq_along(barcode_submitter_mapping$barcode) / chunk_size))
  
  # results_list <- lapply(barcode_chunks, function(chunk) {
  #   process_chunk(chunk, "TCGA-BRCA", "Primary Tumor", "path/to/directory")
  # })
  # Processing only the first chunk
  first_chunk_result <- process_chunk(barcode_chunks[[1]], project, tissue, custom_dir)  # "RangedSummarizedExperiment"   just like result of main.get_gdc_data()
  # processing rna data before merge to demographic data
  # subset: unnec
  first_five_sample_ids_first_chunk_result <- colnames(first_chunk_result)[1:sample_size]  # limit the samples we will retrieve
  data_first_five_sample_ids_first_chunk_result <- first_chunk_result[, first_five_sample_ids_first_chunk_result]
  # assay_tumor_data_subset = process_and_write_csv(data_first_five_sample_ids_first_chunk_result, "first_chunk_result_subset.csv", TRUE)
  assay_tumor_data_subset = process_and_write_csv(data_first_five_sample_ids_first_chunk_result, write_rawcount_path, TRUE)

  
  #  parsing clinical data ? 
  # get the corresponding data from clinical
  submitter_ids_from_gene_exp <- colnames(assay_tumor_data_subset)  # same as first_five_sample_ids_first_chunk_result
  
  # Initialize an empty list to store the clinical data for each submitter_id
  matched_clinical_data <- list()
  
  # Loop through each submitter_id and get the corresponding clinical data
  for(id in submitter_ids_from_gene_exp) {
    # Find the row in the clinical dataset that matches the submitter_id
    matched_row <- clinical[clinical$barcode == id, ]
    
    # Add the matched row to the list
    matched_clinical_data[[id]] <- matched_row
  }
  combined_clinical_data <- do.call(rbind, matched_clinical_data)
  # this contains all of the clinical variables and has the correct Barcode identifier as rowname
  # we need to exl
  # save this as a CSV in the final dir  tcga_clinical.csv
  # 
  tcga_keep_clinic <- c('ajcc_pathologic_stage',  'primary_diagnosis',
                       'ajcc_pathologic_t', 'morphology', 'ajcc_pathologic_n', 'ajcc_pathologic_m',
                       'race',  'ethnicity', 'vital_status', 'gender', 'year_of_birth')
  # drop vital     add 
  # tcga_keep_clinic <- c('ajcc_pathologic_stage',  'primary_diagnosis', 
  #                       'ajcc_pathologic_t', 'morphology', 'ajcc_pathologic_n', 'ajcc_pathologic_m', 
  #                       'race',  'ethnicity')
  
  # Selecting the specified columns from combined_clinical_data
  selected_data <- combined_clinical_data %>% select(all_of(tcga_keep_clinic))
  
  # Save the new dataframe as a CSV file
  write.csv(selected_data, write_clinical_path, row.names = TRUE)
  return(list(assay_tumor_data_subset=assay_tumor_data_subset, clinical=clinical))
}


tumor_ass_clin <- main()  # default settings
assay_tumor_data_subset = tumor_ass_clin$assay_tumor_data_subset
clinical = tumor_ass_clin$clinical


norm_ass_clin <- main(
    out_dir = "/home/john/Desktop/cancer_vis",
    write_rawcount_path = "tcga_synthea_test/tcga_normal_rawcount.csv",
    write_clinical_path = 'tcga_synthea_test/tcga_normal_clinical.csv',
    project = "TCGA-KIRC",
    tissue = "Solid Tissue Normal",
    custom_dir = 'GDCdata_kirc',
    chunk_size = 100,
    sample_size = 2 # less than chunk size
)  

norm_data_subset = norm_ass_clin$assay_tumor_data_subset
norm_clinical = norm_ass_clin$clinical 


# Task: load gene names for pathway
get_mapped_ids_for_pathway <- function(pathway_id) {
  # Get the gene names for the pathway
  pathway_info <- keggGet(pathway_id)
  gene_info <- pathway_info[[1]][["GENE"]]
  genes <- gene_info[seq(2, length(gene_info), by = 2)]
  gene_names <- sapply(strsplit(genes, ";"), function(x) x[1])
  
  # Map gene symbols to Ensembl IDs
  mapped_ids <- mapIds(org.Hs.eg.db, 
                       keys = gene_names, 
                       column = "ENSEMBL", 
                       keytype = "SYMBOL",
                       multiVals = "first")
  
  return(mapped_ids)
}


combine_expression_data_and_get_log_fold_changes <- function(expression_data_matrix, expression_data_label, 
                                                             additional_expression_matrix, additional_expression_label, mapped_ids) {
  # browser()

  # Extract column names (sample identifiers) from both expression data matrices
  expression_names <- colnames(expression_data_matrix)
  additional_names <- colnames(additional_expression_matrix)
  
  # Combine the column names from both datasets to create a unified list of sample names
  sample_names <- c(expression_names, additional_names)
  
  # Create a data frame with condition labels for each sample
  sample_info <- data.frame(
    condition = c(rep(expression_data_label, length(expression_names)), 
                  rep(additional_expression_label, length(additional_names))),
    row.names = sample_names
  )
  
  # Combine the two expression data matrices side-by-side (column binding)
  combined_expression_data <- cbind(expression_data_matrix, additional_expression_matrix)
  
  # Reorder columns in the combined matrix to match the order of sample names
  combined_expression_data <- combined_expression_data[, sample_names]
  
  # Filter for genes of interest based on mapped_ids
  combined_expression_data <- combined_expression_data[unname(mapped_ids), ] # this also changes to as.int so below is redundant but leave this like this
  # comment to get all genes  UNCOMMENT THIS ^--^ IF USING FOLD CHANGES FOR PATHWAY DIAGRAM. 
  gene_ids <- rownames(combined_expression_data)
  combined_expression_data <- apply(combined_expression_data, 2, as.integer)  # unnec now since mapped_ids is all genes in this case
  rownames(combined_expression_data) <- gene_ids  # After the conversion, reassign the stored row names back to combined_expression_data.
  
  
  
  
  # Create a DESeqDataSet object for differential expression analysis
  dds <- DESeqDataSetFromMatrix(
    countData = combined_expression_data,
    colData = sample_info,
    design = ~ condition
  )
  
  # Perform differential expression analysis using DESeq2
  dds <- DESeq(dds)
  # Extract the results of the differential expression analysis
  results <- results(dds)
  
  # Performs its own normalization as part of the differential expression analysis process. It uses size factors to normalize the counts data, which helps to account for differences in sequencing depth and RNA composition between samples. Therefore, if you're using DESeq2, you generally don't need to perform separate normalization on your RNA-seq data before feeding it into this pipeline. The normalization is an integral part of the DESeq2 workflow.
  # If you want to output the normalized expression data from a DESeq2 analysis, you can use the counts() function from the DESeq2 package with the argument normalized=TRUE. This function returns the matrix of normalized count data, where normalization is done using the size factors computed by DESeq2.
  normalized_counts <- counts(dds, normalized = TRUE)
  
  # Extract log2 fold changes from the results for further analysis or visualization
  # log_fold_changes <- results$log2FoldChange
  # names(log_fold_changes) <- rownames(results)
  
  # Return the log2 fold changes
  # return(log_fold_changes)
  
  # Return both log2 fold changes and normalized count data
  return(list(results=results, normalized_counts = normalized_counts))
}


library(KEGGREST)
library(AnnotationDbi)
library(org.Hs.eg.db)

pathway_id = "hsa04010"

# for oom Clustering Analysis
# mapped_ids <- get_mapped_ids_for_pathway(pathway_id)

# for GSEA
mapped_ids <- rownames(assay_tumor_data_subset)


library(DESeq2)

resultsAndNormalized <- combine_expression_data_and_get_log_fold_changes(assay_tumor_data_subset, "melanoma", norm_data_subset, "normal", mapped_ids)
# log_fold_changes <- combine_expression_data_and_get_log_fold_changes(assay_tumor_data_subset, "melanoma", assay_normal_data_subset, "normal", mapped_ids)


results <- resultsAndNormalized$results
log_fold_changes <- results$log2FoldChange
names(log_fold_changes) <- rownames(results)

normalized_counts <- resultsAndNormalized$normalized_counts

# 1. Differential Expression Analysis
# Exporting results 
resOrdered <- results[order(results$padj),]
write.csv(as.data.frame(resOrdered), file="DESeq2_results.csv")  # 5.0 MB  60616 genes

# Get most important genes affected in this cancer
# Ordering results by p-value
orderedResults <- results[order(results$padj),]

# Filter out rows where padj is NA
cleanResults <- orderedResults[!is.na(orderedResults$padj), ]
# Extracting gene identifiers from cleanResults
cleanResultsIdentifiers <- rownames(cleanResults)
# Counting the number of genes
numberOfGenes <- length(rownames(cleanResults))  # 23736

# Now, subset for significant results
sigResults <- cleanResults[cleanResults$padj < 0.01, ]
numberOfsigResults <- length(rownames(sigResults))  # < .05 =2280   < 0.01 = 1325


# Filtering for significant genes based on log2 fold change
topGenes <- sigResults[abs(sigResults$log2FoldChange) > 1, ] #  These genes are considered significant in terms of differential expression.
# numberOftopGenes <- length(rownames(topGenes))

# Extracting gene identifiers
geneIdentifiers <- rownames(topGenes)
geneInfo <- sigResults[geneIdentifiers, ]


library(clusterProfiler)
# Loading required libraries
library(org.Hs.eg.db)

# Converting gene identifiers (assuming they are ENSEMBL IDs) to Entrez IDs
entrezIDs <- bitr(geneIdentifiers, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Performing KEGG enrichment analysis (note: adjust parameters as needed)
# This analysis aims to identify pathways that are significantly enriched in your list of significant genes.
keggResult <- enrichKEGG(gene = entrezIDs$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)  # enriched pathways
# Mapping Genes to Pathways: After getting the enriched pathways (keggResult), the goal is to 
# map these pathways back to your significant genes. However, it's important to note that not all significant genes may be part of the identified enriched pathways, which could be why some genes in topGenes don't have an associated pathway.

# Convert enrichResult object to a dataframe
# keggResult <- as.data.frame(keggResult)
# Convert enrichResult object to a dataframe
keggResult_df <- as.data.frame(keggResult)

# Check the column names to know which ones to use
colnames(keggResult_df)

# Extracting pathway information for each gene
pathways <- as.data.frame(keggResult)[, c("ID", "Description", "GeneRatio")]

# Adding KEGG pathway information to geneInfo
# Extract gene identifiers from topGenes
geneIdentifiers <- rownames(topGenes)

# Assuming you have already converted these to Entrez IDs and have the keggResult
# Initialize a column for pathways in topGenes
topGenes$Pathway <- NA  # Initialize the Pathway column

forTopGenes <- function(topGenes) {
  for (i in 1:nrow(topGenes)) {
    geneEntrezIDs <- entrezIDs$ENTREZID[entrezIDs$ENSEMBL == rownames(topGenes)[i]]
    pathwayDescriptions <- c()
    
    for (geneEntrezID in geneEntrezIDs) {
      if (!is.na(geneEntrezID)) {
        for (j in 1:nrow(keggResult_df)) {
          splitGeneIDs <- strsplit(as.character(keggResult_df$geneID[j]), "/")[[1]]
          if (geneEntrezID %in% splitGeneIDs) {
            pathwayDescriptions <- unique(c(pathwayDescriptions, keggResult_df$Description[j]))
          }
        }
      }
    }
    
    topGenes$Pathway[i] <- if(length(pathwayDescriptions) > 0) paste(pathwayDescriptions, collapse = "; ") else NA
  }
  topGenesdf <- (as.data.frame(topGenes))
  return(topGenesdf) 
}
# If a gene is part of a pathway in keggResult, that pathway's description is added to the gene in topGenes.
topGenesdf <- forTopGenes(topGenes)
# ok then get rid of NA's 
# so we get only the significant genes in the significant pathways

# Filter out rows where the Pathway is NA
topGenesdf_filtered <- topGenesdf[!is.na(topGenesdf$Pathway), ]

# Splitting the pathway information, as one gene can be associated with multiple pathways
topGenesdf_filtered$Pathway <- strsplit(as.character(topGenesdf_filtered$Pathway), "; ")

# Using rownames as the gene identifiers
normalizedPathways <- stack(setNames(topGenesdf_filtered$Pathway, rownames(topGenesdf_filtered)))

# Counting the number of occurrences of each pathway
pathwayCounts <- table(normalizedPathways$values)
pathwayCounts <- as.numeric(table(normalizedPathways$values))

# Unlist the pathways to normalize the data
unlistedPathways <- unlist(topGenesdf_filtered$Pathway)

# Counting the number of occurrences of each pathway
pathwayCounts <- table(unlistedPathways)

# Convert to a named numeric vector
# pathwayCounts <- as.numeric(pathwayCounts)
# names(pathwayCounts) <- names(attr(pathwayCounts, "dimnames")[[1]])


# Assuming keggResult_df contains the pathway IDs in a column named 'ID'
library(KEGGREST)

pathwayIDs <- keggResult_df$ID

# Function to get gene count for a pathway
getGeneCountForPathway <- function(pathwayID) {
  # Retrieve pathway details
  pathwayDetails <- tryCatch(keggGet(pathwayID), error = function(e) return(NA))
  
  # If pathwayDetails is not NA, extract the number of genes
  if (!is.na(pathwayDetails)) {
    genes <- pathwayDetails[[1]]$GENE
    return(length(genes))
  } else {
    return(NA)
  }
}

# Apply the function to each pathway ID
pathwayGeneCounts <- sapply(pathwayIDs, getGeneCountForPathway, USE.NAMES = TRUE)

# Check the gene counts
head(pathwayGeneCounts)

namedTotalGenes <- setNames(pathwayGeneCounts, names(pathwayCounts))

# Calculating proportions
# Ensure that the names in pathwayCounts match those in namedTotalGenes
proportions <- pathwayCounts / namedTotalGenes[names(pathwayCounts)]

# proportions <- pathwayCounts / totalGenesPerPathway

# Creating a dataframe for pathway ranking
pathwayRanking <- data.frame(Pathway = names(proportions), Proportion = proportions)
pathwayRanking <- pathwayRanking[order(-pathwayRanking$Proportion), ]








# 2. Gene Set Enrichment Analysis (GSEA)
# Using clusterProfiler:

# Preparing the input data
# Assuming `geneList` is a vector of genes ranked by log fold changes
geneList <- sort(log_fold_changes, decreasing = TRUE)
entrez_ids <- mapIds(org.Hs.eg.db, 
                     keys = names(geneList), 
                     column = "ENTREZID", 
                     keytype = "ENSEMBL", 
                     multiVals = "first")

# Update the names with Entrez IDs, filtering out NA values
valid_entrez_ids <- entrez_ids[!is.na(entrez_ids)]
# geneList <- geneList[names(geneList) %in% names(valid_entrez_ids)]
names(geneList) <- valid_entrez_ids

# Filtering out NA values and handling duplicates
# valid_entrez_ids <- entrez_ids[!is.na(entrez_ids)]
# unique_entrez_ids <- unique(valid_entrez_ids)
# geneList <- geneList[names(geneList) %in% names(unique_entrez_ids)]
# names(geneList) <- unique_entrez_ids
# Ensure no duplicate names exist
# geneList <- geneList[!duplicated(names(geneList))]

# Filtering out NA values and handling duplicates
# valid_entrez_ids <- entrez_ids[!is.na(entrez_ids)]
unique_entrez_ids <- unique(valid_entrez_ids)

# Create a named vector for subsetting geneList
geneList_subset <- geneList[names(geneList) %in% (unique_entrez_ids)]



# Run GSEA
gseaResult <- gseGO(geneList = geneList_subset, 
                    ont = "BP", 
                    nPerm = 1000, 
                    minGSSize = 10, 
                    maxGSSize = 500, 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 0.05, 
                    verbose = FALSE,
                    OrgDb = org.Hs.eg.db)


# Viewing results
head(gseaResult)

# Assuming pathway IDs are in a column named 'ID' or similar
slotNames(gseaResult)
# [1] "result"      "organism"    "setType"     "geneSets"    "geneList"    "keytype"    
# [7] "permScores"  "params"      "gene2Symbol" "readable"    "termsim"     "method"     
# [13] "dr" 

gseaData <- gseaResult@result
head(gseaData)

# getting significantly affected pathways
# Using Normalized Enrichment Score (NES)
# To select the top 20 pathways based on the absolute value of NES (indicating the strength of enrichment):
# Sort by absolute NES and take the top 20
top_pathways_by_NES <- gseaData[order(abs(gseaData$NES), decreasing = TRUE), ][1:100, ]

# Using Adjusted P-Value
# To select the top 20 pathways based on statistical significance:
# Sort by adjusted p-value and take the top 20
# top_pathways_by_padj <- gseaData[order(gseaData$p.adjust), ][1:100, ]
# this will be identicalt to top_pathways_by_NES because gseaData$p.adjust and gseaData$qvalue the same for each gene for some reason

# Extracting pathway IDs (replace 'ID' with the actual column name for pathway IDs)
top_pathway_ids_by_NES_id <- top_pathways_by_NES$ID
# top_pathway_ids_by_padj_id <- top_pathways_by_padj$ID

# Print the pathway IDs
print(top_pathway_ids_by_NES_id)
# print(top_pathway_ids_by_padj_id)

# Remember that a high NES represents a strong enrichment signal (either up or down), and a low adjusted p-value represents high statistical significance.

# Convert GO IDs to KEGG IDs
# gets all the gene names from that pathway
selectedPath <- top_pathway_ids_by_NES_id[2] # [3] errs Error in .testForValidKeys(x, keys, keytype, fks) : 
  # None of the keys entered are valid keys for 'GO'. Please use the keys method to see a listing of valid arguments.
kegg_ids_NES <- bitr(selectedPath, fromType="GO", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
head(kegg_ids_NES$ENSEMBL)
ensembl_ids <- kegg_ids_NES$ENSEMBL

# Extract the log fold changes for these ENSEMBL IDs
selected_logFC <- log_fold_changes[names(log_fold_changes) %in% ensembl_ids]

# Create a data frame for plotting
plot_data <- data.frame(ensembl = names(selected_logFC), logFC = selected_logFC)

# Visualization using ggplot2
# of log fold changes between tumor and normal gruops
library(ggplot2)
ggplot(plot_data, aes(x = ensembl, y = logFC)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Gene Expression Changes for Selected GO Term", x = "Gene ID", y = "Log Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x labels for readability
# makes nice bar chart for the selected pertinent pathway

# Biomarker id



# Corrlation analysis




# 5. Network Analysis
# Using STRINGdb for network analysis:
library(STRINGdb)

# Assuming `geneList` is your list of genes
string_db <- STRINGdb$new(version="11.0", species=9606, score_threshold=400)
network <- string_db$get_network(geneList)
string_db$plot_network(network)



# 6. Survival Analysis
# Using Cox proportional hazards model:
library(survival)

# need to include for patients in the data to get anything meaningful
# surv_data <- read from the clinical csv for BRCA patients only
# if not x$Death  is not "Alive"  (its datetime of death)  
# also use the duration of survival after diagnosis

# Assuming `surv_data` is your survival data with time, status, and gene expression
cox <- coxph(Surv(time, status) ~ gene_expression, data=surv_data)
summary(cox)

# Plotting survival curves
plot(survfit(cox), xlab="Time", ylab="Survival Probability")

  

# 7 Immune Infiltration Analysis:
  # In cancer research, understanding the tumor microenvironment, especially the immune contexture, can be valuable.
# Tools: CIBERSORT, xCell, or similar tools for assessing immune cell infiltration based on gene expression data.





# Visualization Strategies for Multiple Pathways
#######################################################
#   Heatmaps: Displaying multiple pathways in a heatmap allows you to compare gene expression across different pathways. You can structure the heatmap with genes on one axis and GO pathways on the other.
#     Aggregate Data: Combine the gene lists from multiple pathways. Be aware of any overlap in genes across the selected pathways.
#     Prepare Visualization Data: Depending on the visualization type, you may need to aggregate or transform your data. For instance, for a heatmap, create a matrix with genes as rows and pathways as columns.

#%$
# use this code block for Biomarker identification
# Gene Ontology (GO) Enrichment:
  # Identify overrepresented GO terms in your list of differentially expressed genes to understand the biological processes, molecular functions, and cellular components involved.
#%$

library(pheatmap)

# Assuming 'plot_data' has columns 'ensembl', 'logFC', and 'GO_term' for multiple pathways
selectedPath <- top_pathway_ids_by_NES_id[5]  # group of pathways
kegg_ids_NES <- bitr(selectedPath, fromType="GO", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
head(kegg_ids_NES$ENSEMBL)
ensembl_ids <- kegg_ids_NES$ENSEMBL
unique_ensembl_ids <- unique(ensembl_ids) # Remove Duplicates: Keep only unique ENSEMBL IDs for your analysis.
# it's appropriate to focus on the unique set of ENSEMBL IDs, as the duplicates do not provide additional independent information for your analysis.

# Extract the log fold changes for these ENSEMBL IDs
selected_logFC <- log_fold_changes[names(log_fold_changes) %in% unique_ensembl_ids]

# Filter the kegg_ids_NES to include only unique ENSEMBL IDs
filtered_kegg_ids_NES <- kegg_ids_NES[kegg_ids_NES$ENSEMBL %in% unique_ensembl_ids, ]

# Remove duplicates based on ENSEMBL ID
filtered_kegg_ids_NES_unique <- filtered_kegg_ids_NES[!duplicated(filtered_kegg_ids_NES$ENSEMBL), ]

# Filter to include only ENSEMBL IDs present in selected_logFC
filtered_kegg_ids_NES_corrected <- filtered_kegg_ids_NES_unique[filtered_kegg_ids_NES_unique$ENSEMBL %in% names(selected_logFC), ]

# Ensure that filtered_kegg_ids_NES_corrected has the same number of rows as selected_logFC
if (nrow(filtered_kegg_ids_NES_corrected) != length(selected_logFC)) {
  stop("Mismatch in the number of rows after filtering")
}


# Create plot_data using the aligned and filtered data
plot_data <- data.frame(
  GO = filtered_kegg_ids_NES_corrected$GO, 
  ensembl = names(selected_logFC), 
  logFC = selected_logFC
)

# Proceed with the heatmap creation
heatmap_data <- reshape2::dcast(plot_data, ensembl ~ GO, value.var = "logFC")

# Set ENSEMBL IDs as row names
row.names(heatmap_data) <- heatmap_data$ensembl

# Remove the 'ensembl' column
# heatmap_data <- heatmap_data[, -1]
# Remove the 'ensembl' column but retain the structure as a data frame
heatmap_data <- heatmap_data[, -1, drop = FALSE]

# Check the class again
class(heatmap_data)

# Convert to matrix
heatmap_matrix <- as.matrix(heatmap_data)


# Replace NA values with zeros (or another appropriate value)
heatmap_matrix[is.na(heatmap_matrix)] <- 0
# only one pathway
pheatmap(heatmap_matrix, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE)



# # Create the heatmap
# # pheatmap(heatmap_matrix, scale = "row")
pheatmap(heatmap_matrix, scale = "none")
# 
# # Example of log transformation, adding a small constant to avoid log(0)
# heatmap_matrix_transformed <- log(heatmap_matrix + 1)
# pheatmap(heatmap_matrix_transformed, scale = "row")
# 
# # Example of a transformation that handles negative values
# heatmap_matrix_transformed <- sqrt(heatmap_matrix^2)
# pheatmap(heatmap_matrix_transformed, scale = "row")
# 
# 
# 
# # Check for NA, NaN, or Inf values in heatmap_data
# if (any(is.na(heatmap_data)) || any(is.nan(heatmap_data))) {
#   # Handle problematic values here
#   # For example, you can replace them with zeros or mean values
#   # This is a simple example of replacing them with zero
#   heatmap_data[is.na(heatmap_data) | is.nan(heatmap_data)] <- 0
# }
# 
# row.names(heatmap_data) <- heatmap_data$ensembl
# heatmap_data <- heatmap_data[,-1]
# 
# pheatmap(heatmap_data, scale = "row")

######################################################################################################



#   Circos Plots: If you include multiple pathways, circos plots can be very effective in showing the relationships and overlaps between genes and pathways.

#   Network Diagrams: Visualize the interconnections between pathways and genes. This can be particularly informative for understanding complex biological processes and pathway interactions.

#   GO Bubble Plots: These can represent multiple GO terms, with the size of each bubble corresponding to the number of genes in the term and the color representing the significance or average log fold change.



# other ideas or GSEA Vis


# Interactive Network Visualization: Tools like Cytoscape, which can be integrated with R, 
# allow for interactive network visualizations. You can create networks where nodes represent 
# GO terms or genes, and edges represent relationships or shared genes.

# Combining Multiple Data Types: If you have access to additional data types (like proteomics, metabolomics, etc.), integrating these with your gene expression data can provide a more comprehensive view.

# Advanced Plotting with ggplot2 and Extensions: Utilizing ggplot2 with extensions like ggrepel for better label management or ggforce for more complex layouts can enhance your visualizations.

# Using ReactomePA for Pathway Analysis: Another R package, ReactomePA, is designed for pathway analysis and visualization and might offer more suitable tools for your needs.


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##############################################

# we use statistics to get potentially related pathways but it relly must be manually curated but in our case thats unfeasible
# Notes:
  # Domain Knowledge: This process relies heavily on domain knowledge. If you're not familiar with the pathology, collaborating with an expert is advisable.
  # Dynamic Process: The selection of relevant pathways might evolve as new information becomes available or as your understanding of the pathology deepens.



# 4. Clustering Analysis
# Creating a heatmap:
library(gplots)
# Simplify the Plot: Reduce the complexity of the plot by avoiding the computation of dendrograms:
heatmap.2(as.matrix(normalized_counts), 
          scale="row", 
          Colv=NA, 
          Rowv=NA, 
          dendrogram="none", 
          trace="none", 
          margin=c(10, 6))

library(ComplexHeatmap)
Heatmap(as.matrix(normalized_counts))

# heatmap.2 function can be very memory-intensive due to the computation of dendrograms and the internal handling of data. Here are some strategies to mitigate this issue:
# Reduce Data Size: If your dataset is very large, consider working with a subset of the data, focusing on 
# the most relevant genes or samples. You can select a subset based on variance, expression levels, 
# or biological relevance.
# Assuming `data` is your normalized expression data
heatmap.2(as.matrix(normalized_counts), 
          scale="row", 
          Colv=NA, 
          dendrogram="row", 
          trace="none", 
          margin=c(10, 6))
# Error in plot.new() : figure margins too large
# Error in par(op) : invalid value specified for graphical parameter "pin"

# this works well 
png("heatmap.png", width=2000, height=2000)
heatmap.2(as.matrix(normalized_counts), 
          scale="row", 
          Colv=NA, 
          dendrogram="row", 
          trace="none", 
          margin=c(10, 6))
dev.off()

# this looks visually wrong somehow too many lines not enough genes
heatmap.2(as.matrix(normalized_counts), 
          scale="row", 
          Colv=NA, 
          dendrogram="row", 
          trace="none", 
          margin=c(5, 5))  # Adjust these values as needed



# this contains all of the genes as ensembl and has the correct Barcode identifier as colname

# get log fold changes or whatever normalization

# done proc RNA data -- merge to demopgrahic
# for each of the gene data

# TODO SHINY
# NEW GITHUB REPO
# ADD DEA (DONT ALLOW IF NOT RUNNING LOCALLY  SHOW SAMPLE OUT) AND OTHER ANALYSIS OPTIONS

# ABOUT TAB (AT BOTTOM OF SIDEBAR) THAT DESCRIBES HOW TO USE ETC

# Now we have a directory complete with all of our samples and the integrated Synthea data: each file can be refered to by sample_id
# observations.csv
# reports.csv
# careplans.csv
# patient.csv
# documents.csv
# encounters.csv
# medications.csv
# procedures.csv
# allergies.csv
# condition.csv
# vaccinations.csv
# tcga_clinical.csv
# tcga_rawcount.csv

# tcga_clinical  tooltips 
# ajcc_pathologic_stage: The extent of a cancer, especially whether the disease has spread from the original site to other parts of the body based on AJCC staging criteria.
# primary_diagnosis: Text term used to describe the patient's histologic diagnosis, as described by the World Health Organization's (WHO) International Classification of Diseases for Oncology (ICD-O).
# ajcc_pathologic_t: Code of pathological T (primary tumor) to define the size or contiguous extension of the primary tumor (T), using staging criteria from the American Joint Committee on Cancer (AJCC).
# morphology: The third edition of the International Classification of Diseases for Oncology, published in 2000 used principally in tumor and cancer registries for coding the site (topography) and the histology (morphology) of neoplasms. The study of the structure of the cells and their arrangement to constitute tissues and, finally, the association among these to form organs. In pathology, the microscopic process of identifying normal and abnormal morphologic characteristics in tissues, by employing various cytochemical and immunocytochemical stains. A system of numbered categories for representation of data. 
# ajcc_pathologic_n: 	The codes that represent the stage of cancer based on the nodes present (N stage) according to criteria based on multiple editions of the AJCC's Cancer Staging Manual.
# ajcc_pathologic_m: Code to represent the defined absence or presence of distant spread or metastases (M) to locations via vascular channels or lymphatics beyond the regional lymph nodes, using criteria established by the American Joint Committee on Cancer (AJCC).



# test new analytical methods for the data
# 
# Identifying Potential Biomarkers:
#   
#   Use statistical and bioinformatics methods to identify potential biomarkers in the RNA-seq data. This could involve differential gene expression analysis, survival analysis, and pathway analysis.
# Tools like DESeq2 or edgeR in R can be useful for differential expression analysis.
# For survival analysis, you can use the survival package in R.
# Additional Analyses:
#   
#   Machine Learning: Employ machine learning techniques to predict clinical outcomes or to classify patients based on their gene expression profiles. R packages like caret or Python libraries like scikit-learn can be useful here.
# Pathway Analysis: Use tools like GSEA (Gene Set Enrichment Analysis) to understand the biological pathways involved.
# Clustering and Visualization: Perform cluster analysis (like hierarchical clustering) on gene expression data to identify patterns or subgroups within the data.

# Correlation with Clinical Outcomes:
# genetics leading to primary diagnosis
# Link the genetics data with the corresponding clinical information, such as age at diagnosis, tumor stage, and primary diagnosis. This will allow you to analyze genetic markers in the context of clinical outcomes.

# Demographic Analysis:
#     Explore the distribution of genetic markers across different demographics (like age, gender, race) to uncover patterns or disparities in cancer genetics.

# Predictive Modeling:
# Develop predictive models that use genetic and clinical data to forecast patient outcomes or treatment responses.

# Visual Analytics:
# Create interactive visualizations that allow users to explore the relationships between genetics, clinical characteristics, and outcomes.



# first_five_sample_ids_normal <- colnames(data_all_primary_normal)[1:2]
# data_all_primary_normal_subset <- data_all_primary_normal[, first_five_sample_ids_normal]
# assay_normal_data_subset = process_and_write_csv(data_all_primary_normal_subset, "primary_normal_test_rawcounts.csv", TRUE)




# Notes

# 3. Integrating with Synthea:
  # Synthea generates synthetic patient data. To use TCGA data for setting parameters in Synthea:
  #Parameter Mapping: Map the demographic data from TCGA to Synthea's patient generation parameters. This could involve aligning age groups, race distributions, etc., to ensure the synthetic patients reflect the demographics of your TCGA data.
# Clinical dictionary https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-entity-list&anchor=clinical
## demog https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=demographic
## diag https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=diagnosis&_top=1
# noted relevant to using with Synthea
# relevant_clinic <- c('ajcc_pathologic_stage', 'days_to_last_follow_up', 'primary_diagnosis', 
#                      'ajcc_pathologic_t', 'morphology', 'ajcc_pathologic_n', 'ajcc_pathologic_m', 
#                      'race', 'gender', 'ethnicity', 'vital_status', 'age_at_index', 'days_to_death')

# data poking
# Count the frequency of each gender
gender_distribution <- table(clinical$gender)  # female   male  1098     12 
morphology <- table(clinical$morphology)
primary_diagnosis <- table(clinical$primary_diagnosis)
table_df <- table(clinical$morphology, clinical$primary_diagnosis)
# note all primary_diagnosis each only have 1 type of morphology => all Infiltrating duct carcinoma, NOS 788  have morph 8500/3 788    8520/3 204

ajcc_pathologic_t <- table(clinical$ajcc_pathologic_t)
vital_status <- table(clinical$vital_status)

# Create a bar plot of the gender distribution
barplot(gender_distribution, main="Gender Distribution", xlab="Gender", ylab="Frequency")

# Create a contingency table

# Print the table
print(table_df)

# there is little info on treatments so we just have to work with : have yes/no/NA but might be useful
#   treatments_pharmaceutical_treatment_or_therapy    treatment_or_therapy in GDGdict
#  treatments_radiation_treatment_or_therapy

  # not used but maybe useful: age_at_diagnosis (days) (instead us age_at_index-its in years)
# other?  diagnosis_id    
# considerations 
# can probably catagorize statify the dataset by 
# primary_diagnosis
# ajcc_staging_system_edition  https://chat.openai.com/share/e9e688d6-6c93-47d0-9515-e5c18bbefe21
#   considering the cancer staging system is crucial, ESPECIALLY IF YOU'RE USING STAGING AS A VARIABLE IN YOUR ANALYSIS
#   Changes in staging criteria can reflect advancements in understanding the disease, which might influence treatment decisions or prognostic assessments.
#   Different editions of the AJCC manual may have significant changes in how stages are defined. For example, what was classified as Stage II in one edition might be classified differently in another. This can affect survival analysis, risk stratification, and other outcomes.
# if we are using staging as a variable  which we will for synthea it is probably critical, 
# Create Subgroups: Divide your dataset into subgroups based on the AJCC edition used for each case.
# remove NA's for edition
# Analyze Separately: Conduct your analysis separately for each subgroup. This allows you to assess how staging (and potentially other variables) affects outcomes under each staging system.


# Synthea dictionary 
#   testing Synthea use 
# cd ~/Desktop/cancer_vis/synthea/synthea
# ./run_synthea -p 10 -g F -a 40-50
# does NOT get any local/sub modules ./run_synthea -g F -a 40-60 -c ./src/main/resources/synthea_brca.properties -d ./src/main/resources/modules/brca/
# ./run_synthea -g F -a 40-60 -c ./src/main/resources/synthea_brca.properties
# ./run_synthea -g F -a 40-60 -c ./src/main/resources/synthea_brca.properties -k ../../keep_brca.json

# ./run_synthea -p 10 -k ../../keep.json
# document: this is the final synthea run; trick was to use the Synthea to generate the keep file
# ./run_synthea -p 10 -g F -a  -k ../../keep.json
# get F only

# FIRST THING TO DO IN SHINY 
#   Create dashboard for viewing the patient data for the selected patient only 
#     the various csvs can just be contained in tabs individually
# SECOND THING TO DO IN SHINY
#   analysis tab of all the loaded patients 

# other    we still need a way of getting Noraml tissue patients more easily
# add to git

# from TCGA we MUST replace/retain in the Synthea data  all data we can from the samples: MUST MATCH SYNTHEA DATA: gender, vital_status, 
  # doesnt need to match/will be hard to find matches: age+-3, race
# from TCGA we just need to KEEP: 
  # pathologic, stage data, 
# we CAN keep race, ethnicity

#not  'days_to_last_follow_up', 'gender' (we are only getting F from Synthea anyway), 'age_at_index', 'days_to_death'

# note in documentation exactly what data is true to the gene data and which is artificial from Synthea
# document where to go for definition of factors fr synthea and tcga

# we are keeping all TCGA data and adding to it the 
#data from Synthea: 
# Name, Address, City, State, Postal Code

# TCGA we GET RID OF anything conflicting with PROCEDURE or MED date-times, 
# days_to_last_followup

# WE NEED TO KEEP AL SYNTHEA DATE AND AGES

# documentatino 
#  i want the useres to be able to get their own data similarly so we give them the tools/functions and they can demo other..
#  select ur active condition of interest from Synthea api creator
#   select FHIR for outut 
#   run ./
#   in python script select the dir with the fhir file
#     *generate CSV file for our use where the colnames are barcode, condition, condition_x, etc 
#  get the TCGA data u want 
#  run R script to create the CSV and gene expression CSV

# now you have a dir with  
# from TCGA...
# gene expression csv   colnames are samples in the form of barcodes
# tcga_clinical.csv
# ...then from Synthea
# conditions    colnames are 
# 
# 




# Skip Data Analysis until integrate with Synthea: Once the data is merged and organized as per your requirements, you can proceed with your specific analyses, like differential gene expression analysis, survival analysis, etc.

# goal is to generate a unique set of synthetic data for each case in the TCGA dataset, using fixed records in Synthea would indeed be the more appropriate approach. This method allows for greater precision and control in aligning the synthetic data with specific characteristics of each TCGA case. Here's a detailed plan on how to proceed with this approach:

# Keep is the solution, Fixed records can only control the beginning of the module state
# Selective Filtering: The keep module allows you to define criteria that synthetic patients must meet. This can include specific demographic information (like age, gender, race), as well as clinical criteria (such as cancer staging).
# Flexibility in Criteria: You can set a wide range of criteria in a keep module, which gives you greater flexibility in filtering the synthetic data. This is particularly useful if you are dealing with complex or specific case profiles.

# 3 Generate Synthetic Records:
#   
# Data Validation and Quality Checks:
# Integration with Your R/Shiny Project:
#   
#   Import the generated synthetic data into your R/Shiny application.
# Use this data in conjunction with the TCGA RNASeq data and other relevant information in your interactive cancer pathway visualizer.


# TCGA extras   https://www.costalab.org/wp-content/uploads/2020/11/R_class_D3.html
# We can check all the available projects at TCGA with the command bellow. Since there are many lets look at the first 6 projects using the command head().
GDCprojects = getGDCprojects()

head(GDCprojects[c("project_id", "name")])
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")
# 
query_TCGA = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts")

# To visualize the query results in a more readable way, we can use the command getResults.
lihc_res = getResults(query_TCGA) # make results as table
# head(lihc_res) # data of the first 6 patients.
colnames(lihc_res) # columns present in the table

head(lihc_res$sample_type) # first 6 types of tissue.

summary(factor(lihc_res$sample_type)) # summary of distinct tissues types present in this study
# Metastatic       Primary Tumor   Solid Tissue Normal 
# 14                2222                 226 

query_TCGA = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor"))

GDCdownload(query = query_TCGA, directory = 'GDCdata_BRCA')

tcga_data = GDCprepare(query_TCGA)

dim(tcga_data)

colnames(colData(tcga_data))

# 
# 






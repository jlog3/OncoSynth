# 
# integrating tcga with synthea data for EMR + cancer pathway interactive visualizer
# to extract plausible biomarkers 
#
# 
# 
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

out_dir <- "/home/john/Desktop/cancer_vis"
setwd(out_dir)
project <- "TCGA-BRCA"
tissue <- "Primary Tumor"
custom_dir <- 'GDCdata_BRCA'

clinical <- get_gdc_clinical_data(project)  # can remove download parameter since it will find them if custom dir is same as last downloaded Of the 1183 files for download 1183 already exist.   All samples have been already downloaded
# only identifier for case is submitter_id  TCGA-E9-A1NF  
# head(clinical[[1]][[1]], 2)

# rnaseq_data <- get_gdc_rna_data("TCGA-BRCA", "Primary Tumor", download = TRUE, customdir = 'GDCdata_BRCA')
# 
barcode_submitter_mapping <- get_gdc_rna_data(project, tissue, download = TRUE, customdir = custom_dir)
# Merge with clinical data to add barcodes
clinical <- merge(clinical, barcode_submitter_mapping, by = "submitter_id")

# barcodes <- rna_data$barcodes
# submitter_ids <- rna_data$submitter_ids

chunk_size <- 100  # Adjust based on your memory capacity
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
first_five_sample_ids_first_chunk_result <- colnames(first_chunk_result)[1:2]
data_first_five_sample_ids_first_chunk_result <- first_chunk_result[, first_five_sample_ids_first_chunk_result]
# assay_tumor_data_subset = process_and_write_csv(data_first_five_sample_ids_first_chunk_result, "first_chunk_result_subset.csv", TRUE)
assay_tumor_data_subset = process_and_write_csv(data_first_five_sample_ids_first_chunk_result, "tcga_synthea_test/tcga_rawcount.csv", TRUE)

# this contains all of the genes as ensembl and has the correct Barcode identifier as colname

# get log fold changes or whatever normalization


# 

# done proc RNA data -- merge to demopgrahic
# for each of the gene data
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
library(dplyr)
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
write.csv(selected_data, 'tcga_synthea_test/tcga_clinical.csv', row.names = TRUE)

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




# To Do
# design new Shiny app 
# look at EMR and Shiyn exampes 

# Shiny to do 
# use global.R
# soruce funtinos from otuside server.R to clear up space





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






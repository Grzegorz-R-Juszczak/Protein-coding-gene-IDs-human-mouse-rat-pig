# part A – identification of ambiguous synonyms
mart_export_1 <- read.delim('mart_export.txt', sep=',', header = TRUE) # data import
head(mart_export_1, n=40) # showing first n rows

library(naniar) # package selection
mart_export_1x = mart_export_1 %>%
  naniar::replace_with_na(replace = list(Gene.Synonym = c(""))) # marking missing items in selected columns with “NA”. 
head(mart_export_1x, n=40) # showing first n rows

mart_export_1y = tidyr:: drop_na(mart_export_1x, Gene.Synonym) # removing rows with missing items in the column Gene.Synonym
head(mart_export_1y, n=40) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(mart_export_1y, tolower (`Gene.Synonym`)) ) # grouping data according to items in Gene.Synonym column regardless of lowercase or uppercase letters 
integrated <- purrr::map_dfr(
  .x = grouped$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        
        paste(unique(column), collapse = ", ")
      })
  } )
grouped  <- cbind(grouped, integrated)
grouped $data <- NULL
readr::write_csv(x = grouped, file = "grouped.csv") # data export
head(grouped, n=4) # showing first n rows

grouped_1 <- read.delim('grouped.csv', sep=',', header = TRUE) # data import
head(grouped_1, n=4) # showing first n rows

grouped_1$Ambiguous <- ifelse(grepl(", ", grouped_1$Gene.name),'TRUE','FALSE') # testing for the presence of “,” in items in the column Gene.name. Presence of the “,” indicates that more than one official gene symbol is assigned to the synonym
head(grouped_1, n=4) # showing first n rows

grouped_ambiguous1 = dplyr::filter(grouped_1, Ambiguous == "TRUE") # selecting gene synonyms with assigned more than one official gene symbol
head(grouped_ambiguous1, n=4) # showing first n rows

grouped_ambiguous2 = dplyr::mutate(grouped_ambiguous1, Temporary_column =  'Gene_synonym') # adding a temporary column with words “Gene_Synonym”.
head(grouped_ambiguous2, n=4) # showing first n rows

grouped_ambiguous3 = tidyr::unite(grouped_ambiguous2, 'Anchored_Gene_Synonym', c("Temporary_column","tolower.Gene.Synonym."), sep = "_", remove = FALSE, na.rm = FALSE) # creating a new column containing the words “Gene_Synonym” combined with the gene synonym
head(grouped_ambiguous3, n=4) # showing first n rows

grouped_ambiguous4 = dplyr::mutate(grouped_ambiguous3, Temporary_column = NULL) # removing the temporary column
head(grouped_ambiguous4, n=4) # showing first n rows

readr::write_csv(x = grouped_ambiguous4, file = "Ambiguous_gene_synonyms.csv") # data export

### part B – final editing including column names and additional description of the data

Ambigous_Gene_Synonyms_M1 <- read.delim('Ambiguous_gene_synonyms.csv', sep=',', header = TRUE) # data import
head(Ambigous_Gene_Synonyms_M1, n=4) # showing first n rows

Ambigous_Gene_Synonyms_M2 = dplyr::mutate(Ambigous_Gene_Synonyms_M1, Species = 'Mouse', Category = 'Ambiguous_gene_synonym') # adding new columns with data description (mouse genes)
# Ambigous_Gene_Synonyms_M2 = dplyr::mutate(Ambigous_Gene_Synonyms_M1, Species = 'Rat', Category = 'Ambiguous_gene_synonym') # adding new columns with data description (rat genes)
# Ambigous_Gene_Synonyms_M2 = dplyr::mutate(Ambigous_Gene_Synonyms_M1, Species = 'Human', Category = 'Ambiguous_gene_synonym') # adding new columns with data description (human genes)
# Ambigous_Gene_Synonyms_M2 = dplyr::mutate(Ambigous_Gene_Synonyms_M1, Species = 'Pig', Category = 'Ambiguous_gene_synonym') # adding new columns with data description (pig genes)

head(Ambigous_Gene_Synonyms_M2, n=4) # showing first n rows

Ambigous_Gene_Synonyms_M3  = dplyr::mutate(Ambigous_Gene_Synonyms_M2,  Ambiguous = NULL) # removing unnecessary column
head(Ambigous_Gene_Synonyms_M3, n=4) # showing first n rows

Ambigous_Gene_Synonyms_M4  <- Ambigous_Gene_Synonyms_M3[, c("Species", "Category", "tolower.Gene.Synonym.", "Anchored_Gene_Synonym", "Gene.Synonym", "Gene.name", "Gene.stable.ID")] # reordering columns
head(Ambigous_Gene_Synonyms_M4, n=4) # showing first n rows

colnames(Ambigous_Gene_Synonyms_M4) <- c("Species", "Category", "tolower_Gene_Synonym", "Anchored_Gene_Synonym", "Gene_Synonym_variants", "Assigned_Official_symbols", "Ensembl_stable_ID") # setting the final column names
head(Ambigous_Gene_Synonyms_M4, n=4) # showing first n rows

readr::write_csv(x = Ambigous_Gene_Synonyms_M4, file = "Ambiguous_gene_synonyms_final.csv") # data export




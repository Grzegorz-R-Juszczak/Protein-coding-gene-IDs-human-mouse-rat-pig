mart_export_1 <- read.delim('mart_export.txt', sep=',', header = TRUE) # data import
head(mart_export_1, n=4) # showing first n rows

mart_export_1A = dplyr::mutate(mart_export_1, TemporaryCol =  'Gene') # Adding temporary column with the word “gene”
head(mart_export_1A, n=4) # showing first n rows

mart_export_1B = tidyr::unite(mart_export_1A, 'Anchored_gene_name', c("TemporaryCol","Gene.name"), sep = "_", remove = FALSE, na.rm = FALSE) # adding new column with joint items from two other columns
head(mart_export_1B, n=4) # showing first n rows

mart_export_1C = dplyr::mutate(mart_export_1B,  TemporaryCol = NULL) # removing temporary column
head(mart_export_1C, n=4) # showing first n rows

library(naniar)
mart_export_2a = mart_export_1C %>%
  naniar::replace_with_na(replace = list(Gene.name = c(""))) # marking missing values with NA
head(mart_export_2a, n=15) # showing first n rows

No_gene_symbol = mart_export_2a[is.na(mart_export_2a$Gene.name), ] # selecting genes with missing gene name (rows with NA in column Gene.name)
head(No_gene_symbol, n=4) # showing first n rows
readr::write_csv(x = No_gene_symbol, file = "Genes_without_gene_symbol.csv") # data export

mart_export_2b = tidyr:: drop_na(mart_export_2a, Gene.name) # removing rows with missing values in column Gene.name
head(mart_export_2b, n=15) # showing first n rows

duplicates_from_first = duplicated(mart_export_2b$Gene.name, fromLast = FALSE) # identification of duplicated gene names (forward)
duplicates_from_last = duplicated(mart_export_2b$Gene.name, fromLast = TRUE) # identification of duplicated gene names (reverse)

mart_export_2c  <- cbind(mart_export_2b, duplicates_from_first, duplicates_from_last) # joining the dataset with results of duplication testing
head(mart_export_2c, n=4) # showing first n rows

duplicates = dplyr::filter(mart_export_2c, duplicates_from_first == "TRUE" | duplicates_from_last == "TRUE") # selection of duplicated gene symbols
head(duplicates, n=4) # showing first n rows

List_of_duplicated_gene_symbols = dplyr::arrange(duplicates, Gene.name) # sorting according to the gene name
head(List_of_duplicated_gene_symbols, n=10) # showing first n rows
readr::write_csv(x = List_of_duplicated_gene_symbols, file = "Ensembl_multiplied_gene_symbols.csv") # data export

mart_export_3 <- tidyr::nest( dplyr::group_by(mart_export_2b, `Gene.name`) ) # grouping data according to the items in Gene.name column
mart_export_4 <- purrr::map_dfr(
  .x = mart_export_3$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        
        paste(unique(column), collapse = ", ") 
      })
  } )
mart_export_3  <- cbind(mart_export_3, mart_export_4)
mart_export_3$data <- NULL
head(mart_export_3, n=4) # showing first n rows
readr::write_csv(x = mart_export_3, file = "All_unique_protein_coding.csv") # data export


All_unique_protein_coding <- read.delim('All_unique_protein_coding.csv', sep=',', header = TRUE) # data import
head(All_unique_protein_coding, n=4) # showing first n rows

All_unique_protein_coding_symbols = dplyr::select (All_unique_protein_coding, Gene.name) # selecting “Gene.name” column
readr::write_csv(x = All_unique_protein_coding_symbols, file = "All_unique_protein_coding_symbols.csv") # data export (“Gene.name” column)
head(All_unique_protein_coding_symbols, n=4) # showing first n rows



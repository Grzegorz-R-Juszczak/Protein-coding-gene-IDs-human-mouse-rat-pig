### Part A: identification of ambiguous official gene symbols

library(dplyr) # package selection
mart_export_0 <- read.delim('mart_export.txt', sep=',', header = TRUE) # data import 
head(mart_export_0, n=4) # showing first n rows

mart_export_0A = dplyr::mutate(mart_export_0, TemporaryCol =  'Gene') # adding temporary column with the word “Gene”
head(mart_export_0A, n=4) # showing first n rows

mart_export_0C = tidyr::unite(mart_export_0A, 'Anchored_gene_name', c("TemporaryCol","Gene.name"), sep = "_", remove = FALSE, na.rm = FALSE) # adding new column with joint items from two other columns
head(mart_export_0C, n=4) # showing first n rows

mart_export_0D = dplyr::mutate(mart_export_0C, TemporaryCol = NULL) # removing temporary column
head(mart_export_0D, n=4) # showing first n rows

mart_export_0D$the_same_blind_to_capitalization <- ifelse(
  test = tolower(mart_export_0D[["Gene.name"]]) == tolower(mart_export_0D[["Gene.Synonym"]]),
  yes = T,
  no = F) # testing whether items in column gene.name are the same as in column Gene.Synonym regardless of lowercase and uppercase letters
head(mart_export_0D, n=4) # showing first n rows

mart_export_2 <- dplyr::arrange(mart_export_0D, the_same_blind_to_capitalization) # data sorting according to the the_same_blind_to_capitalization column
head(mart_export_2, n=4) # showing first n rows

mart_export_3 <- dplyr::filter(mart_export_2, the_same_blind_to_capitalization == "FALSE") # selecting rows with different items in columns gene.name and Gene.Synonym
head(mart_export_3, n=4) # showing first n rows

Ambiguous_symbols_protein_other <- dplyr::filter(mart_export_3, Gene.type != "protein_coding") # selecting rows with non protein-coding genes
head(Ambiguous_symbols_protein_other, n=4) # showing first n rows
readr::write_csv(x = Ambiguous_symbols_protein_other, file = "Ambiguous_symbols_protein_other.csv") # data export

mart_export_3b <- dplyr::filter(mart_export_3, Gene.type == "protein_coding") # selecting rows with protein-coding genes
head(mart_export_3b, n=4) # showing first n rows

mart_export_4 = dplyr::mutate(mart_export_3b,  the_same_blind_to_capitalization = NULL) # removing column
head(mart_export_4, n=4) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(mart_export_4, tolower (`Gene.Synonym`)) ) # grouping data according to items in Gene.Synonym column regardless of lowercase or uppercase letters
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

readr::write_csv(x = grouped, file = "mart_export_5.csv") # data export

### Part B: Combining information about input and output data

mart_export_6 <- read.delim('mart_export_5.csv', sep=',', header = TRUE) # data import
head(mart_export_6, n=4) # showing first n rows

All_unique_protein_coding <- read.delim('All_unique_protein_coding.csv', sep=',', header = TRUE) # data import
head(All_unique_protein_coding, n=4) # showing first n rows

All_unique_protein_coding2 = dplyr::mutate(All_unique_protein_coding,  tolower_Gene_name = tolower(Gene.name)) # adding new column with items from Gene.name column written only with lowercase letters
head(All_unique_protein_coding2, n=4) # showing first n rows

mart_export_7  = dplyr:: right_join(All_unique_protein_coding2, mart_export_6, by = c("tolower_Gene_name" = "tolower.Gene.Synonym."), keep = TRUE, na_matches = "never") # joining data according to the items in columns tolower_Gene_name and tolower.Gene.Synonym (all data from the second dataset and matching data from the first dataset).
head(mart_export_7, n=4) # showing first n rows

mart_export_8  = dplyr::mutate(mart_export_7, tolower_Gene_name = NULL) # removing column
head(mart_export_8  , n=4) # showing first n rows

mart_export_9 = mart_export_8 %>% relocate(Anchored_gene_name.y, .after = Gene.name.y) # reordering columns
head(mart_export_9, n=4) # showing first n rows

mart_export_10 = mart_export_9 %>% relocate(Gene.type, .after = Gene.stable.ID.y) # reordering columns
head(mart_export_10, n=4) # showing first n rows

mart_export_11 = mart_export_10 %>% relocate(Gene.stable.ID.x, .after = Anchored_gene_name.x) # reordering columns
head(mart_export_11, n=4) # showing first n rows

readr::write_csv(x = mart_export_11, file = "Ambiguous_symbols_protein_protein.csv") # data export

### Part C: Final editing including changing column names and adding data descriptions

Ambiguous_symbols_protein_protein_M1 <- read.delim('Ambiguous_symbols_protein_protein.csv', sep=',', header = TRUE) # data import
head(Ambiguous_symbols_protein_protein_M1, n=4) # showing first n rows

Ambiguous_symbols_protein_protein_M2 = dplyr::mutate(Ambiguous_symbols_protein_protein_M1, Species = 'Mouse', Category = 'Ambiguous official symbols') # adding new columns with data description (mouse genome)
# Ambiguous_symbols_protein_protein_M2 = dplyr::mutate(Ambiguous_symbols_protein_protein_M1, Species = 'Rat', Category = 'Ambiguous official symbols') # adding new columns with data description (rat genome)
# Ambiguous_symbols_protein_protein_M2 = dplyr::mutate(Ambiguous_symbols_protein_protein_M1, Species = 'Human', Category = 'Ambiguous official symbols') # adding new columns with data description (human genome)
# Ambiguous_symbols_protein_protein_M2 = dplyr::mutate(Ambiguous_symbols_protein_protein_M1, Species = 'Pig', Category = 'Ambiguous official symbols') # adding new columns with data description (pig genome)

head(Ambiguous_symbols_protein_protein_M2, n=4) # showing first n rows

Ambiguous_symbols_protein_protein_M3  = dplyr::mutate(Ambiguous_symbols_protein_protein_M2, Gene.type = NULL) # removing unnecessary column
head(Ambiguous_symbols_protein_protein_M3, n=4) # showing first n rows

Ambiguous_symbols_protein_protein_M4  <- Ambiguous_symbols_protein_protein_M3[, c("Species", "Category", "Gene.name.x", "Anchored_gene_name.x", "Gene.stable.ID.x",
                                                                                  "tolower.Gene.Synonym.", "Gene.Synonym", "Gene.name.y", "Anchored_gene_name.y", "Gene.stable.ID.y")] # column reordering
head(Ambiguous_symbols_protein_protein_M4, n=4) # showing first n rows

colnames(Ambiguous_symbols_protein_protein_M4) <- c("Species", "Category", "Ambiguous_official_symbols", "Anchored_ambiguous_official_symbols", "Ensembl_stable_IDs",
                                                    "tolower_Gene_Synonyms", "Gene_Synonym_variants", "Alternative_official_symbols", "Anchored_Alternative_official_symbols", "Alternative_Ensembl_stable_IDs") # setting final column names
head(Ambiguous_symbols_protein_protein_M4, n=4) # showing first n rows

Ambiguous_symbols_protein_protein_M5 = tidyr:: drop_na(Ambiguous_symbols_protein_protein_M4, Ambiguous_official_symbols) # removing rows with missing items in the column “ambiguous_official_symbol”
head(Ambiguous_symbols_protein_protein_M5, n=4) # showing first n rows

readr::write_csv(x = Ambiguous_symbols_protein_protein_M5, file = "Ambiguous_symbols_protein_protein_final.csv") # data export




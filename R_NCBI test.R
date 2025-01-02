ncbi_dataset = read.delim('ncbi_dataset.tsv', sep='\t', header = TRUE) # data import (tab deleamited)
head(ncbi_dataset, n = 4) # showing first n rows

ncbi_dataset_2  = dplyr::mutate(ncbi_dataset,  Accession = NULL, Begin = NULL, End = NULL, Chromosome = NULL, Orientation = NULL, Name = NULL, Transcripts.accession = NULL, Protein.accession = NULL, Protein.length = NULL, Locus.tag = NULL) # removal of unnecessary columns
head(ncbi_dataset_2, n = 4) # showing first n rows

### missing_IDs ###

missing_IDs = read.delim('missing_IDs.csv', sep=',', header = TRUE) # data import
head(missing_IDs, n = 4) # showing first n rows

missing_IDs_2  = dplyr::mutate(missing_IDs,  HGNC.ID = NULL, any_missing = NULL) # removal of unnecessary columns

head(missing_IDs_2, n = 4) # showing first n rows

missing_IDs_3  = na_rows <- missing_IDs_2[is.na(missing_IDs_2$NCBI.gene..formerly.Entrezgene..ID), ]  # selecting rows with missing items in column NCBI.gene..formerly.Entrezgene..ID
head(missing_IDs_3, n = 4) # showing first n rows

Selected_missing_NCBI = dplyr:: left_join (missing_IDs_3, ncbi_dataset_2, by = c("Gene.name" = "Symbol"), keep = TRUE, na_matches = "never") #joining databases according to the gene symbols, Included are all data from missing_IDs_3 and common data from the ncbi_dataset_2.

head(Selected_missing_NCBI, n = 4) # showing first n rows

Selected_missing_NCBI_2 <- Selected_missing_NCBI[!duplicated(Selected_missing_NCBI), ] # removal of duplicated rows
head(Selected_missing_NCBI_2, n = 4) # showing first n rows

Selected_missing_NCBI_3  = dplyr::mutate(Selected_missing_NCBI_2, Species = 'Mouse', Category = 'IDs_missing_in_Ensembl') # adding new columns with data description (mouse genome)
#Selected_missing_NCBI_3  = dplyr::mutate(Selected_missing_NCBI_2, Species = 'Rat', Category = 'IDs_missing_in_Ensembl') # adding new columns with data description (rat genome)
#Selected_missing_NCBI_3  = dplyr::mutate(Selected_missing_NCBI_2, Species = 'Human', Category = 'IDs_missing_in_Ensembl') # adding new columns with data description (human genome)
# Selected_missing_NCBI_3  = dplyr::mutate(Selected_missing_NCBI_2, Species = 'Pig', Category = 'IDs_missing_in_Ensembl') # adding new columns with data description (pig genome)
head(Selected_missing_NCBI_3, n = 4) # showing first n rows

Selected_missing_NCBI_4  <- Selected_missing_NCBI_3[, c("Species", "Category", "Gene.name", "Gene.stable.ID", "NCBI.gene..formerly.Entrezgene..ID", "Symbol", "Gene.ID", "Gene.Type")] # column reordering
head(Selected_missing_NCBI_4, n=4) # showing first n rows

colnames(Selected_missing_NCBI_4) <- c("Species", "Category", "Gene.symbol_in_Ensembl", "Ensembl_stable_ID", "NCBI_ID_in_Ensembl", "Gene.symbol_in_NCBI", " NCBI_ID_in_NCBI", "Gene.type_in_NCBI") # setting final column names
head(Selected_missing_NCBI_4, n=4) # showing first n rows

readr::write_csv(x = Selected_missing_NCBI_4, file = "Verification_missing_NCBI_IDs.csv") # data export

### Multiplied_NCBI_IDs ###

Duplicated_NCBI_IDs_final = read.delim('Multiplied_NCBI_IDs_final.csv', sep=',', header = TRUE) # data import
head(Duplicated_NCBI_IDs_final, n = 10) # showing first n rows

Duplicated_NCBI_IDs_final_2  = dplyr::mutate(Duplicated_NCBI_IDs_final,  Duplicated_NCBI_IDs_A = NULL, Duplicated_NCBI_IDs_B = NULL) # removal of unnecessary columns
head(Duplicated_NCBI_IDs_final_2, n = 4) # showing first n rows

Selected_duplicated_NCBI_1 = dplyr:: left_join (Duplicated_NCBI_IDs_final_2, ncbi_dataset_2, by = c("NCBI.gene..formerly.Entrezgene..ID" = "Gene.ID"), keep = TRUE, na_matches = "never") #joining databases according to the NCBI IDs. Included are all data from  Duplicated_NCBI_IDs_final_2 and common data from the ncbi_dataset_2.

head(Selected_duplicated_NCBI_1, n = 10) # showing first n rows

Selected_duplicated_NCBI_2 <- Selected_duplicated_NCBI_1[!duplicated(Selected_duplicated_NCBI_1), ] # removal of duplicated rows

head(Selected_duplicated_NCBI_2, n = 10) # showing first n rows


Selected_duplicated_NCBI_2$matched.Symbols <- ifelse(Selected_duplicated_NCBI_2$Gene.name == Selected_duplicated_NCBI_2$Symbol, 'yes', 'no') # testing whether gene symbols from Ensembl (column Gene.name) are the same as in the NCBI data (column Symbol)
head(Selected_duplicated_NCBI_2, n = 10) # showing first n rows


Selected_duplicated_NCBI_3 <- tidyr::nest( dplyr::group_by(Selected_duplicated_NCBI_2, `NCBI.gene..formerly.Entrezgene..ID`) ) # grouping data according to the NCBI IDs in Ensembl data (column “NCBI.gene..formerly.Entrezgene..ID”).
mart_export_3 <- purrr::map_dfr(
  .x = Selected_duplicated_NCBI_3$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        
        paste(unique(column), collapse = ", ") 
      })
  } )
Selected_duplicated_NCBI_3  <- cbind(Selected_duplicated_NCBI_3, mart_export_3)
Selected_duplicated_NCBI_3 $data <- NULL
head(Selected_duplicated_NCBI_3, n=4) # showing first n rows
readr::write_csv(x = Selected_duplicated_NCBI_3, file = "Selected_duplicated_NCBI_3.csv") # data export

Selected_duplicated_NCBI_4 <- read.delim('Selected_duplicated_NCBI_3.csv', sep=',', header = TRUE) # data import
head(Selected_duplicated_NCBI_4, n=4) # showing first n rows

Selected_duplicated_NCBI_4$Ambigous <- ifelse(grepl(", ", Selected_duplicated_NCBI_4$Symbol),'TRUE','FALSE')# Testing for presence of “, ” in column Symbol indicating that more than one gene symbol from NCBI data was linked to the NCBI ID. 
head(Selected_duplicated_NCBI_4, n=4) # showing first n rows

Selected_duplicated_NCBI_5  = dplyr::mutate(Selected_duplicated_NCBI_4, Species = 'Mouse', Category = 'Ambiguous_NCBI_IDs_in_Ensembl') # adding new columns with data description (mouse genome)
#Selected_duplicated_NCBI_5  = dplyr::mutate(Selected_duplicated_NCBI_4, Species = 'Rat', Category = 'Ambiguous_NCBI_IDs_in_Ensembl') # adding new columns with data description (rat genome)
#Selected_duplicated_NCBI_5  = dplyr::mutate(Selected_duplicated_NCBI_4, Species = 'Human', Category = 'Ambiguous_NCBI_IDs_in_Ensembl') # adding new columns with data description (human genome)
# Selected_duplicated_NCBI_5  = dplyr::mutate(Selected_duplicated_NCBI_4, Species = 'Pig', Category = 'Ambiguous_NCBI_IDs_in_Ensembl') # adding new columns with data description (pig genome)
head(Selected_duplicated_NCBI_5, n = 4) # showing first n rows

Selected_duplicated_NCBI_6  <- Selected_duplicated_NCBI_5[, c("Species", "Category", "NCBI.gene..formerly.Entrezgene..ID", "Gene.name", "Gene.stable.ID", "Symbol", "Gene.ID", "Gene.Type", "matched.Symbols", "Ambigous")] # column reordering
head(Selected_duplicated_NCBI_6, n=4) # showing first n rows

colnames(Selected_duplicated_NCBI_6) <- c("Species", "Category",  "NCBI_ID_in_Ensembl", "Gene.symbol_in_Ensembl", "Ensembl_stable_ID", "Gene.symbol_in_NCBI", "NCBI_ID_in_NCBI", "Gene.type_in_NCBI", "Matched.symbols", "Ambigous_in_NCBI") # setting final column names
head(Selected_duplicated_NCBI_6, n=4) # showing first n rows

readr::write_csv(x = Selected_duplicated_NCBI_6, file = "Summary_NCBI_IDs_ambiguous_in_Ensembl.csv") # data export

Selected_duplicated_NCBI_7  = dplyr::filter(Selected_duplicated_NCBI_6,  Ambigous_in_NCBI == "TRUE") # selection of IDs linked with more than one gene symbol in NCBI data
head(Selected_duplicated_NCBI_7, n=4) # showing first n rows
readr::write_csv(x = Selected_duplicated_NCBI_7, file = "NCBI_IDs_ambiguous_in_NCBI.csv") # data export

#### THE END ####


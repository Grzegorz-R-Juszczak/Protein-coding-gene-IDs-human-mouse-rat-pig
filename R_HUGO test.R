hgnc <- read.delim('hgnc-search-1734202353314.txt', sep='\t', header = TRUE) # data import
head(hgnc, n=4) # showing first n rows

hgnc_2 = dplyr::mutate(hgnc,  
                       Name = NULL,
                       URL = NULL) #  removing unnecessary rows
head(hgnc_2, n=4) # showing first n rows

#### missing_IDs ####

missing_IDs = read.delim('missing_IDs.csv', sep=',', header = TRUE) # data import
head(missing_IDs, n = 4) # # showing first n rows

missing_IDs_2  = dplyr::mutate(missing_IDs,  NCBI.gene..formerly.Entrezgene..ID = NULL, any_missing = NULL) #  removing unnecessary rows
head(missing_IDs_2, n = 4) # showing first n rows

missing_IDs_3  = na_rows <- missing_IDs_2[is.na(missing_IDs_2$HGNC.ID), ]  # selecting data with missing items in column HGNC.ID
head(missing_IDs_3, n = 4) # showing first n rows

Selected_missing_HUGO = dplyr:: left_join (missing_IDs_3, hgnc_2, by = c("Gene.name" = "Symbol"), keep = TRUE, na_matches = "never") #joining databases (all data from first dataset and matching data from the second dataset)
head(Selected_missing_HUGO, n = 4) # showing first n rows
Final_editing_1  = dplyr::mutate(Selected_missing_HUGO, Species = 'Human', Category = 'IDs_missing_in_Ensembl') # adding new columns with data description (human genome)
head(Final_editing_1, n = 4) # showing first n rows

Final_editing_2 <- Final_editing_1[, c("Species", "Category", "Gene.name", "Gene.stable.ID", "HGNC.ID", "Type", "Symbol", "ID", "Locus.type", "status")] # column reordering
head(Final_editing_2, n=4) # showing first n rows
colnames(Final_editing_2) <- c("Species", "Category", "Gene.symbol_in_Ensembl", "Ensembl_stable_ID", "HGNC_ID_in_Ensembl", "Type_in_HGNC", "Gene.symbol_in_HGNC", "HGNC_ID_in_HGNC", "Locus.type_in_HGNC", "Status_in_HGNC") # setting final column names
head(Final_editing_2, n=4) # showing first n rows
readr::write_csv(x = Final_editing_2, file = "Verification_missing_HUGO_IDs.csv") # data export

#### The END ####



vgnc <- read.delim('pig_vgnc_gene_set_All.txt', sep='\t', header = FALSE) # data import
head(vgnc, n=4) # showing first n rows
colnames(vgnc) <- c("vgnc_id", "symbol", "name", "locus_group", "locus_type", "status", "location", "location_sortable", "alias_symbol", "alias_name", "prev_symbol", "prev_name", "gene_family", "gene_family_id", "date_approved_reserved", "date_symbol_changed", "date_name_changed", "date_modified", "ncbi_id", "ensembl_gene_id", "uniprot_ids", "pubmed_id", "horde_id", "error_column", "hgnc_orthologs") # setting column names including erroneous column V24
head(vgnc, n=4) # showing first n rows
vgnc2 = vgnc [-1,] # removing first row with column names imported from the file
head(vgnc2, n=4) # showing first n rows

vgnc3 = dplyr::mutate(vgnc2,  
                      name = NULL,
                      location = NULL, 
                      location_sortable = NULL, 
                      alias_symbol = NULL,
                      alias_name = NULL,
                      prev_symbol = NULL,
                      prev_name = NULL,
                      gene_family = NULL,
                      gene_family_id = NULL,
                      date_approved_reserved = NULL,
                      date_symbol_changed = NULL,
                      date_name_changed = NULL,
                      date_modified = NULL,
                      ncbi_id = NULL,
                      ensembl_gene_id = NULL,
                      uniprot_ids = NULL,
                      pubmed_id = NULL,
                      horde_id = NULL,
                      status = NULL,
                      error_column = NULL,
                      hgnc_orthologs = NULL) # removing unnecessary columns
head(vgnc3, n=4) # showing first n rows

#### missing_IDs in Ensembl ####

missing_IDs = read.delim('missing_IDs.csv', sep=',', header = TRUE) # data import
head(missing_IDs, n = 4) #  showing first n rows

missing_IDs_2  = dplyr::mutate(missing_IDs,  NCBI.gene..formerly.Entrezgene..ID = NULL, any_missing = NULL) # removing unnecessary columns
head(missing_IDs_2, n = 4) # showing first n rows

missing_IDs_3  = na_rows <- missing_IDs_2[is.na(missing_IDs_2$VGNC.ID), ]  # selecting rows with missing items in the column VGNC.ID
head(missing_IDs_3, n = 4) # showing first n rows

Selected_missing_VGNC = dplyr:: left_join (missing_IDs_3, vgnc3, by = c("Gene.name" = "symbol"), keep = TRUE, na_matches = "never") # joining databases (all data from the first dataset and matching data from the second dataset
head(Selected_missing_VGNC, n = 4) # showing first n rows

Final_editing_1  = dplyr::mutate(Selected_missing_VGNC, Species = 'Pig', Category = 'IDs_missing_in_Ensembl') # adding new columns with data description (human genome)
head(Final_editing_1, n = 4) # showing first n rows


Final_editing_2 <- Final_editing_1[, c("Species", "Category", "Gene.name", "Gene.stable.ID", "VGNC.ID", "vgnc_id", "symbol", "locus_group", "locus_type")] # column reordering
head(Final_editing_2, n=4) # showing first n rows

colnames(Final_editing_2) <- c("Species", "Category", "Gene.symbol_in_Ensembl", "Ensembl_stable_ID ", " VGNC_ID_in_Ensembl", "VGNC_ID_in_VGNC", "Gene.symbol_in_VGNC", "Locus_group_in_VGNC ", "Locus_type_in_VGNC") # setting final column names
head(Final_editing_2, n=4) # showing first n rows
readr::write_csv(x = Final_editing_2, file = "Verification_missing_VGNC_IDs.csv") # data export

#### The END ####



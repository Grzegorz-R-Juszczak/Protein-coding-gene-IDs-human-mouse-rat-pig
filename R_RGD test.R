########################################################################################################
#### Trial dataset for this R script can be downloaded from                                         ####
#### https://data.mendeley.com/datasets/454s2vw255/1/files/41150316-400d-48ad-ad12-01f89256c1d6     ####
#### The data should be extracted from the compressed folder before running the script.             #### 
########################################################################################################

GENES_RAT = read.delim('GENES_RAT.txt', sep='\t', header = TRUE) # importing data (tab deleamited)
head(GENES_RAT, n = 4) # showing first n rows

GENES_RAT_2  = dplyr::mutate(GENES_RAT, 
                             NAME = NULL, 
                             GENE_DESC = NULL, 
                             CHROMOSOME_CELERA = NULL, 
                             CHROMOSOME_mRatBN7.2 = NULL, 
                             CHROMOSOME_RGSC_v3.4 = NULL, 
                             FISH_BAND = NULL, 
                             START_POS_CELERA = NULL, 
                             STOP_POS_CELERA = NULL, 
                             STRAND_CELERA = NULL, 
                             START_POS_mRatBN7.2 = NULL,
                             STOP_POS_mRatBN7.2 = NULL,
                             STRAND_mRatBN7.2 = NULL,
                             START_POS_RGSC_v3.4 = NULL,
                             STOP_POS_RGSC_v3.4 = NULL,
                             STRAND_RGSC_v3.4 = NULL,
                             CURATED_REF_RGD_ID = NULL,
                             CURATED_REF_PUBMED_ID = NULL,
                             UNCURATED_PUBMED_ID = NULL,
                             NCBI_GENE_ID = NULL,
                             UNIPROT_ID = NULL,
                             GENBANK_NUCLEOTIDE = NULL,
                             GENBANK_PROTEIN = NULL,
                             CANONICAL_PROTEIN = NULL,
                             MARKER_RGD_ID = NULL,
                             MARKER_SYMBOL = NULL,
                             OLD_SYMBOL = NULL,
                             OLD_NAME = NULL,
                             QTL_RGD_ID = NULL,
                             QTL_SYMBOL = NULL,
                             SPLICE_RGD_ID = NULL,
                             SPLICE_SYMBOL = NULL,
                             ENSEMBL_ID = NULL,
                             X.UNUSED. = NULL,
                             CHROMOSOME_Rnor_5.0 = NULL,
                             START_POS_Rnor_5.0 = NULL,
                             STOP_POS_Rnor_5.0 = NULL,
                             STRAND_Rnor_5.0 = NULL,
                             CHROMOSOME_Rnor_6.0 = NULL,
                             START_POS_Rnor_6.0 = NULL,
                             STOP_POS_Rnor_6.0 = NULL,
                             STRAND_Rnor_6.0 = NULL,
                             CHROMOSOME_ENSEMBL = NULL,
                             START_POS_ENSEMBL = NULL,
                             STOP_POS_ENSEMBL = NULL,
                             STRAND_ENSEMBL = NULL,
                             CHROMOSOME_GRCr8 = NULL,
                             START_POS_GRCr8 = NULL,
                             STOP_POS_GRCr8 = NULL,
                             STRAND_GRCr8 = NULL,
                             CHROMOSOME_UTH_SHR = NULL,
                             START_POS_UTH_SHR = NULL,
                             STOP_POS_UTH_SHR = NULL,
                             STRAND_UTH_SHR = NULL,
                             CHROMOSOME_UTH_SHRSP = NULL,
                             START_POS_UTH_SHRSP = NULL,
                             STOP_POS_UTH_SHRSP = NULL,
                             STRAND_UTH_SHRSP = NULL,
                             CHROMOSOME_UTH_WKY = NULL,
                             START_POS_UTH_WKY = NULL,
                             STOP_POS_UTH_WKY = NULL,
                             TIGR_ID = NULL,
                             STRAND_UTH_WKY = NULL) # removing unnecessary columns
head(GENES_RAT_2, n = 4) # showing first n rows


### missing_IDs in Ensembl ###

missing_IDs = read.delim('missing_IDs.csv', sep=',', header = TRUE) # data import
head(missing_IDs, n = 4) # showing first n rows

missing_IDs_2  = dplyr::mutate(missing_IDs, NCBI.gene..formerly.Entrezgene..ID = NULL, any_missing = NULL) # removal of unnecessary columns
head(missing_IDs_2, n = 4) # showing first n rows

missing_IDs_3  = na_rows <- missing_IDs_2[is.na(missing_IDs_2$RGD.ID), ]  # selecting rows with missing items in column RGD.ID
head(missing_IDs_3, n = 4) # showing first n rows

#### Multiplied_RGD_IDs in Ensembl ####

Duplicated_Committee_IDs_final = read.delim('Multiplied_Committee_IDs_final.csv', sep=',', header = TRUE) #wczytywanie pliku
head(Duplicated_Committee_IDs_final, n = 10) # showing first n rows

Duplicated_Committee_IDs_final_2  = dplyr::mutate(Duplicated_Committee_IDs_final, Duplicated_Committee_IDs_A = NULL, Duplicated_Committee_IDs_B = NULL) # wywalanie zbednych kolumn
head(Duplicated_Committee_IDs_final_2 , n = 4) # showing first n rows

Selected_duplicated_rgd_1 = dplyr:: left_join (Duplicated_Committee_IDs_final_2, GENES_RAT_2, by = c("RGD.ID" = "GENE_RGD_ID"), keep = TRUE, na_matches = "never") #joining databases according to the RGD IDs (all data from the first dataset and matching data from the second dataset).
head(Selected_duplicated_rgd_1, n = 10) # showing first n rows

Selected_duplicated_rgd_1$matched.Symbol <- ifelse(Selected_duplicated_rgd_1$Gene.name == Selected_duplicated_rgd_1$SYMBOL, 'yes', 'no') # testing whether gene symbols from Ensembl are the same as in RGD data 
head(Selected_duplicated_rgd_1, n = 10) # showing first n rows

Selected_duplicated_rgd_2  = na_rows <- Selected_duplicated_rgd_1[is.na(Selected_duplicated_rgd_1$GENE_RGD_ID), ]  # selecting rows with missing items in column GENE_RGD_ID
head(Selected_duplicated_rgd_2, n = 10) # showing first n rows
readr::write_csv(x = Selected_duplicated_rgd_2, file = "Ensembl_RGD_IDs_missing_in_RGD.csv") # data export


Selected_duplicated_rgd_3 <- tidyr::nest( dplyr::group_by(Selected_duplicated_rgd_1, `RGD.ID`) ) # grouping data according to the RGD IDs in Ensembl data (column “RGD.ID”).
mart_export_3 <- purrr::map_dfr(
  .x = Selected_duplicated_rgd_3$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        
        paste(unique(column), collapse = ", ") 
      })
  } )
Selected_duplicated_rgd_3 <- cbind(Selected_duplicated_rgd_3, mart_export_3)
Selected_duplicated_rgd_3$data <- NULL
head(Selected_duplicated_rgd_3, n=4) # showing first n rows
readr::write_csv(x = Selected_duplicated_rgd_3, file = " Selected_duplicated_rgd_3.csv") # data export

Selected_duplicated_rgd_4 <- read.delim(' Selected_duplicated_rgd_3.csv', sep=',', header = TRUE) # data import
head(Selected_duplicated_rgd_4, n=4) # showing first n rows

Selected_duplicated_rgd_4$Ambigous <- ifelse(grepl(", ", Selected_duplicated_rgd_4$SYMBOL),'TRUE','FALSE')# Testing for presence of “, ” in column SYMBOL indicating that more than one gene symbol from RGD data was linked to the RGD ID. 
head(Selected_duplicated_rgd_4, n=4) # showing first n rows

Final_editing_1  = dplyr::mutate(Selected_duplicated_rgd_4, Species = 'Rat', Category = 'Ambiguous_RGD_IDs_in_Ensembl') # adding new columns with data description (rat genome)
head(Final_editing_1, n = 4) # showing first n rows

Final_editing_2 <- Final_editing_1[, c("Species", "Category", "RGD.ID", "Gene.name", "Gene.stable.ID", "GENE_RGD_ID", "SYMBOL", "GENE_REFSEQ_STATUS", "NOMENCLATURE_STATUS", "GENE_TYPE","matched.Symbol", "Ambigous")] # column reordering
head(Final_editing_2, n=4) # showing first n rows


colnames(Final_editing_2) <- c("Species", "Category", "RGD_ID_in_Ensembl", "Gene.symbol_in_Ensembl", "Ensembl_stable_ID", "RGD_ID_in_RGD", "Gene.symbol_in_RGD", "GENE_REFSEQ_STATUS_in_RGD", "NOMENCLATURE_STATUS_in_RGD", "GENE_TYPE_in_RGD","Matched.symbol", "Ambigous_in_RGD") # setting final column names
head(Final_editing_2, n=4) # showing first n rows

Final_editing_3  = dplyr::mutate(Final_editing_2, RGD_ID_in_RGD = NULL, GENE_REFSEQ_STATUS_in_RGD = NULL, NOMENCLATURE_STATUS_in_RGD = NULL) # removal of unnecessary columns
head(Final_editing_3, n = 4) # showing first n rows

readr::write_csv(x = Final_editing_3, file = "Summary_RGD_IDs_ambiguous_in_Ensembl.csv") # data export

Final_editing_4  = dplyr::filter(Final_editing_3, Ambigous_in_RGD == "TRUE") # selection of IDs linked to more than one gene symbol in RGD data
head(Final_editing_4, n=4) # showing first n rows

readr::write_csv(x = Final_editing_4, file = "RGD_IDs_ambiguous_in_RGD.csv") # data export

##### THE END #####

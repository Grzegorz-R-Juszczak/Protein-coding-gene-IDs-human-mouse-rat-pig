########################################################################################################
#### Trial dataset for this R script can be downloaded from                                         ####
#### https://data.mendeley.com/datasets/454s2vw255/1/files/79a080ef-3dfb-45e4-88d0-c9f1274d9b9c     ####
#### The data should be extracted from the compressed folder before running the script.             #### 
########################################################################################################

mart_export_1 <- read.delim('mart_export.txt', sep=',', header = TRUE) # data importing
head(mart_export_1, n=20) # showing first n rows

library(naniar) # package selection
mart_export_1b = naniar::replace_with_na_all(data = mart_export_1,
                                             condition = ~.x == "") # marking missing values with “NA” in all columns
head(mart_export_1b, n=20) # showing first n rows

mart_export_2 <- tidyr::nest( dplyr::group_by(mart_export_1b, `Gene.name`) ) # grouping data according to the items in the column “Gene.name”.
mart_export_3 <- purrr::map_dfr(
  .x = mart_export_2$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        
        paste(unique(column), collapse = ", ") 
      })
  } )
mart_export_2  <- cbind(mart_export_2, mart_export_3)
mart_export_2$data <- NULL
head(mart_export_2, n=4) # showing first n rows
readr::write_csv(x = mart_export_2, file = "mart_export_2.csv") # data export

mart_export_4 <- read.delim('mart_export_2.csv', sep=',', header = TRUE) # data importing
head(mart_export_4, n=4) # showing first n rows

mart_export_6 = mart_export_4  %>% add_label_missings() # Adding a column describing if there are any missing items in each row
head(mart_export_6, n=20) # showing first n rows

mart_export_7  = dplyr::filter(mart_export_6,  any_missing == "Missing") # selection of rows with at least one missing item
head(mart_export_7, n=20) # showing first n rows
readr::write_csv(x = mart_export_7, file = "Missing_IDs.csv") # exporting to text file

###### Multiplied Ensembl IDs ######

Ensembl_IDs = dplyr::select (mart_export_6, c(Gene.name, Gene.stable.ID)) # selecting columns with gene symbols and Ensembl IDs (Gene.stable.ID)
head(Ensembl_IDs, n=20) # showing first n rows

library(tidyr) # package selection
Ensembl_IDs_2  = Ensembl_IDs %>% separate_longer_delim(Gene.stable.ID, delim = ', ') # Splitting Ensembl stable IDs into rows
head(Ensembl_IDs_2, n=20) # showing first n rows

Duplicated_Ensembl_IDs_A = duplicated(Ensembl_IDs_2$Gene.stable.ID, incomparables = c(NA, 'NA'), fromLast = FALSE) # Identification of duplicated Ensembl stable IDs (forward)
head(Duplicated_Ensembl_IDs_A, n=20) # showing first n rows
Duplicated_Ensembl_IDs_B = duplicated(Ensembl_IDs_2$Gene.stable.ID, incomparables = c(NA, 'NA'), fromLast = TRUE) # Identification of duplicated Ensembl stable IDs (reverse)
head(Duplicated_Ensembl_IDs_B, n=20) # showing first n rows

Duplicated_Ensembl_IDs_2  <- cbind(Ensembl_IDs_2, Duplicated_Ensembl_IDs_A, Duplicated_Ensembl_IDs_B) # assembling the main dataset with identifiers of duplicated items
head(Duplicated_Ensembl_IDs_2, n=4) # showing first n rows

Duplicated_Ensembl_IDs_3  = dplyr::filter(Duplicated_Ensembl_IDs_2, Duplicated_Ensembl_IDs_A == "TRUE" | Duplicated_Ensembl_IDs_B == "TRUE") # selection of duplicated Ensembl IDs
head(Duplicated_Ensembl_IDs_3, n=10) # showing first n rows

Duplicated_Ensembl_IDs_4 <- dplyr::arrange(Duplicated_Ensembl_IDs_3, Gene.stable.ID) # data sorting
head(Duplicated_Ensembl_IDs_4, n=10) # showing first n rows

readr::write_csv(x = Duplicated_Ensembl_IDs_4, file = "Multiplied_Ensembl_IDs_final.csv") # exporting to text file
###### Multiplied NCBI IDs #######

NCBI_IDs = dplyr::select (mart_export_6, c(Gene.name, Gene.stable.ID, NCBI.gene..formerly.Entrezgene..ID)) # selecting columns with gene symbols, Ensembl IDs and NCBI IDs
head(NCBI_IDs, n=20) # showing first n rows

NCBI_IDs_2  = NCBI_IDs %>% separate_longer_delim(NCBI.gene..formerly.Entrezgene..ID, delim = ', ') # Splitting NCBI IDs into rows
head(NCBI_IDs_2, n=20) # showing first n rows
Duplicated_NCBI_IDs_A = duplicated(NCBI_IDs_2$NCBI.gene..formerly.Entrezgene..ID, incomparables = c(NA, 'NA'), fromLast = FALSE) # Identification of duplicated NCBI stable IDs (forward)
head(Duplicated_NCBI_IDs_A, n=20) # showing first n rows
Duplicated_NCBI_IDs_B = duplicated(NCBI_IDs_2$NCBI.gene..formerly.Entrezgene..ID, incomparables = c(NA, 'NA'), fromLast = TRUE) # Identification of duplicated NCBI stable IDs (reverse)
head(Duplicated_NCBI_IDs_B, n=20) # showing first n rows

Duplicated_NCBI_IDs_2  <- cbind(NCBI_IDs_2, Duplicated_NCBI_IDs_A, Duplicated_NCBI_IDs_B) # assembling the main dataset with identifiers of duplicated items
head(Duplicated_NCBI_IDs_2, n=20) # showing first n rows

Duplicated_NCBI_IDs_3  = dplyr::filter(Duplicated_NCBI_IDs_2, Duplicated_NCBI_IDs_A == "TRUE" | Duplicated_NCBI_IDs_B == "TRUE") # selection of duplicated NCBI IDs
head(Duplicated_NCBI_IDs_3, n=10) # showing first n rows

Duplicated_NCBI_IDs_4 <- dplyr::arrange(Duplicated_NCBI_IDs_3, NCBI.gene..formerly.Entrezgene..ID) # data sorting
head(Duplicated_NCBI_IDs_4, n=10) # showing first n rows

readr::write_csv(x = Duplicated_NCBI_IDs_4, file = "Multiplied_NCBI_IDs_final.csv") # exporting to text file

###### Multiplied committee IDs #######

Committee_IDs = dplyr::select (mart_export_6, c(Gene.name, Gene.stable.ID, HGNC.ID)) # selecting columns with gene symbols, Ensembl IDs and HGNC IDs. 
head(Committee_IDs, n=20) # showing first n rows

Committee_IDs_2  = Committee_IDs %>% separate_longer_delim(HGNC.ID, delim = ', ') # Splitting HGNC IDs into rows
head(Committee_IDs_2, n=20) # showing first n rows
Duplicated_Committee_IDs_A = duplicated(Committee_IDs_2$HGNC.ID, incomparables = c(NA, 'NA'), fromLast = FALSE) # Identification of duplicated HGNC IDs (forward)
head(Duplicated_Committee_IDs_A, n=20) # showing first n rows
Duplicated_Committee_IDs_B = duplicated(Committee_IDs_2$HGNC.ID, incomparables = c(NA, 'NA'), fromLast = TRUE) # Identification of duplicated HGNC IDs (reverse)
head(Duplicated_Committee_IDs_B, n=20) # showing first n rows

Duplicated_Committee_IDs_2  <- cbind(Committee_IDs_2, Duplicated_Committee_IDs_A, Duplicated_Committee_IDs_B) # assembling the main dataset with identifiers of duplicated items
head(Duplicated_Committee_IDs_2, n=20) # showing first n rows

Duplicated_Committee_IDs_3  = dplyr::filter(Duplicated_Committee_IDs_2, Duplicated_Committee_IDs_A == "TRUE" | Duplicated_Committee_IDs_B == "TRUE") # selection of duplicated committee IDs
head(Duplicated_Committee_IDs_3, n=10) # showing first n rows

Duplicated_Committee_IDs_4 <- dplyr::arrange(Duplicated_Committee_IDs_3, HGNC.ID) # data sorting
head(Duplicated_Committee_IDs_4, n=10) # showing first n rows

readr::write_csv(x = Duplicated_Committee_IDs_4, file = "Multiplied_Committee_IDs_final.csv") # exporting to text file

#### THE END ####

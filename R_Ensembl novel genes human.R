########################################################################################################
#### Trial dataset for this R script can be downloaded from                                         ####
#### https://data.mendeley.com/datasets/454s2vw255/1/files/1dbebafa-3a45-4982-9682-cbdf2477bc60     ####
#### The data should be extracted from the compressed folder before running the script.             #### 
########################################################################################################

#### Ensembl data import ####

Ensembl <- read.delim('mart_export.txt', sep=',', header = TRUE) # data importing
head(mart_export_1, n=5) # showing first n rows

#### HGNC data import and initial processing ####

HGNC = read.delim('results.txt', sep='\t', header = TRUE) # data import (tab deleamited)
head(HGNC, n = 4) # wyświetlanie 4 rzędów

HGNC$Ensembl_ID_as_Gene.symbol <- ifelse(grepl("ENSG00", HGNC$Approved.symbol), 'yes', 'no') # Testing the presence of Ensembl IDs in the column with gene symbols
head(HGNC, n=4) # showing first n rows

Ensembl_ID_testHGNC  = dplyr::filter(HGNC, Ensembl_ID_as_Gene.symbol == "yes") # selection of rows with Ensembl IDs used as a gene symbols. 
head(Ensembl_ID_testHGNC, n=4) # showing first n rows

HGNC  = dplyr::filter(HGNC, Ensembl_ID_as_Gene.symbol == "no") # removal of rows with Ensembl IDs used as a gene symbols. 
head(HGNC, n=4) # showing first n rows

HGNC = dplyr::mutate(HGNC,  
                     Ensembl_ID_as_Gene.symbol = NULL) # removing unnecessary columns
head(HGNC, n = 4) # showing first n rows

#### Joining datasets step 1 ####
Ensembl_HGNC_V1 = dplyr:: left_join (Ensembl, HGNC, by = c("HGNC.ID" = "HGNC.ID"), keep = TRUE, na_matches = "never") # joining databases (all data from the first dataset and matching data from the second dataset
head(Ensembl_HGNC_V1, n = 4) # showing first n rows
Ensembl_HGNC_V1_2 = tidyr:: drop_na(Ensembl_HGNC_V1, HGNC.ID.y) # removal of rows with missing HGNC IDs in HGNC data (column “HGNC.ID.y”)
head(Ensembl_HGNC_V1_2, n = 4) # showing first n rows

#### Joining datasets step 2 ####

Ensembl_HGNC_V2  = na_rows <- Ensembl_HGNC_V1[is.na(Ensembl_HGNC_V1$HGNC.ID.y), ]  # selection of rows with missing HGNC IDs in HGNC data (column “HGNC.ID.y”)
head(Ensembl_HGNC_V2, n = 4) # showing first n rows

Ensembl_HGNC_V2_2  = dplyr::mutate(Ensembl_HGNC_V2,  
                                   HGNC.ID.y = NULL,
                                   Status = NULL, 
                                   Approved.symbol = NULL, 
                                   Approved.name = NULL,
                                   Ensembl.gene.ID = NULL) # removing columns with missing data from HGNC
head(Ensembl_HGNC_V2_2, n = 4) # showing first n rows

Ensembl_HGNC_V2_3 = dplyr:: left_join (Ensembl_HGNC_V2_2, HGNC, by = c("Gene.stable.ID" = "Ensembl.gene.ID"), keep = TRUE, na_matches = "never") # joining databases (all data from the first dataset and matching data from the second dataset
head(Ensembl_HGNC_V2_3, n = 4) # showing first n rows

Ensembl_HGNC_V2_4 = dplyr:: rename(Ensembl_HGNC_V2_3, HGNC.ID.y = HGNC.ID) # changing column name to enable proper binding of data from step 1 and 2
head(Ensembl_HGNC_V2_4, n = 4) # showing first n rows


#### Joining data from step 1 and 2 ####

Ensembl_HGNC = dplyr:: bind_rows(Ensembl_HGNC_V1_2, Ensembl_HGNC_V2_4) # row binding 
head(Ensembl_HGNC, n = 4) # showing first n rows

readr::write_csv(x = Ensembl_HGNC, file = "Ensembl_HGNC.csv") # data export


#### NCBI data import and editing ####

NCBI = read.delim('ncbi_dataset.tsv', sep='\t', header = TRUE) # data import (tab deleamited)
head(NCBI, n = 4) # showing first n rows

NCBI_2  = dplyr::mutate(NCBI,  
                        Taxonomic.Name = NULL, 
                        Common.Name = NULL, 
                        Transcripts = NULL,
                        Gene.Group.Identifier = NULL,
                        Gene.Group.Method = NULL, 
                        Chromosomes = NULL, 
                        Nomenclature.ID = NULL,
                        Annotation.Genomic.Range.Accession = NULL,
                        Annotation.Genomic.Range.Start = NULL, 
                        Annotation.Genomic.Range.Stop = NULL, 
                        OMIM.IDs = NULL,
                        Orientation = NULL,
                        Proteins = NULL, 
                        Synonyms = NULL, 
                        Taxonomic.ID = NULL,
                        SwissProt.Accessions = NULL) # removing unnecessary columns
head(NCBI_2, n = 4) # showing first n rows

library(naniar)
NCBI_2  = NCBI_2  %>%
  naniar::replace_with_na(replace = list(Ensembl.GeneIDs = c(""))) # marking missing items in selected column with „NA”.

head(NCBI_2, n = 4) # showing first n rows

NCBI_2$Ensembl_ID_as_Gene.symbol <- ifelse(grepl("ENSG00", NCBI_2$Symbol), 'yes', 'no') # Testing the presence of Ensembl IDs in the column with gene symbols
head(NCBI_2, n=4) # showing first n rows

Ensembl_ID_testNCBI  = dplyr::filter(NCBI_2, Ensembl_ID_as_Gene.symbol == "yes") # selection of rows with Ensembl IDs used as a gene symbols. 
head(Ensembl_ID_testNCBI, n=4) # showing first n rows

NCBI_2  = dplyr::filter(NCBI_2, Ensembl_ID_as_Gene.symbol == "no") # removal of rows with Ensembl IDs used as a gene symbols. 
head(NCBI_2, n=4) # showing first n rows

NCBI_2 = dplyr::mutate(NCBI_2,  
                       Ensembl_ID_as_Gene.symbol = NULL) # removing unnecessary columns
head(NCBI_2, n = 4) # showing first n rows

#### Joining datasets step 1 ####

Ensembl_HGNC_NCBI_V1 = dplyr:: left_join (Ensembl_HGNC, NCBI_2, by = c("NCBI.gene..formerly.Entrezgene..ID" = "NCBI.GeneID"), keep = TRUE, na_matches = "never") # joining databases (all data from the first dataset and matching data from the second dataset
head(Ensembl_HGNC_NCBI_V1, n = 4) # showing first n rows

Ensembl_HGNC_NCBI_V1_2 = tidyr:: drop_na(Ensembl_HGNC_NCBI_V1, NCBI.GeneID) # removal of rows with missing NCBI IDs in NCBI data (column “NCBI.GeneID”)
head(Ensembl_HGNC_NCBI_V1_2, n = 4) # showing first n rows

#### Joining datasets step 2 ####

Ensembl_HGNC_NCBI_V2  = na_rows <- Ensembl_HGNC_NCBI_V1[is.na(Ensembl_HGNC_NCBI_V1$NCBI.GeneID), ]  # selection of rows with missing NCBI IDs in NCBI data (column “NCBI.GeneID”)
head(Ensembl_HGNC_NCBI_V2, n = 4) # showing first n rows

Ensembl_HGNC_NCBI_V2_2  = dplyr::mutate(Ensembl_HGNC_NCBI_V2,  
                                        NCBI.GeneID = NULL,
                                        Symbol = NULL,
                                        Description = NULL,
                                        Gene.Type = NULL, 
                                        Ensembl.GeneIDs = NULL) # removing columns with missing data from NCBI
head(Ensembl_HGNC_NCBI_V2_2, n = 4) # showing first n rows

Ensembl_HGNC_NCBI_V2_3 = dplyr:: left_join (Ensembl_HGNC_NCBI_V2_2, NCBI_2, by = c("Gene.stable.ID" = "Ensembl.GeneIDs"), keep = TRUE, na_matches = "never") # joining databases (all data from the first dataset and matching data from the second dataset
head(Ensembl_HGNC_NCBI_V2_3, n = 4) # showing first n rows

#### Joining data from step 1 and 2 ####

Ensembl_HGNC_NCBI = dplyr:: bind_rows(Ensembl_HGNC_NCBI_V1_2, Ensembl_HGNC_NCBI_V2_3) # row binding 
head(Ensembl_HGNC_NCBI, n = 4) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(Ensembl_HGNC_NCBI, `Gene.stable.ID`) ) # grouping data according to the items in Gene.stable.ID column (Ensembl IDs from Ensembl)
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
grouped$data <- NULL

head(grouped, n = 4) # showing first n rows
readr::write_csv(x = grouped, file = "Ensembl_HGNC_NCBI_grouped.csv") # data export

Ensembl_HGNC_NCBI_grouped = read.delim('Ensembl_HGNC_NCBI_grouped.csv', sep=',', header = TRUE) # data import
head(Ensembl_HGNC_NCBI_grouped, n = 4) # showing first n rows

colnames(Ensembl_HGNC_NCBI_grouped) <- c("Ensembl_ID_in_Ensembl", "Gene.symbol_in_Ensembl", 
                                         "HGNC_ID_in_Ensembl", 
                                         "NCBI_ID_in_Ensembl",
                                         "HGNC_ID_in_HGNC",
                                         "Status_in_HGNC",
                                         "Gene.symbol_in_HGNC",
                                         "Approved.name_in_HGNC",
                                         "Ensembl_ID_in_HGNC",
                                         "NCBI_ID_in_NCBI", 
                                         "Gene.symbol_in_NCBI", 
                                         "Description_in_NCBI",
                                         "Gene.Type_in_NCBI", 
                                         "Ensembl_ID_in_NCBI") # setting final column names
head(Ensembl_HGNC_NCBI_grouped, n=4) # showing first n rows

Ensembl_HGNC_NCBI_grouped$AmbiguousHGNC<- ifelse(grepl(", ", Ensembl_HGNC_NCBI_grouped$Gene.symbol_in_HGNC),'TRUE','FALSE') # testing for the presence of “,” in items in the column Gene.symbol_in_HGNC. Presence of the “,” indicates that more than one gene symbol in HGNC is assigned to Ensembl novel gene
head(Ensembl_HGNC_NCBI_grouped, n=4) # showing first n rows

Ensembl_HGNC_NCBI_grouped$AmbiguousNCBI<- ifelse(grepl(", ", Ensembl_HGNC_NCBI_grouped$Gene.symbol_in_NCBI),'TRUE','FALSE') # testing for the presence of “,” in items in the column Gene.symbol_in_NCBI. Presence of the “,” indicates that more than one gene symbol in NCBI is assigned to Ensembl novel gene
head(Ensembl_HGNC_NCBI_grouped, n=4) # showing first n rows

readr::write_csv(x = Ensembl_HGNC_NCBI_grouped, file = "Ensembl_HGNC_NCBI_grouped.csv") # data export

#### The End ####


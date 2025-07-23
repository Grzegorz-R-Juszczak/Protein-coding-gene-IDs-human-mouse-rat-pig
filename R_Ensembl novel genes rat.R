########################################################################################################
#### Trial dataset for this R script can be downloaded from                                         ####
#### https://data.mendeley.com/datasets/454s2vw255/1/files/84daa43d-5ebc-463d-b006-f2ed347374d4     ####
#### The data should be extracted from the compressed folder before running the script.             #### 
########################################################################################################

#### Ensembl data import ####

mart_export_1 <- read.delim('mart_export.txt', sep=',', header = TRUE) # data importing
head(mart_export_1, n=5) # showing first n rows

#### RGD data import and initial processing ####

RGD = read.delim('GENES_RAT.txt', sep='\t', header = TRUE) # importing data (tab deleamited)
head(RGD, n = 4) # showing first n rows

RGD_2  = dplyr::mutate(RGD, 
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
                       STRAND_UTH_WKY = NULL,
                       GENE_REFSEQ_STATUS = NULL,
                       NOMENCLATURE_STATUS = NULL) # removing unnecessary columns
head(RGD_2, n = 4) # showing first n rows

library(tidyr) # package selection
RGD_3 = RGD_2 %>% separate_longer_delim(ENSEMBL_ID, delim = ';') # Splitting Ensembl symbols into separate rows
head(RGD_3, n=4) # showing first n rows

RGD_3$Ensembl_ID_as_Gene.symbol <- ifelse(grepl("ENSRNOG00", RGD_3$SYMBOL), 'yes', 'no') # Testing the presence of Ensembl IDs in the column with gene symbols
head(RGD_3, n=4) # showing first n rows

Ensembl_ID_testRGD  = dplyr::filter(RGD_3, Ensembl_ID_as_Gene.symbol == "yes") # selection of rows with Ensembl IDs used as a gene symbols. 

head(Ensembl_ID_testRGD, n=4) # showing first n rows


RGD_3  = dplyr::filter(RGD_3, Ensembl_ID_as_Gene.symbol == "no") # removal of rows with Ensembl IDs used as a gene symbols. 
head(RGD_3, n=4) # showing first n rows
RGD_3 = dplyr::mutate(RGD_3,  
                      Ensembl_ID_as_Gene.symbol = NULL) # removing unnecessary columns
head(RGD_3, n = 4) # showing first n rows

#### Joining datasets step 1 ####
Ensembl_RGD_V1 = dplyr:: left_join (mart_export_1, RGD_3, by = c("RGD.ID" = "GENE_RGD_ID"), keep = TRUE, na_matches = "never") # joining databases (all data from the first dataset and matching data from the second dataset
head(Ensembl_RGD_V1, n = 4) # showing first n rows

Ensembl_RGD_V1_2 = tidyr:: drop_na(Ensembl_RGD_V1, GENE_RGD_ID) # removal of rows with missing RGD IDs in RGD data (column “GENE_RGD_ID”)
head(Ensembl_RGD_V1_2, n = 4) # showing first n rows

#### Joining datasets step 2 ####

Ensembl_RGD_V2  = na_rows <- Ensembl_RGD_V1[is.na(Ensembl_RGD_V1$GENE_RGD_ID), ]  # selection of rows with missing RGD IDs in RGD data (column “GENE_RGD_ID”)
head(Ensembl_RGD_V2, n = 4) # showing first n rows

Ensembl_RGD_V2_2  = dplyr::mutate(Ensembl_RGD_V2,  
                                  GENE_RGD_ID = NULL,
                                  SYMBOL = NULL,
                                  GENE_TYPE= NULL,
                                  ENSEMBL_ID = NULL) # removing columns with missing data from RGD
head(Ensembl_RGD_V2_2, n = 4) # showing first n rows

Ensembl_RGD_V2_3 = dplyr:: left_join (Ensembl_RGD_V2_2, RGD_3, by = c("Gene.stable.ID" = "ENSEMBL_ID"), keep = TRUE, na_matches = "never") # joining databases (all data from the first dataset and matching data from the second dataset
head(Ensembl_RGD_V2_3, n = 4) # showing first n rows

#### Joining data from step 1 and 2 ####

Ensembl_RGD = dplyr:: bind_rows(Ensembl_RGD_V1_2, Ensembl_RGD_V2_3) # row binding 
head(Ensembl_RGD, n = 4) # showing first n rows

readr::write_csv(x = Ensembl_RGD, file = "Ensembl_RGD.csv") # data export

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


NCBI_2$Ensembl_ID_as_Gene.symbol <- ifelse(grepl("ENSRNOG00", NCBI_2$Symbol), 'yes', 'no') # Testing the presence of Ensembl IDs in the column with gene symbols
head(NCBI_2, n=4) # showing first n rows

Ensembl_ID_testNCBI  = dplyr::filter(NCBI_2, Ensembl_ID_as_Gene.symbol == "yes") # selection of rows with Ensembl IDs used as a gene symbols. 
head(Ensembl_ID_testNCBI, n=4) # showing first n rows

NCBI_2  = dplyr::filter(NCBI_2, Ensembl_ID_as_Gene.symbol == "no") # removal of rows with Ensembl IDs used as a gene symbols. 
head(NCBI_2, n=4) # showing first n rows

NCBI_2 = dplyr::mutate(NCBI_2,  
                       Ensembl_ID_as_Gene.symbol = NULL) # removing unnecessary columns
head(NCBI_2, n = 4) # showing first n rows


#### Joining datasets step 1 ####

Ensembl_RGD_NCBI_V1 = dplyr:: left_join (Ensembl_RGD, NCBI_2, by = c("NCBI.gene..formerly.Entrezgene..ID" = "NCBI.GeneID"), keep = TRUE, na_matches = "never") # joining databases (all data from the first dataset and matching data from the second dataset
head(Ensembl_RGD_NCBI_V1, n = 4) # showing first n rows

Ensembl_RGD_NCBI_V1_2 = tidyr:: drop_na(Ensembl_RGD_NCBI_V1, NCBI.GeneID) # removal of rows with missing NCBI IDs in NCBI data (column “NCBI.GeneID”)
head(Ensembl_RGD_NCBI_V1_2, n = 4) # showing first n rows

#### Joining datasets step 2 ####

Ensembl_RGD_NCBI_V2  = na_rows <- Ensembl_RGD_NCBI_V1[is.na(Ensembl_RGD_NCBI_V1$NCBI.GeneID), ]  # selection of rows with missing NCBI IDs in NCBI data (column “NCBI.GeneID”)
head(Ensembl_RGD_NCBI_V2, n = 4) # showing first n rows

Ensembl_RGD_NCBI_V2_2  = dplyr::mutate(Ensembl_RGD_NCBI_V2,  
                                       NCBI.GeneID = NULL,
                                       Symbol = NULL,
                                       Description = NULL,
                                       Gene.Type = NULL, 
                                       Ensembl.GeneIDs = NULL) # removing columns with missing data from NCBI
head(Ensembl_RGD_NCBI_V2_2, n = 4) # showing first n rows

Ensembl_RGD_NCBI_V2_3 = dplyr:: left_join (Ensembl_RGD_NCBI_V2_2, NCBI_2, by = c("Gene.stable.ID" = "Ensembl.GeneIDs"), keep = TRUE, na_matches = "never") # joining databases (all data from the first dataset and matching data from the second dataset
head(Ensembl_RGD_NCBI_V2_3, n = 4) # showing first n rows

#### Joining data from step 1 and 2 ####

Ensembl_RGD_NCBI = dplyr:: bind_rows(Ensembl_RGD_NCBI_V1_2, Ensembl_RGD_NCBI_V2_3) # row binding 
head(Ensembl_RGD_NCBI, n = 4) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(Ensembl_RGD_NCBI, `Gene.stable.ID`) ) # grouping data according to the items in Gene.stable.ID column (Ensembl IDs from Ensembl)
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
readr::write_csv(x = grouped, file = "Ensembl_RGD_NCBI_grouped.csv") # data export

Ensembl_RGD_NCBI_grouped = read.delim('Ensembl_RGD_NCBI_grouped.csv', sep=',', header = TRUE) # data import
head(Ensembl_RGD_NCBI_grouped, n = 4) # showing n first rows

colnames(Ensembl_RGD_NCBI_grouped) <- c("Ensembl_ID_in_Ensembl", "Gene.symbol_in_Ensembl", 
                                        "RGD_ID_in_Ensembl", 
                                        "NCBI_ID_in_Ensembl",
                                        "RGD_ID_in_RGD",
                                        "Gene.symbol_in_RGD",
                                        "Gene.Type_in_RGD",
                                        "Ensembl_ID_in_RGD",
                                        "NCBI_ID_in_NCBI", 
                                        "Gene.symbol_in_NCBI", 
                                        "Description_in_NCBI",
                                        "Gene.Type_in_NCBI", 
                                        "Ensembl_ID_in_NCBI") # setting final column names
head(Ensembl_RGD_NCBI_grouped, n=4) # showing first n rows

Ensembl_RGD_NCBI_grouped$AmbiguousRGD<- ifelse(grepl(", ", Ensembl_RGD_NCBI_grouped$Gene.symbol_in_RGD),'TRUE','FALSE') # testing for the presence of “,” in items in the column Gene.symbol_in_RGD. Presence of the “,” indicates that more than one gene symbol in RGD is assigned to Ensembl novel gene
head(Ensembl_RGD_NCBI_grouped, n=4) # showing first n rows

Ensembl_RGD_NCBI_grouped$AmbiguousNCBI<- ifelse(grepl(", ", Ensembl_RGD_NCBI_grouped$Gene.symbol_in_NCBI),'TRUE','FALSE') # testing for the presence of “,” in items in the column Gene.symbol_in_NCBI. Presence of the “,” indicates that more than one gene symbol in NCBI is assigned to Ensembl novel gene
head(Ensembl_RGD_NCBI_grouped, n=4) # showing first n rows

readr::write_csv(x = Ensembl_RGD_NCBI_grouped, file = "Ensembl_RGD_NCBI_grouped.csv") # data export

#### The End ####




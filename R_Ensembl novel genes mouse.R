#### Ensembl data import ####

mart_export_1 <- read.delim('mart_export.txt', sep=',', header = TRUE) # data importing
head(mart_export_1, n=5) # showing first n rows

#### MGI data import and initial processing ####

MRK_List2 <- read.delim('MRK_List2.rpt', sep='\t', header = TRUE) # data import (the gz file has to be decompressed before importing)
head(MRK_List2, n=4) # showing first n rows

MRK_List2_2 = dplyr::mutate(MRK_List2,  
                            Chr = NULL,
                            cM.Position = NULL, 
                            genome.coordinate.start = NULL, 
                            genome.coordinate.end = NULL, 
                            strand = NULL,
                            Marker.Synonyms..pipe.separated. = NULL,
                            Status = NULL,
                            Marker.Name = NULL,
                            Marker.Type = NULL) # removing unnecessary columns
head(MRK_List2_2, n=4) # showing first n rows

MRK_ENSEMBL_1 <- read.delim('MRK_ENSEMBL.rpt', sep='\t', header = FALSE) # data import
head(MRK_ENSEMBL_1, n=4) # showing first n rows

MRK_ENSEMBL_2 = dplyr::mutate(MRK_ENSEMBL_1,  
                              V2 = NULL,
                              V3 = NULL, 
                              V4 = NULL, 
                              V5 = NULL, 
                              V7 = NULL,
                              V8 = NULL,
                              V9 = NULL,
                              V10 = NULL,
                              V11 = NULL,
                              V12 = NULL,
                              V13 = NULL) # removing unnecessary columns
head(MRK_ENSEMBL_2, n=4) # showing first n rows

MGI = dplyr:: full_join(MRK_List2_2, MRK_ENSEMBL_2, by = c("MGI.Accession.ID" = "V1"), keep = TRUE, na_matches = "never") # combining datasets
head(MGI, n = 4) # showing first n rows


MGI_2 = dplyr::mutate(MGI,  
                      V1 = NULL) # removing unnecessary columns
head(MGI_2, n=4) # showing first n rows

MGI_3 = dplyr:: rename(MGI_2, Ensembl_ID_in_MGI = V6) # renaming column
head(MGI_3, n=4) # showing first n rows

MGI_3$Ensembl_ID_as_Gene.symbol <- ifelse(grepl("ENSMUSG00", MGI_3$Marker.Symbol), 'yes', 'no') # Testing the presence of Ensembl IDs in the column with gene symbols
head(MGI_3, n=4) # showing first n rows

Ensembl_ID_testMGI  = dplyr::filter(MGI_3, Ensembl_ID_as_Gene.symbol == "yes") # selection of rows with Ensembl IDs used as a gene symbols. 
head(Ensembl_ID_testMGI, n=4) # showing first n rows


MGI_3  = dplyr::filter(MGI_3, Ensembl_ID_as_Gene.symbol == "no") # removal of rows with Ensembl IDs used as a gene symbols. 
head(MGI_3, n=4) # showing first n rows

MGI_3 = dplyr::mutate(MGI_3,  
                      Ensembl_ID_as_Gene.symbol = NULL) # removing unnecessary columns
head(MGI_3, n = 4) # showing first n rows


#### Joining datasets step 1 ####
Ensembl_MGI_V1 = dplyr:: left_join (mart_export_1, MGI_3, by = c("MGI.ID" = "MGI.Accession.ID"), keep = TRUE, na_matches = "never") # joining databases (all data from the first dataset and matching data from the second dataset
head(Ensembl_MGI_V1, n = 4) # showing first n rows

Ensembl_MGI_V1_2 = tidyr:: drop_na(Ensembl_MGI_V1, MGI.Accession.ID) # removal of rows with missing MGI IDs in MGI data (column “MGI.Accession.ID”)
head(Ensembl_MGI_V1_2, n = 4) # showing first n rows

#### Joining datasets step 2 ####

Ensembl_MGI_V2  = na_rows <- Ensembl_MGI_V1[is.na(Ensembl_MGI_V1$MGI.Accession.ID), ]  # selection of rows with missing MGI IDs in MGI data (column “MGI.Accession.ID”)
head(Ensembl_MGI_V2, n = 4) # showing first n rows

Ensembl_MGI_V2_2  = dplyr::mutate(Ensembl_MGI_V2,  
                                  MGI.Accession.ID = NULL,
                                  Marker.Symbol = NULL, 
                                  Feature.Type = NULL, 
                                  Ensembl_ID_in_MGI = NULL) # removing columns with missing data from MGI
head(Ensembl_MGI_V2_2, n = 4) # showing first n rows

Ensembl_MGI_V2_3 = dplyr:: left_join (Ensembl_MGI_V2_2, MGI_3, by = c("Gene.stable.ID" = "Ensembl_ID_in_MGI"), keep = TRUE, na_matches = "never") # joining databases (all data from the first dataset and matching data from the second dataset
head(Ensembl_MGI_V2_3, n = 4) # showing first n rows

#### Joining data from step 1 and 2 ####

Ensembl_MGI = dplyr:: bind_rows(Ensembl_MGI_V1_2, Ensembl_MGI_V2_3) # row binding 
head(Ensembl_MGI, n = 4) # showing first n rows

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

NCBI_2$Ensembl_ID_as_Gene.symbol <- ifelse(grepl("ENSMUSG00", NCBI_2$Symbol), 'yes', 'no') # Testing the presence of Ensembl IDs in the column with gene symbols
head(NCBI_2, n=4) # showing first n rows

Ensembl_ID_testNCBI  = dplyr::filter(NCBI_2, Ensembl_ID_as_Gene.symbol == "yes") # selection of rows with Ensembl IDs used as a gene symbols. 
head(Ensembl_ID_testNCBI, n=4) # showing first n rows

NCBI_2  = dplyr::filter(NCBI_2, Ensembl_ID_as_Gene.symbol == "no") # removal of rows with Ensembl IDs used as a gene symbols. 
head(NCBI_2, n=4) # showing first n rows

NCBI_2 = dplyr::mutate(NCBI_2,  
                       Ensembl_ID_as_Gene.symbol = NULL) # removing unnecessary columns
head(NCBI_2, n = 4) # showing first n rows

#### Joining datasets step 1 ####

Ensembl_MGI_NCBI_V1 = dplyr:: left_join (Ensembl_MGI, NCBI_2, by = c("NCBI.gene..formerly.Entrezgene..ID" = "NCBI.GeneID"), keep = TRUE, na_matches = "never") # joining databases (all data from the first dataset and matching data from the second dataset
head(Ensembl_MGI_NCBI_V1, n = 4) # showing first n rows

Ensembl_MGI_NCBI_V1_2 = tidyr:: drop_na(Ensembl_MGI_NCBI_V1, NCBI.GeneID) # removal of rows with missing NCBI IDs in NCBI data (column “NCBI.GeneID”)
head(Ensembl_MGI_NCBI_V1_2, n = 4) # showing first n rows

#### Joining datasets step 2 ####

Ensembl_MGI_NCBI_V2  = na_rows <- Ensembl_MGI_NCBI_V1[is.na(Ensembl_MGI_NCBI_V1$NCBI.GeneID), ]  # selection of rows with missing NCBI IDs in NCBI data (column “NCBI.GeneID”)
head(Ensembl_MGI_NCBI_V2, n = 4) # showing first n rows

Ensembl_MGI_NCBI_V2_2  = dplyr::mutate(Ensembl_MGI_NCBI_V2,  
                                       NCBI.GeneID = NULL,
                                       Symbol = NULL,
                                       Description = NULL,
                                       Gene.Type = NULL, 
                                       Ensembl.GeneIDs = NULL) # removing columns with missing data from NCBI
head(Ensembl_MGI_NCBI_V2_2, n = 4) # showing first n rows

Ensembl_MGI_NCBI_V2_3 = dplyr:: left_join (Ensembl_MGI_NCBI_V2_2, NCBI_2, by = c("Gene.stable.ID" = "Ensembl.GeneIDs"), keep = TRUE, na_matches = "never") # joining databases (all data from the first dataset and matching data from the second dataset
head(Ensembl_MGI_NCBI_V2_3, n = 4) # showing first n rows

#### Joining data from step 1 and 2 ####

Ensembl_MGI_NCBI = dplyr:: bind_rows(Ensembl_MGI_NCBI_V1_2, Ensembl_MGI_NCBI_V2_3) # row binding 
head(Ensembl_MGI_NCBI, n = 4) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(Ensembl_MGI_NCBI, `Gene.stable.ID`) ) # grouping data according to the items in Gene.stable.ID column (Ensembl IDs from Ensembl)
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
readr::write_csv(x = grouped, file = "Ensembl_MGI_NCBI_grouped.csv") # data export

Ensembl_MGI_NCBI_grouped = read.delim('Ensembl_MGI_NCBI_grouped.csv', sep=',', header = TRUE) #wczytywanie pliku
head(Ensembl_MGI_NCBI_grouped, n = 4) # showing first n rows

colnames(Ensembl_MGI_NCBI_grouped) <- c("Ensembl_ID_in_Ensembl", "Gene.symbol_in_Ensembl",
                                        "MGI_ID_in_Ensembl", 
                                        "NCBI_ID_in_Ensembl",
                                        "MGI_ID_in_MGI",
                                        "Gene.symbol_in_MGI",
                                        "Gene.Type_in_MGI",
                                        "Ensembl_ID_in_MGI",
                                        "NCBI_ID_in_NCBI", 
                                        "Gene.symbol_in_NCBI", 
                                        "Description_in_NCBI",
                                        "Gene.Type_in_NCBI", 
                                        "Ensembl_ID_in_NCBI") # setting final column names
head(Ensembl_MGI_NCBI_grouped, n=4) # showing first n rows

Ensembl_MGI_NCBI_grouped$AmbiguousMGI<- ifelse(grepl(", ", Ensembl_MGI_NCBI_grouped$Gene.symbol_in_MGI),'TRUE','FALSE') # testing for the presence of “,” in items in the column Gene.symbol_in_MGI. Presence of the “,” indicates that more than one gene symbol in MGI is assigned to Ensembl novel gene
head(Ensembl_MGI_NCBI_grouped, n=4) # showing first n rows

Ensembl_MGI_NCBI_grouped$AmbiguousNCBI<- ifelse(grepl(", ", Ensembl_MGI_NCBI_grouped$Gene.symbol_in_NCBI),'TRUE','FALSE') # testing for the presence of “,” in items in the column Gene.symbol_in_NCBI. Presence of the “,” indicates that more than one gene symbol in NCBI is assigned to Ensembl novel gene
head(Ensembl_MGI_NCBI_grouped, n=4) # showing first n rows

readr::write_csv(x = Ensembl_MGI_NCBI_grouped, file = "Ensembl_MGI_NCBI_grouped.csv") # data export

#### The End ####




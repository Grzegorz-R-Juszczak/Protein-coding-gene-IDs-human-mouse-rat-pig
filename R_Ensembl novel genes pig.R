########################################################################################################
#### Trial dataset for this R script can be downloaded from                                         ####
#### https://data.mendeley.com/datasets/454s2vw255/1/files/b5ca5b96-c7a0-48f1-aee9-43565b5c322a    ####
#### The data should be extracted from the compressed folder before running the script.             #### 
########################################################################################################

#### Ensembl data import ####

mart_export_1 <- read.delim('mart_export.txt', sep=',', header = TRUE) # data importing
head(mart_export_1, n=5) # showing first n rows

library(naniar)
mart_export_1 = mart_export_1  %>%
  naniar::replace_with_na(replace = list(VGNC.ID = c(""))) # marking missing items in selected column with „NA”.
head(mart_export_1, n=5) # showing first n rows

#### VGNC data import and initial processing ####

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
                      locus_type = NULL,
                      uniprot_ids = NULL,
                      pubmed_id = NULL,
                      horde_id = NULL,
                      status = NULL,
                      error_column = NULL,
                      hgnc_orthologs = NULL) # removing unnecessary columns
head(vgnc3, n=4) # showing first n rows

vgnc3$Ensembl_ID_as_Gene.symbol <- ifelse(grepl("ENSSSCG00", vgnc3$symbol), 'yes', 'no') # Testing the presence of Ensembl IDs in the column with gene symbols
head(vgnc3, n=4) # showing first n rows

Ensembl_ID_testVGNC  = dplyr::filter(vgnc3, Ensembl_ID_as_Gene.symbol == "yes") # selection of rows with Ensembl IDs used as a gene symbols. 
head(Ensembl_ID_testVGNC, n=4) # showing first n rows


vgnc3  = dplyr::filter(vgnc3, Ensembl_ID_as_Gene.symbol == "no") # removal of rows with Ensembl IDs used as a gene symbols. 
head(vgnc3, n=4) # showing first n rows

vgnc3 = dplyr::mutate(vgnc3,  
                      Ensembl_ID_as_Gene.symbol = NULL) # removing unnecessary columns
head(vgnc3, n = 4) # showing first n rows

#### Joining datasets step 1 ####
Ensembl_VGNC_V1 = dplyr:: left_join (mart_export_1, vgnc3, by = c("VGNC.ID" = "vgnc_id"), keep = TRUE, na_matches = "never") # joining databases (all data from the first dataset and matching data from the second dataset
head(Ensembl_VGNC_V1, n = 4) # showing first n rows

Ensembl_VGNC_V1_2 = tidyr:: drop_na(Ensembl_VGNC_V1, vgnc_id) # removal of rows with missing VGNC IDs in VGNC data (column “vgnc_id”)
head(Ensembl_VGNC_V1_2, n = 4) # showing first n rows

#### Joining datasets step 2 ####

Ensembl_VGNC_V2  = na_rows <- Ensembl_VGNC_V1[is.na(Ensembl_VGNC_V1$vgnc_id), ]  # selection of rows with missing VGNC IDs in VGNC data (column “vgnc_id”)
head(Ensembl_VGNC_V2, n = 4) # showing first n rows

Ensembl_VGNC_V2_2  = dplyr::mutate(Ensembl_VGNC_V2,  
                                   vgnc_id = NULL,
                                   symbol = NULL, 
                                   locus_group = NULL, 
                                   ensembl_gene_id = NULL) # removing columns with missing data from VGNC
head(Ensembl_VGNC_V2_2, n = 4) # showing first n rows

Ensembl_VGNC_V2_3 = dplyr:: left_join (Ensembl_VGNC_V2_2, vgnc3, by = c("Gene.stable.ID" = "ensembl_gene_id"), keep = TRUE, na_matches = "never") # joining databases (all data from the first dataset and matching data from the second dataset
head(Ensembl_VGNC_V2_3, n = 4) # showing first n rows

#### Joining data from step 1 and 2 ####

Ensembl_VGNC = dplyr:: bind_rows(Ensembl_VGNC_V1_2, Ensembl_VGNC_V2_3) # row binding 
head(Ensembl_VGNC, n = 4) # showing first n rows

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

NCBI_2$Ensembl_ID_as_Gene.symbol <- ifelse(grepl("ENSSSCG00", NCBI_2$Symbol), 'yes', 'no') # Testing the presence of Ensembl IDs in the column with gene symbols
head(NCBI_2, n=4) # showing first n rows

Ensembl_ID_testNCBI  = dplyr::filter(NCBI_2, Ensembl_ID_as_Gene.symbol == "yes") # selection of rows with Ensembl IDs used as a gene symbols. 
head(Ensembl_ID_testNCBI, n=4) # showing first n rows

NCBI_2  = dplyr::filter(NCBI_2, Ensembl_ID_as_Gene.symbol == "no") # removal of rows with Ensembl IDs used as a gene symbols. 
head(NCBI_2, n=4) # showing first n rows

NCBI_2 = dplyr::mutate(NCBI_2,  
                       Ensembl_ID_as_Gene.symbol = NULL) # removing unnecessary columns
head(NCBI_2, n = 4) # showing first n rows

#### Joining datasets step 1 ####

Ensembl_VGNC_NCBI_V1 = dplyr:: left_join (Ensembl_VGNC, NCBI_2, by = c("NCBI.gene..formerly.Entrezgene..ID" = "NCBI.GeneID"), keep = TRUE, na_matches = "never") # joining databases (all data from the first dataset and matching data from the second dataset
head(Ensembl_VGNC_NCBI_V1, n = 4) # showing first n rows

Ensembl_VGNC_NCBI_V1_2 = tidyr:: drop_na(Ensembl_VGNC_NCBI_V1, NCBI.GeneID) # removal of rows with missing NCBI IDs in NCBI data (column “NCBI.GeneID”)
head(Ensembl_VGNC_NCBI_V1_2, n = 4) # showing first n rows

#### Joining datasets step 2 ####

Ensembl_VGNC_NCBI_V2  = na_rows <- Ensembl_VGNC_NCBI_V1[is.na(Ensembl_VGNC_NCBI_V1$NCBI.GeneID), ]  # selection of rows with missing NCBI IDs in NCBI data (column “NCBI.GeneID”)
head(Ensembl_VGNC_NCBI_V2, n = 4) # showing first n rows

Ensembl_VGNC_NCBI_V2_2  = dplyr::mutate(Ensembl_VGNC_NCBI_V2,  
                                        NCBI.GeneID = NULL,
                                        Symbol = NULL,
                                        Description = NULL,
                                        Gene.Type = NULL, 
                                        Ensembl.GeneIDs = NULL) # removing columns with missing data from NCBI
head(Ensembl_VGNC_NCBI_V2_2, n = 4) # showing first n rows

Ensembl_VGNC_NCBI_V2_3 = dplyr:: left_join (Ensembl_VGNC_NCBI_V2_2, NCBI_2, by = c("Gene.stable.ID" = "Ensembl.GeneIDs"), keep = TRUE, na_matches = "never") # joining databases (all data from the first dataset and matching data from the second dataset
head(Ensembl_VGNC_NCBI_V2_3, n = 4) # showing first n rows

#### Joining data from step 1 and 2 ####

Ensembl_VGNC_NCBI = dplyr:: bind_rows(Ensembl_VGNC_NCBI_V1_2, Ensembl_VGNC_NCBI_V2_3) # row binding 
head(Ensembl_VGNC_NCBI, n = 4) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(Ensembl_VGNC_NCBI, `Gene.stable.ID`) ) # grouping data according to the items in Gene.stable.ID column (Ensembl IDs from Ensembl)
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
readr::write_csv(x = grouped, file = "Ensembl_VGNC_NCBI_grouped.csv") # data export

Ensembl_VGNC_NCBI_grouped = read.delim('Ensembl_VGNC_NCBI_grouped.csv', sep=',', header = TRUE) # data import
head(Ensembl_VGNC_NCBI_grouped, n = 4) # showing first n rows

colnames(Ensembl_VGNC_NCBI_grouped) <- c("Ensembl_ID_in_Ensembl", "Gene.symbol_in_Ensembl", 
                                         "VGNC_ID_in_Ensembl", 
                                         "NCBI_ID_in_Ensembl",
                                         "VGNC_ID_in_VGNC",
                                         "Gene.symbol_in_VGNC",
                                         "Gene.Type_in_VGNC",
                                         "Ensembl_ID_in_VGNC",
                                         "NCBI_ID_in_NCBI", 
                                         "Gene.symbol_in_NCBI", 
                                         "Description_in_NCBI",
                                         "Gene.Type_in_NCBI", 
                                         "Ensembl_ID_in_NCBI") # setting final column names
head(Ensembl_VGNC_NCBI_grouped, n=4) # showing first n rows

Ensembl_VGNC_NCBI_grouped$AmbiguousVGNC<- ifelse(grepl(", ", Ensembl_VGNC_NCBI_grouped$Gene.symbol_in_VGNC),'TRUE','FALSE') # testing for the presence of “,” in items in the column Gene.symbol_in_VGNC. Presence of the “,” indicates that more than one gene symbol in VGNC is assigned to Ensembl novel gene
head(Ensembl_VGNC_NCBI_grouped, n=4) # showing first n rows

Ensembl_VGNC_NCBI_grouped$AmbiguousNCBI<- ifelse(grepl(", ", Ensembl_VGNC_NCBI_grouped$Gene.symbol_in_NCBI),'TRUE','FALSE') # testing for the presence of “,” in items in the column Gene.symbol_in_NCBI. Presence of the “,” indicates that more than one gene symbol in NCBI is assigned to Ensembl novel gene
head(Ensembl_VGNC_NCBI_grouped, n=4) # showing first n rows

readr::write_csv(x = Ensembl_VGNC_NCBI_grouped, file = "Ensembl_VGNC_NCBI_grouped.csv") # data export

#### The End ####




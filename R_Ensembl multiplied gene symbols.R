mart_export_1 <- read.delim('mart_export.txt', sep=',', header = TRUE) # data import
head(mart_export_1, n=4) # showing first n rows

mart_export_1$Scaffold <- ifelse(grepl("[_.]", mart_export_1$Chromosome.scaffold.name),'TRUE','FALSE') # Testing the presence of “.” or “_”in the items in the column “Chromosome.scaffold.name”. Presence of a dot or _ indicates that the Ensemble gene ID is assigned to the scaffold instead of chromosome. 
head(mart_export_1, n=20) # showing first n rows
library(dplyr) # package selection

mart_export_1_Scaffold = dplyr::filter(mart_export_1, Scaffold == "TRUE") # selection of gene IDs assigned to scaffolds. 
head(mart_export_1_Scaffold, n=10) # showing first n rows
readr::write_csv(x = mart_export_1_Scaffold, file = "Genes_with_scaffold_locus.csv") # data export

mart_export_1_Chromosome = dplyr::filter(mart_export_1, Scaffold == "FALSE") # selection of gene IDs assigned to chromosomes.
head(mart_export_1_Chromosome, n=10) # showing first n rows

mart_export_2_Chromosome = dplyr::mutate(mart_export_1_Chromosome, Scaffold = NULL) # removing unnecessary column 
head(mart_export_2_Chromosome, n=4) # showing first n rows

mart_export_3 <- tidyr::nest( dplyr::group_by(mart_export_2_Chromosome, `Gene.name`) ) # grouping data according to gene symbols in the column “Gene.name”.
mart_export_4 <- purrr::map_dfr(
  .x = mart_export_3$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        
        paste(unique(column), collapse = ", ") 
      })
  } )
mart_export_3  <- cbind(mart_export_3, mart_export_4) 
mart_export_3$data <- NULL
head(mart_export_3, n=4) # showing first n rows
readr::write_csv(x = mart_export_3, file = "symbol_grouped.csv") # data export

symbol_grouped <- read.delim('symbol_grouped.csv', sep=',', header = TRUE) # data import
head(symbol_grouped, n=4) # showing first n rows

symbol_grouped$Duplicated <- ifelse(grepl(", ", symbol_grouped$Gene.stable.ID),'TRUE','FALSE') # testing presence of “,” in items in the column “Gene.stable.ID” indicating that more than one Ensembl ID is assigned to single gene symbol in the dataset restricted to genes assigned only to chromosomes (scaffolds removed). 
head(symbol_grouped, n=20) # showing first n rows

symbol_grouped_Duplicated = dplyr::filter(symbol_grouped, Duplicated == "TRUE") # selecting gene symbols with assigned more than one Ensembl gene ID attributed to chromosomes. 
head(symbol_grouped_Duplicated, n=10) # showing first n rows

symbol_grouped_Duplicated$DifferentChromosemes <- ifelse(grepl(", ", symbol_grouped_Duplicated$Chromosome.scaffold.name),'TRUE','FALSE') # Testing the presence of “,” in items present in column “Chromosome.scaffold.name”. Presence of the  “,” indicates that the multiplied gene symbol is assigned to different chromosomes 
head(symbol_grouped_Duplicated, n=20) # showing first n rows

symbol_grouped_Duplicated_DifferentChromosemes = dplyr::filter(symbol_grouped_Duplicated, DifferentChromosemes == "TRUE") # selection of multiplied gene symbols assigned to different chromosomes.
head(symbol_grouped_Duplicated_DifferentChromosemes, n=10) # showing first n rows

symbol_grouped_Duplicated_DifferentChromosemes$Somatic <- ifelse(grepl("[1-9]", symbol_grouped_Duplicated_DifferentChromosemes$Chromosome.scaffold.name),'TRUE','FALSE') # testing for presence of numbers in chromosome names indicating assignment to somatic chromosomes.
head(symbol_grouped_Duplicated_DifferentChromosemes, n=20) # showing first n rows

symbol_grouped_Duplicated_DifferentChromosemes_Somatic = dplyr::filter(symbol_grouped_Duplicated_DifferentChromosemes, Somatic == "TRUE") # selection of gene symbols assigned to somatic chromosomes
head(symbol_grouped_Duplicated_DifferentChromosemes_Somatic, n=10) # showing first n rows
readr::write_csv(x = symbol_grouped_Duplicated_DifferentChromosemes_Somatic, file = "Multiplied_symbols_DifferentChromosemes_Autosomal.csv")

symbol_grouped_Duplicated_DifferentChromosemes_Sex = dplyr::filter(symbol_grouped_Duplicated_DifferentChromosemes, Somatic == "FALSE") # selection of gene symbols assigned only to sex chromosomes
head(symbol_grouped_Duplicated_DifferentChromosemes_Sex, n=10) # showing first n rows
readr::write_csv(x = symbol_grouped_Duplicated_DifferentChromosemes_Sex, file = "Multiplied_symbols_DifferentChromosemes_Sex.csv") # data export

symbol_grouped_Duplicated_TheSameChromosome = dplyr::filter(symbol_grouped_Duplicated, DifferentChromosemes == "FALSE") #  selection of gene symbols assigned to single chromosome
head(symbol_grouped_Duplicated_TheSameChromosome, n=10) # showing first n rows

symbol_grouped_Duplicated_TheSameChromosome$different_strand <- ifelse(grepl(",", symbol_grouped_Duplicated_TheSameChromosome$Strand),'TRUE','FALSE') # testing for presence of “,” in items from the column “strand”. Presence of “,” indicates that multiplied gene symbols are assigned to different strands. 
head(symbol_grouped_Duplicated_TheSameChromosome, n=20) # showing first n rows

symbol_grouped_Duplicated_TheSameChromosome_different_strand = dplyr::filter(symbol_grouped_Duplicated_TheSameChromosome, different_strand == "TRUE") # selection of multiplied gene symbols assigned to different strands
head(symbol_grouped_Duplicated_TheSameChromosome_different_strand, n=10) # showing first n rows
readr::write_csv(x = symbol_grouped_Duplicated_TheSameChromosome_different_strand, file = "Multiplied_symbols_TheSameChromosome_Different_strand.csv") # data export

symbol_grouped_Duplicated_TheSameChromosome_TheSame_strand = dplyr::filter(symbol_grouped_Duplicated_TheSameChromosome, different_strand == "FALSE") # selection of multiplied gene symbols assigned to the same strand strands
head(symbol_grouped_Duplicated_TheSameChromosome_TheSame_strand, n=10) # showing first n rows

symbol_grouped_Duplicated_TheSameChromosome_TheSame_strand2 = dplyr::mutate(symbol_grouped_Duplicated_TheSameChromosome_TheSame_strand, Duplicated
                                                                            = NULL, DifferentChromosemes = NULL, different_strand = NULL) # removing unnecessary columns
head(symbol_grouped_Duplicated_TheSameChromosome_TheSame_strand2, n=4) # showing first n rows

symbol_grouped_Duplicated_TheSameChromosome_TheSame_strand2$Start <- ifelse(grepl(", ", symbol_grouped_Duplicated_TheSameChromosome_TheSame_strand2$Gene.start..bp.),'TRUE','FALSE') # testing for the presence of “,” in items in the “Gene.start..bp.” column. Presence of “,”  indicates different start points of multiplied gene symbols. 
head(symbol_grouped_Duplicated_TheSameChromosome_TheSame_strand2, n=20) # showing first n rows

symbol_grouped_Duplicated_TheSameChromosome_The_same_start = dplyr::filter(symbol_grouped_Duplicated_TheSameChromosome_TheSame_strand2, Start == "FALSE") # selection of multiplied gene symbols with the same starting point
head(symbol_grouped_Duplicated_TheSameChromosome_The_same_start, n=20) # showing first n rows
readr::write_csv(x = symbol_grouped_Duplicated_TheSameChromosome_The_same_start, file = "Multiplied_symbols_TheSameChromosome_The_same_start.csv") # data export

The_same_chromosome_Different_start = dplyr::filter(symbol_grouped_Duplicated_TheSameChromosome_TheSame_strand2, Start == "TRUE") # selection of multiplied gene symbols with the different starting point
head(The_same_chromosome_Different_start, n=20) # showing first n rows

The_same_chromosome_Different_start = dplyr::mutate(The_same_chromosome_Different_start, Start = NULL) # removing unnecessary columns
head(The_same_chromosome_Different_start, n=4) # showing first n rows

library(tidyr) # package selection
The_same_chromosome_Different_start2 = The_same_chromosome_Different_start %>% separate_wider_delim(Gene.start..bp., delim = ', ', names_sep = "", too_few = "align_start") # splitting grouped starting points into separate columns
head(The_same_chromosome_Different_start2, n=4) # showing first n rows

The_same_chromosome_Different_start3 = The_same_chromosome_Different_start2 %>% separate_wider_delim(Gene.end..bp., delim = ', ', names_sep = "", too_few = "align_start") # splitting grouped end points into separate columns
head(The_same_chromosome_Different_start3, n=4) # showing first n rows

readr::write_csv(x = The_same_chromosome_Different_start3, file = "Multiplied_symbols_TheSameChromosome_Different_start.csv") # data export

Multiplied_symbols_TheSameChromosome_Different_start <- read.delim('Multiplied_symbols_TheSameChromosome_Different_start.csv', sep=',', header = TRUE) # data import
head(Multiplied_symbols_TheSameChromosome_Different_start, n=4) # showing first n rows

Multiplied_symbols_TheSameChromosome_Different_start$Overlap1v2 <- ifelse(
  Multiplied_symbols_TheSameChromosome_Different_start$Gene.end..bp.1 < Multiplied_symbols_TheSameChromosome_Different_start$Gene.start..bp.2,
  'no overlap',
  'X') # testing whether the end position of a gene with assigned one Ensembl ID is smaller than the start position of a gene with another Ensembl ID to identify non overlapping genes
head(Multiplied_symbols_TheSameChromosome_Different_start, n=4) # showing first n rows

Multiplied_symbols_TheSameChromosome_Different_start$Overlap2v1 <- ifelse(
  Multiplied_symbols_TheSameChromosome_Different_start$Gene.end..bp.2 < Multiplied_symbols_TheSameChromosome_Different_start$Gene.start..bp.1,
  'no overlap',
  'X') # testing whether the end position of a gene with assigned one Ensembl ID is smaller than the start position of a gene with another Ensembl ID to identify non overlapping genes
head(Multiplied_symbols_TheSameChromosome_Different_start, n=4) # showing first n rows

Multiplied_symbols_TheSameChromosome_Different_start$Overlap1v3 <- ifelse(
  Multiplied_symbols_TheSameChromosome_Different_start$Gene.end..bp.1 < Multiplied_symbols_TheSameChromosome_Different_start$Gene.start..bp.3,
  'no overlap',
  'X') # testing whether the end position of a gene with assigned one Ensembl ID is smaller than the start position of a gene with another Ensembl ID to identify non overlapping genes
head(Multiplied_symbols_TheSameChromosome_Different_start, n=100) # showing first n rows

Multiplied_symbols_TheSameChromosome_Different_start$Overlap3v1 <- ifelse(
  Multiplied_symbols_TheSameChromosome_Different_start$Gene.end..bp.3 < Multiplied_symbols_TheSameChromosome_Different_start$Gene.start..bp.1,
  'no overlap',
  'X') # testing whether the end position of a gene with assigned one Ensembl ID is smaller than the start position of a gene with another Ensembl ID to identify non overlapping genes
head(Multiplied_symbols_TheSameChromosome_Different_start, n=100) # showing first n rows

Multiplied_symbols_TheSameChromosome_Different_start$Overlap2v3 <- ifelse(
  Multiplied_symbols_TheSameChromosome_Different_start$Gene.end..bp.2 < Multiplied_symbols_TheSameChromosome_Different_start$Gene.start..bp.3,
  'no overlap',
  'X') # testing whether the end position of a gene with assigned one Ensembl ID is smaller than the start position of a gene with another Ensembl ID to identify non overlapping genes
head(Multiplied_symbols_TheSameChromosome_Different_start, n=100) # showing first n rows

Multiplied_symbols_TheSameChromosome_Different_start$Overlap3v2 <- ifelse(
  Multiplied_symbols_TheSameChromosome_Different_start$Gene.end..bp.3 < Multiplied_symbols_TheSameChromosome_Different_start$Gene.start..bp.2,
  'no overlap',
  'X') # testing whether the end position of a gene with assigned one Ensembl ID is smaller than the start position of a gene with another Ensembl ID to identify non overlapping genes
head(Multiplied_symbols_TheSameChromosome_Different_start, n=100) # showing first n rows

readr::write_csv(x = Multiplied_symbols_TheSameChromosome_Different_start, file = "Multiplied_symbols_TheSameChromosome_Different_start_final.csv") # data export


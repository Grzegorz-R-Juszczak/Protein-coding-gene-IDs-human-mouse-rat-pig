AGS1 = read.delim('Ambiguous_gene_synonyms_final.csv', sep=',', header = TRUE) # importing data
head(AGS1, n=4) # showing first n rows

AGS2 = dplyr::select (AGS1, Assigned_Official_symbols) # selecting a column with assigned official symbols
head(AGS2, n=4) # showing first n rows

library(tidyr) # package selection
AGS3 = AGS2 %>% separate_longer_delim(Assigned_Official_symbols, delim = ', ') # Splitting official symbols into separate rows
head(AGS3, n=4) # showing first n rows
colnames(AGS3) <- c("Gene.name") # changing column name to enable merging the data from two input files
head(AGS3, n=7) # showing first n rows

AS1 = read.delim('Ambiguous_symbols_protein_protein_final.csv', sep=',', header = TRUE) # data import
head(AS1, n=4) # showing first n rows

AS2A = dplyr::select (AS1, Alternative_official_symbols) # selecting column with alternative official symbols
head(AS2A, n=7) # showing first n rows

library(tidyr) # package selection
AS3A = AS2A %>% separate_longer_delim(Alternative_official_symbols, delim = ', ') # Splitting official symbols into separate rows
head(AS3A, n=7) # showing first n rows
colnames(AS3A) <- c("Gene.name") # changing column name to enable merging the data from two input files
head(AS3A, n=7) # showing first n rows

AS2B = dplyr::select (AS1, Ambiguous_official_symbols) # selecting column with ambiguous official symbols
head(AS2B, n=7) # showing first n rows
colnames(AS2B) <- c("Gene.name") # changing column name to enable merging the data from two input files
head(AS2B, n=7) # showing first n rows

all = dplyr:: bind_rows(AGS3, AS3A, AS2B) # data binding
head(all, n=7) # showing first n rows

all2 = dplyr::mutate(all,  tolower_Gene_name = tolower(Gene.name)) # adding a new column with gene symbols written only in lower case letters
head(all2, n=7) # showing first n rows

all2uniq = unique(all2$tolower_Gene_name) # selecting unique gene symbols
head(all2uniq, n=7) # showing first n rows

all2uniq2 = data.frame(all2uniq) # changing string into data frame
head(all2uniq2, n=7) # showing first n rows
readr::write_csv(x = all2uniq2, file = "All_genes_linked_to_ambiguous_symbols.csv", col_names = FALSE) # data export



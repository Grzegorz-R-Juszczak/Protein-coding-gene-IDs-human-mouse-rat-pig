########################################################################################################
#### Trial dataset for this R script can be downloaded from                                         ####
#### https://data.mendeley.com/datasets/454s2vw255/1/files/43562d2c-c52d-42ff-ba93-da869fb83453     ####
#### The data should be extracted from the compressed folder before running the script.             #### 
########################################################################################################

mart_export_1 <- read.delim('mart_export.txt', sep=',', header = TRUE) # importing data
head(mart_export_1, n=4) # showing first n rows

mart_export_2 = dplyr::mutate(mart_export_1,  tolower_Gene_Synonym = tolower(Gene.Synonym)) # adding new column with items from Gene.Synonym column written only with lowercase letters
head(mart_export_2, n=4) # showing first n rows

mart_export_3 = dplyr::mutate(mart_export_2,  tolower_Gene_name = tolower(Gene.name)) # adding new column with items from Gene.name column written only with lowercase letters
head(mart_export_3, n=4) # showing first n rows

library(naniar) # package selection
mart_export_4 = mart_export_3 %>%
  naniar::replace_with_na(replace = list(Gene.name = c(""), Gene.Synonym = c(""), Gene.stable.ID = c(""))) # marking missing items in selected columns with “NA”
head(mart_export_4, n=4) # showing first n rows

mart_export_5 = tidyr:: drop_na(mart_export_4, Gene.name) # removing rows with missing gene symbols (“NA” in column Gene.name).
head(mart_export_5, n=4) # showing first n rows

Unique_gene_symbols = unique(mart_export_5$tolower_Gene_name) # list of unique gene symbols
head(Unique_gene_symbols, n=4) # showing first n rows

write.table(x = Unique_gene_symbols, file="Unique_gene_symbols.csv", row.names=F, col.names=F) # data export

mart_export_6 = tidyr:: drop_na(mart_export_4, Gene.Synonym) # removing rows without synonyms („NA” in column Gene.Synonym). 
head(mart_export_6, n=4) # showing first n rows

Unique_gene_synonyms = unique(mart_export_6$tolower_Gene_Synonym) # list o unique synonyms
head(Unique_gene_synonyms, n=4) # showing first n rows
write.table(x = Unique_gene_synonyms, file="Unique_gene_synonyms.csv", row.names=F, col.names=F) # data export

VennDiagram::venn.diagram(
  x = list(Unique_gene_symbols, Unique_gene_synonyms), # selection of data for the venn diagram
  category.names = c("Unique_gene_symbols" , "Unique_gene_synonyms"), # setting labels for data in the venn diagram
  
  filename = 'venn.png', # output file with venn diagram
  output=TRUE,
  imagetype="png" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 300,
  lwd = 3, # line width of the circles' circumferences
  col=c('red', 'green'),  #  colours of the circles' circumferences
  cex = 1.0, # size of the areas' labels
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = .6, # size of the category names
  cat.default.pos = "outer", # default location of category names 
  cat.dist = 0.03, # the distances of the category names from the edges of the circles
  cat.pos = -180, # the positions (in degrees) of the category names along the circles, with 0 (default) at the 12 o'clock location
  ext.pos = -180, # the positions (in degrees) of the external area labels along the circles, with 0 (default) at 12 o'clock (line indicating small overlapping areas)
  ext.dist = -0.19, # how far to place the external area labels relative to its anchor point (setting shorter than the default) 
) # venn diagram


Mouse <- read.csv("Mice_All_unique_protein_coding_symbols.csv", header = TRUE) # bez header = FALSE # data import
head(Mouse, n=5) # showing first n rows

Mouse_2 = dplyr::mutate(Mouse,  tolower_Gene.name = tolower(Gene.name)) # adding new column with gene symbols written in lower case letters
head(Mouse_2, n=5) # showing first n rows

Mouse_3 <- Mouse_2 [,2] # creating numeric vector with lowercase symbols without column name
head(Mouse_3, n=20) # showing first n rows

Rat <- read.csv("Rat_All_unique_protein_coding_symbols.csv", header = TRUE) # data import
head(Rat, n=5) # showing first n rows

Rat_2 = dplyr::mutate(Rat,  tolower_Gene.name = tolower(Gene.name)) # adding new column with gene symbols written in lower case letters
head(Rat_2, n=5) # showing first n rows

Rat_3 <- Rat_2 [,2] # creating numeric vector with lowercase symbols without column name
head(Rat_3, n=20) # showing first n rows

Human <- read.csv("Human_All_unique_protein_coding_symbols.csv", header = TRUE) # data import
head(Human, n=5) # showing first n rows

Human_2 = dplyr::mutate(Human,  tolower_Gene.name = tolower(Gene.name)) # adding new column with gene symbols written in lower case letters
head(Human_2, n=5) # showing first n rows

Human_3 <- Human_2 [,2] # creating numeric vector with lowercase symbols without column name
head(Human_3, n=20) # showing first n rows


Pig = read.csv("Pig_All_unique_protein_coding_symbols.csv", header = TRUE) # data import
head(Pig, n=5) # showing first n rows

Pig_2 = dplyr::mutate(Pig,  tolower_Gene.name = tolower(Gene.name)) # adding new column with gene symbols written in lower case letters
head(Pig_2, n=5) # showing first n rows

Pig_3 <- Pig_2 [,2] # creating numeric vector with lowercase symbols without column name
head(Pig_3, n=20) # showing first n rows

library(VennDiagram) # package selection
venn.diagram(
  x = list(Rat_3, Human_3, Mouse_3, Pig_3), # # selection of data for the venn diagram
  category.names = c("Rat", "Human", "Mouse", "Pig"), # setting labels for data in the venn diagram
  
  filename = 'venn.png', # output file with venn diagram
  output=TRUE,
  imagetype="png" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 300,
  lwd = 3, # line width of the circles' circumferences
  col=c('blue', 'green', 'red', 'black'),  #  colours of the circles' circumferences
  cex = 1.0, # size of the areas' labels
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = .6, # size of the category names
  cat.default.pos = "outer", # default location of category names 
  cat.dist = .09, # the distances of the category names from the edges of the circles
)

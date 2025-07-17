# Protein-coding-gene-IDs-human-mouse-rat-pig

Repository Protein-coding-gene-IDs-human-mouse-rat-pig

Detailed information is available in manuscript “Protein-coding genes in humans and model mammals (mouse, rat and pig): gene identifiers and disambiguation of gene nomenclature retrieved from the Ensembl genome browser” submitted to BMC Genomics (Juszczak et al., 2025). Below we provide short description of repository and R programs. 

Application
The repository contains programs enabling extraction of information about gene identifiers from Ensembl / Biomart genome datasets.

Background
Gene nomenclature changed over time due to new discoveries and between-species standardization efforts. Therefore, it contains both current official symbols and synonyms (aliases) that were used in the past. The issue is further complicated by the fact that each gene has a different history of scientific investigation and associated assignment of gradually refined names. Therefore, genetic datasets contain both gene symbols that did not changed since the time of publication (official symbols) and some number of obsolete symbols (synonyms) that were updated after the publication of the data. This means that comparison between genetic datasets obtained from different sources requires standardization of gene nomenclature. 
 
##################################################################################################################################################################################################

Programs R_Ensembl gene symbol search (REgeness) adjusted for individual species:
R_Ensembl gene symbol search mouse.R
R_ensembl gene symbol search human.R
R_ensembl gene symbol search pig.R
R_ensembl gene symbol search rat.R

Application
Updating gene nomenclature in case when only gene symbols are available and checking for ambiguous assignment between gene symbols.  Symbol ambiguity occurs in case when different genes share the same symbol created from different gene names. 

Description
The program imports list of gene symbols, performs a double Ensembl search with input symbols treated as official symbols and synonyms and integrates data to provide a single list of updated symbols with assigned additional IDs and annotation about unique or ambiguous association between input and output gene symbols. Running the script requires installation of the biomaRt, tidyr, dplyr, purrr and readr R packages. 

Input file format
Data requirement is a list of gene symbols in the csv file named InputData with one column named Input_gene_symbols 

Example of input data (data arrangement)

Input_gene_symbols

Ttr

Kl

Arg 

Output data
The updated official symbols together with information about their ambiguity, gene description and additional IDs is exported to the Final_search_results.csv file. Additionally, the script imports and saves biomaRt codes for selecting Ensembl databases (file EnsemblDatabases.csv), datasets (EnsemblDatasets.csv), filters (file EnsemblFilters.csv) and output data called attributes in Ensembl (EnsemblAttributes.csv). These data provide information about current code including the version of Ensembl used for data downloading (file EnsemblDatabases.csv) and enable code modifications including species (EnsemblDatasets.csv) and output data associated with updated genes (file EnsemblAttributes.csv).

##################################################################################################################################################################################################

Program
R_protein-coding unique.R

Application
The program extracts lists of unique gene symbols, gene symbols with more than one stable Ensembl ID and novel genes with assigned Ensembl IDs but without a gene symbol. The output files are used by other programs deposited in this repository. 

Description
The program extracts information about gene symbols. Additionally, it adds a column with gene symbols combined with the word “Gene” (for example Gene_Usp9y) to prevent an unintentional transformation of symbols by Excel. These composite identifiers can be used to check whether there are any unintentional transformation of gene symbols after uploading the data to Excel. The script used for data processing requires installation of naniar, tidyr, purrr, readr and dplyr R packages.  

Input file format
Csv file with a default name “mart_export” downloaded from 
https://www.ensembl.org/biomart/martview/  using Ensembl Genes database containing genetic datasets for various species. The selected Ensembl filter setting should be  “Gene / Gene type / protein_coding” while the attributes settings (data included in the result file) should be “Gene stable ID” and “Gene name”. 

Example of input data (data arrangement)

Gene stable ID,Gene name
ENSMUSG00000064341,mt-Nd1
ENSMUSG00000064345,mt-Nd2
ENSMUSG00000064351,mt-Co1
ENSMUSG00000064354,mt-Co2

Output data
Results are exported to following files:
All_unique_protein_coding.csv (list of unique gene symbols together with Ensembl stable IDs)
All_unique_protein_coding_symbols.csv (lists of unique gene symbols), 
Ensembl_multiplied_gene_symbols.csv (lists of gene symbols with more than one stable Ensembl ID) 
Genes_without_gene_symbol.csv (lists of novel genes with assigned Ensembl IDs but without a gene symbol)

##################################################################################################################################################################################################

Program
R_venn between species comparison.R

Application
Creation of the Venn diagram showing a number of overlapping gene symbols in different species. 

Description
Creation of the Venn diagram showing a number of overlapping gene symbols in different species. It is required to install the dplyr and VennDiagram packages to run the script. 

Input file format
Files All_unique_protein_coding_symbols.csv after adding  prefix Human_, Mice_, Pig_ or Rat_ to the file name depending on the analyzed species (for example Mice_All_unique_protein_coding_symbols.csv).  The file All_unique_protein_coding_symbols.csv  is generated by the script R_protein-coding unique. 

Example of input data (data arrangement)	

Gene.name
mt-Nd1
mt-Nd2
mt-Co1
mt-Co2

Output data
Venn.png file

##################################################################################################################################################################################################

Program
R_Ensembl multiplied gene symbols.R

Application 
Extraction of information concerning genomic localization of genes with assigned more than one Ensembl stable ID.

Description
The analysis includes separation of Ensembl IDs with chromosomal localization from IDs with scaffold localization.  The data restricted to genes with Ensembl IDs assigned only to chromosomes are again tested for the presence of symbols with more than one Ensembl ID. The new list containing genes with more than one Ensembl ID assigned to chromosome is analyzed to identify cases when alternative Ensembl IDs are assigned to different types of chromosomal localization. More details are provided in description of output files. Additionally, alternative Ensembl IDs assigned to a different start position at the same strand are tested for an overlap by comparison between end and start positions of alternative IDs. The number of comparisons should be adjusted for each species based on the maximum number of  alternative Ensembl IDs assigned to a single gene symbol. The script used for data processing requires installation of dplyr, tidyr, readr, dplyr and purrr R packages

Input file format
Csv file with a default name “mart_export” downloaded from 
https://www.ensembl.org/biomart/martview/  using Ensembl Genes database with genetic datasets for various species. Selected Ensembl filter setting should be Gene /  Input external references ID list / Gene stable ID(s). The Ensembl gene stable IDs (e.g. ENSMUSG00000020807) for external references list is available in file Ensembl_multiplied_gene_symbols.csv generated by the script R_protein-coding unique.R.  The output Ensembl attributes  should be Gene stable ID, Gene name, Chromosome/scaffold name, Gene start (bp), Gene end (bp) and Strand.


Example of input data (data arrangement)

Gene stable ID,Gene name,Chromosome/scaffold name,Gene start (bp),Gene end (bp),Strand
ENSMUSG00000006378,Gcat,15,78915074,78922553,1
ENSMUSG00000006471,Ndor1,2,25134833,25146034,-1
ENSMUSG00000020807,4933427D14Rik,11,72044755,72098285,-1

Output data
Genes_with_scaffold_locus.csv (list of IDs assigned to scaffolds) 
Multiplied_symbols_DifferentChromosemes_Somatic.csv (list of alternative Ensembl IDs that are assigned to different chromosomes including at least one somatic chromosome) 
Multiplied_symbols_DifferentChromosemes_Sex.csv (list of alternative Ensembl IDs that are assigned only to sex chromosomes)
Multiplied_symbols_TheSameChromosome_Different_strand.csv (list of alternative Ensembl IDs that are assigned to the same chromosome but opposite strands)
Multiplied_symbols_TheSameChromosome_The_same_start.csv  (list of alternative Ensembl IDs that are assigned to the same chromosome and start position)
Multiplied_symbols_TheSameChromosome_Different_start_final.csv (list of alternative Ensembl IDs that are assigned to the same chromosome but different start position)

##################################################################################################################################################################################################

Program
R_summary gene synonyms and official symbols.R

Application
The program creates a lists of unique gene synonyms and official symbols and a Venn diagram comparing these two types of symbols. 

Description
The first step of data processing is the creation of new columns containing gene synonyms and official gene symbols written only in lowercase letters to avoid confusion caused by a diverse usage of the uppercase and lowercase letters in gene symbols. Missing synonyms and official symbols are marked with an “NA” symbol to enable removal of such entries, and the lists of unique synonyms and official symbols are exported to csv files. Finally, the lists of genes are used to create the Venn diagram showing the number of gene synonyms, official symbols and the overlap between them. The script used for this analysis requires installation of naniar, tidyr, dplyr and VennDiagram R packages. 

Input file format
Csv file with a default name “mart_export” downloaded from 
https://www.ensembl.org/biomart/martview/  using Ensembl Genes database with genetic datasets for various species. Selected Ensembl filter setting should be “Gene / Gene type / protein_coding” while the attributes settings should be “Gene Synonym”, “Gene name” and “Gene stable ID”. 

Example of input data (data arrangement)

Gene Synonym,Gene name,Gene stable ID
0610040O15Rik,mt-Nd1,ENSMUSG00000064341
ND1,mt-Nd1,ENSMUSG00000064341
protein 1,mt-Nd1,ENSMUSG00000064341

Output data
Unique_gene_synonyms.csv (lists of unique synonyms)
Unique_gene_symbols.csv files (lists of official symbols)
Venn.png (Venn diagram showing the number of gene synonyms, official symbols and the overlap between them)

##################################################################################################################################################################################################

Program
R_ambiguous official symbols.R

Application
Identification of ambiguous symbols that are used both as a current official symbol of one gene and as a synonym for another gene.

Description
The program identifies ambiguous symbols that are used both as a current official symbol of one gene and as a synonym for another gene. The data retrieved from the Ensembl Biomart genome browser are not sufficient for identification of truly ambiguous symbols because in some cases there is the same symbol in the gene synonym and gene name columns suggesting that a symbol that was a synonym in the past became later an official symbol during the process of nomenclature standardization. Such cases are removed with the R script. The comparison between synonyms and official symbols is  done without recognizing the uppercase and lowercase letters. Next, the data restricted to protein-coding genes are grouped by gene synonym to create the final list of ambiguous gene symbols. The program removes also cases when ambiguous symbol refers to non protein-coding genes and such cases are exported to separate file. The final list of ambiguous gene symbols is merged with the data in the All_unique_protein_coding.csv file  that is generated by the script R_protein-coding unique and contains input symbols used to search the Ensembl Biomart database together with Ensembl stable IDs. This step allows for bringing together a complete set of additional identifiers that are necessary for disambiguation of symbols. The final stage of data processing includes editing the names of the columns and providing additional information enabling easy identification of data type and species. The part C of the script also contains the command removing unpaired input (Ambiguous_official_symbols) and output data (tolower_Gene_Synonyms). The lack of congruence between input and output data is a rare error caused by the occurrence of erroneous space inside of the gene symbol. The script is ready for analyzing mouse data while commands specific for rat, human and pig data are silenced with #.
The script used for this analysis requires installation of tidyr, purrr, readr and dplyr R packages.

Input file format
Csv file with a default name “mart_export” downloaded from 
https://www.ensembl.org/biomart/martview/  using Ensembl Genes database with genetic datasets for various species. Selected Ensembl filter setting should be  “Gene / Input external references ID list / Gene Synonym(s)”.  Gene Synonyms for external references list are available in file All_unique_protein_coding_symbols.csv generated by the script R_protein-coding unique.  The attributes settings should be “Gene Synonym, Gene name, Gene type and Gene stable ID”. 

Example of input data (data arrangement)

Gene Synonym,Gene name,Gene type,Gene stable ID
ACAT2,Soat2,protein_coding,ENSMUSG00000023045
ADAM12,Adam12,protein_coding,ENSMUSG00000054555
ADAM3,Adam3,protein_coding,ENSMUSG00000031553

Output data
Ambiguous_symbols_protein_and_other.csv (cases of symbol ambiguity resulting from association with non protein-coding genes)
mart_export_5.csv (list of ambiguous official symbols / synonyms grouped by gene synonym)
Ambiguous_symbols_protein_protein_final.csv file (list of ambiguous official symbols / synonyms after final editing)

##################################################################################################################################################################################################

Program
R_ambiguous synonyms.R

Application
Identification of ambiguous synonyms that can be assigned to at least two different official gene symbols in Ensembl genome data.

Description
The first step of data processing is to remove all genes without synonym symbols. Next, the dataset is grouped by gene synonym to bring together all official symbols linked to each gene synonym. The grouping is done without recognizing the uppercase and lowercase letters in gene synonyms. Next, the program tests whether there is more than one official symbol assigned to each synonym. The final stage of the data processing is performed with the part B of the script to edit the names of the columns and to provide additional information enabling easy identification of the data type and species. The R script also creates an additional column with gene symbols anchored to words “Gene synonym” (for example Gene_synonym_nd1) to prevent unintentional transformations of symbols by Excel. The script is ready for analyzing mouse data while commands specific for rat, human and pig data are silenced with #. It is required to install the naniar, tidyr, purrr, readr and dplyr R packages to run the script. 

Input file format
Csv file with a default name “mart_export” downloaded from 
https://www.ensembl.org/biomart/martview/  using Ensembl Genes database with genetic datasets for various species. Selected Ensembl filter setting should be “Gene / Gene type / protein_coding” while the attributes settings should be “Gene Synonym”, “Gene name” and “Gene stable ID”. 

Example of input data (data arrangement)

Gene Synonym,Gene name,Gene stable ID
0610040O15Rik,mt-Nd1,ENSMUSG00000064341
ND1,mt-Nd1,ENSMUSG00000064341
protein 1,mt-Nd1,ENSMUSG00000064341

Output data

Ambiguous_gene_synonyms.csv file (list of ambiguous synonyms) 
Ambiguous_gene_synonyms_final.csv (list of ambiguous synonyms after final editing) 

##################################################################################################################################################################################################

Programs
R_Ensembl novel genes human.R
R_Ensembl novel genes mouse.R
R_Ensembl novel genes pig.R
R_Ensembl novel genes rat.R

Application
Joining Ensembl data with data from other databases to verify missing Ensembl gene symbols.

Description
The program is adjusted for each species and joins Ensembl data with data from other databases to verify Ensembl genes without assigned symbol. The data from alternative resources are additionally screened for Ensembl stable IDs used as gene symbols and such entries are removed. Finally, the data retrieved from the Ensembl are merged with the data collected from the other sites according to external IDs available in the Ensembl and Ensembl IDs provided in alternative databases. It is required to install the naniar, dplyr, tidyr, purrr and readr R packages to run the script. 

Input file format
1) Csv file with a default name “mart_export” downloaded from 
https://www.ensembl.org/biomart/martview/ using Ensembl Genes database with genetic datasets for various species. Selected Ensembl filter setting should be “Gene / Input external references ID list / Gene stable ID(s)”. The Ensembl gene stable IDs (e.g. ENSMUSG00000020807) for external references list is available in file Genes_without_gene_symbol.csv generated by the script R_protein-coding unique.R.  
The attributes should be “Gene name”, “Gene stable ID”, “NCBI gene (formerly Entrezgene) ID” and one of the IDs provided by the nomenclature committee that is mouse “MGI ID”, rat “RGD ID”, human “HGNC ID” and pig “VGNC ID” (selected in external references). 

Example of input data (data arrangement)

Gene stable ID,Gene name,MGI ID,NCBI gene (formerly Entrezgene) ID
ENSMUSG00000074720,,,
ENSMUSG00000079190,,,100041057
ENSMUSG00000079192,,,

2) MGI	 
Files MRK_List2.rpt and MRK_ENSEMBL.rpt downloaded from  https://www.informatics.jax.org/downloads/reports/index.html. Compressed files should be unzipped before importing to the R studio.                     

3) NCBI data. ncbi_dataset.tsv file downloaded from 	
https://www.ncbi.nlm.nih.gov/datasets/gene/  with selected all columns.

4) RGD rat genes.  File GENES_RAT.txt downloaded from https://download.rgd.mcw.edu/data_release/RAT/

5) HGNC human genes. File results.txt downloaded from https://biomart.genenames.org/ with following attributes: HGNC ID, Status, Approved symbol, Approved name and  Ensembl gene ID. 

6) VGNC data pig genes. File pig_vgnc_gene_set_All.txt downloaded from 
https://vertebrate.genenames.org/download/statistics-and-files/. 
Dada selection:
Filter statistics and download files by species: Pig
Filter statistics and download files by chromosome: All 
Statistics: Total Approved Symbols

Output data
Ensembl_MGI_NCBI_grouped.csv (Mouse Ensembl data integrated with MGI and NCBI data)
Ensembl_RGD_NCBI_grouped.csv (Rat Ensembl data integrated with RGD and NCBI data)
Ensembl_HGNC_NCBI_grouped.csv (Human Ensembl data integrated with HGNC and NCBI data)
Ensembl_VGNC_NCBI_grouped.csv (Pig Ensembl data integrated with VGNC and NCBI data)

##################################################################################################################################################################################################

Program
R_genes linked to ambiguous symbols.R

Application
Extraction of complete list of all official symbols that can be affected by the symbol ambiguity

Description
The program uses result from other analyses to extract a complete list of all official symbols that can be affected by the symbol ambiguity. The script requires installation of 
tidyr, dplyr and readr R packages. 

Input file format
Ambiguous_gene_synonyms_final.csv file generated by the script  R_ambiguous synonyms.R

Ambiguous_symbols_protein_protein_final.csv file generated by the scipt R_ambiguous official symbols.R

Output data
All_genes_linked_to_ambiguous_symbols.csv. (list of all genes linked to ambiguous symbols)

##################################################################################################################################################################################################

Programs
R_test of additional IDs in Ensembl Human.R
R_test of additional IDs in Ensembl Mouse.R
R_test of additional IDs in Ensembl Pig.R
R_test of additional IDs in Ensembl Rat.R

Application
The scripts identify missing additional IDs and IDs assigned to more than one official gene symbol in data retrieved from Ensembl database.

Description
The scripts identify missing additional IDs and IDs assigned to more than one official gene symbol in data retrieved from Ensembl database. It is required to install the naniar, tidyr, dplyr and readr R packages to run the script. 

Input file format
Csv file with a default name “mart_export” downloaded from 
https://www.ensembl.org/biomart/martview/ using Ensembl Genes database with genetic datasets for various species.  Selected Ensembl filter setting should be “Gene / Input external references ID list / Gene Name(s)”. The list of gene symbols (e.g. avpr2) for external references list is available in file All_genes_linked_to_ambiguous_symbols.csv generated by the script R_genes linked to ambiguous symbols.R.  The attributes in Gene/Ensembl and External/External references directories should be “Gene name”, “Gene stable ID”, “NCBI gene (formerly Entrezgene) ID” and one of the IDs provided by the nomenclature committee that is mouse “MGI ID”, rat “RGD ID”, human “HGNC ID” or pig “VGNC ID”. 

Example of input data (data arrangement)

Gene name,Gene stable ID,NCBI gene (formerly Entrezgene) ID,MGI ID
mt-Nd1,ENSMUSG00000064341,17716,MGI:101787
mt-Co1,ENSMUSG00000064351,17708,MGI:102504
Phax,ENSMUSG00000008301,56698,MGI:1891839

Output data
Missing_IDs.csv (list of missing gene IDs in Ensembl data)
Multiplied_Ensembl_IDs_final.csv (list of Ensembl gene IDs assigned to more than one official gene symbol in Ensembl database)
Multiplied_NCBI_IDs_final.csv (list of NCBI gene IDs assigned to more than one official gene symbol in Ensembl database)
Multiplied_Committee_IDs_final.csv (list of nomenclature committee (MGi, RGD, HGNC or VGNC)  gene IDs assigned to more than one official gene symbol in Ensembl database)

##################################################################################################################################################################################################

Program
R_HUGO test.R

Application
1) Checking whether missing  HGNC (Hugo Gene Nomenclature Committee) IDs in Ensembl are available in  Hugo Gene Nomenclature Committee database. 

Description
1) The program restricts Ensembl data to cases with missing HGNC IDs and resulting list merges with genome data retrieved from HGNC database. The datasets are joined based on official gene symbols. The data exported to the csv file enable verification whether missing HGNC IDs in Ensembl are also missing in HGNC database. It is required to install dplyr and readr packages to run the script. 

Input file format
1) Missing_IDs.csv generated by script R_test of additional IDs in Ensembl Human.R
2) File hgnc-search-1734202353314.txt (adjust file name to enable proper data import into R Studio) downloaded from https://www.genenames.org/tools/search/#!/?query=&rows=20&start=0&filter=document_type:gene with filter by type selection: Gene. 

Output data
Verification_missing_HUGO_IDs.csv (list of Ensembl genes with missing  HGNC IDs merged with gene identifiers retrieved directly from HGNC database)

##################################################################################################################################################################################################

Program
R_NCBI test.R

Application
1) Checking whether missing NCBI IDs in Ensembl are available in NCBI database 
2) Checking whether NCBI IDs assigned in Ensembl database to more than one official gene symbol are also assigned to more than one gene in NCBI database. 

Description
1) The program restricts Ensembl data to cases with missing NCBI IDs and resulting list merges with genome data retrieved from NCBI. The datasets are joined based on official gene symbols. The data exported to the csv file enable verification whether missing NCBI IDs in Ensembl are also missing in NCBI database. 
2) The program merges list of duplicated NCBI IDs in Ensembl genome data (the same NCBI ID assigned to more than one official gene symbol in Ensembl) with genome data retrieved directly from NCBI database. The datasets are joined based on NCBI IDs. The data exported to the csv file enable verification of assignment of NCBI IDs to more than 1 gene. 
The script is ready for analyzing mouse data while commands specific for rat, human and pig data are silenced with #.
It is required to install dplyr, tidyr, purrr and readr packages to run the script. 

Input file format
1) Missing_IDs.csv (generated by scripts: 
R_test of additional IDs in Ensembl Human.R
R_test of additional IDs in Ensembl Mouse.R 
R_test of additional IDs in Ensembl Rat.R)
2) Multiplied_NCBI_IDs_final.csv (generated by scripts: 
R_test of additional IDs in Ensembl Human.R
R_test of additional IDs in Ensembl Mouse.R
R_test of additional IDs in Ensembl Rat.R)
3) File ncbi_dataset.tsv downloaded without any filter selection from:
https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_000001635.27/  (mouse reference genome GRCm39) with all selected columns
https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_015227675.2/  (Rattus norvegicus genome mRatBN7.2) with all selected columns
https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_000001405.40/ (human reference genome GRCh38.p14) with all selected columns

Output data
Verification_missing_NCBI_IDs.csv (list of Ensembl genes with missing  NCBI IDs merged with gene identifiers retrieved directly from NCBI database)
Summary_NCBI_IDs_ambiguous_in_Ensembl.csv (list of NCBI IDs assigned in Ensembl database to more than one official gene symbol together with gene identifiers retrieved directly from NCBI database)

##################################################################################################################################################################################################

Program
R_RGD test.R

Application 
Checking whether RGD IDs assigned in Ensembl database to more than one official gene symbol are also assigned to more than one gene in RGD database. 

Description
The program merges list of duplicated RGD IDs in Ensembl genome data (the same RGD ID assigned to more than one official gene symbol in Ensembl) with genome data retrieved directly from RGD database. The datasets are joined based on RGD IDs. The data exported to the csv file enable verification of assignment of RGD IDs to more than 1 gene. It is required to install dplyr, tidyr, purrr and readr packages to run the script. 

Input file format
1) Missing_IDs.csv (generated by script R_test of additional IDs in Ensembl Rat.R)
2) Multiplied_NCBI_IDs_final.csv (generated by script R_test of additional IDs in Ensembl Rat.R)
3) File GENES_RAT.txt  downloaded from https://download.rgd.mcw.edu/data_release/RAT/

Output data
Summary_RGD_IDs_ambiguous_in_Ensembl.csv (list of RGD IDs assigned in Ensembl database to more than one official gene symbol together with gene identifiers retrieved directly from RGD database)


##################################################################################################################################################################################################

Program
R_vgnac test.R

Application
Checking whether missing VGNC IDs in Ensembl are available in VGNC database.

Description
The program restricts Ensembl data to cases with missing VGNC IDs and resulting list merges with genome data retrieved from VGNC database. The datasets are joined based on official gene symbols. The import of VGNC data is followed by correction of column name assignment. The data exported to the csv file enable verification whether missing VGNC IDs in Ensembl are also missing in VGNC database. 
It is required to install dplyr and readr packages to run the script. 

Input file format
1) Missing_IDs.csv (generated by script R_test of additional IDs in Ensembl Pig.R)
2) File pig_vgnc_gene_set_All.txt downloaded from 
https://vertebrate.genenames.org/download/statistics-and-files/. 
Dada selection
Filter statistics and download files by species: Pig
Filter statistics and download files by chromosome: All 
Statistics: Total Approved Symbols

Output data
Verification_missing_VGNC_IDs.csv (list of Ensembl genes with missing  VGNC IDs merged with gene identifiers retrieved directly from VGNC database)

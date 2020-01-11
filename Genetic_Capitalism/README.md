# Genetic Capitalism and Stabilizing Selection of Antimicrobial Resistance Genotypes in _Escherichia coli_
## Supplemental Materials

|                                              Authors                                              |    Paper    |   Citation  | Figures |
|-------------------------------------------------------------------------------------------------|:-----------:|:-----------:|:--------:|
| Colby T. Ford, Ph.D.<br>Kevin Smith<br>Gabriel Lopez Zenarosa, Ph.D.<br>David Brown<br>John Williams<br>Daniel Janies, Ph.D. | [Coming Soon - Cladistics](#) | ```Coming Soon``` | [Tableau Viz](https://public.tableau.com/profile/cford38#!/vizhome/E_coliGenotypeSetsViz/GeneticCapitalism) |

-------------------------------------

|     | File Type        | File Name                                                                                                | Description                                                                                                                                                          |
|-----|------------------|----------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1.  | Script           | [clustering.R](clustering.R)                                                                             | Script for performing cluster analysis of AMR genotype gains and losses.                                                                                             |
| 2.  | Script           | [distribution_tests.R](distribution_tests.R)                                                             | Script for performing distribution tests by AMR group & resistance mechanism.                                                                                        |
| 3.  | Script           | [TNT/TNTFileGenerator/TNT_file_generator.R](TNT/TNTFileGenerator/TNT_file_generator.R)                   | Script for generating the data file for the TNT-based phylogenetic analysis.                                                                                         |
| 4.  | Reference Data   | [TNT/PDG000000004.1024.reference_target.tree.newick](TNT/PDG000000004.1024.reference_target.tree.newick) | Reference Tree from NCBI. (Current PDG trees for _E. coli_ can be found [here](https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Escherichia_coli_Shigella/).)                |
| 5.  | Reference Method | [NCBI_TreeMethods.txt](NCBI_TreeMethods.txt)                                                             | Information of the NCBI processing pipeline for generating the SNP trees. (Can also be found on NCBI's site [here](https://ftp.ncbi.nlm.nih.gov/pathogen/Methods.txt).) |
| 6.  | Reference Data   | [IsolationData.xlsx](IsolationData.xlsx)                                                                 | Roll-Up of Isolation Metadata (Source country, organism, etc.)                                                                                                       |
| 7.  | Reference Data   | [AMR_FunctionalMechanisms.csv](AMR_FunctionalMechanisms.csv)                                             | Functional mechanisms of the 409 AMR genotypes in this study. (From [CARD](https://card.mcmaster.ca/))                                                               |
| 8.  | Results          | [clusters.csv](clusters.csv)                                                                             | Resulting AMR genotype clusters from the cluster analysis of gains and losses.                                                                                       |

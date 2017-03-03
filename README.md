# Targeted sequencing reveals EPHA genes as markers of paclitaxel-induced peripheral neuropathy
***
### Authors: 
Maria Apellániz-Ruiz1, Héctor Tejero2, Lucía Inglada-Perez1,3, Lara Sánchez-Barroso1, Gerardo Gutiérrez-Gutiérrez4, Isabel Calvo5,6, Beatriz Castelo7, Andrés Redondo7, Jesus García-Donas8, Nuria Romero8, Maria Sereno9, María Merino9, Maria Currás1, Cristina Montero-Conde1, Veronika Mancikova1, Elisabeth Åvall-Lundqvist10, Henrik Green11,12, Fátima Al-Shahrour2, Alberto Cascón1,3, Mercedes Robledo1,3, Cristina Rodríguez-Antona1,3
***
## Introduction 

This is the repository of the code used for the analysis carried out in the article __Targeted sequencing reveals EPHA genes as markers of paclitaxel-induced peripheral neuropathy__. 

The objective of the article is to find rare variants associated to paclitaxel-induced neuropathy in patients with ovarian cancer. An extreme phenotype design has been used. 

The study used 228 samples. 196 were processed using the TruSeq Custom Amplicon Kit (Illumina) covering the coding plus 25 bp intronic flanking region of 39 genes described in the paper. The other 32 samples come from exome sequencing. 

## Navigation 

We have divided the code in three folders, each containing a different step of the analysis. 

### Exome analysis

This folder contains the code and templates for the alignment of the exomes using RUbioSeq. It also contains the code in order to get extract the gene panel sequencies from the exomes. 

### Filtering and annotation 

This folder contains the code used to filter and annotate the variants obtained according to quality criterions and the type of mutation. 

### Statistical Analysis using SKAT

This folder contains the R scripts used to carry out the association analysis using SKAT. 





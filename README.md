# Cryptobiota and Climate Change
This is a repository for code and data used in the analysis of the manuscript entitled: "Biodiversity of coral reef cryptobiota shuffles but does not decline under the combined stressors of ocean warming and acidification".

## Fastq Sequences
Fastq files were submited to NCBI through Genome -https://geome-db.org/workbench/project-overview. The NCBI Sequence Read Archive (SRA) numbers for this project are SRS7105074-SRS7105095 and can be obtained here: https://www.ncbi.nlm.nih.gov/sra. 

The demultiplexed fastq files are named ArmsMeso1 thru 22. There are two files, one forward and one reverse for each treatment: ArmsMeso1.1 = forward and ArmsMeso1.2 = reverse. The treatment names that correspond to the ArmsMeso numbers are in the MetaData.csv file located in the Data folder.

Sequences were processed through the JAMP pipeline - https://github.com/VascoElbrecht/JAMP. See manuscript supplement for details.

The JAMP pipeline output is the 3_Raw_MOTU_Table.csv file located in the Data Folder. 

## Data and R Code

All data and R code used to obtain the results are presented in the Data and R Code folders. R packages needed to process the data are provided in the respective heading of each R file. See manuscript supplement for the details on sequence annotations and conducted statistical analyses. All Figures and accompanying Tables presented in the manuscript are highlighted within the R files and are found here:

- Figure 1 - EnvironmentalData.R
- Figure 2 - Diversity.R 
- Figure 3 - Phyla_Family_Composition.R
- Fig. S3 - MOTU_Processing.R
- Fig. S4 and Fig. S5 Diversity.R
- Fig. S6 - Phyla_Family_Composition.R
- Tables S1:S4 - MOTU_Processing
- Tables S5:S10 - Diversity.R
- Tables S11:S14 - Phyla_Family_Composition.R

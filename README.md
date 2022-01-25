# REGENIE Association Analysis Results UKB 450K (Gene Burden & Single Variant)

Created by Saleh Shekari (ss1173@exeter.ac.uk), Amy Dawes (a.dawes@exeter.ac.uk) and Gareth Hawkes (g.hawkes2@exeter.ac.uk)

Please contact [Saleh](mailto:ss1173@exeter.ac.uk?subject=[GitHub]%20SAIGE-GENE%20Association%20Analysis%20Summery) to report bugs or suggest and features. Pull requests welcome and encouraged.

## Contents
- [Introduction](#Introduction)
- [Quick-start](#Quick-start-instructions)




## Introduction
Briefly, when run using the instructions below, this code will achieve the following from UKB-450K-WES REGENIE-GENE trait association results files as input files:
* For each variant class (e.g. LoF, Missense etc.) extract each genes transcript burden across all chromosomes (1-24) from REGENIE-GENE association output.
* Generate summary tables for each variant class and MAF cut-off filter (MAF<0.1%, and MAF<0.01%).
* Generate master tables for each MAF cut-off filter across all variant classed combined.
* Generate top hits tables for each MAF cut-off filter with variant class labelled.
* Generate top hits tables for each variant class with MAF<0.1%.
* Generate directory for masks that test failed with list of gene/transcript/masks which failed.
* Generate a master file for single variant association results.
* Generate top hits table for single variant association results.
* Generate top hits table for rare single variant association results.
* Generate Q-Q plots for each variant class and MAF cut-off filter as .pdf format
* Generate Q-Q plots for your traits GWAS nearset genes
* Generate Manhattan plots for MAF cut-off filter (MAF<0.1%, and MAF<0.01%)

__Before start you need to have a .txt file for your GWAS SNPs detailes (first column "Chr"; second column "Position (b37)").__
__In GWAS SNPs list use X for chromosome 23, this is required for the script to work.__

__Please run this code from one of your own working directories by following the quick start instructions below.__


# Quick-start-instructions

**1. Log into Slade**

**2. Navigate to your working directory**

```
cd /full/path/to/<your_working_directory>
```

**3. Submit script to Slade job nodes**
* Please make sure in this step that you input the full/absolute filepath to your working directory, as this is required for the script to work.

```
nohup bash /slade/projects/UKBB/UKBB_general_working_area/exWAS_script_Sal_Amy_Gareth/exWAS_UKB_450K_REGENIE_summary_script.sh </full/path/to/<directory_you_Want_to_Archive> <your_Interest_Trait> <Date_REGENIE_File_Made> </full/path/to/GWAS_SNPs_List> &
```
* Double Enter

**4. Check job status**
```
jobs
```

**5. Check job submission's overview**
```
less nohup.out 
```

https://github.com/SalShk/REGENIE-GENE_Association_Analysis_Results_UKB_450K/issues/1#issue-1114191185

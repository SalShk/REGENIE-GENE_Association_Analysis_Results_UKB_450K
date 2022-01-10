## Introduction
Briefly, when run using the instructions below, this code will achieve the following from UKB-450K-WES REGENIE-GENE trait association results files as input files:
* For each variant class (e.g. LoF, Missense etc.) extract each genes transcript burden across all chromosomes (1-24) from REGENIE-GENE association output.
* Generate summary tables for each variant class and MAF cut-off filter (MAF<0.1%, and MAF<0.01%)
* Generate master tables for each MAF cut-off filter across all variant classed combined.
* Generate top hits tables for each MAF cut-off filter with variant class labelled
* Generate directory for masks that test failed with list of gene/transcript/masks which failed
* Generate Q-Q plots for each variant class and MAF cut-off filter as .pdf format
* Generate Q-Q plots for your traits GWAS nearset genes
* Generate Manhattan plots for MAF cut-off filter (MAF<0.1%, and MAF<0.01%)

__Before start you need to have a .txt file for your GWAS SNPs detailes (first column "Chr"; second column "Position (b37)").__

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
```

**4. Check job status**
```
squeue -u <user_ID>
```

**5. Check job submission's overview**
```
less slurm-<job_ID>.out 
```

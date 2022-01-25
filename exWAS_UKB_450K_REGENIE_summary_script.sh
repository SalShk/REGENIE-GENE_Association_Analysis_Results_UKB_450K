
#!/bin/bash


DIRECTORY=${1}
TRAIT=${2}
Date=${3}
gwasfile=${4}

 cd ${DIRECTORY}

# If it doesn't exist already, a new directory will be created:

 mkdir ${TRAIT}_regenie_UKB_450K

 cd ${TRAIT}_regenie_UKB_450K



#### STEP 1. COMBINE ALL CHROMOSOMES FILES FOR EACH VARIANT CLASS AND MAF FILTER (0.0001 AND 0.001) -  There will be 2 files produced (one for MAF<0.01% and one for MAF<0.1%) for each variant subclass:
#(LoF, LoF_SNVs, Missense, Missenes CADD>30 ,LoF+Missenes CADD>30, Synonymous)


  cat /slade/projects/UKBB/DNA_Nexus/${TRAIT}_regenie_burden_dnanexus_${Date}/${TRAIT}_Step2_Chr*_${TRAIT}.regenie > ${TRAIT}_transcript_gene_master_file_450k.txt

  grep ".0.001\|singleton" ${TRAIT}_transcript_gene_master_file_450k.txt > ${TRAIT}_transcript_gene.0.001_450k.txt

  grep ".0.0001\|singleton" ${TRAIT}_transcript_gene_master_file_450k.txt > ${TRAIT}_transcript_gene.0.0001_450k.txt


  sed -e 's/ /\t/g' ${TRAIT}_transcript_gene.0.001_450k.txt > ${TRAIT}_transcript_gene.0.001_450k1.txt
  rm ${TRAIT}_transcript_gene.0.001_450k.txt
  echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene.0.001_450k1.txt > ${TRAIT}_transcripts_genes.0.001_450k.txt
  grep "TEST_FAIL" ${TRAIT}_transcripts_genes.0.001_450k.txt > ${TRAIT}_TEST_FAIL_transcripts_genes.0.001_450k.txt
  sed '/TEST_FAIL/d' -i ${TRAIT}_transcripts_genes.0.001_450k.txt
  awk '{if($12>5) {print}}' ${TRAIT}_transcripts_genes.0.001_450k.txt > ${TRAIT}_transcripts_genes.0.001_450k_tophits.txt
  rm ${TRAIT}_transcript_gene.0.001_450k1.txt
  grep synonymous ${TRAIT}_transcripts_genes.0.001_450k.txt > ${TRAIT}_transcripts_genes.0.001_synonymous_450k.txt
  sed '/synonymous/d' -i ${TRAIT}_transcripts_genes.0.001_450k.txt

  sed -e 's/ /\t/g' ${TRAIT}_transcript_gene.0.0001_450k.txt > ${TRAIT}_transcript_gene.0.0001_450k1.txt
  rm ${TRAIT}_transcript_gene.0.0001_450k.txt
  echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene.0.0001_450k1.txt > ${TRAIT}_transcripts_genes.0.0001_450k.txt
  grep "TEST_FAIL" ${TRAIT}_transcripts_genes.0.0001_450k.txt > ${TRAIT}_TEST_FAIL_transcripts_genes.0.0001_450k.txt
  sed '/TEST_FAIL/d' -i ${TRAIT}_transcripts_genes.0.0001_450k.txt
  awk '{if($12>5) {print}}' ${TRAIT}_transcripts_genes.0.0001_450k.txt > ${TRAIT}_transcripts_genes.0.0001_450k_tophits.txt
  rm ${TRAIT}_transcript_gene.0.0001_450k1.txt
  grep synonymous ${TRAIT}_transcripts_genes.0.0001_450k.txt > ${TRAIT}_transcripts_genes.0.0001_synonymous_450k.txt
  sed '/synonymous/d' -i ${TRAIT}_transcripts_genes.0.0001_450k.txt


  echo "All masks for each MAF cut-off filed successfully"

#######################============================================================ More classification if you want (more fun) ==========================================================================##########
##                                                                                                                                                                                                               ##
##===============================================================================================================================================================================================================##



  grep "LoF_HC.singleton" ${TRAIT}_transcript_gene_master_file_450k.txt > ${TRAIT}_transcript_gene_LoF_HC.singleton_450k.txt
  sed -e 's/ /\t/g' ${TRAIT}_transcript_gene_LoF_HC.singleton_450k.txt > ${TRAIT}_transcript_gene_LoF_HC.singleton_450k1.txt
  rm ${TRAIT}_transcript_gene_LoF_HC.singleton_450k.txt
  echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_LoF_HC.singleton_450k1.txt > ${TRAIT}_transcript_gene_LoF_HC.singleton_450k.txt
  sed '/TEST_FAIL/d' -i ${TRAIT}_transcript_gene_LoF_HC.singleton_450k.txt
  rm ${TRAIT}_transcript_gene_LoF_HC.singleton_450k1.txt
  awk '{if($12>5) {print}}' ${TRAIT}_transcript_gene_LoF_HC.singleton_450k.txt > ${TRAIT}_transcript_gene_LoF_HC.singleton_450k_tophits.txt
  #echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_LoF_HC.singleton_450k1_tophits.txt > ${TRAIT}_transcript_gene_LoF_HC.singleton_450k_tophits.txt
  #rm ${TRAIT}_transcript_gene_LoF_HC.singleton_450k1_tophits.txt

  echo "LoF_HC.singleton filed successfully"


  grep "LoF_HC.0.001" ${TRAIT}_transcript_gene_master_file_450k.txt > ${TRAIT}_transcript_gene_LoF_HC.0.001_450k.txt
  sed -e 's/ /\t/g' ${TRAIT}_transcript_gene_LoF_HC.0.001_450k.txt >  ${TRAIT}_transcript_gene_LoF_HC.0.001_450k1.txt
  rm  ${TRAIT}_transcript_gene_LoF_HC.0.001_450k.txt
  echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_LoF_HC.0.001_450k1.txt > ${TRAIT}_transcript_gene_LoF_HC.0.001_450k.txt
  rm ${TRAIT}_transcript_gene_LoF_HC.0.001_450k1.txt
  sed '/TEST_FAIL/d' -i ${TRAIT}_transcript_gene_LoF_HC.0.001_450k.txt
  awk '{if($12>5) {print}}' ${TRAIT}_transcript_gene_LoF_HC.0.001_450k.txt > ${TRAIT}_transcript_gene_LoF_HC.0.001_450k_tophits.txt
  #echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_LoF_HC.0.001_450k1_tophits.txt > ${TRAIT}_transcript_gene_LoF_HC.0.001_450k_tophits.txt
  #rm  ${TRAIT}_transcript_gene_LoF_HC.0.001_450k1_tophits.txt

  echo "LoF_HC.0.001 filed successfully"

  grep "LoF_SNV.singleton" ${TRAIT}_transcript_gene_master_file_450k.txt > ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k.txt
  sed -e 's/ /\t/g' ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k.txt >  ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k1.txt
  rm  ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k.txt
  echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k1.txt > ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k.txt
  rm ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k1.txt
  sed '/TEST_FAIL/d' -i ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k.txt
  awk '{if($12>5) {print}}' ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k.txt > ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k_tophits.txt
  #echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k1_tophits.txt > ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k_tophits.txt
  #rm  ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k1_tophits.txt

  echo "LoF_SNV.singleton filed successfully"

  grep "LoF_SNV.0.001" ${TRAIT}_transcript_gene_master_file_450k.txt > ${TRAIT}_transcript_gene_LoF_SNV.0.001_450k.txt
  sed -e 's/ /\t/g' ${TRAIT}_transcript_gene_LoF_SNV.0.001_450k.txt > ${TRAIT}_transcript_gene_LoF_SNV.0.001_450k1.txt
  rm  ${TRAIT}_transcript_gene_LoF_SNV.0.001_450k.txt
  echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_LoF_SNV.0.001_450k1.txt > ${TRAIT}_transcript_gene_LoF_SNV.0.001_450k.txt
  rm ${TRAIT}_transcript_gene_LoF_SNV.0.001_450k1.txt
  sed '/TEST_FAIL/d' -i ${TRAIT}_transcript_gene_LoF_SNV.0.001_450k.txt
  awk '{if($12>5) {print}}' ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k.txt > ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k_tophits.txt
  #echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k_tophits1.txt > ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k_tophits.txt
  #rm  ${TRAIT}_transcript_gene_LoF_SNV.singleton_450k_tophits1.txt

  echo "LoF_SNV.0.001 filed successfully"

  grep "missense.singleton" ${TRAIT}_transcript_gene_master_file_450k.txt > ${TRAIT}_transcript_gene_missense.singleton_450k.txt
  sed -e 's/ /\t/g' ${TRAIT}_transcript_gene_missense.singleton_450k.txt > ${TRAIT}_transcript_gene_missense.singleton_450k1.txt
  rm  ${TRAIT}_transcript_gene_missense.singleton_450k.txt
  echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_missense.singleton_450k1.txt > ${TRAIT}_transcript_gene_missense.singleton_450k.txt
  rm ${TRAIT}_transcript_gene_missense.singleton_450k1.txt
  sed '/TEST_FAIL/d' -i ${TRAIT}_transcript_gene_missense.singleton_450k.txt
  awk '{if($12>5) {print}}' ${TRAIT}_transcript_gene_missense.singleton_450k.txt > ${TRAIT}_transcript_gene_missense.singleton_450k_tophits.txt
  #echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_missense.singleton_450k_tophits1.txt > ${TRAIT}_transcript_gene_missense.singleton_450k_tophits.txt
  #rm  ${TRAIT}_transcript_gene_missense.singleton_450k_tophits1.txt

  echo "missense.singleton filed successfully"

  grep "missense.0.001" ${TRAIT}_transcript_gene_master_file_450k.txt > ${TRAIT}_transcript_gene_missense.0.001_450k.txt
  sed -e 's/ /\t/g' ${TRAIT}_transcript_gene_missense.0.001_450k.txt > ${TRAIT}_transcript_gene_missense.0.001_450k1.txt
  rm  ${TRAIT}_transcript_gene_missense.0.001_450k.txt
  echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_missense.0.001_450k1.txt > ${TRAIT}_transcript_gene_missense.0.001_450k.txt
  rm ${TRAIT}_transcript_gene_missense.0.001_450k1.txt
  sed '/TEST_FAIL/d' -i ${TRAIT}_transcript_gene_missense.0.001_450k.txt
  awk '{if($12>5) {print}}' ${TRAIT}_transcript_gene_missense.0.001_450k.txt > ${TRAIT}_transcript_gene_missense.0.001_450k_tophits.txt
  #echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_missense.0.001_450k_tophits1.txt > ${TRAIT}_transcript_gene_missense.0.001_450k_tophits.txt
  #rm  ${TRAIT}_transcript_gene_missense.0.001_450k_tophits1.txt

  echo "missense 0.001 filed successfully"

  grep "missense_cadd30.singleton" ${TRAIT}_transcript_gene_master_file_450k.txt > ${TRAIT}_transcript_gene_missense_cadd30.singleton_450k.txt
  sed -e 's/ /\t/g' ${TRAIT}_transcript_gene_missense_cadd30.singleton_450k.txt > ${TRAIT}_transcript_gene_missense_cadd30.singleton_450k1.txt
  rm  ${TRAIT}_transcript_gene_missense_cadd30.singleton_450k.txt
  echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_missense_cadd30.singleton_450k1.txt > ${TRAIT}_transcript_gene_missense_cadd30.singleton_450k.txt
  rm ${TRAIT}_transcript_gene_missense_cadd30.singleton_450k1.txt
  sed '/TEST_FAIL/d' -i ${TRAIT}_transcript_gene_missense_cadd30.singleton_450k.txt
  awk '{if($12>5) {print}}' ${TRAIT}_transcript_gene_missense_cadd30.singleton_450k.txt > ${TRAIT}_transcript_gene_missense_cadd30.singleton_450k_tophits.txt
  #echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_missense_cadd30.singleton_450k_tophits1.txt > ${TRAIT}_transcript_gene_missense_cadd30.singleton_450k_tophits.txt
  #rm  ${TRAIT}_transcript_gene_missense_cadd30.singleton_450k_tophits1.txt

  echo "missense_cadd30.singleton filed successfully"

  grep "missense_cadd30.0.001" ${TRAIT}_transcript_gene_master_file_450k.txt > ${TRAIT}_transcript_gene_missense_cadd30.0.001_450k.txt
  sed -e 's/ /\t/g' ${TRAIT}_transcript_gene_missense_cadd30.0.001_450k.txt > ${TRAIT}_transcript_gene_missense_cadd30.0.001_450k1.txt
  rm  ${TRAIT}_transcript_gene_missense_cadd30.0.001_450k.txt
  echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_missense_cadd30.0.001_450k1.txt > ${TRAIT}_transcript_gene_missense_cadd30.0.001_450k.txt
  rm ${TRAIT}_transcript_gene_missense_cadd30.0.001_450k1.txt
  sed '/TEST_FAIL/d' -i ${TRAIT}_transcript_gene_missense_cadd30.0.001_450k.txt
  awk '{if($12>5) {print}}' ${TRAIT}_transcript_gene_missense_cadd30.0.001_450k.txt > ${TRAIT}_transcript_gene_missense_cadd30.0.001_450k_tophits.txt
  #echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_missense_cadd30.0.001_450k_tophits1.txt > ${TRAIT}_transcript_gene_missense_0.001_450k_tophits1.txt
  #rm  ${TRAIT}_transcript_gene_missense_cadd30.0.001_450k_tophits1.txt

  echo "missense_cadd30.0.001 filed successfully"


  grep "lof_or_missensecadd30.singleton" ${TRAIT}_transcript_gene_master_file_450k.txt > ${TRAIT}_transcript_gene_lof_or_missensecadd30.singleton_450k.txt
  sed -e 's/ /\t/g' ${TRAIT}_transcript_gene_lof_or_missensecadd30.singleton_450k.txt > ${TRAIT}_transcript_gene_lof_or_missensecadd30.singleton_450k1.txt
  rm  ${TRAIT}_transcript_gene_lof_or_missensecadd30.singleton_450k.txt
  echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_lof_or_missensecadd30.singleton_450k1.txt > ${TRAIT}_transcript_gene_lof_or_missensecadd30.singleton_450k.txt
  rm ${TRAIT}_transcript_gene_lof_or_missensecadd30.singleton_450k1.txt
  sed '/TEST_FAIL/d' -i ${TRAIT}_transcript_gene_lof_or_missensecadd30.singleton_450k.txt
  awk '{if($12>5) {print}}' ${TRAIT}_transcript_gene_lof_or_missensecadd30.singleton_450k.txt > ${TRAIT}_transcript_gene_lof_or_missensecadd30.singleton_450k_tophits.txt
  #echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_lof_or_missensecadd30.singleton_450k_tophits1.txt > ${TRAIT}_transcript_gene_lof_or_missensecadd30.singleton_450k_tophits.txt
  #rm  ${TRAIT}_transcript_gene_lof_or_missensecadd30.singleton_450k_tophits1.txt

  echo "lof_or_missensecadd30.singleton filed successfully"


  grep "lof_or_missensecadd30.0.001" ${TRAIT}_transcript_gene_master_file_450k.txt > ${TRAIT}_transcript_gene_lof_or_missensecadd30.0.001_450k.txt
  sed -e 's/ /\t/g' ${TRAIT}_transcript_gene_lof_or_missensecadd30.0.001_450k.txt > ${TRAIT}_transcript_gene_lof_or_missensecadd30.0.001_450k1.txt
  rm  ${TRAIT}_transcript_gene_lof_or_missensecadd30.0.001_450k.txt
  echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_lof_or_missensecadd30.0.001_450k1.txt > ${TRAIT}_transcript_gene_lof_or_missensecadd30.0.001_450k.txt
  rm ${TRAIT}_transcript_gene_lof_or_missensecadd30.0.001_450k1.txt
  sed '/TEST_FAIL/d' -i ${TRAIT}_transcript_gene_lof_or_missensecadd30.0.001_450k.txt
  awk '{if($12>5) {print}}' ${TRAIT}_transcript_gene_lof_or_missensecadd30.0.001_450k.txt > ${TRAIT}_transcript_gene_lof_or_missensecadd30.0.001_450k_tophits.txt
  #echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_lof_or_missensecadd30.0.001_450k_tophits1.txt > ${TRAIT}_transcript_gene_lof_or_missensecadd30.0.001_450k_tophits.txt
  #rm  ${TRAIT}_transcript_gene_lof_or_missensecadd30.0.001_450k_tophits1.txt

  echo "lof_or_missensecadd30.0.001 filed successfully"

  grep "synonymous.0.001" ${TRAIT}_transcript_gene_master_file_450k.txt > ${TRAIT}_transcript_gene_synonymous.0.001_450k.txt
  sed -e 's/ /\t/g' ${TRAIT}_transcript_gene_synonymous.0.001_450k.txt > ${TRAIT}_transcript_gene_synonymous.0.001_450k1.txt
  rm  ${TRAIT}_transcript_gene_synonymous.0.001_450k.txt
  echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_synonymous.0.001_450k1.txt > ${TRAIT}_transcript_gene_synonymous.0.001_450k.txt
  rm ${TRAIT}_transcript_gene_synonymous.0.001_450k1.txt
  sed '/TEST_FAIL/d' -i ${TRAIT}_transcript_gene_synonymous.0.001_450k.txt
  awk '{if($12>5) {print}}' ${TRAIT}_transcript_gene_synonymous.0.001_450k.txt > ${TRAIT}_transcript_gene_synonymous.0.001_450k_tophits.txt
  #echo -e "CHROM\tPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_transcript_gene_synonymous.0.001_450k_tophits1.txt > ${TRAIT}_transcript_gene_synonymous.0.001_450k_tophits.txt
  #rm  ${TRAIT}_transcript_gene_synonymous.0.001_450k_tophits1.txt

  echo "synonymous.0.001 filed successfully"

  echo "All masks for each MAF cut-off filed successfully"

#============================================================================ end of working progress section ==============================================================================================||



#===========================================================================================================================================================================================================||
#=================== Created by Saleh Shekari Jan 2022. This is a script to generate Single Variant files                           ========================================================================||
#===========================================================================================================================================================================================================||


  source /slade/local/UKBB/UKBB_general_working_area/exWAS_script_Sal_Amy_Gareth/Single_variants_archive_script.sh

  echo -e "CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" | cat - ${TRAIT}_Single_Variant_master_file_450k1.txt > ${TRAIT}_Single_Variant_master_file_450k.txt
  rm ${TRAIT}_Single_Variant_master_file_450k1.txt
  grep "TEST_FAIL" ${TRAIT}_Single_Variant_master_file_450k.txt > ${TRAIT}_TEST_FAIL_Single_Variant_master_file_450k.txt
  sed '/TEST_FAIL/d' -i ${TRAIT}_Single_Variant_master_file_450k.txt
  awk '{if($12>7) {print}}' ${TRAIT}_Single_Variant_master_file_450k.txt > ${TRAIT}_Single_Variant_450k_tophits.txt
  awk -F "\t" '{if(($6<0.01) && ($12>7)) {print}}' ${TRAIT}_Single_Variant_master_file_450k.txt > ${TRAIT}_Single_Variant_450k_tophits_rare.txt


  cut -f3 ${TRAIT}_Single_Variant_450k_tophits_rare.txt > ${TRAIT}_Single_Variant_rare_tophits_ids.txt


  mkdir ${TRAIT}_transcript_gene_450k_tophits
  mv *450k_tophits.txt ${TRAIT}_transcript_gene_450k_tophits.txt ${TRAIT}_Single_Variant_450k_tophits_rare.txt
  mkdir ${TRAIT}_TEST_FAIL
  mv ${TRAIT}_TEST_FAIL_* ${TRAIT}_TEST_FAIL

  echo "All tophits & TEST_FAIL filed successfully"



#============================================================================ end of working progress section ============================================================================#



#=======================================================================================================================================================================================================#
#===================         Created by Gareth  Jan 2022. This is a script to generate high pli and GWAS signal plotting file===========================================================================#
#=======================================================================================================================================================================================================#

module load SciPy-bundle/2019.03-foss-2019a
module load matplotlib/3.2.1-foss-2020a-Python-3.8.2

python "/slade/local/UKBB/UKBB_general_working_area/exWAS_script_Sal_Amy_Gareth/pli_gene_script_plotting_v3.py" --trait ${TRAIT} --date ${Date} --gwasfile ${gwasfile}

echo "Q-Q plots and nearest genes files generate successfully"


#============================================================================ end of working progress section ===========================================================================================#



#=======================================================================================================================================================================================================||
#============== Created by Amy Jan 2022.In this section we will genarate manhattan plots for minor allele frequency (MAF) cut of of 0.0001 and 0.001 ===================================================||
#=======================================================================================================================================================================================================||

####  LOAD R AND EXECUTE RSCRIPT(S)

module purge;module load R-bundle-Packages;
echo "R loaded successfully "
Rscript /slade/local/UKBB/UKBB_general_working_area/exWAS_script_Sal_Amy_Gareth/manhattan_plot.R ${DIRECTORY} ${TRAIT}

  rm test.eps


echo "Manhattans plots generate successfully"

#============================================================================ end of working progress section ============================================================================#





echo "Done"

#!/bin/bash

touch chr1_22.txt

     for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
     do
     grep -v "CHROM" /slade/projects/UKBB/DNA_Nexus/${TRAIT}_regenie_burden_dnanexus_${Date}/Single_Variant_${TRAIT}_Step2_Chr${i}_${TRAIT}.regenie >> chr1_22.txt
     done

     sed -i 's/ /\t/g' chr1_22.txt
     awk '{print $1,$2,$3,$4,$5,$6,$8,$9,$10,$11,$12,$13,$14}' chr1_22.txt > chr_1_22.txt
     rm chr1_22.txt
     sed -i 's/ /\t/g' chr_1_22.txt

touch chrX_Y.txt

      for i in X Y
      do
      grep -v "CHROM" /slade/projects/UKBB/DNA_Nexus/${TRAIT}_regenie_burden_dnanexus_${Date}/Single_Variant_${TRAIT}_Step2_Chr${i}_${TRAIT}.regenie >> chrX_Y.txt
      done

     sed -i 's/ /\t/g' chrX_Y.txt
     cat chr*.txt >> ${TRAIT}_Single_Variant_master_file_450k1.txt
     rm chrX_Y.txt chr_1_22.txt

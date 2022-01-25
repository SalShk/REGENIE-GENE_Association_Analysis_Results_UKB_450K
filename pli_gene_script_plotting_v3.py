import pandas as pd
#from qmplot import qqplot
import numpy as np
import matplotlib.pyplot as plt
import csv
import math 
from scipy.stats import beta
import argparse
import sys

parser = argparse.ArgumentParser(description='A script for collapsing P-values from Burden Tests on SAIGE-GENE output. Email Gareth Hawkes (g.hawkes2@exeter.ac.uk) to report a bug')

parser.add_argument('--trait',  dest='trait', required=True, help="Name of phenotype variable used")
parser.add_argument('--date',    dest='date',   required=True, help="Tab delimited input file of genes and chromosomes")
parser.add_argument('--gwasfile',    dest='gwas',   required=True, help="Tab delimited input file of genes and chromosomes")
args = parser.parse_args()


exwaspath = "/slade/local/UKBB/DNA_Nexus/" + args.trait + "_regenie_burden_dnanexus_" + args.date + "/" + args.trait  +"_Step2_Chr"

data_genes =  pd.read_csv("/slade/local/UKBB/DNA_Nexus/qq_plot_filter_script/gnomad.v2.1.1.lof_metrics.by_transcript.txt", header=0, sep ="\t")
#data_genes = data_genes[data_genes["canonical"]==True]
# data_genes = data_genes.drop_duplicates('gene').sort_index()
# data_genes = data_genes.sort_values('pLI', ascending=False).drop_duplicates('gene').sort_index()
data_genes["ID"] = data_genes["gene"] +"_" +  data_genes["transcript"]
data_genes["GENE"] = data_genes["ID"].str.split("_").str[0]
# data_gwas = pd.read_csv('Copy of BW_SNPs_for_gareth.txt', header=0, sep ="\t")
# data_gwas = pd.read_csv('bmi_locke_snps.txt', header=0, sep ="\t")
data_gwas = pd.read_csv(args.gwas, header=0, sep ="\t")
data_pli = data_genes[["ID", "pLI"]]

genes_tmp_300kb =[]
betas_tmp_300kb =[]
se_tmp_300kb =[]
p_tmp_300kb =[]
mask_tmp_300kb = []
pli_tmp_300kb = []
transcript_tmp_300kb = []


genes_tmp_near =[]
betas_tmp_near =[]
se_tmp_near =[]
p_tmp_near=[]
mask_tmp_near = []
pli_tmp_near = []
transcript_tmp_near = []

# data_exwas_master = pd.read_csv('./ldl_corr_raw_sin_regenie_burden/regenie_out_chr' + str(1) +  '_ldl_corr_raw_sin.regenie', header=1, sep =" ")
data_exwas_master = pd.read_csv(exwaspath + "1_" + args.trait + ".regenie", header=1, sep =" ")
data_exwas_master['ALLELE1'] =data_exwas_master['ALLELE1'].astype("string")
# data_exwas_master = data_exwas_master.loc[(data_exwas_master["ALLELE1"] == "LoF.0.001") | (data_exwas_master["ALLELE1"] == "LoF.0.0001") | (data_exwas_master["ALLELE1"] == "missense.0.001") |  (data_exwas_master["ALLELE1"] == "missense.0.0001")  ]
data_exwas_master = data_exwas_master.loc[(data_exwas_master["ALLELE1"] == "LoF_HC.0.001") | (data_exwas_master["ALLELE1"] == "LoF_HC.0.0001") | (data_exwas_master["ALLELE1"] == "missense.0.001") |  (data_exwas_master["ALLELE1"] == "missense.0.0001") |   (data_exwas_master["ALLELE1"] == "lof_or_missensecadd30.0.001")  |   (data_exwas_master["ALLELE1"] == "lof_or_missensecadd30.0.0001") |   (data_exwas_master["ALLELE1"] == "synonymous.0.0001") |   (data_exwas_master["ALLELE1"] == "synonymous.0.001")]
data_exwas_master["GENE_ID"] = data_exwas_master["ID"].str.split(".").str[0]
data_exwas_master["GENE_NAME"] = data_exwas_master["GENE_ID"].str.split("_").str[0]
data_exwas_master = data_exwas_master.merge(data_pli, how = "left", left_on = "GENE_ID", right_on="ID")
data_exwas_master.dropna(subset = ["ID_y"], inplace = True)
for chr, chr_data in data_gwas.groupby("Chr"):
    # print(chr)
    # if str(chr)!="16":
    #     continue
    # if float(chr) > 1:
    #     break
    data_genes_chr = data_genes.loc[data_genes["chromosome"]==str(chr)]
    # data_exwas = pd.read_csv('./ldl_corr_raw_sin_regenie_burden/regenie_out_chr' + str(chr) +  '_ldl_corr_raw_sin.regenie', header=1, sep =" ")
    data_exwas = pd.read_csv(exwaspath + str(chr) + '_' + args.trait + ".regenie", header=1, sep =" ")
    
    
    data_exwas['ALLELE1'] =data_exwas['ALLELE1'].astype("string")
    
    # data_exwas = data_exwas.loc[(data_exwas["ALLELE1"] == "LoF.0.001") | (data_exwas["ALLELE1"] == "LoF.0.0001") | (data_exwas["ALLELE1"] == "missense.0.001") |  (data_exwas["ALLELE1"] == "missense.0.0001")  ]
    
    data_exwas = data_exwas.loc[(data_exwas["ALLELE1"] == "LoF_HC.0.001") | (data_exwas["ALLELE1"] == "LoF_HC.0.0001") | (data_exwas["ALLELE1"] == "missense.0.001") |  (data_exwas["ALLELE1"] == "missense.0.0001") |   (data_exwas["ALLELE1"] == "lof_or_missensecadd30.0.001")  |   (data_exwas["ALLELE1"] == "lof_or_missensecadd30.0.0001") |   (data_exwas["ALLELE1"] == "synonymous.0.0001") |   (data_exwas["ALLELE1"] == "synonymous.0.001") ]

    data_exwas["GENE_ID"] = data_exwas["ID"].str.split(".").str[0]
    data_exwas["GENE_ID"] = data_exwas["GENE_ID"].astype(str)
    data_exwas["GENE_NAME"] = data_exwas["GENE_ID"].str.split("_").str[0]
    #concat_tmp = [data_exwas_master, data_exwas]
    

    data_exwas = data_exwas.merge(data_pli, how = "left", left_on = "GENE_ID", right_on="ID")
    data_exwas.dropna(subset = ["ID_y"], inplace = True)
    concat_tmp = [data_exwas_master, data_exwas]
    data_exwas_master =pd.concat(concat_tmp, ignore_index=False)
    
    for index, row in chr_data.iterrows():
       data_genes_chr["distup"] = abs(data_genes_chr["end_position"] - row["Position (b37)"])
       data_genes_chr["distdwn"]  =abs (data_genes_chr["start_position"] - row["Position (b37)"])
       # if row["Position (b37)"]!=25928179:
            # continue
       data_genes_tmp = data_genes_chr[(abs(data_genes_chr["start_position"] - row["Position (b37)"])  < 300000) | (abs(data_genes_chr["end_position"] - row["Position (b37)"])  < 300000) | ((data_genes_chr["start_position"] < row["Position (b37)"]) & (data_genes_chr["end_position"] > row["Position (b37)"]))]
       data_genes_tmp["distup"] = abs(data_genes_tmp["end_position"] - row["Position (b37)"])
       data_genes_tmp["distdwn"]  =abs (data_genes_tmp["start_position"] - row["Position (b37)"])
       data_genes_nearest_down = data_genes_chr[(abs(data_genes_chr["start_position"] - row["Position (b37)"])) == min(abs(data_genes_chr["start_position"] - row["Position (b37)"]))]
       data_genes_nearest_up = data_genes_chr[(abs(data_genes_chr["end_position"] - row["Position (b37)"])) == min(abs(data_genes_chr["end_position"] - row["Position (b37)"]))]
       data_genes_nearest_inside = data_genes_chr[((data_genes_chr["start_position"] < row["Position (b37)"]) & (data_genes_chr["end_position"] > row["Position (b37)"]))]
       if len(data_genes_nearest_inside) ==1:
           data_exwas_near = data_exwas.loc[data_exwas["GENE_ID"] == data_genes_nearest_inside["ID"].tolist()[0]]
           
       elif np.min(data_genes_nearest_up["distup"]) < np.min(data_genes_nearest_up["distup"]):
            data_exwas_near = data_exwas[data_exwas["GENE_ID"] ==data_genes_nearest_up["ID"].tolist()[0]]
       else:
           data_exwas_near = data_exwas.loc[data_exwas["GENE_ID"] ==data_genes_nearest_down["ID"].tolist()[0]]
       for mask, mask_data in data_exwas_near.groupby("ALLELE1"):
            for index, row2 in mask_data.iterrows():
                genes_tmp_near.append(row2["GENE_NAME"])
                betas_tmp_near.append(row2["BETA"])
                se_tmp_near.append(row2["SE"])
                transcript_tmp_near.append(row2["ID_x"])
                p_tmp_near.append("{:e}".format((10**-float(row2["LOG10P"]))))
                mask_tmp_near.append(mask)
                pli_tmp_near.append(row2["pLI"])
       
       #print()
       
       data_exwas_tmp = data_exwas[data_exwas["GENE_ID"].isin(data_genes_tmp["ID"])]
       
       for mask, mask_data in data_exwas_tmp.groupby("ALLELE1"):
           for index, row2 in mask_data.iterrows():
               genes_tmp_300kb.append(row2["GENE_NAME"])
               betas_tmp_300kb.append(row2["BETA"])
               transcript_tmp_300kb.append(row2["ID_x"])
               se_tmp_300kb.append(row2["SE"])
               p_tmp_300kb.append("{:e}".format((10**-float(row2["LOG10P"]))))
               mask_tmp_300kb.append(mask)
               pli_tmp_300kb.append(row2["pLI"])
               
               
               
bp_300kbp_frame = pd.DataFrame(data = [genes_tmp_300kb,transcript_tmp_300kb,  betas_tmp_300kb, se_tmp_300kb, p_tmp_300kb, mask_tmp_300kb, pli_tmp_300kb]).transpose()
bp_300kbp_frame.drop_duplicates(inplace = True)
bp_300kbp_frame.columns = ["GENE","TRANSCRIPT", "BETA", "SE", "P", "MASK","pLI"]
              
near_frame = pd.DataFrame(data = [genes_tmp_near,transcript_tmp_near, betas_tmp_near, se_tmp_near,p_tmp_near,mask_tmp_near,pli_tmp_near]).transpose()
near_frame.drop_duplicates(inplace = True)
near_frame.columns = ["GENE","TRANSCRIPT", "BETA", "SE", "P", "MASK","pLI"]

for mask, mask_data in near_frame.groupby("MASK"):

    bp_300_data = bp_300kbp_frame[bp_300kbp_frame["MASK"]==mask]
    
    
    fig, (ax1, ax2, ax3,ax4,ax5) = plt.subplots(1, 5, figsize=(36, 16), dpi=80)
    
    data_exwas_tmp = data_exwas_master.loc[(data_exwas_master["ALLELE1"] == mask)]
    data_exwas_tmp["P"] = 10**-data_exwas_tmp["LOG10P"]
    data_exwas_tmp.dropna(subset = ["GENE_NAME"], inplace = True)
    data_exwas_tmp.reset_index(inplace = True)
    
    mask_data["P"] = mask_data["P"].astype(str).astype(float)
    bp_300_data["P"] = bp_300_data["P"].astype(str).astype(float)
    
    data_exwas_tmp.sort_values(by = "P", inplace = True)
    data_exwas_tmp.reset_index(inplace=True)
    index = []
    # index2
    upperlist =-np.log10(beta.ppf(0.975, [i for i in range(1, len(data_exwas_tmp)+1)], [(len(data_exwas_tmp)-j) for j in range(1, len(data_exwas_tmp)+1)]))
    lowerlist =-np.log10(beta.ppf(0.025, [i for i in range(1, len(data_exwas_tmp)+1)], [(len(data_exwas_tmp)-j) for j in range(1, len(data_exwas_tmp)+1)]))
    upper_upperlist =-np.log10(beta.ppf(0.00001, [i for i in range(1, len(data_exwas_tmp)+1)], [(len(data_exwas_tmp)-j) for j in range(1, len(data_exwas_tmp)+1)]))
    #     logupper=-log10(qbeta(0.975, seq(1,nrow(data)), (nrow(data)-seq(1,nrow(data))+1)))
    # loglower=-log10(qbeta(0.025, seq(1,nrow(data)), (nrow(data)-seq(1,nrow(data))+1)))
    for i in range(1, len(data_exwas_tmp)+1):
        index.append(-math.log10(i/len(data_exwas_tmp)))
    ax1.plot(index,index, color= "black")
    ax1.plot(index, upperlist, color = "red")
    ax1.plot(index, lowerlist, color = "red")
    ax1.scatter(index,(data_exwas_tmp["LOG10P"]), color = "blue")
    ax1.set_xlabel("Expected (-log10(p))", fontsize =16)
    ax1.set_ylabel("Observed (-log10(p))", fontsize =16)
    ax1.set_title(mask + " (All Genes)", fontsize =16)
    # # upper = dir =1
    dir =-1
    for i, txt in enumerate(data_exwas_tmp["GENE_NAME"]):

        if (i <= 10) & (data_exwas_tmp["LOG10P"][i]>lowerlist[i]):
            ax1.annotate(txt, (index[i]+0.7*dir, data_exwas_tmp["LOG10P"][i]))
        dir=-dir

    
    
    mask_data["LOG10P"] = -np.log10(mask_data["P"].astype("float"))
    mask_data.sort_values(by = "P", inplace = True)
    mask_data.reset_index(inplace=True)
    data_exwas_tmp = mask_data
    index = []
    upperlist =-np.log10(beta.ppf(0.975, [i for i in range(1, len(data_exwas_tmp)+1)], [(len(data_exwas_tmp)-j) for j in range(1, len(data_exwas_tmp)+1)]))
    lowerlist =-np.log10(beta.ppf(0.025, [i for i in range(1, len(data_exwas_tmp)+1)], [(len(data_exwas_tmp)-j) for j in range(1, len(data_exwas_tmp)+1)]))
    upper_upperlist =-np.log10(beta.ppf(0.00001, [i for i in range(1, len(data_exwas_tmp)+1)], [(len(data_exwas_tmp)-j) for j in range(1, len(data_exwas_tmp)+1)]))
    for i in range(1, len(data_exwas_tmp)+1):
        index.append(-math.log10(i/len(data_exwas_tmp)))
    ax2.plot(index,index, color= "black")
    ax2.plot(index, upperlist, color = "red")
    ax2.plot(index, lowerlist, color = "red")
    ax2.scatter(index,(data_exwas_tmp["LOG10P"]), color = "blue")
    ax2.set_xlabel("Expected (-log10(p))", fontsize =16)
    ax2.set_ylabel("Observed (-log10(p))", fontsize =16)
    ax2.set_title(mask + " (Nearest Gene)", fontsize =16)
    # # upper = dir =1
    dir =-1
    for i, txt in enumerate(data_exwas_tmp["GENE"]):

        if (i <= 10) & (data_exwas_tmp["LOG10P"][i]>lowerlist[i]):
            ax2.annotate(txt, (index[i]+0.7*dir, data_exwas_tmp["LOG10P"][i]))
        dir=-dir
        
    mask_data = mask_data[mask_data["pLI"] > 0.9]
    mask_data.sort_values(by = "P", inplace = True)
    mask_data.reset_index(inplace=True)
    data_exwas_tmp = mask_data
    index = []
    upperlist =-np.log10(beta.ppf(0.975, [i for i in range(1, len(data_exwas_tmp)+1)], [(len(data_exwas_tmp)-j) for j in range(1, len(data_exwas_tmp)+1)]))
    lowerlist =-np.log10(beta.ppf(0.025, [i for i in range(1, len(data_exwas_tmp)+1)], [(len(data_exwas_tmp)-j) for j in range(1, len(data_exwas_tmp)+1)]))
    upper_upperlist =-np.log10(beta.ppf(0.00001, [i for i in range(1, len(data_exwas_tmp)+1)], [(len(data_exwas_tmp)-j) for j in range(1, len(data_exwas_tmp)+1)]))
    for i in range(1, len(data_exwas_tmp)+1):
        index.append(-math.log10(i/len(data_exwas_tmp)))
    ax4.plot(index,index, color= "black")
    ax4.plot(index, upperlist, color = "red")
    ax4.plot(index, lowerlist, color = "red")
    ax4.scatter(index,(data_exwas_tmp["LOG10P"]), color = "blue")
    ax4.set_xlabel("Expected (-log10(p))", fontsize =16)
    ax4.set_ylabel("Observed (-log10(p))", fontsize =16)
    ax4.set_title(mask + " (Nearest Gene, pLI > 0.9)", fontsize =16)
    # # upper = dir =1
    dir =-1
    for i, txt in enumerate(data_exwas_tmp["GENE"]):

        if (i <= 10) & (data_exwas_tmp["LOG10P"][i]>lowerlist[i]):
            ax4.annotate(txt, (index[i]+0.7*dir, data_exwas_tmp["LOG10P"][i]))
        dir=-dir
    
    
    bp_300_data["LOG10P"] = -np.log10(bp_300_data["P"].astype("float"))
    bp_300_data.sort_values(by = "P", inplace = True)
    bp_300_data.reset_index(inplace=True)
    data_exwas_tmp = bp_300_data
    index = []
    upperlist =-np.log10(beta.ppf(0.975, [i for i in range(1, len(data_exwas_tmp)+1)], [(len(data_exwas_tmp)-j) for j in range(1, len(data_exwas_tmp)+1)]))
    lowerlist =-np.log10(beta.ppf(0.025, [i for i in range(1, len(data_exwas_tmp)+1)], [(len(data_exwas_tmp)-j) for j in range(1, len(data_exwas_tmp)+1)]))
    upper_upperlist =-np.log10(beta.ppf(0.00001, [i for i in range(1, len(data_exwas_tmp)+1)], [(len(data_exwas_tmp)-j) for j in range(1, len(data_exwas_tmp)+1)]))
    for i in range(1, len(data_exwas_tmp)+1):
        index.append(-math.log10(i/len(data_exwas_tmp)))
    ax3.plot(index,index, color= "black")
    ax3.plot(index, upperlist, color = "red")
    ax3.plot(index, lowerlist, color = "red")
    ax3.scatter(index,(data_exwas_tmp["LOG10P"]), color = "blue")
    ax3.set_xlabel("Expected (-log10(p))", fontsize =16)
    ax3.set_ylabel("Observed (-log10(p))", fontsize =16)
    ax3.set_title(mask + " (300kbp Genes)", fontsize =16)
    # # upper = dir =1
    dir =-1
    for i, txt in enumerate(data_exwas_tmp["GENE"]):

        if (i <= 10) & (data_exwas_tmp["LOG10P"][i]>lowerlist[i]):
            ax3.annotate(txt, (index[i]+0.7*dir, data_exwas_tmp["LOG10P"][i]))
        dir=-dir
        
    bp_300_data = bp_300_data[bp_300_data["pLI"]>0.9]
    bp_300_data.sort_values(by = "P", inplace = True)
    bp_300_data.reset_index(inplace=True)
    data_exwas_tmp = bp_300_data
    index = []
    upperlist =-np.log10(beta.ppf(0.975, [i for i in range(1, len(data_exwas_tmp)+1)], [(len(data_exwas_tmp)-j) for j in range(1, len(data_exwas_tmp)+1)]))
    lowerlist =-np.log10(beta.ppf(0.025, [i for i in range(1, len(data_exwas_tmp)+1)], [(len(data_exwas_tmp)-j) for j in range(1, len(data_exwas_tmp)+1)]))
    upper_upperlist =-np.log10(beta.ppf(0.00001, [i for i in range(1, len(data_exwas_tmp)+1)], [(len(data_exwas_tmp)-j) for j in range(1, len(data_exwas_tmp)+1)]))
    for i in range(1, len(data_exwas_tmp)+1):
        index.append(-math.log10(i/len(data_exwas_tmp)))
    ax5.plot(index,index, color= "black")
    ax5.plot(index, upperlist, color = "red")
    ax5.plot(index, lowerlist, color = "red")
    ax5.scatter(index,(data_exwas_tmp["LOG10P"]), color = "blue")
    ax5.set_xlabel("Expected (-log10(p))", fontsize =16)
    ax5.set_ylabel("Observed (-log10(p))", fontsize =16)
    ax5.set_title(mask + " (300kbp Genes, pLI > 0.9)", fontsize =16)
    # # upper = dir =1
    dir =-1
    for i, txt in enumerate(data_exwas_tmp["GENE"]):

        if (i <= 10) & (data_exwas_tmp["LOG10P"][i]>lowerlist[i]):
            ax5.annotate(txt, (index[i]+0.7*dir, data_exwas_tmp["LOG10P"][i]))
        dir=-dir
    
    #
    plt.savefig(fname = "qq_plots_"  + args.trait + "_" + args.date + "_snps_" + str(mask) + ".png")
near_frame.to_csv("nearest_genes_" + args.trait +"_" + args.date + "_.txt", sep="\t", quoting=csv.QUOTE_NONE, index = False)
bp_300kbp_frame.to_csv("genes_300kb__" + args.trait + "_" + args.date+"_.txt", sep="\t", quoting=csv.QUOTE_NONE, index = False)

        #   qqplot(p_tmp)
  #  break
               

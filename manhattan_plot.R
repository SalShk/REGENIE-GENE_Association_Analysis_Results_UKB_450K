#options(repos = "https://www.stats.bris.ac.uk/R/")
# read in the argument(s) that was/were given on command line from the previous step, trailing from the execution of this Rscript from the bash script file:
args <- commandArgs(trailingOnly = TRUE)
# set working directory to that specified from the first argument provided on the command line, and then print working directory to command line - check that the script has loaded and working directory set correctly.
setwd(paste(paste(args[1], args[2], sep="/"), "_regenie_UKB_450K/", sep = ""))
getwd()
#load qqman library, and set x11 device to cairo 

library(qqman)
#options(devide='x11')
#x11(type="cairo")
library(Cairo)



# load in file to plot
maf0.001_file = read.table(paste(paste(paste(paste(args[1], args[2], sep="/"), "_regenie_UKB_450K/", sep = ""), args[2], sep = ""), "_transcripts_genes.0.001_450k.txt", sep=""),  header = T, stringsAsFactors = F)
print("maf_0.001_file loaded")
maf0.0001_file = read.table(paste(paste(paste(paste(args[1], args[2], sep="/"), "_regenie_UKB_450K/", sep = ""),  args[2], sep = ""), "_transcripts_genes.0.0001_450k.txt", sep=""), header = T, stringsAsFactors = F)
print("maf_0.0001_file loaded")

## Create the plot
dev.list()	# list the active graphic devices

CairoPNG(file = paste( args[2], "manhattan_all_transcripts_maf_0.001_450k.png", sep=""), width = 1800, height = 800) #call the pdf command to start the plot
#pdf(file = paste(paste(paste(args[1], args[2], sep = "/"), args[2], sep = "/"), "_transcripts_genes_0.001_450k.pdf")  #Call the pdf command to start the plot
maf0.001_file$LOG10P = as.numeric(maf0.001_file$LOG10P)
maf0.001_file = transform(maf0.001_file, P = 10^(-1*LOG10P))
maf0.001_file$P = as.numeric(maf0.001_file$P)
manhattan(maf0.001_file, chr="CHROM", bp="POS", snp="ID", p="P", logp = TRUE, genomewideline = FALSE, suggestiveline = FALSE, annotatePval = 0.00001, 
annotateTop =TRUE, main = paste(args[2], "exWAS UKB 450K MAF<0.1%", sep= " "), cex.axis = 0.9)

dev.off() # turn off plot devices
setEPS()
postscript('test.eps')

## Create the plot
dev.list()	# list the active graphic devices

CairoPNG(file = paste( args[2], "manhattan_all_transcripts_maf_0.0001_450k.png", sep=""), width = 1800, height = 800) #call the pdf command to start the plot
#pdf(file = paste(paste(paste(args[1], args[2], sep = "/"), args[2], sep = "/"), "_transcripts_genes_0.001_450k.pdf")  #Call the pdf command to start the plot
maf0.0001_file$LOG10P = as.numeric(maf0.0001_file$LOG10P)
maf0.0001_file = transform(maf0.0001_file, P = 10^(-1*LOG10P))
maf0.0001_file$P = as.numeric(maf0.0001_file$P)
manhattan(maf0.0001_file, chr="CHROM", bp="POS", snp="ID", p="P", logp = TRUE, genomewideline = FALSE, suggestiveline = FALSE, annotatePval = 0.00001,
annotateTop =TRUE, main = paste(args[2], "exWAS UKB 450K MAF<0.01%", sep= " "),  cex.axis = 0.9)

dev.off() # turn off plot devices
setEPS()
postscript('test.eps')





#png(file = paste( args[2], "manhattan_all_transcripts_maf_0.0001_450k.png", sep="")) #call the pdf command to start the plot
#pdf(file = paste(paste(paste(args[1], args[2], sep = "/"), args[2], sep = "/"), "_transcripts_genes_0.0001_450k.pdf")  #Call the pdf command to start the plot
#maf0.0001file$LOG10P = as.numeric(maf0.0001_file$LOG10P)
#maf0.0001_file = transform(maf0.0001_file, P = 10^(-1*LOG10P))
#maf0.0001_file$P = as.numeric(maf0.0001_file$P)
#manhattan(maf0.0001_file, chr="CHROM", bp="POS", snp="ID", p="P", logp = TRUE, genomewideline = FALSE, suggestiveline = FALSE,
#annotatePval = 0.0001, annotateTop =TRUE, main = paste(args[2], "exWAS UKB 450K MAF<0.01%", sep= " "), cex.axis = 0.9)

#dev.off() # turn off plot devices


## PLOT 2 


# List of SNPs to highlight are in the snpsOfInterest object
# We will use ggrepel for the annotation
#library(ggrepel)
#library(ggplot2)
#library(dplyr)
# Prepare the dataset
#don <- Exwas_Results_0.001  %>% 
  
  # Compute chromosome size
#  group_by(CHROM) %>% 
#  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
#  mutate(tot=cumsum(chr_len)-chr_len) %>%
#  select(-chr_len) %>%
  
  # Add this info to the initial dataset
#  left_join(Exwas_Results_0.001, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
#  arrange(CHROM, POS) %>%
#  mutate( POScum=POS+tot) %>%

  # Add highlight and annotation information
 ## #mutate( is_highlight=ifelse(ID %in% snpsOfInterest, "yes", "no")) %>%
#  mutate( is_annotate=ifelse(-log10(LOG10P)>4, "yes", "no")) 

# Prepare X axis
#axisdf <- don %>% group_by(CHROM) %>% summarize(center=( max(POScum) + min(POScum) ) / 2 )

# Make the plot
#ggplot(don, aes(x=POScum, y=-log10(LOG10P))) +
    
    # Show all points
 #   geom_point( aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
 #   scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    
    # custom X axis:
 #   scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center ) +
 #   scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis

    # Add highlighted points
 #   geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  
    # Add label using ggrepel to avoid overlapping
 #   geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=ID), size=2) +

    # Custom the theme:
 #   theme_bw() +
 #   theme( 
 #     legend.position="none",
 #     panel.border = element_blank(),
 #     panel.grid.major.x = element_blank(),
 #     panel.grid.minor.x = element_blank()
 #   )

#ggsave("customised_manhattan.pdf")

##**********************************************************##
## Global Screening Array Illumina (GSA) Post-imputation QC ##
##**********************************************************##
## Tim Hasenbein
## PhD. Guillermo G. Torres
## Last modification 03.2020
## Creation: 01.2019
## Aim: Perform post-imputation QC on imputed data set
library(data.table)
library(tidyverse)


##### Set envorionment #####
mypath <- if(class(try(rstudioapi::getSourceEditorContext()$path,silent=TRUE))=="try-error") {paste0(dirname(substring(argv[grep("--file=", argv)], 8)),'/')} else paste(head(strsplit(rstudioapi::getSourceEditorContext()$path,'/')[[1]],-1),collapse='/')
setwd(mypath)
    

##### Extract files #####
system("for i in {1..22}; do 7z e chr_$i.zip -pSpiOZ7xgZab9Z0; done", show.output.on.console = TRUE) 


##### Prepare input files #####
system("sudo ./vcfparse.pl -d /home/tim/Desktop/Master/post_imp_with_ExChip/post_imputation -o post_imputation -g")


##### Check quality of imputation #####
system("7z e chr1.dose.vcf.cut.gz")
chr <- fread("chr1.dose.vcf.cut", header=T)
chr <- separate(chr, INFO, into = c("V8","V9","V10","V11"), sep = ";")
chr <- chr%>%rename("R2" = "V10")%>%select("POS", "R2")
chr <- separate(chr, R2, into = c("del","R2"), sep = "=")%>%select(-"del")  
chr_plot <- ggplot(data = chr, mapping =  aes(x = as.numeric(POS), y = as.numeric(R2))) + 
geom_point(alpha = 0.05, size = 0.01, colour = "black") +
  labs(title = "Chromosome 1",x = "Position (bp)", y = "R2-value") + 
  ylim(0,1) + 
  xlim(13380, 249239762)
chr_plot
ggsave('chr1_plot.png', width=8,height=6)  


##### Apply r2 threshold  #####
system("bcftools view -i 'R2>.75' -Oz chr20.dose.vcf.gz > chr20.filtered.vcf.gz") # for one file
system("for i in {1..3}; do bcftools view -i 'R2>.75' -Oz chr$i.dose.vcf.gz > chr$i.filtered.vcf.gz; done") # for multiple files

                                                              
##### Convert filtered chr to plink  #####
recode <- function(filename){
  for (chr in chrs){
    system(paste0("plink --vcf ",chr," --double-id  --make-bed --out ",chr))
  }
}
chrs <- Sys.glob(paste0(mypath,"/chr*.filtered.vcf.gz")) 
recode(filename = chrs)


##### Merge files  #####
all_chrs <- gsub("\\.bed","",Sys.glob("chr*.bed"))
write.table(all_chrs, "all_chrs.txt", row.names = F, col.names = F, sep = " ", quote = F) 
system("plink --bfile chr1 --merge-list list.txt --make-bed --out GSAdata_imp")
# 13323318 variants


##### Step 7: Get Case/ctrl + sex status back  #####
# take the .fam file before imputation and rename it to GSAdata_imp.fam -> old one is GSAdata_imp_backup.fam


##### Update bim file variant name #####
GSAdata_imp <- fread("GSAdata_imp.bim", header = F)
head(GSAdata_imp)
update <- GSAdata_imp%>%mutate(V2=paste0("chr",GSAdata_imp$V1,":",GSAdata_imp$V4))
head(update)
write.table(update, "GSAdata_imp.bim", row.names = F, col.names = F, sep = " ", quote = F)


############ PER-Marker QC ############


###### Identification of SNPs showing a significant deviation from HWE  
system(paste0("plink --bfile ",results,"GSA_impQC --hwe 0.000000001 --make-bed --out ",results,"GSA_impQC_temp")) 
tempstep(paste0(results,"GSA_impQC_temp"),paste0(results,"GSA_impQC")) # 13089455 variants pass


###### Identification of SNPs with significantly different missing genotype rates btw cases/ controls ######
system(paste0("plink --bfile ",results,"GSA_impQC --test-missing --out ",results,"GSA_impQC"))  
inform_missing <- fread(paste0(results,'GSA_impQC.missing'), header=T)
fail_diffmiss_test <- data.frame(inform_missing$SNP[inform_missing$P < 0.00001]) 
write.table(fail_diffmiss_test, paste0(results,"fail_diffmiss_test"), row.names = F, col.names = F, sep = " ", quote = F) 
system(paste0("plink --bfile ",results,"GSA_impQC --exclude ",results,"fail_diffmiss_test --make-bed --out ",results,"GSA_impQC_temp")) # 13085738 variants pass
tempstep(paste0(results,"GSA_impQC_temp"),paste0(results,"GSA_impQC"))


###### emoval of all markers with low MAF ######
system("plink --bfile GSA_impQC --maf 0.01 --make-bed --out GSA_impQC_3103") # 7363285 variants pass





  


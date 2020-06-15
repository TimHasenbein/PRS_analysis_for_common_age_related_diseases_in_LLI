##************************************************************##
## Global Screening Array Illumina (GSA) Association analysis ##
##************************************************************##
## Tim Hasenbein
## PhD. Guillermo G. Torres
## Last modification 03.2020
## Creation: 01.2020
## Aim: Perform single-variant association analysis with unimputed/imputed data set 
library(qqman)
library(data.table)
library(tidyverse)
require(GWASTools)


##### set up env. var ######
mypath <- if(class(try(rstudioapi::getSourceEditorContext()$path,silent=TRUE))=="try-error") {paste0(dirname(substring(argv[grep("--file=", argv)], 8)),'/')} else paste(head(strsplit(rstudioapi::getSourceEditorContext()$path,'/')[[1]],-1),collapse='/')
setwd(mypath)
genoQCed <- "F:/QC_1603/results/"
impQCed <- "F:/Imputation/Post_imp/filtered_recoded/results/"


############### Association analysis: Unimputed data ###############


##### Association testing #####
system(paste0("plink --bfile ",genoQCed,"GSAdata_QCed_3003_assoc  --logistic --covar GSAhapmap.eigenvec --covar-number 1-5 --out GSAunimp"))
GSAassoc <- fread("GSAunimp.assoc.logistic", data.table = F)
GSAassoc_add <- GSAassoc[GSAassoc$TEST == "ADD", ] # 447285 variants
sigSNPs <- GSAassoc_add %>% arrange(P) %>% filter(P <= 5e-08)
SNPlistraw <- fread("GSAdata_QCed_3003_assoc.bim", data.table = F)
sigSNPs <- sigSNPs %>% left_join(SNPlistraw, by=c('BP' = 'V4')) %>% select(-"V1", -"V3",-"V5", -"V6")
SNPshighlight <- sigSNPs %>% select("SNP")
write.table(SNPshighlight,"SNPshighlight", col.names = T, row.names = F, sep = " ", quote = F)
system("plink --bfile GSAdata_QCed_3003_assoc  --assoc --adjust --out GSAassoc") 


##### Data visualization #####
qqPlot(GSAassoc_add$P, thinThreshold = 2, ci = T, ylim = c(0,20), cex = .5) 
ggsave("qqplot_unimp.pdf")
GSAassoc <- GSAassoc_add[ , c(1,2,3,9)]  
GSAassoc <- GSAassoc[!is.na(GSAassoc$P), ]
write.table(GSAassoc,"GSAassoc", col.names = T, row.names = F, sep = " ", quote = F)
chr19 <- GSAassoc_add %>% filter(CHR == 19)
apoe <- chr19 %>% filter(BP %in% 45000000:45800000)
apoesnps <- apoe$SNP
manhattan(GSAassoc, chr="CHR",bp="BP",p="P",snp="SNP", cex = 0.3, ylim = c(0, 25), highlight = apoesnps) 
manhattan(subset(GSAassoc, CHR == 19), highlight = apoesnps, cex = 0.3) 

  
############### Association analysis: Imputed data ###############


##### Association testing #####
system(paste0("plink --bfile ",impQCed,"GSA_imp_3103_QCed  --logistic --covar GSAhapmap.eigenvec --covar-number 1-5 --out GSAimp_3103"))
GSAassocImp <- fread(paste0(impQCed,"GSAimp_3103.assoc.logistic"), data.table = F)
GSAassocImp_add <- GSAassocImp[GSAassocImp$TEST == "ADD", ] # 7363285 variants
sigSNPsImp <- GSAassocImp_add %>% arrange(P) %>% filter(P <= 5e-08)
SNPlistrawImp <- fread(paste0(impQCed,"GSA_imp_3103_QCed.bim"), data.table = F)
sigSNPsImp <- sigSNPsImp %>% left_join(SNPlistrawImp, by=c('BP' = 'V4')) %>% select(-"V1", -"V3",-"V5", -"V6")
SNPshighlightImp <- sigSNPsImp %>% select("SNP")
write.table(SNPshighlightImp,"SNPshighlightImp", col.names = T, row.names = F, sep = " ", quote = F)
system(paste0("plink --bfile ",impQCed,"GSA_imp_3103_QCed  --assoc --adjust --out GSAImpassoc"))


##### Data visualization #####
qqPlot(GSAassocImp_add$P, thinThreshold = 2, ci = T, ylim = c(0,20), cex = .5)
GSAassocImp <- GSAassocImp_add[ , c(1,2,3,9)]  
GSAassocImp <- GSAassocImp[!is.na(GSAassocImp$P), ]
write.table(GSAassocImp,"GSAassocImp", col.names = T, row.names = F, sep = " ", quote = F)
chr19Imp <- GSAassocImp_add %>% filter(CHR == 19)
apoeImp <- chr19Imp %>% filter(BP %in% 45000000:45800000)
apoesnpsImp <- apoeImp$SNP
sig <- c("chr8:24754946", "chr5:126086995", "chr5:126087063","chr6:61940050")
highlight <- c(apoesnpsImp, sig)
manhattan(GSAassocImp, chr="CHR",bp="BP",p="P",snp="SNP", cex = 0.3, ylim = c(0, 25), highlight = highlight) 
manhattan(subset(GSAassocImp, CHR == 19), highlight = apoesnpsImp, cex = 0.3,) 



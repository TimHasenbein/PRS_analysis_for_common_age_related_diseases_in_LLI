##***********************************************************************##
## Global Screening Array Illumina (GSA) Pre-imputation data preparation ##
##***********************************************************************##
## Tim Hasenbein
## PhD. Guillermo G. Torres
## Last modification 03.2020
## Creation: 01.2020
## Aim: Prepare data for MIS imputation
library(data.table)
library(tidyverse)


###### Set up env. var ######
mypath <- if(class(try(rstudioapi::getSourceEditorContext()$path,silent=TRUE))=="try-error") {paste0(dirname(substring(argv[grep("--file=", argv)], 8)),'/')} else paste(head(strsplit(rstudioapi::getSourceEditorContext()$path,'/')[[1]],-1),collapse='/')
setwd(mypath)
datafolder <- mypath
GSAdata <- "F:/QC_1603/results/"
results <- paste0(datafolder, "/results/")
tempstep <- function(infile,outfile){system(paste0("plink --bfile ",infile," --make-bed --out ",outfile))} 

  
##### Merge SNPs (cases) with ExChip cases #####
# select exome cases (ctrls are not the same)
Exome_fam <- read.table(paste0(datafolder,'/Ageing_ExomeChip_QCed.fam'),quote='',header=F,as.is=T)%>%as_tibble()
Exome_cases <- Exome_fam%>%filter(str_detect(V2,"AGI+"))
write.table(Exome_cases, paste0(results,"Exome_cases.txt"), row.names = F, col.names = F, sep = " ", quote = F) 
system(paste0("plink --bfile ",datafolder,"/Ageing_ExomeChip_QCed --keep ",results,"Exome_cases.txt --make-bed --out ",results,"ExChipQCed_cases")) 
# merge exome and GSA
Exome_fam <- read.table(paste0(results,'ExChipQCed_cases.fam'),quote='',header=F,as.is=T)%>%as_tibble()
Exome_merge <- Exome_fam%>%
  mutate("V7"= V2)%>%
  separate(V7,into = c("V7","V8","ID"), sep = "_")%>%
  select(-"V7", -"V8")
GSA_fam <- read.table(paste0(GSAdata,'GSAdata_QCed_1603.fam'),quote='',header=F,as.is=T)%>%as_tibble()
GSA_cases <- GSA_fam%>%filter(str_detect(V2,"AGE+"))
GSA_merge <- GSA_cases%>%
  mutate("V7"= V2)%>%
  separate(V7,into = c("V7","ID","V8"), sep = "_")%>%
  select("V2", "ID")
GSA_exome <- Exome_merge %>%
  left_join(GSA_merge, by = "ID")%>%
  filter(!is.na(V2.y))
GSA_exome_ind <- GSA_exome%>%rename("V2"="V2.x")%>%select(-"ID", -"V2.y") 
write.table(GSA_exome_ind, paste0(results,"GSA_exome_ind.txt"), row.names = F, col.names = F, sep = " ", quote = F) 
GSA_bim <- read.table(paste0(GSAdata,'GSAdata_QCed_1603.bim'),quote='',header=F,as.is=T)%>%as_tibble()
Exome_bim <- read.table(paste0(results,'ExChipQCed_cases.bim'),quote='',header=F,as.is=T)%>%as_tibble()
ebim <- Exome_bim%>%mutate(ID=paste0('chr',V1,':',V4))
write.table(ebim%>%select(V2,ID), paste0(results,"ExChip_newNames.txt"), row.names = F, col.names = F, sep = " ", quote = F) 
Exome_SNPnoinGSA<- anti_join(ebim,GSA_bim,by=c('ID'='V2'))
write.table(Exome_SNPnoinGSA%>%select(ID), paste0(results,"ExChip_SNPsUnique.txt"), row.names = F, col.names = F, sep = " ", quote = F) 
system(paste0("plink --bfile ",datafolder,"/Ageing_ExomeChip_QCed --keep ",results,"GSA_exome_ind.txt --make-bed --out ",results,"ExChipQCed_cases")) 
system(paste0("plink --bfile ",results,"ExChipQCed_cases --update-map ",results,"ExChip_newNames.txt --update-name --make-bed --out ",results,"ExChipQCed_cases_temp --allow-no-sex"))
system(paste0("plink --bfile ",results,"ExChipQCed_cases_temp --extract ",results,"ExChip_SNPsUnique.txt --make-bed --out ",results,"ExChipQCed_cases --allow-no-sex")) 
write.table(GSA_exome%>%select(V1,V2.x,V2.y)%>%mutate(Y=V1)%>%select(V1,V2.x,Y,V2.y), paste0(results,"ExChipQ_cases_UpdateIDs.txt"), row.names = F, col.names = F, sep = " ", quote = F)
system(paste0("plink --bfile ",results,"ExChipQCed_cases --update-ids ",results,"ExChipQ_cases_UpdateIDs.txt --make-bed --out ",results,"ExChipQCed_cases_temp --allow-no-sex")) 
tempstep(paste0(results,"ExChipQCed_cases_temp"),paste0(results,"ExChipQCed_cases"))
system(paste0("plink --bfile ",GSAdata,"GSAdata_QCed_1603 --bmerge ",results,"ExChipQCed_cases --make-bed --out ",results,"GSA_ExChip"))
system(paste0("plink --bfile ",results,"ExChipQCed_cases --exclude ",results,"GSA_ExChip-merge.missnp --make-bed --out ",results,"ExChipQCed_cases_GSAX")) 
system(paste0("plink --bfile ", GSAdata, "GSAdata_QCed_1603 --exclude ",results,"GSA_ExChip-merge.missnp  --make-bed --out ",results,"GSAdata_EchipX")) # 634196 variants 
system(paste0("plink --bfile ",results,"GSAdata_EchipX --bmerge ",results,"ExChipQCed_cases_GSAX --make-bed --out ",results,"GSA_ExChip")) # 824489 variants and 5596 people pass 


##### Pre-imputation check for HRC imputation #####
system("wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.7.zip")
system("wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz")
system(paste0("plink --bfile ",subsets,"GSA_ExChip --freq --out ",imputation,"GSA_ExChip"))
sytem("perl HRC-1000G-check-bim.pl -b GSA_ExChip.bim -f GSA_ExChip.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h")
system("sh Run-plink.sh") # 712837 variants


##### Data preparation according to MIS ####
# recode
imputation_MichiganServer <- function(out_folder,filename){
  for (chr in chrs){
    system(paste0("plink --bfile ",chr," --recode vcf --out ",chr))
  }
}
chrs <- gsub("\\.bed","",Sys.glob(paste0(imputation,"GSA_ExChip-updated-chr*.bed"))) 
imputation_MichiganServer(out_folder=imputation,filename=chrs)
# sort chr and compress files
for (chr in chrs){
  system(paste0("bcftools sort ",chr,".vcf -Oz -o ",chr,".vcf.gz"))
}


##### Upload files at MIS #####
# Genotype Imputation (Minimac4) 1.2.4
# Reference Panel: HRC r1.1 2016 (GRCh37/hg19)
# Array Build: GRCh37/hg19
# rsq Filter: No
# Phasing:  Eagle4.0
# Population: EUR
# Mode: Quality Control & Imputation



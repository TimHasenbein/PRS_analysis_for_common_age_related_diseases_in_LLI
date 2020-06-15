##*********************************************************##
## Global Screening Array Illumina (GSA) Quality control   ##
##*********************************************************##
## Tim Hasenbein
## PhD. Guillermo G. Torres 
## Last modification 03.2020
## Creation: 11.2019
## QC on GSA raw genotype data


###### Get access to toolbox with functions ######
argv <- commandArgs(trailingOnly = FALSE)
mypath <- if(class(try(rstudioapi::getSourceEditorContext()$path,silent=TRUE))=="try-error") {paste0(dirname(substring(argv[grep("--file=", argv)], 8)),'/')} else paste(head(strsplit(rstudioapi::getSourceEditorContext()$path,'/')[[1]],-1),collapse='/')
toolbox <- paste0(mypath,'/GSA_toolbox.R')
source(toolbox)
packages(c('tidyverse','openxlsx','data.table'))
dircreations <- function(dir_path){
  if(dir.exists(dir_path)){message(paste0(dir_path,' Folder already exist, overwriting...'))
  }else dir.create(dir_path)
}


###### Set up env. var ######
rootfolder <- mypath
datafolder <- paste0(rootfolder, "/rawdata/")
hapmap <-  paste0(datafolder, "hapmap/")
results <- paste0(rootfolder, "/results/")
figures <- paste0(results, "figures/")
subsets <- paste0(results, "subsets/")
IBD <- paste0(results, "IBD/")
stratification <- paste0(results, "stratification/")
dircreations(datafolder) ;dircreations(hapmap);dircreations(results);dircreations(figures);dircreations(subsets);dircreations(IBD);dircreations(stratification)
tempstep <- function(infile,outfile){system(paste0("plink --bfile ",infile," --make-bed --out ",outfile))} # for windows to delete files (input cannot be output)


#############  Pre-quality control data preparation  ############# 

###### General Statistics and MAF Trimming ######
system(paste0("plink --bfile ",datafolder,"GSAdata_raw --freq --out ",results,"GSAdata_raw")) 
system(paste0("plink --bfile ",datafolder,"GSAdata_raw --freq case-control --out ",results,"GSAdata_raw")) 
system(paste0("plink --bfile ",datafolder,"GSAdata_raw --freqx --out ",results,"GSAdata_raw"))
# 700078 variants, 6031 people (3085 males, 2946 females)


###### Update SNPs with lost physical information ######
## where CHR=0 POS=0 IDS=... 
## chr24=Y; chr23=X; chr26=M (mitochondria)
bim <- read.table(paste0(datafolder,"GSAdata_raw.bim"),quote='',header=F,as.is=T)%>%as_tibble()
snp0 <-  bim%>%filter(V4==0)
Y <- snp0%>%filter(str_detect(V2,"Y-[[:digit:]]+"))
X <- snp0%>%filter(str_detect(V2,"X-[[:digit:]]+"))
GSAchrpos <- snp0%>%filter(str_detect(V2,"^GSA-[[:digit:]]+:[[:digit:]]+"))
rss <- snp0%>%filter(str_detect(V2,"rs[[:digit:]]+"))
rsstable <- tibble(ID=rss%>%pull(V2), RS_ID=rss%>%dplyr::pull(V2)%>%str_extract("rs[[:alnum:]]+"))
rss_chpos <- getSNPsPosH19(rsstable$RS_ID) 
DF <- rsstable%>%left_join(rss_chpos,by=c('RS_ID'))
GSAdf <- chrposDF(GSAchrpos,prefix="GSA-",splited=':')
Ydf <- chrposDF(Y,prefix="[[:digit:]]+[[:punct:]][[:digit:]]{1,2}[[:punct:]]",splited='-')
Xdf <- chrposDF(X,prefix="IDS-chr",splited='-')
DF <- bind_rows(DF,GSAdf,Ydf,Xdf)
snpswithchrpos <- bim%>%filter(V4!=0)
snpswithchrposdf <- chrposDFv2(snpswithchrpos)
DFtot <- bind_rows(DF,snpswithchrposdf) # 7378 variants had lost informations 
write_tsv(DFtot,paste0(results,"Total_GSA_SNPset.txt")) 
write_tsv(DFtot%>%filter(Chr!='0')%>%.[,'ID'],paste0(results,"includeSNPs.txt"),col_names=F) 
write_tsv(DFtot%>%filter(Chr!='0')%>%dplyr::select(ID,Chr_Pos),paste0(results,"updateSNPsName.txt"),col_names=F) 
write_tsv(DFtot%>%filter(Chr!='0')%>%dplyr::select(ID,Chr),paste0(results,"updateSNPsChr.txt"),col_names=F) 
write_tsv(DFtot%>%filter(Chr!='0')%>%dplyr::select(ID,Pos),paste0(results,"updateSNPsPos.txt"),col_names=F)
## Remove markers with still Chr=0 pos=0   
system(paste0("plink --bfile ",datafolder,"GSAdata_raw --extract ",results,"includeSNPs.txt --make-bed --out ",results,"GSAdata")) # 699940 variants and 6031 people (3085 males, 2946 females)
## Updating Chromosome and physical position
system(paste0("plink --bfile ",results,"GSAdata --update-map ",results,"updateSNPsPos.txt --make-bed --out ",results,"GSAdata_temp")) # Position
tempstep(paste0(results,"GSAdata_temp"),paste0(results,"GSAdata"))
system(paste0("plink --bfile ",results,"GSAdata --update-map ",results,"updateSNPsChr.txt --update-chr --make-bed --out ",results,"GSAdata_temp")) # Chromosome
tempstep(paste0(results,"GSAdata_temp"),paste0(results,"GSAdata"))


###### Exclude duplicated markers ######
system(paste0("plink --bfile ",results,"GSAdata --list-duplicate-vars require-same-ref --missing --out ",results,"GSAdata")) 
snpdup <- read_tsv(paste0(results,'GSAdata.dupvar'),quote='',col_names=T)
lmiss <- read.table(paste0(results, "GSAdata.lmiss"),quote='',header=T,as.is=T)%>%as.tibble()
## Separate duplicated SNPs and make table with corresponding F_MISS
keepSNP <- rep(NA,nrow(snpdup))
snpdupfmis <- snpdup%>%separate_rows(IDS,sep=" ")%>%left_join(lmiss%>%dplyr::select(SNP,F_MISS),by=c("IDS"="SNP"))
pb <- txtProgressBar(min = 0, max = nrow(snpdup), style = 3)
for(r in 1:nrow(snpdup)){
  rowfmis <- snpdupfmis%>%filter(CHR==snpdup[r,]$CHR & POS==snpdup[r,]$POS)
  keepSNP[r] <- rowfmis%>%filter(F_MISS==min(F_MISS))%>%pull(IDS)%>%.[1]
  setTxtProgressBar(pb, r)
}
snps_remove <- snpdupfmis%>%filter(!IDS%in%keepSNP)%>%dplyr::select(IDS) # 1094 variants
write.table(snps_remove, paste0(results,"snps_remove.txt"), row.names = F, col.names = F, sep = " ", quote = F)
system(paste0("plink --bfile ",results,"GSAdata --exclude ",results,"snps_remove.txt --make-bed --out ",results,"GSAdata_temp --allow-no-sex")) # 698846 variants and 6031 people pass
tempstep(paste0(results,"GSAdata_temp"),paste0(results,"GSAdata"))
## Update variant ids to chr_:bp 
system(paste0("plink --bfile ",results,"GSAdata --update-map ",results,"updateSNPsName.txt --update-name --make-bed --out ",results,"GSAdata_temp")) 
tempstep(paste0(results,"GSAdata_temp"),paste0(results,"GSAdata"))


#############  QUALITY CONTROL  #############
############ PER-INDIVIDUAL QC ############ 


###### Identification of individuals with discordant sex information  ###### 
system(paste0("plink --bfile ",results,"GSAdata --check-sex --out ",results,"GSAdata"))
sexcheck <- read.table(paste0(results,'GSAdata.sexcheck'),quote='',header=T,as.is=T)%>%as_tibble()
genderPlot(sexcheck,outplot=paste0(figures,"01_Gender.pdf"),outfiles=results) 
system(paste0("plink --bfile ",results,"GSAdata --remove ",results,"Indiv.sexproblem --make-bed --out ",results,"GSAdata_temp")) # 698846 variants and 6012 people pass (19 indiviudals removed, 1347 cases, 4665 controls)
tempstep(paste0(results,"GSAdata_temp"),paste0(results,"GSAdata"))


##### Identification of individuals with outlying missing genotype/heterozygosity rate ######
system(paste0("plink --bfile ",results,"GSAdata --missing --het --out ",results,"GSAdata --allow-no-sex"))
imiss <- read.table(paste0(results,'GSAdata.imiss'),quote='',header=T,as.is=T)%>%as.tibble()
het <- read.table(paste0(results,'GSAdata.het'),quote='',header=T,as.is=T)%>%as.tibble()
ihfail <- imiss_hetPlot(imiss,het,outplot=paste0(figures,"02_Indiv_missingInfo.pdf"),outfiles=results,pmiss=0.081,timesSD=4)## 71 individuals failed heterozygosity and missigness 
system(paste0("plink --bfile ",results,"GSAdata --remove ",results,"Indiv.imisshetfail --make-bed --out ",results,"GSAdata_temp")) # 5941  people pass (3051 males, 2890 females)
tempstep(paste0(results,"GSAdata_temp"),paste0(results,"GSAdata"))


###### Indentification of individuals with divergent ancestry ######
system(paste0("plink --bfile ",results,"GSAdata --exclude ", hapmap,"high-LD-regions.txt --range --indep-pairwise 50 5 0.2 --make-bed --out ",stratification,"GSAdataIBD_LDexclude --allow-no-sex")) # 677965 variants and 5941 people pass
system(paste0("plink --bfile ",stratification,"GSAdataIBD_LDexclude --extract ",stratification,"GSAdataIBD_LDexclude.prune.in --make-bed --out ",stratification,"GSAdataIBD_LDexclude_pruned")) # 279860 variants and 5941 people pass
## Create SNPlist with variants present in ancestry + GSA data set and extract these
hapmapsnplist <- data.frame(read.table(paste0(hapmap,"hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps_hg19.snplist")))%>%as_tibble()
gsasnplist <- read.table(paste0(stratification,"GSAdataIBD_LDexclude_pruned.bim"),quote='',header=F,as.is=T)%>%as_tibble()
snplist_ancestry <- inner_join(gsasnplist,hapmapsnplist,by=c('V2'='V1'))
write.table(snplist_ancestry%>%pull(V2),paste0(stratification,"common_snps.txt"), row.names=F, col.names = F, sep = " ", quote = F)
system(paste0("plink --bfile ",stratification,"GSAdataIBD_LDexclude_pruned --extract ",stratification,"common_snps.txt --make-bed --out ",stratification,"GSA_hapmapsnps")) # 53909 variants and 5941 people 
system(paste0("plink --bfile ",hapmap,"hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps_hg19 --extract ",stratification,"common_snps.txt --make-bed --out ",stratification,"hapmap_commonsnps")) # 1012421 variants loaded, 53896 variants and 395 people pass
## Merge GSA_hapmapsnps now with hapmap.founders and remove 3+ allelic variants
system(paste0("plink --bfile ",stratification,"GSA_hapmapsnps --bmerge ", stratification,"hapmap_commonsnps --make-bed --out ",stratification,"GSAhapmap")) 
system(paste0("plink --bfile ",stratification,"GSA_hapmapsnps --exclude ",stratification,"GSAhapmap-merge.missnp --make-bed --out ",stratification,"GSA_excl")) # 53856 variants and 5941 people 
system(paste0("plink --bfile ", stratification,"hapmap_commonsnps --exclude ",stratification,"GSAhapmap-merge.missnp  --make-bed --out ",stratification,"hapmap_commonsnps_excl")) # 53855 variants and 395 people 
system(paste0("plink --bfile ",stratification,"GSA_excl --bmerge ",stratification,"hapmap_commonsnps_excl --make-bed --out ",stratification,"GSAhapmap")) # 53855 variants and 6336  people
## Perform PCA with merged data set in low LD
system(paste0("plink --bfile ",stratification,"GSAhapmap --pca header tabs --make-rel --out ",stratification,"GSAhapmap"))
## Get relationships and eigenvec
relationships_w_pops <- read.table(paste0(hapmap,"/relationships_w_pops_041510_hapmap.txt"),quote='',header=T,as.is=T)%>%as.tibble()
eigenvecs <- read.table(paste0(stratification,"GSAhapmap.eigenvec"),quote='',header=T,as.is=T)%>%as.tibble()
eigenvec <- eigenvecs %>% 
  left_join(relationships_w_pops, by = c("IID")) %>%
  select(-FID.y, -dad, -mom, -sex, -pheno)
eigenvec[is.na(eigenvec)] <- "Germans"
eigenvec <- eigenvec%>%mutate(CaseControl = ifelse(str_detect(IID,"AGE"),"Case","Control"))
## Create plot and remove outliers
pcaPlot(eigenvec,pcs=c('PC1-PC2','PC1-PC3','PC2-PC3'),outplot=paste0(figures,"03_Strat_hapmap"),outfiles=results,getOutilers=T,lofth=2.1,pop='Germans')
system(paste0("plink --bfile ",results,"GSAdata --remove ",results,"Indiv.PopStratFail --make-bed --out ",results,"GSAdata_temp --allow-no-sex")) #  698846 variants and 5796 people pass
tempstep(paste0(results,"GSAdata_temp"),paste0(results,"GSAdata"))
file.remove(list.files(stratification,'GSAdataIBD_LDexclude|GSA_excl|hapmap_commonsnps',full.names = T),recursive=TRUE)


###### Identification of duplicated or related individuals ###### 
system(paste0("plink --bfile ",results,"GSAdata --exclude ", hapmap,"high-LD-regions.txt --range --indep-pairwise 50 5 0.2 --make-bed --out ",IBD,"GSAdataIBD_LDexclude --allow-no-sex")) 
system(paste0("plink --bfile ",IBD,"GSAdataIBD_LDexclude --extract ",IBD,"GSAdataIBD_LDexclude.prune.in --make-bed --out ",IBD,"GSAdataIBD_LDexclude_pruned")) # 273397 variants and 5796 people pass
system(paste0("plink --bfile ",IBD,"GSAdataIBD_LDexclude_pruned --genome --out ",IBD,"GSAdata_IBD --allow-no-sex"))
ibd <- fread(paste0(IBD,"GSAdata_IBD.genome"), data.table = F)
ibd_pi <- ibd%>%filter(PI_HAT>0.185)%>%as_tibble()
keepIID <- c()
RMIID <- c()
rowkeep <- c()
ibd_pi_imis <- ibd_pi%>%dplyr::select(IID1,IID2)%>%gather(c('IID1','IID2'),key='Indiv',value='IDS')%>%left_join(imiss%>%dplyr::select(IID,F_MISS),by=c("IDS"="IID"))%>%distinct(IDS, .keep_all = TRUE)
pb <- txtProgressBar(min = 0, max = nrow(ibd_pi), style = 3)
c=1
for(r in 1:nrow(ibd_pi)){
  right <- ibd_pi[r,]$IID1;  left <- ibd_pi[r,]$IID2
  if(right%>%str_detect("AGE") & left%>%str_detect("AGE")) {
    rowimis <- ibd_pi_imis%>%filter(IDS%in%c(right,left))
    min <- rowimis%>%filter(F_MISS==min(F_MISS))%>%pull(IDS)%>%.[1]
    max <- rowimis%>%filter(F_MISS==max(F_MISS))%>%pull(IDS)%>%.[1]
    if(!min%in%keepIID){keepIID <- c(keepIID,min)};if(!max%in%RMIID){RMIID <- c(RMIID,max)}
    c <- c+1
  }else if((right%>%str_detect("AGE") & !left%>%str_detect("AGE")) | (!right%>%str_detect("AGE") & left%>%str_detect("AGE")) ){
    age <- which(str_detect(c(right,left),"AGE"));noage <- which(!str_detect(c(right,left),"AGE"))
    if(!c(right,left)[age]%in%keepIID){keepIID <- c(keepIID,c(right,left)[age])};if(!c(right,left)[noage]%in%RMIID){RMIID <- c(RMIID,c(right,left)[noage])}
    c <- c+1
  }else{rowkeep <- c(rowkeep,r)}
  setTxtProgressBar(pb, c)
}
ibd_pi2 <- ibd_pi[rowkeep,]
for(r in 1:nrow(ibd_pi2)){
  right <- ibd_pi[r,]$IID1
  left <- ibd_pi[r,]$IID2
  if (left%in%RMIID){
    if(!right%in%c(RMIID,keepIID)){keepIID <- c(keepIID,right);c <- c+1
    }else {c <- c+1}
  }else if (!left%in%c(RMIID,keepIID)){
    if(!right%in%c(RMIID,keepIID)){
      rowimis <- ibd_pi_imis%>%filter(IDS%in%c(right,left))
      min <- rowimis%>%filter(F_MISS==min(F_MISS))%>%pull(IDS)%>%.[1]
      max <- rowimis%>%filter(F_MISS==max(F_MISS))%>%pull(IDS)%>%.[1]
      if(!min%in%keepIID){keepIID <- c(keepIID,min)};if(!max%in%RMIID){RMIID <- c(RMIID,max)};c <- c+1
    }else if (right%in%c(RMIID)){
      keepIID <- c(keepIID,left);c <- c+1
    }else if (right%in%c(keepIID)){
      rowimis <- ibd_pi_imis%>%filter(IDS%in%c(right,left))
      min <- rowimis%>%filter(F_MISS==min(F_MISS))%>%pull(IDS)%>%.[1]
      max <- rowimis%>%filter(F_MISS==max(F_MISS))%>%pull(IDS)%>%.[1]
      if(length(which(keepIID==max))!=0) {keepIID <- keepIID[-which(keepIID==max)]}
      if(!min%in%keepIID){keepIID <- c(keepIID,min)};if(!max%in%RMIID){RMIID <- c(RMIID,max)};c <- c+1
    }
  }else if (left%in%c(keepIID)){
    if(!right%in%c(RMIID,keepIID)){
      rrowimis <- ibd_pi_imis%>%filter(IDS%in%c(right,left))
      min <- rowimis%>%filter(F_MISS==min(F_MISS))%>%pull(IDS)%>%.[1]
      max <- rowimis%>%filter(F_MISS==max(F_MISS))%>%pull(IDS)%>%.[1]
      if(length(which(keepIID==max))!=0) {keepIID <- keepIID[-which(keepIID==max)]}
      if(!min%in%keepIID){keepIID <- c(keepIID,min)};if(!max%in%RMIID){RMIID <- c(RMIID,max)};c <- c+1
    }else if (right%in%c(RMIID)){c <- c+1
    }else if (right%in%c(keepIID)){
      rowimis <- ibd_pi_imis%>%filter(IDS%in%c(right,left))
      min <- rowimis%>%filter(F_MISS==min(F_MISS))%>%pull(IDS)%>%.[1]
      max <- rowimis%>%filter(F_MISS==max(F_MISS))%>%pull(IDS)%>%.[1]
      if(length(which(keepIID==max))!=0) {keepIID <- keepIID[-which(keepIID==max)]}
      if(!min%in%keepIID){keepIID <- c(keepIID,min)};if(!max%in%RMIID){RMIID <- c(RMIID,max)};c <- c+1
    }
  }
  setTxtProgressBar(pb, c)
}
table(duplicated(keepIID)) 
table(duplicated(RMIID))
Ind_to_remove <- ibd_pi_imis%>%filter(IDS%in%RMIID)%>%dplyr::select(IDS)%>%left_join(imiss%>%dplyr::select(FID,IID),by=c("IDS"="IID"))%>%dplyr::select(FID,IDS)
write.table(Ind_to_remove, paste0(results,"ibd_test_fail_ID"), row.names = F, col.names = F, sep = " ", quote = F)
file.remove(list.files(IBD,'GSAdataIBD_LDexclude',full.names = T),recursive=TRUE)
system(paste0("plink --bfile ",results,"GSAdata --remove ",results,"ibd_test_fail_ID --make-bed --out ",results,"GSAdata_temp --allow-no-sex")) # 698846 variants and 5596 people pass
tempstep(paste0(results,"GSAdata_temp"),paste0(results,"GSAdata"))


#########  TOTAL IDIVIDUAL QCed  ##########
#* 698846 Markers and 5596 people pass    #
#* 5596 people (2840 males, 2756 females) # 
#* 1295 are cases and 4301 are controls   #
###########################################


###### Split the data set ######
## According to Plink, the autosomes should be coded 1 through 22. The following other codes can be used to specify other chromosome types:
# X    X chromosome                    -> 23
# Y    Y chromosome                    -> 24
# XY   Pseudo-autosomal region of X    -> 25
# MT   Mitochondrial                   -> 26
## 1. Subsetting only mt markers
system(paste0("plink --bfile ",results,"GSAdata --chr 26 --make-bed --out ",subsets,"GSAdata_mt")) # 139 variants and 5596 people pass 
## 2. Including only sex chr (nonmale Y chr and .hh file markers)
system(paste0("plink --bfile ",results,"GSAdata --chr 23 24 --make-bed --out ",subsets,"GSAdata_sex")) # 19975 variants and 5596 people pass 


############ PER-MARKER QC ############


###### Identification of SNPs with an excessive missing genotype ######
system(paste0("plink --bfile ",results,"GSAdata --missing --out ",results,"GSAdata"))
lmiss <- read.table(paste0(results,"GSAdata.lmiss"),quote='',header=T,as.is=T)%>%as.tibble()
lmissfails <- lmissPlot(lmiss,th=0.05,outplot=paste0(figures,"05_SNPmissingnes.pdf"),outfiles=results) # Number of SNPs with high missigness (>0.05): 62503
system(paste0("plink --bfile ",results,"GSAdata --geno 0.05  --make-bed --out ",results,"GSAdata_temp")) # 636343 variants and 5596 people pass
tempstep(paste0(results,"GSAdata_temp"),paste0(results,"GSAdata"))


###### Identification of SNPs showing a significant deviation from HWE  ######
system(paste0("plink --bfile ",results,"GSAdata --hardy --filter-controls --out ",results,"GSAdata"))
hwe <- read.table(paste0(results,'GSAdata.hwe'),quote='',header=T,as.is=T)%>%as.tibble()
fam <- read.table(paste0(results,'GSAdata.fam'),quote='',header=F,as.is=T)%>%as.tibble()
hweFails <- hwePlot(hwe,pvalth=10^-5,outplot=paste0(figures,"06_SNP_HWEfail_ctrl.png"),outfiles=results,numof.cases=fam%>%dplyr::filter(V6==2)%>%dim()%>%.[1],numof.controls=fam%>%dplyr::filter(V6==1)%>%dim()%>%.[1]) # 1910 SNPs failed HWE in controls
system(paste0("plink --bfile ",results,"GSAdata --exclude ",results,"SNP.hweCRTLfail --make-bed --out ",results,"GSAdata_temp")) # 634433 variants and 5596 people pass 
tempstep(paste0(results,"GSAdata_temp"),paste0(results,"GSAdata"))


###### Identification of SNPs with significantly different missing genotype rates btw cases/controls ######
system(paste0("plink --bfile ",results,"GSAdata_QC_1603 --test-missing --out ",results,"GSAdata"))  
inform_missing <- read.table(paste0(results,'GSAdata.missing'),quote='',header=T,as.is=T)%>%as.tibble()
fail_diffmiss_test <- data.frame(inform_missing$SNP[inform_missing$P < 0.00001]) 
write.table(fail_diffmiss_test, paste0(results,"fail_diffmiss_test"), row.names = F, col.names = F, sep = " ", quote = F) 
system(paste0("plink --bfile ",results,"GSAdata --exclude ",results,"fail_diffmiss_test --make-bed --out ",results,"GSAdata_temp")) # 626658 variants and 5596 people pass
tempstep(paste0(results,"GSAdata_temp"),paste0(results,"GSAdata_QCed_3003"))


##### Removal of sex chromosomal markers ######
system(paste0("plink --bfile ",results,"GSAdata_QCed_3003 --autosome --make-bed --out ",results,"GSAdata_QCed_3003nosex")) # 610941 variants and 5596 people pass


###########  TOTAL MARKER QCed  ###########
#* GSAdata_QCed_1603                      #            
#* 610941 variants                        #
#* 5596 people (2840 males, 2756 females) #          
#* 1295 are cases and 4301 are controls   #     
###########################################


###### For association: Removal of all markers with a very low MAF ######
system(paste0("plink --bfile ",results,"GSAdata_QCed_3003nosex --maf 0.01 --make-bed --out ",results,"GSAdata_QCed_3003_assoc")) # 447285  variants

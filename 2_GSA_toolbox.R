##*************************************************##
## Global Screening Array Illumina (GSA) Toolbox ####
##*************************************************##

#==================================================#
########---- load/install packages ----#############
#==================================================#
packages <- function(requirements,quiet=FALSE){
  has   <- requirements %in% rownames(installed.packages())
  if(any(!has)){
    message("Installing packages...")
    setRepositories(ind=c(1:5))
    #options(install.packages.check.source = "no")
    install.packages(requirements[!has],repos="https://cran.uni-muenster.de/")
  }
  if(quiet){
    for(r in requirements){suppressMessages(require(r,character.only=TRUE))}
  }else for(r in requirements){message(paste(r,suppressMessages(require(r,character.only=TRUE)),sep=': '))}
}

#=====================================================#
########---- SNP annotation functions ----#############
#=====================================================#

#library(BSgenome.Hsapiens.UCSC.hg19);library(BSgenome.Hsapiens.UCSC.hg38)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)#,'biomaRt','tidyverse','TxDb.Hsapiens.UCSC.hg19.knownGene','liftOver','SNPlocs.Hsapiens.dbSNP151.GRCh38'
snpdf <- function(GRobject,prefix=c('hg18','hg19','hg38')){
  if(prefix=='hg19')prefix=''
  if (class(GRobject)=="UnstitchedGPos"){
    df <- tibble(RS_ID=mcols(GRobject)$RefSNP_id,Chr=seqnames(GRobject)%>%as.numeric(),Pos=start(GRobject)%>%as.numeric())%>%mutate(Chr_Pos=paste0(prefix,'chr',Chr,':',Pos))
    return(df)
  }
}
getSNPsPosH19 <- function(rsnums){
  packages(c('tidyverse','SNPlocs.Hsapiens.dbSNP151.GRCh38','SNPlocs.Hsapiens.dbSNP144.GRCh37','SNPlocs.Hsapiens.dbSNP.20120608'))
  #  check if the rs are in the genome version hg19
  snps37 <- SNPlocs.Hsapiens.dbSNP144.GRCh37 # snp info of hg19
  rsnums <- rsnums[!duplicated(rsnums)]
  snps <- snpsById(snps37,rsnums,ifnotfound="drop")
  snpsdf <- snpdf(snps,prefix='hg19')
  XO <- setdiff(rsnums,snpsdf$RS_ID)
  #if not found are in genome version 18 ?:
  snps18 <- SNPlocs.Hsapiens.dbSNP.20120608 # snp info of hg18
  snpsx <- snpsById(snps18,XO,ifnotfound="drop")
  snpsdf <- bind_rows(snpsdf,snpdf(snpsx,prefix='hg18'))
  XO <- setdiff(rsnums,snpsdf$RS_ID)
  #if not found are in genome version 38 ?:
  snps38 <- SNPlocs.Hsapiens.dbSNP151.GRCh38 # snp info of hg38
  snpsxx <- snpsById(snps18,XO,ifnotfound="drop")
  snpsdf <- bind_rows(snpsdf,snpdf(snpsxx,prefix='hg38'))
  XO <- setdiff(rsnums,snpsdf$RS_ID)
  #if not found are in genome version 38 ?:
  if (length(XO)!=0){
    try({
    snpdbmanual <- read_tsv('docs/snpsEnsembl37.txt')%>%filter(!str_detect(`Chromosome name`,"^H.+"),progress=F)
    snpdbmanual <- snpdbmanual%>%mutate(`Chromosome name`=lapply(snpdbmanual[,3],function(x) ifelse(!str_detect(x,"[XYM]"),x,ifelse(x=='Y',24,ifelse(x=='X',23,ifelse(x=='M',26,0)))))%>%unlist()%>%as.numeric())
    dbsnp <- tibble(RS_ID=snpdbmanual[,1,drop=T],Chr=snpdbmanual[,3,drop=T],Pos=snpdbmanual[,4,drop=T])
    dbsnp <- dbsnp%>%mutate(Chr_Pos=paste0('chr',Chr,':',Pos))%>%filter(RS_ID%in%XO)
    snpsdf <- bind_rows(snpsdf,dbsnp)
    XO <- setdiff(rsnums,snpsdf$RS_ID)
    },silent=T)
  }
  # Rs with chr 0 pos 0 are still reported for manual curation
  XOdf <- tibble(RS_ID=XO,Chr=0,Pos=0)%>%mutate(Chr_Pos=paste0('chr',Chr,':',Pos))
  snpsdf <- bind_rows(snpsdf,XOdf)
  return(snpsdf)
}
chrposDF <- function(famslice,prefix="",splited=':'){
  #chr24=Y; chr23=X; chr26=M (mitochondria)
  v2 <- famslice%>%pull(V2)
  chr <- lapply(v2, function(x){
    c <- str_remove(x,prefix)%>%str_split(splited)%>%unlist()%>%.[1]
    ifelse(c=='X',c <- 23,ifelse(c=='Y',c <- 24,ifelse(c=='M',c <- 26,c <- as.numeric(c) )))
    })%>%unlist
  pos <- lapply(v2, function(x) str_remove(x,prefix)%>%str_split(splited)%>%unlist()%>%.[2]%>%as.numeric)%>%unlist
  cp <- paste0('chr',lapply(v2, function(x) str_remove(x,prefix))%>%unlist)%>%str_replace(splited,':')
  df <- tibble(ID=v2,RS_ID='0',Chr=chr,Pos=pos,Chr_Pos=cp)
  return(df)
}
chrposDFv2 <- function(famslice){
  v2 <- famslice%>%pull(V2)
  rs <- lapply(v2, function(x) str_extract(x,"rs[[:digit:]]+"))%>%unlist
  rs[is.na(rs)] <- '0'
  cp <- paste0('chr',famslice%>%pull(V1),':',famslice%>%pull(V4))
  df <- tibble(ID=v2,RS_ID=rs,Chr=famslice%>%pull(V1)%>%as.numeric,Pos=famslice%>%pull(V4)%>%as.numeric,Chr_Pos=cp)
}
#snpDF <- DFtot
snpAnnotation <- function(snpDF){
  packages(c('AnnotationHub','TxDb.Hsapiens.UCSC.hg19.knownGene'))
  hub <- AnnotationHub()
  unique(hub$dataprovider)
  q <- query(hub,c('GRanges','Homo sapiens','hg19','UCSC track'))
  q$sourcetype
  q$title
  res <- q[[103]]
  
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  methods(class='TxDb')
  txs <- transcripts(txdb)
  head(seqlevels(txdb))
  gr <- GRangesFilter(GenomicRanges::GRanges("chr1:44000000-55000000"))
  transcripts(txdb19, filter=~(symbol %startsWith% "SNORD") | symbol == "ADA")
}
#=========================================#
########---- QC functions ----#############
#=========================================#
ancestryMarkers <- function(){
  ancestry <- read.table(paste0(rootfolder,'/data/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps_hg19.bim'),quote='',header=F,as.is=T)%>%as.tibble()
  ExomeChip_Annotation <- read.table(paste0(rootfolder,'/bin/ExomeChip_SNPsInfo_ann.txt'),quote='',header=T,as.is=T,sep="\t")%>%as.tibble()
  aMarkers <- ExomeChip_Annotation%>%dplyr::select(exmID,Chr.position)%>%inner_join(ancestry,by=c('Chr.position'='V2'))%>%dplyr::select(exmID)
  write.table(aMarkers,paste0(rootfolder,'/data/ancestryMarkers_ExomeChip.txt'),quote=F,sep="\t",col.names=F,row.names=F)
}


genderPlot <- function(sexcheck,outplot,outfiles){
  packages('tidyverse')
  if (length(which(sexcheck$PEDSEX==1 | sexcheck$PEDSEX==2))<=3) {
    cat(paste("There are not enough male(1) or female(2) in PEDSEX of file ",dataFile,"; Will not generate the figure\n",sep=""))
  } else {
    
    Fth <- sexcheck%>%dplyr::filter(PEDSEX==2 & F>=0 &STATUS!='OK')%>%summarise(th=min(F))%>%pull(th)
    Mth <- sexcheck%>%dplyr::filter(PEDSEX==1 & STATUS!='OK')%>%summarise(th=max(F))%>%pull(th)
    hf <- hist(sexcheck%>%dplyr::filter(PEDSEX==2)%>%pull(F), breaks=100,main="Female",xlab="Chr X inbreeding estimate")
    fh = cut(hf$breaks, c(-Inf,Fth, Inf))
    hm <- hist(sexcheck%>%dplyr::filter(PEDSEX==1)%>%pull(F),breaks =100,main="Male",xlab="Chr X inbreeding estimate")
    mh = cut(hm$breaks, c(-Inf,Fth, Inf))
    pdf(outplot,useDingbats=F,height=7,width=10,bg='transparent')
    par(mfcol=c(1,2), pty="m")
    plot(hf,main="Female",xlab="Chr X inbreeding estimate",col=c("grey95","red")[fh])
    abline(v=Fth,col='red',lty=3)
    plot(hm,main="Male",xlab="Chr X inbreeding estimate",col=c("red","grey95")[mh])
    abline(v=Mth,col='red',lty=3)
    dev.off()
    sexproblem <- sexcheck %>% filter(STATUS == "PROBLEM")
    write.table(sexproblem, paste0(outfiles,"Indiv.sexproblem"), row.names = F, col.names = F, sep = " ", quote = F) 
    cat(paste(outplot," was successfully generated.\n",sep=""))
  }
}
#
imiss_hetPlot <- function(imiss,het,outplot,outfiles,pmiss=0.01,timesSD=2){
  if(!is.numeric(pmiss)) stop(paste0('pmiss (portion of missing rates) is not numeric'))
  if(!is.numeric(timesSD)) stop(paste0('timesSD (times of sd threshold) is not numeric'))
  # timesSD = int or float: times of sd threshold default 3.2 #
  # pmiss = float: portion of missing rates default 0.01 #
  library(tidyverse)
  het <- het%>%dplyr::mutate(hetRate=round((N.NM.-O.HOM.)/N.NM., 5))
  imiss <- imiss%>%dplyr::mutate(F_MISSlog=log10(F_MISS),F_MISScol=log10(1+F_MISS))
  colors  <- densCols(imiss%>%pull(F_MISScol),het$hetRate)
  colors2  <- densCols(imiss%>%pull(F_MISScol),het$hetRate,colramp=colorRampPalette(c("red","white")))
  i <- imiss%>%dplyr::filter(F_MISSlog>=log10(pmiss))%>%pull(IID)
  h <- het%>%dplyr::filter(hetRate>=(mean(hetRate)+(timesSD*sd(hetRate))) | hetRate<=(mean(hetRate)-(timesSD*sd(hetRate))))
  # Low heterozygosity may indicate inbreeding, which can lead to reduce fitness of a population with many homozygous genotypes, and high heterozygosity may indicate contamination. 
  failed <- het%>%dplyr::filter(IID%in%c(h$IID,i))%>%dplyr::inner_join(imiss%>%dplyr::select(IID,F_MISS),by=c('IID'='IID'))
 if(imiss%>%dplyr::filter(N_MISS==0)%>%nrow!=0){maxSNVs <- imiss%>%dplyr::filter(N_MISS==0)%>%pull(N_GENO)%>%max}
  
  pdf(outplot, width=10, height=5,useDingbats=F,bg='transparent')
  par(mfcol=c(1,1), pty="m")
  plot(imiss%>%pull(F_MISSlog),het%>%pull(hetRate),col=ifelse(het$hetRate %in% h$hetRate,colors2,ifelse(imiss%>%pull(IID) %in% i,colors2,colors)),pch=20,xlab="Proportion of missing genotypes", ylab="Heterozygosity rate",xaxt='n',frame=FALSE,xlim=c(-3,0),ylim=c(mean(het$hetRate)-(10*sd(het$hetRate)),mean(het$hetRate)+(10*sd(het$hetRate))))#,axes=F
  axis(1,at=c(-3,-2,-1,0),labels=c(0.001,0.01,0.1,1))
  #axis(2,at=c(0.22,0.24,0.22,0.26,0.28,0.3,0.32),tick=T)
  #axis(1,at=c(-3,-2,-1,0),labels=c(0.001,0.01,0.1,1))
  #axis(2,at=c(0.05,0.06,0.07,0.08,0.09,0.1),tick=T)
  abline(h=mean(het$hetRate)-(timesSD*sd(het$hetRate)),col="RED",lty=2)
  abline(h=mean(het$hetRate)+(timesSD*sd(het$hetRate)),col="RED",lty=2)
  abline(v=log10(pmiss), col="RED", lty=2)
  result<-dev.off()
  
  write.table(failed,paste0(outfiles,'fail.imisshet'), sep="\t", quote=F, row.names =F)
  write.table(failed[,c(1,2)],paste0(outfiles,"Indiv.imisshetfail"),sep="\t", quote=F, row.names =F, col.names=F)
  cat(paste(outplot," and failed list was successfully generated.\n",NROW(failed)," individuals failed heterozygosity and missigness\n",sep=""))
  return(failed)
}
#eigenvec <- eigenvecGSA
pcaPlot <- function(eigenvec,pcs=c('PC1-PC2'),outplot,outfiles,lofth=1.7,getOutilers=F,pop='Germans'){
  library(RColorBrewer)
  library(tidyverse)
  sampleColors = colorRampPalette(brewer.pal(8, "Set1"))(length(levels(as.factor(eigenvec$population))))
  for (i in pcs){
    j <- i%>%str_split('-')%>%unlist();x_pc <- j[1];y_pc <- j[2]
    ggplot(eigenvec,aes(x=eigenvec[[x_pc]],y=eigenvec[[y_pc]],color=population))+
      geom_point(aes(shape=CaseControl, color=population))+
      scale_shape_manual(values=c(17,16,18,19,20,21),name='Case/Control')+
      scale_color_manual(name="Population",values=sampleColors)+
      xlab(paste0(x_pc))+ylab(paste0(y_pc))+ggtitle("Population stratification")+
      geom_hline(aes(yintercept=0),linetype="dashed",alpha=1/4)+geom_vline(aes(xintercept=0),linetype="dashed",alpha=1/4)+
      theme_minimal()#+geom_vline(xintercept=mean(evecs$PC1)+6*sd(evecs$PC1))+geom_hline(yintercept=mean(evecs$PC2)+6*sd(evecs$PC2))
    ggsave(paste0(outplot,i,'.pdf'),width=8,height=6,device="pdf",useDingbats=FALSE)
    eigenvec2 <- eigenvec%>%filter(population%in%pop)
    sampleColors2 = brewer.pal(8, "Set1")[1:(length(levels(as.factor(eigenvec$CaseControl))))]
    ggplot(eigenvec2,aes(x=eigenvec2[[x_pc]],y=eigenvec2[[y_pc]],color=CaseControl))+
      geom_point(aes(shape=population, color=CaseControl))+
      scale_shape_manual(values=c(16,17,18,19,20,21),name='Population')+
      scale_color_manual(name="Case/Control",values=sampleColors2)+
      xlab(paste0(x_pc))+ylab(paste0(y_pc))+ggtitle("Population stratification")+
      geom_hline(aes(yintercept=0),linetype="dashed",alpha=1/4)+geom_vline(aes(xintercept=0),linetype="dashed",alpha=1/4)+
      theme_minimal()#+geom_vline(xintercept=mean(evecs$PC1)+6*sd(evecs$PC1))+geom_hline(yintercept=mean(evecs$PC2)+6*sd(evecs$PC2))
    ggsave(paste0(outplot,i,'_',pop,'.pdf'),width=8,height=6,device="pdf",useDingbats=FALSE)
  }
  if(getOutilers){
    library(Rlof)
    evecs <- eigenvec2%>%dplyr::select(-FID.x,-IID,-population)%>%.[,1:3]
    lof_factors <- round(lof(evecs,nrow(eigenvec)/2,cores=3),1)
    outliers <- which(lof_factors>=lofth)
    cols <- rep("grey",length(lof_factors))
    cols[outliers]='red'
    pdf(paste0(outplot,'_',pop,'_LocalOutliers.pdf'), width=10, height=5,useDingbats=F,bg='transparent')
    plot(lof_factors,xlab='Sample',ylab='Local outlier factor',pch=19,col=cols,frame=FALSE)
    abline(h=lofth,col='red',lty=2)
    dev.off()
    evecs_plot <- eigenvec2%>%dplyr::select(c(2:12),population,CaseControl)%>%mutate(lof=lof_factors,lof_col=ifelse(lof_factors>=lofth,'Outlier','None'))
    for (i in pcs){
      j <- i%>%str_split('-')%>%unlist();x_pc <- j[1];y_pc <- j[2]
      ggplot(evecs_plot,aes(x=evecs_plot[[x_pc]],y=evecs_plot[[y_pc]]))+
        geom_point(aes(shape=CaseControl, color=lof_col))+
        scale_shape_manual(values=c(17,16,18,19,20,21),name='Case/Control')+
        scale_color_manual(name="Outilers",values=c('grey','red'))+
        xlab(paste0(x_pc))+ylab(paste0(y_pc))+ggtitle("Population stratification")+
        geom_hline(aes(yintercept=0),linetype="dashed",alpha=1/4)+geom_vline(aes(xintercept=0),linetype="dashed",alpha=1/4)+
        theme_minimal()#+geom_vline(xintercept=mean(evecs$PC1)+6*sd(evecs$PC1))+geom_hline(yintercept=mean(evecs$PC2)+6*sd(evecs$PC2))
      ggsave(paste0(outplot,i,'_',pop,'_lofPCA.pdf'),width=8,height=6,device="pdf",useDingbats=FALSE)
    }
    fails <- eigenvec[outliers,]
    evecs_plot2 <- eigenvec%>%filter(population%in%c(pop,'CEU') & !IID%in%fails$IID)
    for (i in pcs){
      j <- i%>%str_split('-')%>%unlist();x_pc <- j[1];y_pc <- j[2]
      sampleColors2 = brewer.pal(8, "Set1")[1:(length(levels(as.factor(evecs_plot2$population))))]
      ggplot(evecs_plot2,aes(x=evecs_plot2[[x_pc]],y=evecs_plot2[[y_pc]],color=population))+
        geom_point(aes(shape=CaseControl, color=population))+
        scale_shape_manual(values=c(17,16,18,19,20,21),name='Case/Control')+
        scale_color_manual(name="Population",values=sampleColors)+
        xlab(paste0(x_pc))+ylab(paste0(y_pc))+ggtitle("Population stratification")+
        geom_hline(aes(yintercept=0),linetype="dashed",alpha=1/4)+geom_vline(aes(xintercept=0),linetype="dashed",alpha=1/4)+
        theme_minimal()#+geom_vline(xintercept=mean(evecs$PC1)+6*sd(evecs$PC1))+geom_hline(yintercept=mean(evecs$PC2)+6*sd(evecs$PC2))
      ggsave(paste0(outplot,i,'_',pop,'_lofPCA_outrmPop.pdf'),width=8,height=6,device="pdf",useDingbats=FALSE)
      sampleColors2 = brewer.pal(8, "Set1")[1:(length(levels(as.factor(evecs_plot2$CaseControl))))]
      ggplot(evecs_plot2,aes(x=evecs_plot2[[x_pc]],y=evecs_plot2[[y_pc]],color=CaseControl))+
        geom_point(aes(shape=population, color=CaseControl))+
        scale_shape_manual(values=c(17,16,18,19,20,21),name='Population')+
        scale_color_manual(name="Case/Control",values=sampleColors2)+
        xlab(paste0(x_pc))+ylab(paste0(y_pc))+ggtitle("Population stratification")+
        geom_hline(aes(yintercept=0),linetype="dashed",alpha=1/4)+geom_vline(aes(xintercept=0),linetype="dashed",alpha=1/4)+
        theme_minimal()#+geom_vline(xintercept=mean(evecs$PC1)+6*sd(evecs$PC1))+geom_hline(yintercept=mean(evecs$PC2)+6*sd(evecs$PC2))
      ggsave(paste0(outplot,i,'_',pop,'_lofPCA_outrm_CC.pdf'),width=8,height=6,device="pdf",useDingbats=FALSE)
    }
    message(paste0('Number of population stratification failures: ',nrow(fails)))
    write.table(fails[,c(1,2)],paste0(outfiles,"Indiv.PopStratFail"),sep="\t", quote=F, row.names =F, col.names=F)
  }
}
#
lmissPlot <- function(lmiss,th=0.03,outplot,outfiles){
  library(tidyverse)
  h <- hist(lmiss%>%dplyr::filter(F_MISS>=th-0.02)%>%pull(F_MISS),xlim=c(th-0.02,round(max(lmiss$F_MISS),1)),breaks=100,main=paste0("SNVs missing rate > ",th),ylab="Number of SNVs",xlab="Fraction of missing data")
  hcuts = cut(h$breaks, c(-Inf,th,Inf))
  pdf(outplot, width=10, height=5,useDingbats=F,bg='transparent')
  plot(h,main=paste0("SNPs missing rate > ",th),ylab="Number of SNPs",xlab="Fraction of missing data",xlim=c(th-0.02,1),col=c("grey95","red")[hcuts])
  #abline(v=th+0.01,col="RED",lty=2)
  dev.off()
  ggplot()
  fail <- lmiss%>%dplyr::filter(F_MISS>=th)
  message(paste0('Number of SNP with high missigness (>',th,'): ',nrow(fail)))
  write.table(fail,paste0(outfiles,'fail.lmiss'), sep="\t", quote=F, row.names =F)
  write.table(fail$SNP,paste0(outfiles,"SNP.lmissfail"),sep="\t", quote=F, row.names =F, col.names=F)
  cat(paste(outplot," and failed list was successfully generated.\n",NROW(fail)," SNPs failed missigness\n",sep=""))
  return(fail)
}
#
hwePlot <- function(hwe,pvalth=10^-5,outplot,outfiles,numof.cases,numof.controls){
  library(tidyverse)
  #conservative P-value threshold like P < 1x10-5
  hwe_cont <- hwe%>%dplyr::filter(TEST=='UNAFF')%>%dplyr::mutate(E_Pval=ppoints(P))
  #pvalth <- round(5/dim(hwe_cont)[1],5)
  hwe_cont_fail <- hwe_cont%>%dplyr::filter(P<pvalth)
  hwe_fail <- hwe_cont_fail%>%
    dplyr::full_join(hwe%>%dplyr::filter(TEST=='AFF' & SNP%in%hwe_cont_fail$SNP)%>%
                       dplyr::select(SNP,TEST,GENO,O.HET., E.HET.,P)%>%dplyr::mutate(failed=ifelse(P<pvalth,'y','n')),by=c('SNP'),suffix=c(".U",".A"))
  
  hwe_cont_filtered <- hwe_cont%>%dplyr::filter(!SNP%in%hwe_fail$SNP)
  # Inflation coeficient
  q <- qchisq(hwe_cont_filtered$P, df=1, ncp=0, lower.tail = FALSE, log.p = FALSE)
  lambda = round(median(q)/qchisq(0.5,1),2) # The median of a chi-squared distribution with one degree of freedom is 0.4549364
  lambda_1000 <- round(1.0 + (lambda - 1.0) * ( (1/numof.controls)/(1/1000) ),2)
  par(mfcol=c(1,1),mar= c(5, 4, 4, 2) + 0.1, pty="m")
#  graph <- ggplot(hwe_cont_filtered,aes(x=hwe_cont_filtered%>%pull(E_Pval)%>%-log10(.)%>%sort,y=hwe_cont_filtered%>%pull(P)%>%-log10(.)%>%sort))+geom_point()+
#    ylab(expression(Observed~~-log[10](italic(p))))+xlab(expression(Expected~~-log[10](italic(p))))+theme_minimal()+
#    ggtitle(paste0("HWE-test QQplot - Controls; Lambda: ",lambda,"; lambda_1000: ",lambda_1000,'; No.Pval: ',dim(hwe_cont_filtered)[1]))+
#    geom_abline(aes(intercept=0,slope=1),linetype="dashed",alpha=0.75,colour="red")#+geom_vline(aes(xintercept=0),linetype="dashed",alpha=1/4)+
#  ggsave(graph,paste0(outplot),width=8,height=6,device="png")
  
  write.table(hwe_fail,paste0(outfiles,'SNPhweCRTLfail.txt'), sep="\t", quote=F, row.names =F)
  write.table(hwe_fail$SNP,paste0(outfiles,"SNP.hweCRTLfail"),sep="\t", quote=F, row.names =F, col.names=F)
  cat(paste(outplot," and failed list was successfully generated.\n",NROW(hwe_fail)," SNPs failed HWE in controls\n",sep=""))
  return(hwe_fail)
  #png(outplot, width=10, height=5,res=150)
  #plot(hwe_cont_filtered%>%pull(E_Pval)%>%-log10(.)%>%sort,hwe_cont_filtered%>%pull(P)%>%-log10(.)%>%sort,ylab=expression(Observed~~-log[10](italic(p))),
  #     xlab=expression(Expected~~-log[10](italic(p))),
  #     main=paste0("HWE-test QQplot - Controls; Lambda: ",lambda,"; lambda_1000: ",lambda_1000,'; No.Pval: ',dim(hwe_cont_filtered)[1]),pch=20,frame=FALSE,
  #     xlim=c(0,round(hwe_cont_filtered%>%pull(E_Pval)%>%-log10(.)%>%max)),ylim=c(0,round(hwe_cont_filtered%>%pull(P)%>%-log10(.)%>%max)))
  #with(hwe_cont_filtered%>%dplyr::filter(P<=round(5/dim(hwe_cont)[1],5)),text(P~E_Pval, labels=SNP, cex=0.6, pos=1, col="red"))
  #abline(a=0,b=1,lty=2,col='red')
  #abline(h=(-log10(pvalth)),col='blue',lty=2)
  #dev.off()
  
}


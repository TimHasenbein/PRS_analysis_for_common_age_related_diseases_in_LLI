##*****************************************************##
## Global Screening Array Illumina (GSA) PRS analysis  ##
##*****************************************************##
## Tim Hasenbein
## PhD. Guillermo G. Torres
## Last modification 06.2020
## PRS calculation using published scores and PLINK: Here for coronary artery disease
library(data.table)
library(tidyverse)


###### set up env. var ######
mypath <- if(class(try(rstudioapi::getSourceEditorContext()$path,silent=TRUE))=="try-error") {paste0(dirname(substring(argv[grep("--file=", argv)], 8)),'/')} else paste(head(strsplit(rstudioapi::getSourceEditorContext()$path,'/')[[1]],-1),collapse='/')
setwd(mypath)
imp_data <- "F:/Imputation/Post_imp/filtered_recoded/results/"


##### Step 1: Preprare basedata #####
base_data <- fread("CoronaryArteryDisease_PRS_LDpred_rho0.001_v3.txt", header = T)
head(base_data) # 6630150 variants in base data
base_data$variant <- paste0("chr",base_data$chr,":",base_data$position_hg19) 
base_data <- base_data%>%select("variant", "effect_allele", "effect_weight") 
write.table(base_data, "basedata_cad", row.names = F, col.names = F, sep = " ", quote = F)


##### Step 2: Calculate PRS #####
system(paste0("plink --bfile ",imp_data,"GSA_imp_PRS --score basedata_cad sum --out cad")) # 5920526 variants


##### Step 3: Standardize polygenic scores #####
data <-fread("cad.profile", header=T)
head(data)
data$SCORESUM=(data$SCORESUM-mean(data$SCORESUM))/sd(data$SCORESUM, na.rm=T)
data <- data%>%rename("PRS" = "SCORESUM")
LLI <- data%>%filter(str_detect(IID,"AGE"))
Ctrl <- data%>%filter(str_detect(IID,"AGE", negate = T))
LLI$PHENO <- 'LLI'
Ctrl$PHENO <- 'Ctrl'
data <- rbind(LLI, Ctrl)


##### Step 4: Visualize data #####
## Histogram
ggplot(data, aes(x=PRS, color=PHENO)) +
  geom_histogram(fill="white", alpha=0.5, position="identity") +
  labs(title="Coronary Artery Disease") + 
  xlab("Polygenic Risk Score")
ggsave(paste0(mypath, '/cad_histo.pdf'))
mu <- data%>%group_by(PHENO)%>%summarise(M=mean(PRS),SD=sd(PRS),Med=median(PRS))
## Boxplot
ggplot(data%>%mutate(PRS=(PRS-mean(PRS))), aes(x=factor(PHENO),y=PRS,fill=factor(PHENO)))+
  stat_boxplot(geom='errorbar', linetype=1, width=0.5)+ 
  geom_boxplot(outlier.shape=1)+    
  stat_summary(fun.y=mean, geom="point", size=2) 
ggsave(paste0(mypath, '/cad_boxplot.pdf'))
## Density plot  
ggplot(data,aes(x=PRS,fill=factor(PHENO),color=factor(PHENO)))+geom_density(alpha=0.1)+
  geom_vline(data=mu,aes(xintercept=M, color=factor(PHENO)),linetype="dashed")+theme_classic()+xlab('Estimated CAD-PRS')+
  scale_color_discrete(name = "", labels = c("Young controls", "LLI"))+scale_fill_discrete(name = "", labels = c("Young controls", "LLI"))
ggsave(paste0(mypath,'/cad_densityplot.pdf'))
## Point plot 
ggplot(data%>%arrange(PRS)%>%rownames_to_column('Pos')%>%mutate(Pos=ifelse(PHENO==2,as.numeric(Pos)+50,as.numeric(Pos))),aes(x=Pos,y=PRS,color=factor(PHENO)))+geom_point(alpha=0.3)+theme_classic()+ylab('Estimated CAD-PRS')+xlab('Individuals')+scale_color_discrete(name = "", labels = c("Young controls", "LLI"))
ggsave(paste0(mypath,'/cad_points.pdf'))


##### Step 5: Statistical analysis #####
# Shapiro-Wilk normality test
shapiro.test(LLI$PRS) 
shapiro.test(Ctrl$PRS) 
# Barlett's test 
bartlett.test(PRS~PHENO, data) 
# T-test for signifant differences between LLI and younger PRS mean
t.test(PRS~PHENO, data = data, var.equal = TRUE) 



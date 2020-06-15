#**********************************************************************************##
## Global Screening Array Illumina (GSA) longevity-PGS using published statistics  ##
##*********************************************************************************##
## Tim Hasenbein
## PhD. Guillermo G. Torres
## Last modification 04.2020
## Generate a PGS for longevity using the summary statistics by Deelen et al. 2019: A meta-analysis of genome-wide association studies identifies multiple longevity genes
library(data.table)
library(tidyverse)
library(pROC)
library(rcompanion)


###### set up env. var ######
mypath <- if(class(try(rstudioapi::getSourceEditorContext()$path,silent=TRUE))=="try-error") {paste0(dirname(substring(argv[grep("--file=", argv)], 8)),'/')} else paste(head(strsplit(rstudioapi::getSourceEditorContext()$path,'/')[[1]],-1),collapse='/')
setwd(mypath)
imp_data <- "F:/Imputation/Post_imp/filtered_recoded/results/"


########## PRS calculation using PRSice-2 ##########


##### Step 1: Obtain GWAS summary stats + preprare base data set #####
base_data <- fread("Results_99th_percentile.txt", header = T)
head(base_data)
base_data$SNP <- paste0("chr",base_data$Chr,":",base_data$Position)
base_data_prsice <- base_data%>%select("SNP", "EA", "Beta", "P-value")
head(base_data_prsice)
write.table(base_data_prsice, "basedata_PRSice", row.names = F, col.names = T, sep = " ", quote = F)


##### Step 2: Run PRS calculation #####
system("Rscript PRSice.R --dir . \
--prsice ./PRSice_linux \
--base basedata_PRSice \
--target GSA_imp_PRS \
--snp SNP \
--A1 EA \
--stat Beta \
--pvalue P-value \
--binary-target T \
--bar-levels 5e-08,5e-07,5e-06,5e-05,5e-04,5e-03,5e-02,5e-01 \
--fastscore \
--no-regress \
--out PGS_final") 


##### Step 3: Standardize polygenic scores #####
data_prsice <-fread("PGS_final.all.score", header=T) 
head(data_prsice)
data_prsice <- data_prsice%>%rename("PRS" = "0.05")
data_prsice <- data_prsice%>%select("FID", "IID", "PRS")
data_prsice$PRS=(data_prsice$PRS-mean(data_prsice$PRS))/sd(data_prsice$PRS, na.rm=T)
LLI <- data_prsice%>%filter(str_detect(IID,"AGE"))
Ctrl <- data_prsice%>%filter(str_detect(IID,"AGE", negate = T))
LLI$PHENO <- 'LLI'
Ctrl$PHENO <- 'Ctrl'
data_prsice <- rbind(LLI, Ctrl)
head(data_prsice)


##### Step 4 Visualize data #####
# Histogram
ggplot(data_prsice, aes(x=PRS, color=PHENO)) +
  geom_histogram(fill="white", alpha=0.5, position="identity") +
  labs(title="PGS for longevity  (PRSice-2)", x = "PGS")
ggsave(paste0('final_PRSicePaper_histoplot.pdf'))
mu <- data_prsice%>%group_by(PHENO)%>%summarise(M=mean(PRS),SD=sd(PRS),Med=median(PRS))
# Boxplot
ggplot(data_prsice%>%mutate(PRS=(PRS-mean(PRS))), aes(x=factor(PHENO),y=PRS,fill=factor(PHENO)))+
  stat_boxplot(geom='errorbar', linetype=1, width=0.5)+ 
  geom_boxplot(outlier.shape=1)+    
  stat_summary(fun.y=mean, geom="point", size=2) +
  labs(y = "PGS", x = "Phenotype")
ggsave(paste0('final_PRSicePaper_boxplot.pdf'))
# Density plot  
ggplot(data_prsice,aes(x=PRS,fill=factor(PHENO),color=factor(PHENO)))+geom_density(alpha=0.1)+
  geom_vline(data=mu,aes(xintercept=M, color=factor(PHENO)),linetype="dashed")+theme_classic()+xlab('PGS')+
  scale_color_discrete(name = "", labels = c("Young controls", "LLI"))+scale_fill_discrete(name = "", labels = c("Young controls", "LLI"))
ggsave(paste0('final_PRSicePaper_densityplot.pdf'))
# Linear plot
ggplot(data_prsice%>%arrange(PRS)%>%rownames_to_column('Pos')%>%mutate(Pos=ifelse(PHENO==2,as.numeric(Pos)+50,as.numeric(Pos))),aes(x=Pos,y=PRS,color=factor(PHENO)))+geom_point(alpha=0.3)+theme_classic()+ylab('Estimated Longevity-PRS (PRSice-2')+xlab('Individuals')+scale_color_discrete(name = "", labels = c("Young controls", "LLI"))
ggsave(paste0(mypath,'final_PRSicePaper_points.pdf'))


##### Step 5: Statistics #####
# Shapiro-Wilk normality test 
shapiro.test(LLI$PRS)  
shapiro.test(Ctrl$PRS) 
# Barlett's test 
bartlett.test(PRS~PHENO, data_prsice) 
# Two Sample t-test for signifant differences between LLI and younger PRS mean 
wilcox.test(PRS~PHENO, data = data_prsice, var.equal = TRUE) 


##### Step 6: Logistic regression model #####
mod1<-glm(as.factor(PHENO)~PRS, data=data_prsice, family = binomial)
summary(mod1)
nagelkerke(mod1) 
auc(data_prsice$PHENO, data_prsice$PRS) 
roc_obj <- roc(data_prsice$PHENO, data_prsice$PRS)
plot.roc(roc_obj, auc.polygon=TRUE, col = "blue", print.auc=TRUE, header = "E", main = "Area under the curve (PRSice-2)") 
ggsave(paste0('fina_PRSice_AUC.pdf'))


########## PRS calculation using LDpred ##########


##### Step 1 Obtain GWAS summary stats + preprare base data set #####
#base_data_ldpred <- base_data%>%select("CHR", " BP", "SNP", "A1", "A2", "BETA", "P")
#write.table(base_data, "basedata_LDpred", row.names = F, col.names = T, sep = " ", quote = F)


##### Step 2 Run also with PRSice-2 #####
# 2.1 Coordinate base and target file
#system("ldpred coord \
#          --gf=GSA_impQC_maf_nosex \
#          --ssf=basedata_jansenAD \
#          --ssf-format=CUSTOM \
#          --chr CHR \
#          --pos BP \
#          --rs SNP \
#          --A1 A1 \
#          --A2 A2 \
#          --eff BETA \
#          --pval P \
#          --N 455258 \
#          --out AD_ldpred_coord")

# 2.2 Recalculate effect weights
#system("ldpred gibbs \
#          --cf AD_ldpred_coord \
#          --ldr 1474 \
#          --ldf AD_ldf \
#          --f 1 \
#          --N 455258 \
#          --ldf AD_ldf \
#          --out AD_ldpred_gibbs")

## 2.3 Calculate PRS using new effect weights
#system("Ldpred score \
#          --gf=GSA_impQC_maf_nosex \
#          --rf AD_ldpred_gibbs \
#          --out AD_ldpred_score")



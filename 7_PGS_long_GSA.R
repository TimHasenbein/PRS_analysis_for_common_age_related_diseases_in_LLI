##****************************************************************************##
## Global Screening Array Illumina (GSA) longevity-PGS GSA summary statistics ##
##****************************************************************************##
## Tim Hasenbein
## PhD. Guillermo G. Torres 
## Last modification 06.2020
## Use PRSice and calculate PRS for longevity. Base data set is from own results of association analysis. Therefore the data set has to be splitted --> Data set for association testing: 70% of cases + controls, for PRS calculation: 30 % of cases and controls
library(tidyverse)
library(data.table)
require(GWASTools)
library(qqman)
library(pROC)
library(rcompanion)


##### set up env. var ######
setwd("/home/tim/Desktop/PGS_longevity/")
genoQCed <- "F:/1_QC/results/"


##### Split controls in 70% and 30% data sets  #####
fam <- fread(paste0(genoQCed,"GSA_imp_PRS.fam")) # 5596 individuals
cases <- fam%>%filter(str_detect(V2, "AGE")) # 1295 cases  907 (70%), 388 (30%)
ctrl <- fam%>%filter(str_detect(V2,"AGE", negate = T)) # 4301 -> 3011 (70%), 1290 (30%)
## Generate data set for association analysis: cases 70%, ctrl 70%
assoc_cases <- sample_n(cases, 907)
assoc_ctrls <- sample_n(ctrl, 3011)
assoc_data <- assoc_cases%>%bind_rows(assoc_ctrls) # 3918 Ind.
assoc_ind <- assoc_data %>% select(V1, V2)
write.table(assoc_ind, "assoc_ind.txt", row.names = F, col.names = F, sep = " ", quote = F)
system(paste0("plink --bfile ",genoQCed,"GSA_imp_PRS --keep assoc_ind.txt --make-bed --out assoc_data")) 
# 7,363,285 variants, 3918 people (907 cases, 3011 controls)

## Generate data set for PRS calculation: cases 30%, ctrl 30%
system(paste0("plink --bfile ",genoQCed,"GSA_imp_PRS --remove assoc_ind.txt --make-bed --out prs_data")) 
# 7,363,285 variants, 1678 people (388  are cases and 1290 are controls)


############### Association analysis: Unimputed data ###############


##### Association testing with data set assoc_data #####
system("plink --bfile assoc_data  --logistic --covar GSAhapmap.eigenvec --covar-number 1-5 --out GSA_longevity")
GSAassoc <- fread("GSA_longevity.assoc.logistic", data.table = F)
GSAassoc_add <- GSAassoc[GSAassoc$TEST == "ADD", ] # 447285 variants
sigSNPs <- GSAassoc_add %>% arrange(P) %>% filter(P <= 5e-08)
SNPlistraw <- fread("assoc_data.bim", data.table = F)
sigSNPs <- sigSNPs %>% left_join(SNPlistraw, by=c('BP' = 'V4')) %>% select(-"V1", -"V3",-"V5", -"V6")
SNPshighlight <- sigSNPs %>% select("SNP")
write.table(SNPshighlight,"SNPshighlight", col.names = T, row.names = F, sep = " ", quote = F)
system("plink --bfile assoc_data  --assoc --adjust --out GSA_longevityassoc") # 1.15608 lambda


##### Data visualization #####
qqPlot(GSAassoc_add$P, thinThreshold = 2, ci = T, ylim = c(0,20), cex = .5) 
GSAassoc <- GSAassoc_add[ , c(1,2,3,9)]  
GSAassoc <- GSAassoc[!is.na(GSAassoc$P), ]
write.table(GSAassoc,"GSAassoc", col.names = T, row.names = F, sep = " ", quote = F)
chr19 <- GSAassoc_add %>% filter(CHR == 19)
apoe <- chr19 %>% filter(BP %in% 45000000:45800000)
apoesnps <- apoe$SNP
manhattan(GSAassoc, chr="CHR",bp="BP",p="P",snp="SNP", cex = 0.3, ylim = c(0, 25), highlight = apoesnps) 


########## PRS calculation using PRSice-2 ##########


##### Step 1: Obtain GWAS summary stats + preprare base data set #####
head(GSAassoc_add)
GSAassoc_base <-  GSAassoc_add%>%select("SNP", "A1", "OR", "P")
write.table(GSAassoc_base, "GSAassoc_base", row.names = F, col.names = T, sep = " ", quote = F)


##### Step 2: Run PRS calculation #####
system("Rscript PRSice.R --dir . \
--prsice ./PRSice_linux \
--base GSAassoc_base \
--target prs_data \
--snp SNP \
--A1 A1 \
--stat OR \
--pvalue P \
--binary-target T \
--bar-levels 5e-08,5e-07,5e-06,5e-05,5e-04,5e-03,5e-02,5e-01 \
--fastscore \
--extract PRS_owngwas.valid \
--no-regress \
--out PRS_owngwas") 


##### Step 3: Standardize polygenic scores #####
data_prsice <-fread("PRS_owngwas.all.score", header=T) 
head(data_prsice)
data_prsice <- data_prsice%>%rename("PRS" = "1")
data_prsice <- data_prsice%>%select("FID", "IID", "PRS")
data_prsice$PRS=(data_prsice$PRS-mean(data_prsice$PRS))/sd(data_prsice$PRS, na.rm=T)
LLI <- data_prsice%>%filter(str_detect(IID,"AGE"))
Ctrl <- data_prsice%>%filter(str_detect(IID,"AGE", negate = T))
LLI$PHENO <- 'LLI'
Ctrl$PHENO <- 'Ctrl'
data_prsice <- rbind(LLI, Ctrl)
head(data_prsice)


##### Step 4: Visualize data #####
# Histogram
ggplot(data_prsice, aes(x=PRS, color=PHENO)) +
  geom_histogram(fill="white", alpha=0.5, position="identity") +
  labs(title="PGS for longevity  (PRSice-2)", x = "PGS")
ggsave(paste0('own_PRSice_histoplot.pdf'))
mu <- data_prsice%>%group_by(PHENO)%>%summarise(M=mean(PRS),SD=sd(PRS),Med=median(PRS))
# Boxplot
ggplot(data_prsice%>%mutate(PRS=(PRS-mean(PRS))), aes(x=factor(PHENO),y=PRS,fill=factor(PHENO)))+
  stat_boxplot(geom='errorbar', linetype=1, width=0.5)+ 
  geom_boxplot(outlier.shape=1)+    
  stat_summary(fun.y=mean, geom="point", size=2) +
  labs(y = "PGS", x = "Phenotype")
ggsave(paste0('own_PRSice_boxplot.pdf'))
# Density plot  
ggplot(data_prsice,aes(x=PRS,fill=factor(PHENO),color=factor(PHENO)))+geom_density(alpha=0.1)+
  geom_vline(data=mu,aes(xintercept=M, color=factor(PHENO)),linetype="dashed")+theme_classic()+xlab('PGS')+
  scale_color_discrete(name = "", labels = c("Young controls", "LLI"))+scale_fill_discrete(name = "", labels = c("Young controls", "LLI"))
ggsave(paste0('own_PRSice_densityplot.pdf'))
# Linear plot
ggplot(data_prsice%>%arrange(PRS)%>%rownames_to_column('Pos')%>%mutate(Pos=ifelse(PHENO==2,as.numeric(Pos)+50,as.numeric(Pos))),aes(x=Pos,y=PRS,color=factor(PHENO)))+geom_point(alpha=0.3)+theme_classic()+ylab('Estimated Longevity-PRS (PRSice-2')+xlab('Individuals')+scale_color_discrete(name = "", labels = c("Young controls", "LLI"))
ggsave(paste0(mypath,'own_PRSice_points.pdf'))


##### Step 5: Statistics #####
# Shapiro-Wilk normality test 
shapiro.test(LLI$PRS)  
shapiro.test(Ctrl$PRS) 
# Barlett's test 
bartlett.test(PRS~PHENO, data_prsice) 
# Two Sample t-test for signifant differences between LLI and younger PRS mean 
wilcox.test(PRS~PHENO, data = data_prsice, var.equal = T) 


##### Step 6: Logistic regression model #####
mod1<-glm(as.factor(PHENO)~PRS, data=data_prsice, family = binomial)
summary(mod1)
nagelkerke(mod1) 
auc(data_prsice$PHENO, data_prsice$PRS) 
roc_obj <- roc(data_prsice$PHENO, data_prsice$PRS)
plot.roc(roc_obj, auc.polygon=TRUE, col = "blue", print.auc=TRUE, header = "E", main = "Area under the curve (PRSice-2)") 



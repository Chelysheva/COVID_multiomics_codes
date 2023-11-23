#### Import libraries ####
#R version 4.0.5
library(edgeR) 
packageVersion("edgeR")
#[1] ‘3.32.1’
library(limma)
packageVersion("limma")
#[1] ‘3.46.0’
library(biomaRt)
packageVersion("biomaRt")
#[1] ‘2.46.3’
library(dplyr)
packageVersion("dplyr")
#[1] ‘1.1.2’
library(ggplot2)
packageVersion("ggplot2")
#[1] ‘3.4.2’
library(ggpubr)
packageVersion("ggpubr")
#[1] ‘0.4.0’
library(EnhancedVolcano)
packageVersion("EnhancedVolcano")
#[1] ‘1.8.0’

#### Import data ####
mycounts=read.table("ONT_stage1_transcript_counts.txt",header = T)

#### Filter transcripts ####
mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
my_genes <- getBM(attributes=c('ensembl_transcript_id',
                               'hgnc_symbol'),
                  mart = mart)
mycounts$Reference=sub('\\..*', '', mycounts$Reference)
allmycounts=merge(my_genes,mycounts,by.x='ensembl_transcript_id',by.y='Reference')
allmycounts<-allmycounts[rowSums(allmycounts[3:41] >= 3) > 4, ]
allmycounts<-allmycounts[!grepl("HB", allmycounts$hgnc_symbol),]
allmycounts<-allmycounts[!grepl("MT", allmycounts$hgnc_symbol),]
allmycounts<-allmycounts[!grepl("RPL", allmycounts$hgnc_symbol),]
allmycounts<-allmycounts[!grepl("RPS", allmycounts$hgnc_symbol),]
allmycounts<-allmycounts[!grepl("RNA", allmycounts$hgnc_symbol),]
allmycounts<-allmycounts[!allmycounts$hgnc_symbol=="",]


allmycounts<-allmycounts[,order(colnames(allmycounts))]
####Import and prepare metadata####
meta=read.table("ONT_metadata.txt",header = T)

meta$cond<-paste(meta$TP,meta$vaccine,meta$NAAT,sep="_")
meta$cond<-gsub("D.*", "D0", meta$cond) #group all D0 together
meta$cond<-gsub("CT_ChAdOx1_0", "CT_0", meta$cond) #group all NAAT-ve together
meta$cond<-gsub("CT_MenACWY_0", "CT_0", meta$cond) #group all NAAT-ve together
table(meta$TP)
#CT D0 
#29 10 
table(meta$vaccine)
#ChAdOx1 MenACWY 
#17      22 
table(meta$NAAT)
#0  1 
#20 19 
table(meta$cond)
#CT_0 CT_ChAdOx1_1 CT_MenACWY_1           D0 
#13            7            9           10 


#Create DGE object
allmycounts=DGEList(as.matrix(allmycounts[,3:41]),genes = allmycounts[,c(1,2)],remove.zeros = T,group=meta$cond)

#Density plot 
allmycounts_cpm <-log(cpm(allmycounts$counts))  # make cpm for the density plot

df_long <- melt(allmycounts_cpm) #transform to long format
ggplot(df_long, aes(x = value, color = Var2)) +
  geom_density(alpha = 0.6) +
  labs(title = "Filtered data", x="Log-cpm", y="Density")+
  xlim(-5,12) +
  theme_pubr(legend="none")

#DGE
groups <- allmycounts$samples$group

design <- model.matrix(~0+groups)
cont.1<-makeContrasts(CT_ChAdOx1_1-CT_MenACWY_1,levels=groups)
cont.2<-makeContrasts(D0-CT_MenACWY_1,levels=groups)
cont.3<-makeContrasts(D0-CT_ChAdOx1_1,levels=groups)
cont.4<-makeContrasts(CT_0-CT_MenACWY_1,levels=groups)
cont.5<-makeContrasts(CT_0-CT_ChAdOx1_1,levels=groups)

cont.matrix<-cbind(cont.1,cont.2,cont.3,cont.4,cont.5)

vm <- voom(allmycounts, design=design)
fit <- lmFit(vm, design=design)
fit <- contrasts.fit(fit, cont.matrix)
fit <- eBayes(fit)
dt <- decideTests(fit)

summary(dt)

#### Save DGE tables ####
g_transcripts_vaccines1 <- topTable(fit, n=Inf, adjust="fdr",p.value = 1,coef = 1)
write.table(g_transcripts_vaccines1, "COVID_ONT_DGE_ChAdOx1_placebo_NAATpos.txt", quote = F, sep= "\t", row.names = F)
g_transcripts_vaccines2 <- topTable(fit, n=Inf, adjust="fdr",p.value = 1,coef = 2)
write.table(g_transcripts_vaccines2, "COVID_ONT_DGE_D0_placebo_NAATpos.txt", quote = F, sep= "\t", row.names = F)
g_transcripts_vaccines3 <- topTable(fit, n=Inf, adjust="fdr",p.value = 1,coef = 3)
write.table(g_transcripts_vaccines3, "COVID_ONT_DGE_D0_ChAdOx1_NAATpos.txt", quote = F, sep= "\t", row.names = F)
g_transcripts_vaccines4 <- topTable(fit, n=Inf, adjust="fdr",p.value = 1,coef = 4)
write.table(g_transcripts_vaccines4, "COVID_ONT_DGE_CT_placebo_NAATpos_NAATneg.txt", quote = F, sep= "\t", row.names = F)
g_transcripts_vaccines5 <- topTable(fit, n=Inf, adjust="fdr",p.value = 1,coef = 5)
write.table(g_transcripts_vaccines5, "COVID_ONT_DGE_CT_ChAdOx1_NAATpos_NAATneg.txt", quote = F, sep= "\t", row.names = F)


#### Volcano plots ####
g_vaccines <- topTable(fit, n=Inf, adjust="fdr",p.value = 1, coef = 1)
volcanodata <- g_vaccines
colnames(volcanodata)
volcanodata<-volcanodata[,c(2,3,6,7)]
# Set up edgeR compatible column names
colnames(volcanodata) <- c("Gene","logFC","pvalue","FDR")
# Set transcript names as row names
row.names(volcanodata)<-make.names(volcanodata[,1], unique=TRUE)

keyvals.colour <- ifelse(
  volcanodata$logFC > 0 & volcanodata$FDR<0.05, 'red',
  ifelse(volcanodata$logFC < 0 & volcanodata$FDR<0.05, 'blue',
         'black'))
keyvals.colour[is.na(keyvals.colour)] <- 'grey'
names(keyvals.colour)[keyvals.colour == 'red'] <- 'high'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'mid'
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'low'
top_genes_covid<- topTable(fit, n=Inf, adjust="fdr",p.value = 0.05, coef = 1)
EnhancedVolcano(volcanodata, 
                title = "",
                subtitle = "",
                legendPosition = 'none',
                lab = rownames(volcanodata),
                selectLab = top_genes_covid$hgnc_symbol,
                x = 'logFC',
                y = 'pvalue', 
                #FCcutoff = 0,
                labSize = 4,
                pointSize = 2.0, 
                xlim = c(-7,7),
                ylim=c(0,8),
                colCustom = keyvals.colour,
                cutoffLineType = 'blank',gridlines.major = F,gridlines.minor = F,vline = 0)

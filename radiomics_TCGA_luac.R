library(NMF)
library(ConsensusClusterPlus)
library(Cairo)
library(Biocomb)
library(pROC)

########### Assign Radiomics data ####################################################
CT <- read.delim('TCGA_result_romove_second_mass.txt')
## Select larger mass if CT shows two masses #
CT$ID <- gsub('_done', '',CT$Original.Num)
PET <- read.delim('TCGA_result_PET.txt')
PET$ID <- gsub('_done', '',PET$Original.Num)

########### Assign TCGA mutation data ###########################################

mutation <- read.delim('./broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf')
duplicated(gsub('-...-...-....-..', '', unique(mutation$Tumor_Sample_Barcode)))

mutation_ID <- gsub('-...-...-....-..', '', mutation$Tumor_Sample_Barcode)
mutation$ID <- mutation_ID

data_mut <- unstack(mutation, ID ~ Hugo_Symbol)

summary(data_mut)

write.table(data_rm_NA$ID, file = 'ID.txt', quote = FALSE, sep = ' ',
            row.names = FALSE,
            col.names = FALSE) 

################ CNV ###########################
cnv <- read.delim('./cnv/gdac.broadinstitute.org_LUAD-TP.CopyNumber_Gistic2.Level_4.2015082100.0.0/all_thresholded.by_genes.txt')
ID_cnv <- gsub('\\....\\....\\.....\\...', '', colnames(cnv)[-1:-3])
stk11_cnv <- data.frame(ID = gsub('\\.', '-', ID_cnv), cnv = as.numeric(cnv[cnv$Gene.Symbol == 'STK11',-1:-3]))
sum(stk11_cnv$cnv < -1)

dataCT$ID %in% mutation_ID
### subset having mutation data ########

dataCT <- CT[CT$ID %in% unique(mutation_ID),]
dataCT$EGFR <- dataCT$ID %in% data_mut$EGFR
dataCT$TP53 <- dataCT$ID %in% data_mut$TP53
dataCT$KRAS <- dataCT$ID %in% data_mut$KRAS
dataCT$BRAF <- dataCT$ID %in% data_mut$BRAF
dataCT$KEAP1 <- dataCT$ID %in% data_mut$KEAP1
dataCT$NF1 <- dataCT$ID %in% data_mut$NF1
dataCT$SETD2 <- dataCT$ID %in% data_mut$SETD2

####### STK11  mutation only one case #####################
############## CNV data not complete about half cases #####

clinical <- read.delim('./nationwidechildrens.org_clinical_patient_luad.txt') 
summary(clinical)
colnames(clinical)

data_ct <- merge(clinical, CT, by.x = 'bcr_patient_barcode', by.y = 'ID', all.y = TRUE)
summary(data_ct)

############## Fusion data from http://54.84.12.177/PanCanFusV2/ ######
fusion <- read.csv('translocation.csv')

dataCT$ID[(dataCT$ID %in% gsub('-01A', '', fusion$ID))]
fusion[fusion$ID == 'TCGA-50-8460-01A',]

dataCT$fusion <- rep('absence', dim(dataCT)[1])
dataCT$fusion [dataCT$ID == 'TCGA-50-8460'] <- 'EML4_ALK'


################# TCGA expression subtype ##############################
exp_subtype <- read.csv('tcga.luad.gene.expression.subtypes.20121025.csv')

dataCT <- merge(dataCT, exp_subtype, by.x = 'ID', by.y = 'sampleId', all.x = TRUE)

######################################################
############### STK11 ################################
##### From TCGA data portal ##########################




library(NMF)

d <- as.matrix(dataCT[,5:165])
d <- scale(d, center = TRUE, scale = TRUE)
d <- sweep(d, 1, apply(d, 1, median, na.rm = TRUE))
d <- d[complete.cases(d),]


############################################
########## consensusclusterplus#############

library(ConsensusClusterPlus)

results = ConsensusClusterPlus(d,maxK=10,reps=5000,pItem=0.8,pFeature=1,
                               title='consensus',
                               clusterAlg="hc",
                               innerLinkage = "ward.D2",
                               finalLinkage = "ward.D2",
                               distance="euclidean",
                               plot="pdf")

results[[4]]$consensusClass

icl = calcICL(results,title='consensus',plot="png")

icl[["clusterConsensus"]]

icl[["itemConsensus"]][1:500,]

data_rm_NA <- dataCT[complete.cases(dataCT[,5:165]),] 

ann <- data.frame(EGFR = factor(dataCT$EGFR)[complete.cases(dataCT[,5:165])],
                  TP53 = factor(dataCT$TP53)[complete.cases(dataCT[,5:165])],
                  KRAS = factor(dataCT$KRAS)[complete.cases(dataCT[,5:165])],
                  BRAF = factor(dataCT$BRAF)[complete.cases(dataCT[,5:165])],
                  KEAP1 = factor(dataCT$KEAP1)[complete.cases(dataCT[,5:165])],
                  #STK11 = factor(dataCT$STK11)[complete.cases(dataCT[,5:165])],
                  NF1 = factor(dataCT$NF1)[complete.cases(dataCT[,5:165])],
                  SETD2 = factor(dataCT$SETD2)[complete.cases(dataCT[,5:165])]
                  #Fusion = factor(dataCT$fusion)[complete.cases(dataCT[,5:165])],
)

#'ward', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'


library(Cairo)
color <- colorRampPalette(brewer.pal(9,'Blues'))(100)

CairoPDF(file = './heatmaps/mcquitty.pdf',
         width =7.5, height = 7.5, pointsize = 16)

ah<- aheatmap(d, 
              distfun = "euclidean",
              hclustfun = "mcquitty",
              annRow = ann)
dev.off()

CairoPDF(file = './heatmaps/complete.pdf',
         width =7.5, height = 7.5, pointsize = 16)

ah<- aheatmap(d, 
              distfun = "euclidean",
              hclustfun = "complete",
              annRow = ann)
dev.off()

CairoPDF(file = './heatmaps/ward_euclidean.pdf',
         width =7.5, height = 7.5, pointsize = 16)

ah<- aheatmap(d, 
              distfun = "euclidean",
              color = color,
              hclustfun = "ward",
              annRow = ann,
              Colv = results[[5]]$consensusTree)
dev.off()

CairoPDF(file = './heatmaps/average.pdf',
         width =7.5, height = 7.5, pointsize = 16)

ah<- aheatmap(d, 
              distfun = "euclidean",
              hclustfun = "average",
              annRow = ann)
dev.off()

CairoPDF(file = './heatmaps/median.pdf',
         width =7.5, height = 7.5, pointsize = 16)

ah<- aheatmap(d, 
              distfun = "euclidean",
              hclustfun = "median",
              annRow = ann)
dev.off()

CairoPDF(file = './heatmaps/centroid.pdf',
         width =7.5, height = 7.5, pointsize = 16)

ah<- aheatmap(d, 
              distfun = "euclidean",
              hclustfun = "centroid",
              annRow = ann)
dev.off()

CairoPDF(file = './heatmaps/single.pdf',
         width =7.5, height = 7.5, pointsize = 16)

ah<- aheatmap(d, 
              distfun = "euclidean",
              hclustfun = "single",
              annRow = ann)
dev.off()


hc <- hclust(dist(d, method = "euclidean"), method = 'ward.D2')
ct_4 <- cutree(hc, 4)

ann$EGFR[ct_4 == 1]
colnames(d)[ct_4 == 2]
colnames(d)[ct_4 == 3]
colnames(d)[ct_4 == 4]


table(ct_4 ==1, ann$EGFR)
fisher.test(ct_4 ==1, ann$EGFR)

#######################################
### How clustering in test.set ?#########
### Rand statistic #####################





##### AUC ##############################
library(pROC)

auc_EGFR <- c()
for (i in 5:165){
  a <- roc( data_rm_NA$EGFR, data_rm_NA[,i],
       direction="auto", auc=TRUE)
  
      auc_EGFR <- rbind(auc_EGFR, 0.5 + abs(a$auc-0.5))
}


auc_TP53 <- c()
for (i in 5:165){
  a <- auc(data_rm_NA$TP53, data_rm_NA[,i])
  auc_TP53 <- rbind(auc_TP53, 0.5 + abs(a -0.5))
}

auc_KRAS <- c()
for (i in 5:165){
  a <- auc(data_rm_NA$KRAS, data_rm_NA[,i])
  auc_KRAS <- rbind(auc_KRAS, 0.5 + abs(a -0.5))
}

auc_BRAF <- c()
for (i in 5:165){
  a <- auc(data_rm_NA$BRAF, data_rm_NA[,i])
  auc_BRAF <- rbind(auc_BRAF, 0.5 + abs(a -0.5))
}

auc_KEAP1 <- c()
for (i in 5:165){
  a <- auc(data_rm_NA$KEAP1, data_rm_NA[,i])
  auc_KEAP1 <- rbind(auc_KEAP1, 0.5 + abs(a -0.5))
}

auc_NF1 <- c()
for (i in 5:165){
  a <- auc(data_rm_NA$NF1, data_rm_NA[,i])
  auc_NF1 <- rbind(auc_NF1, 0.5 + abs(a -0.5))
}

auc_SETD2 <- c()
for (i in 5:165){
  a <- auc(data_rm_NA$SETD2, data_rm_NA[,i])
  auc_SETD2 <- rbind(auc_SETD2, 0.5 + abs(a -0.5))
}

auc_data <- t(cbind(auc_EGFR, auc_TP53, auc_KRAS, auc_BRAF, auc_KEAP1, auc_NF1, auc_SETD2))


CairoPDF(file = 'auc.pdf',
         width =7.5, height = 7.5, pointsize = 16)
aheatmap(auc_data,
         color = color,
         Rowv = NA,
         Colv = results[[4]]$consensusTree)

dev.off()




####### K =5 ############################
CairoPDF(file = './heatmaps/ward_k5.pdf',
         width =7.5, height = 7.5, pointsize = 16)

ah<- aheatmap(d, 
              distfun = "euclidean",
              color = color,
              hclustfun = "ward",
              annRow = ann,
              Colv = results[[5]]$consensusTree)
dev.off()
CairoPDF(file = 'auc_k5.pdf',
         width =7.5, height = 7.5, pointsize = 16)
aheatmap(auc_data,
         color = color,
         Rowv = NA,
         Colv = results[[5]]$consensusTree)

dev.off()


tapply(auc_data[1,], results[[5]]$consensusClass, mean)
tapply(auc_data[2,], results[[5]]$consensusClass, mean)
tapply(auc_data[3,], results[[5]]$consensusClass, mean)
tapply(auc_data[4,], results[[5]]$consensusClass, mean)
tapply(auc_data[5,], results[[5]]$consensusClass, mean)
tapply(auc_data[6,], results[[5]]$consensusClass, mean)
tapply(auc_data[7,], results[[5]]$consensusClass, mean)

tapply(auc_data[1,], results[[5]]$consensusClass, sd)
tapply(auc_data[2,], results[[5]]$consensusClass, sd)
tapply(auc_data[3,], results[[5]]$consensusClass, sd)
tapply(auc_data[4,], results[[5]]$consensusClass, sd)
tapply(auc_data[5,], results[[5]]$consensusClass, sd)
tapply(auc_data[6,], results[[5]]$consensusClass, sd)
tapply(auc_data[7,], results[[5]]$consensusClass, sd)


################# AUC significance ####################
########### http://www.inside-r.org/node/322637 #########
 
compute.auc.permutation(auc_data[1,],results[[5]]$consensusClass,repetitions=1000)

# example
data(data_test)

# class label must be factor
data_test[,ncol(data_test)]<-as.factor(data_test[,ncol(data_test)])

auc.val=compute.aucs(dattable=data_test)
vauc<-auc.val[,"AUC"]
rep.num<-20

p.values=compute.auc.permutation(aucs=vauc,dattable=data_test,rep.num)

######################################################
######### Pairwise correlation for medoids############

cor_table_cluster_1 <- cor(d[,results[[5]]$consensusClass == 1])
mean_pairwise_cor_1 <- apply(cor_table_cluster_1, 1, mean)
medoid_1 <- colnames(d[,results[[5]]$consensusClass == 1])[mean_pairwise_cor == max(mean_pairwise_cor)]
mean_cluster_cor_1 <- mean(cor_table_cluster_1)
sd_cluster_cor_1 <- sd(cor_table_cluster_1)

cor_table_cluster_2 <- cor(d[,results[[5]]$consensusClass == 2])
mean_pairwise_cor_2 <- apply(cor_table_cluster_2, 1, mean)
medoid_2 <- colnames(d[,results[[5]]$consensusClass == 2])[mean_pairwise_cor == max(mean_pairwise_cor)]
mean_cluster_cor_2 <- mean(cor_table_cluster_2)
sd_cluster_cor_2 <- sd(cor_table_cluster_2)

cor_table_cluster_3 <- cor(d[,results[[5]]$consensusClass == 3])
mean_pairwise_cor_3 <- apply(cor_table_cluster_3, 1, mean)
medoid_3 <- colnames(d[,results[[5]]$consensusClass == 3])[mean_pairwise_cor == max(mean_pairwise_cor)]
mean_cluster_cor_3 <- mean(cor_table_cluster_3)
sd_cluster_cor_3 <- sd(cor_table_cluster_3)

cor_table_cluster_4 <- cor(d[,results[[5]]$consensusClass == 4])
mean_pairwise_cor_4 <- apply(cor_table_cluster_4, 1, mean)
medoid_4 <- colnames(d[,results[[5]]$consensusClass == 4])[mean_pairwise_cor == max(mean_pairwise_cor)]
mean_cluster_cor_4 <- mean(cor_table_cluster_4)
sd_cluster_cor_4 <- sd(cor_table_cluster_4)

cor_table_cluster_5 <- cor(d[,results[[5]]$consensusClass == 5])
mean_pairwise_cor_5 <- apply(cor_table_cluster_5, 1, mean)
medoid_5 <- colnames(d[,results[[5]]$consensusClass == 5])[mean_pairwise_cor == max(mean_pairwise_cor)]
mean_cluster_cor_5 <- mean(cor_table_cluster_5)
sd_cluster_cor_5 <- sd(cor_table_cluster_5)




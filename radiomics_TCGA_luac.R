
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

### subset having mutation data ########

dataCT <- CT[CT$ID %in% unique(mutation_ID),]
dataCT$EGFR <- dataCT$ID %in% data_mut$EGFR
dataCT$TP53 <- dataCT$ID %in% data_mut$TP53
dataCT$KRAS <- dataCT$ID %in% data_mut$KRAS


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



library(NMF)

d <- as.matrix(dataCT[,5:164])
d <- scale(d, center = TRUE, scale = TRUE)
d <- sweep(d, 1, apply(d, 1, median, na.rm = TRUE))
d <- d[complete.cases(d),]

ann <- data.frame(EGFR = factor(dataCT$EGFR)[complete.cases(dataCT[,5:165])],
                  TP53 = factor(dataCT$TP53)[complete.cases(dataCT[,5:165])],
                  KRAS = factor(dataCT$KRAS)[complete.cases(dataCT[,5:165])],
                  #Fusion = factor(dataCT$fusion)[complete.cases(dataCT[,5:165])],
                  exp_subtype = factor(dataCT$expression_subtype)[complete.cases(dataCT[,5:165])])


#'ward', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'


library(Cairo)

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

CairoPDF(file = './heatmaps/ward.pdf',
         width =7.5, height = 7.5, pointsize = 16)

ah<- aheatmap(d, 
              distfun = "euclidean",
              hclustfun = "ward",
              annRow = ann)
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


#######################################
### How clustering in test.set ?#########


p_value_KRAS <- c()

for (i in 5:165) {
  a <- wilcox.test(dataCT[,i] ~ dataCT$KRAS)
  p_value_KRAS <- append(p_value_KRAS, a$p.value)
}

t.test()
p_value_KRAS <0.01
colnames(dataCT)[17]

hist(log(dataCT$Energy_val))
boxplot(Energy_val ~ KRAS, dataCT)

p_value_EGFR <- c()

for (i in 5:165) {
  a <- wilcox.test(dataCT[,i] ~ dataCT$EGFR)
  p_value_EGFR <- append(p_value_EGFR, a$p.value)
}

sum(p_value_EGFR <0.01)

p_value_TP53 <- c()
for (i in 5:165) {
  a <- wilcox.test(dataCT[,i] ~ dataCT$TP53)
  p_value_TP53 <- append(p_value_TP53, a$p.value)
}

sum(p_value_TP53 <0.01)

par(mfrow = c(4,4))
for (i in 4:19) {
  hist(dataCT[,i])
}

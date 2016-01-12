######## Download from TCGA website #######################################
url_RPPA <- 'https://tcga-data.nci.nih.gov/docs/publications/luad_2014/LUAD_2014.MDA_RPPA_Core.Level_3/mdanderson.org_LUAD.MDA_RPPA_Core.Level_3.1.5.0.tar.gz'
url_GISTIC <- 'https://tcga-data.nci.nih.gov/docs/integration/adfs/vendor/TCGA.DCC.GenomeWideSNP6.marker.na31.lst.gz'
url_MAF <- 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/other/publications/luad_2014/AN_TCGA_LUAD_PAIR_capture_freeze_FINAL_230.aggregated.capture.tcga.uuid.curated.somatic.maf'

download.file(url_RPPA, './RPPA.tar.gz', "auto", quiet = FALSE, mode = "w",
              cacheOK = TRUE,
              extra = getOption("download.file.extra"))


download.file(url_GISTIC, './GISTIC.lst.gz', "auto", quiet = FALSE, mode = "w",
              cacheOK = TRUE,
              extra = getOption("download.file.extra"))

download.file(url_MAF, './mutation.maf', "auto", quiet = FALSE, mode = "w",
              cacheOK = TRUE,
              extra = getOption("download.file.extra"))

untar('RPPA.tar.gz', exdir = '.')
untar('gdac.broadinstitute.org_LUAD.Clinical_Pick_Tier1.Level_4.2015082100.1.0.tar.gz',
      exdir = '.')
untar('clinical_TCGA_portal.tar', exdir = '.')
untar('clinical_17_sample.tar', exdir = './clinical_sample')

untar('gdac.broadinstitute.org_LUAD.Mutation_Packager_Calls.Level_3.2015082100.0.0.tar.gz',
      exdir = '.')
untar('mutation_TCGA_portal.tar', exdir = '.')


########### Assign Radiomics data ####################################################
CT <- read.delim('TCGA_result_romove_second_mass.txt')
## Select larger mass if CT shows two masses #
CT$ID <- gsub('_done', '',CT$Original.Num)
PET <- read.delim('TCGA_result_PET.txt')
PET$ID <- gsub('_done', '',PET$Original.Num)

########### Assign TCGA data ###########################################

mutation <- read.delim('./mutation/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf')
duplicated(gsub('-...-...-....-..', '', unique(mutation$Tumor_Sample_Barcode)))

mutation_ID <- gsub('-...-...-....-..', '', mutation$Tumor_Sample_Barcode)
mutation$ID <- mutation_ID
library(dplyr)
data_mut <- unstack(mutation, ID ~ Hugo_Symbol)
summary(data_mut)

dataCT <- CT[CT$ID %in% unique(mutation_ID),]
dataCT$EGFR <- dataCT$ID %in% data_mut$EGFR
dataCT$TP53 <- dataCT$ID %in% data_mut$TP53
dataCT$KRAS <- dataCT$ID %in% data_mut$KRAS


clinical <- read.delim('./clinical_TCGA_portal/nationwidechildrens.org_clinical_patient_luad.txt') 
summary(clinical)

data_ct <- merge(clinical, CT, by.x = 'bcr_patient_barcode', by.y = 'ID', all.y = TRUE)
summary(data_ct)

con <- file('mutation.maf')
readLines(con, n=4)
mutation <- read.delim('./mutation/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf')


mutation_ID <- read.delim('./gdac.broadinstitute.org_LUAD.Mutation_Packager_Calls.Level_3.2015082100.0.0/MANIFEST.txt',sep = ' ')[,2]
mutation_ID <- gsub('-01.maf.txt', '', mutation_ID)
read.delim('./gdac.broadinstitute.org_LUAD.Mutation_Packager_Calls.Level_3.2015082100.0.0/MANIFEST.txt')

sum(CT$ID %in% mutation_ID)
sum(CT$ID %in% clinical$bcr_patient_barcode)

fusion <- read.csv('tcga.luad.fusions.20121029.csv')

########### CT data #######################################
library(gplots)
heatmap(scale(as.matrix(CT[,4:164]),center = TRUE, scale = TRUE))


library(NMF)
hc = hclust(dist(data_hc_s), 'complete')


memb_3 <- factor(cutree(hc, k = 3), levels = 1:3,
                 labels = c('A', 'B', 'C'))

data_serous$group <- memb_3

ann_s <- data.frame(Subtype_h = data_serous$histological_type,
                    cluster = data_serous$group,                  
                    Subtype_m = data_serous$SUBTYPE,
                    XIST = log(data_serous$XIST),
                    BRCA1 = data_serous$BRCA1,
                    BRCA1_M = factor(data_serous$BRCA1_meth > 0.4))

ann_colors_s <- list(Subtype_h = c('yellow', 'blue'),
                     cluster = c('darkolivegreen','deepskyblue', 'darkslategray1'),
                     Subtype_m = c('lightblue', 'orange', 'purple', 'green'),
                     BRCA1 = c('greenyellow', 'lightskyblue'),
                     BRCA1_M = c('seagreen', 'mediumpurple')
)

library(NMF)
library(Cairo)
CairoPDF(file = '/home/jun/XIST_UT/Figures/Figure1.pdf',
         width =7.5, height = 7.5, pointsize = 16)
layout(matrix(c(1,1,2,2), ncol = 2, byrow = TRUE),
       widths = c(1,1),
       heights = c(276,95)) 
#########################################################################

ann <- data.frame(EGFR = factor(dataCT$EGFR),
                  TP53 = factor(dataCT$TP53),
                  KRAS = factor(dataCT$KRAS))

ann_colors <- list(Subtype_h = c('yellow', 'blue'),
                     cluster = c('chartreuse', 'gold', 'cyan'),
                     Subtype_m = c('lightblue', 'orange', 'purple', 'green'))

ah<- aheatmap(scale(as.matrix(dataCT[,4:164]),center = TRUE, scale = TRUE), 
              hclustfun=function(d) hclust(d, method="average"),
              #Colv = colv,
              annRow = ann,
              #annCol = ann_col,
              #annColors = ann_colors_e,
              #cex = 2,
              #labRow = rep('',dim(data_hc_e)[1]),
              #labCol = rep('',dim(data_hc_e)[2])
              #fontsize = 12,
              #cexCol = ,
              #fontsize = 16,
              #labCol = rep('',dim(data_hc_e)[2]),
              #Colv = colv
              #reorderfun = function(d, w) reorder(d, 10)
)



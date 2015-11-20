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


########### Assign Radiomics data ####################################################
CT <- read.delim('TCGA_result_romove_second_mass.txt')
PET <- read.delim('TCGA_result_PET.txt')




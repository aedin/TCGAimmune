# Script to get updated panCancer (published) files

#Data are available at https://gdc.cancer.gov/about-data/publications/pancanatlas
# The manifest is at https://gdc.cancer.gov/files/public/file/PanCan-General_Open_GDC-Manifest_0.txt
# All of these files can be downloaded using the gdc-client

# on the command line
# gdc-client download -m https://gdc.cancer.gov/files/public/file/PanCan-General_Open_GDC-Manifest_0.txt

# Using R
library(GenomicDataCommons)
ids<-read.table("https://gdc.cancer.gov/files/public/file/PanCan-General_Open_GDC-Manifest_0.txt", header=TRUE, as.is=TRUE)

datadir ="./data/TCGA"
gdc_TGCA = gdcdata(ids$id, destination_dir=datadir)

#file information is in  extdata
fileinfo<-read.csv("./inst/extdata/TGCA_DataFiles_README.csv", as.is=TRUE)
fileinfo<-read.csv(system.file("extdata", "TGCA_DataFiles_README.csv", package = "curatedTCGAPanImmune"))



#require(maftools)
#maf<-read.maf(file.path(datadir,fileInfo["Mutations",), isTCGA = TRUE)






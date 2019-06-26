# Download datasets

.downloadTIMER<-function() {
  library(S4Vectors)
  urlPath= "https://cistrome.shinyapps.io/timer/_w_b38e86ce/immuneEstimation.txt"
  TIMER<-read.delim(urlPath,as.is = TRUE)
  TIMER <-DataFrame(TIMER)
  save(TIMER, file="./data/TIMER.rda")
}


.downloadTCGA<-function(){library(openxlsx)
  library(S4Vectors)
  SupplUrls<- readr::read_delim("inst/extdata/PanImmune_Manuscript_Suppl_Manifest.txt", delim="\t")
  urls<-SupplUrls[grep("xlsx", SupplUrls$URL),]
  suppl<- lapply(urls$URL, read.xlsx)
  suppl<-lapply(suppl, DataFrame)
  names(suppl) = urls$Shortname
  suppl$sources= urls
  TCGAimmune<-suppl
  save(TCGAimmune, file="./data/TCGApanCancerImmune_suppl.rda")
}


.downloadTamborero<-function(){
  library(openxlsx)
  library(S4Vectors)
  SupplUrls<- readr::read_delim("inst/extdata/Tamborero_Manuscript_Suppl_Manifest.txt", delim="\t")
  urls<-SupplUrls[grep("xlsx", SupplUrls$URL),]
  suppl<- lapply(urls$URL, read.xlsx)
  suppl<-lapply(suppl, DataFrame)
  names(suppl) = urls$Shortname
  suppl$sources= urls
  Tamborero<-suppl
  names(Tamborero)[2]<- "ImmuneSubtypes"
  Tamborero$sources$Shortname[2] =  names(Tamborero)[2]
  save(Tamborero, file="./data/TamboreroImmune_suppl.rda")
  moGSA<-readxl::read_excel("./data/cl16_clust.summary.xls", sheet=3)
  save(moGSA, file="./data/moGSA.rda")
}

ImmuneMAE<-function(){
  data("TIMER")
  data("TamboreroImmune_suppl")
  data("TCGApanCancerImmune_suppl")
  TIMER$TCGA.Participant.Barcode<-substr(as.character(TIMER$barcode),1,12)

  # Make MAE
  require(MultiAssayExperiment)
  ilist <- list(TIMER = TIMER, Tamborero = Tamborero$ImmuneSubtypes,
                 TCGA = TCGAimmune$features)
  iSampleMap <- listToMap(ilist)
  iSampleMap$primary= c(TIMER$TCGA.Participant.Barcode, TCGAimmune$features$TCGA.Participant.Barcode, Tamborero$ImmuneSubtypes$Patient.TCGA.ID)
  iSampleMap$colname= c(as.character(TIMER$barcode), TCGAimmune$features$TCGA.Participant.Barcode, Tamborero$ImmuneSubtypes$Patient.TCGA.ID)

  tDataFrame<-function(x){
    nam = as.character(x[,1])
    x= x[,-1]
    x<-as.data.frame(x)
    x= t(x)
    colnames(x) = nam
    return(x)
  }

  ilist<-lapply(ilist,tDataFrame)



   colDat = TCGAimmune$features[,c("TCGA.Participant.Barcode", "TCGA.Study",
                                   "OS",  "OS.Time","PFI","PFI.Time")]
  rownames(colDat) = colDat[,1]
  TCGAimmuneMAE <- MultiAssayExperiment(ilist,colData = colDat, sampleMap = iSampleMap)

}

#subsetByColumn(TCGAimmuneMAE[,,"TIMER"],TCGAimmuneMAE$TCGA.Study=="KIRC")


ImmuneList<-function(subset=NULL){

  ilist <- list(TIMER = TIMER, Tamborero = Tamborero$ImmuneSubtypes,
                TCGA = TCGAimmune$features)

}
}

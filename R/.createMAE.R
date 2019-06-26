# Merge Immune Data into MAE
library(tibble)
library(magrittr)
library(MultiAssayExperiment)
library(curatedTCGAPanImmune)
library(curatedTCGAData)
library(ComplexHeatmap)
library(magrittr)
library(SummarizedExperiment)

#Tamborero Immune Subtypes in Kidney Cancer
# Changed name on Tamborero so need to rebuild
load("~aculhane/Dropbox (Partners HealthCare)/TCGA_immune/curatedTCGAPanImmune/data/Tamborero.rda")
#data(Tamborero)

#moGSA
data(moGSA)
moGSA<-t(moGSA%>% as.data.frame %>% column_to_rownames("barcode"))

#TIMER
data(TIMER)
TIMER<-t(TIMER %>%as.data.frame %>% column_to_rownames("barcode"))

#data(Tamborero)
Tamborero= Tamborero$ImmuneSubtypes
Tamborero<-t(Tamborero %>%as.data.frame %>% column_to_rownames("Patient.TCGA.ID"))

#Features
data(TCGApanCancerImmune_suppl)
Thorsson <-t(TCGAimmune$features %>%as.data.frame %>% column_to_rownames("TCGA.Participant.Barcode"))
pheno<-as.data.frame(t(Thorsson[c(1,32:35),]))
pheno$barcode = rownames(pheno)
Thorsson<- Thorsson[-c(1,32:35),]


TCGAImmune<-list(Tamborero=Tamborero, TIMER=TIMER, Thorsson=Thorsson,moGSA=moGSA )


# Create sample map for each datasets
# The sample map needs assay, primary id, colname of the assay
library(MultiAssayExperiment)
sampleMap<- DataFrame(
  assay=rep(names(sapply(TCGAImmune, ncol)), sapply(TCGAImmune, ncol)),
  primary= substr(unlist(sapply(TCGAImmune, colnames)),1,12),
  colname= unlist(sapply(TCGAImmune, colnames))
)


TCGAimmuneMAE <- MultiAssayExperiment(TCGAImmune, pheno, sampleMap )
save(TCGAimmuneMAE, file="~/Dropbox/TCGA_immune/curatedTCGAPanImmune/data/TCGAimmuneMAE.rda")
upsetSamples(TCGAimmuneMAE)


# Add immune subtypes to the colData for each assay type
library(SummarizedExperiment)
TCGAimmuneMAE[[1]] <- SummarizedExperiment(TCGAimmuneMAE[[1]])
TCGAimmuneMAE[[2]] <- SummarizedExperiment(TCGAimmuneMAE[[2]])
TCGAimmuneMAE[[3]] <- SummarizedExperiment(TCGAimmuneMAE[[3]])
TCGAimmuneMAE[[4]] <- SummarizedExperiment(TCGAimmuneMAE[[4]])

# Create colData for each with the Immune Subtypes

colData(TCGAimmuneMAE[[1]]) =DataFrame(t(assay(TCGAimmuneMAE[[1]])))
colData(TCGAimmuneMAE[[2]]) =DataFrame(t(assay(TCGAimmuneMAE[[2]])))

colData(TCGAimmuneMAE[[3]]) =DataFrame(t(assay(TCGAimmuneMAE[[3]])))

colData(TCGAimmuneMAE[[4]]) =DataFrame(t(assay(TCGAimmuneMAE[[4]])))

ind<-!apply(assay(TCGAimmuneMAE[[4]]), 2, function(x) all(is.na(x)))
cInd<-grep("^L", rownames(TCGAimmuneMAE[[4]]))

# Thorsson
TCGAimmuneMAE[[3]]$Thorsson_Subtypes = TCGAimmuneMAE[[3]]$Immune.Subtype
identical(rownames(colData(TCGAimmuneMAE)),colnames(TCGAimmuneMAE[[3]]))
TCGAimmuneMAE$Thorsson_Subtypes=TCGAimmuneMAE[[3]]$Thorsson_Subtypes

#moGSA
cl<-apply(assay(TCGAimmuneMAE[[4]])[cInd, ind],2, which.max)
cl<-rownames(TCGAimmuneMAE[[4]])[cInd][cl]
cl[cl=="NULL"] <-NA
# ignore C9
cl[cl=="L9"] <-NA
TCGAimmuneMAE[[4]]$moGSA_Subtypes = cl
assay(TCGAimmuneMAE[[4]])<-rbind(assay(TCGAimmuneMAE[[4]]), moGSA_Subtypes=cl)
x3<-wideFormat(TCGAimmuneMAE[1,,4], colDataCols = "Thorsson_Subtypes")
x1<-as.data.frame(colData(TCGAimmuneMAE[[4]]))
x1$id= rownames(colData(TCGAimmuneMAE[[4]]))
x1$barcode = substr(x1$id,1,12)
xM<-merge(as.data.frame(x1), as.data.frame(x3)[,c(1,3)], by.x="barcode", by.y="primary", all=TRUE)
rownames(xM) = xM$id
identical(rownames(xM), colnames(TCGAimmuneMAE[[4]]))
xM<-xM[colnames(TCGAimmuneMAE[[4]]),]
colData(TCGAimmuneMAE[[4]])<-DataFrame(xM)

# TIMER
TCGAimmuneMAE[[2]]$TIMER_kmeans_Subtypes<-kmeans(t(assay(TCGAimmuneMAE[[2]])),6)$cluster
TCGAimmuneMAE[[2]]$TIMER_mclust_Subtypes<-mclust:::Mclust(t(assay(TCGAimmuneMAE[[2]])))$classification
TCGAimmuneMAE[[2]]$barcode= substr(colnames(TCGAimmuneMAE[[2]]),1,12)
x2<-wideFormat(TCGAimmuneMAE[2,,2], colDataCols = "Thorsson_Subtypes")[,1:2]
x1<-as.data.frame(colData(TCGAimmuneMAE[[2]]))
x1$id= rownames(colData(TCGAimmuneMAE[[2]]))
xM<-merge(as.data.frame(x1), as.data.frame(x2), by.x="barcode", by.y="primary", all=TRUE)
rownames(xM) = xM$id
xM<-xM[rownames(colData(TCGAimmuneMAE[[2]])),]
identical(rownames(xM), colnames(TCGAimmuneMAE[[2]]))
colData(TCGAimmuneMAE[[2]])<-DataFrame(xM)


#Tamborero
TCGAimmuneMAE[[1]]$Tamborero_Subtypes = TCGAimmuneMAE[[1]]$Immune.phenotype
x2<-wideFormat(TCGAimmuneMAE[2,,1], colDataCols = "Thorsson_Subtypes")
x2<-cbind(colData(TCGAimmuneMAE[[1]])[x2$primary,], x2[,c(1,3)])
colData(TCGAimmuneMAE[[1]])<-x2[rownames(colData(TCGAimmuneMAE[[1]])),]


save(TCGAimmuneMAE, file="~/Dropbox/TCGA_immune/curatedTCGAPanImmune/data/TCGAimmuneMAE_withColData.rda")
upsetSamples(TCGAimmuneMAE)
# Lets merge all colDatas


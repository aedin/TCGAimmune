

# DLG
# dgibbs@systemsbiology.org
# March 29, 2017

# In several cases, more coarse grained cell content estimates are useful.
# Since each individual sums to 1. We can sum sub-catagories by sample.

###################################################################################
library(stringr)
library(readr)
library(dplyr)
mannot <- read_tsv("data/merged_sample_quality_annotations.tsv")
mannot_selected <- mannot %>% filter(Do_not_use == 'False' & AWG_excluded_because_of_pathology == 0)
mannot_unique <- mannot_selected %>% select(patient_barcode, `cancer type`) %>% unique()

ciber <- read.table("data/correlates/TCGA.Kallisto.fullIDs.cibersort.relative.tsv", stringsAsFactors=F, header=T)
ciber$pan.samplesID  <- str_sub(ciber$SampleID, start=1, end=12)
ciber$pan.samplesID <- str_replace_all(ciber$pan.samplesID, pattern="\\.", replacement="-")


###################################################################################

#Aggragate 1


newCiber <- data.frame(
  patient_barcode=ciber$pan.samplesID,
  Lymphocytes=c(ciber$B.cells.naive+ciber$B.cells.memory+ciber$T.cells.CD4.naive+ciber$T.cells.CD4.memory.resting+ciber$T.cells.CD4.memory.activated+ciber$T.cells.follicular.helper+
      ciber$T.cells.regulatory..Tregs+ciber$T.cells.gamma.delta+ciber$T.cells.CD8+ciber$NK.cells.resting+ciber$NK.cells.activated+ciber$Plasma.cells),
  Neutrophils=ciber$Neutrophils,
  Eosinophils=ciber$Eosinophils,
  Mast.cells=(ciber$Mast.cells.resting + ciber$Mast.cells.activated),
  Dendritic.cells=(ciber$Dendritic.cells.resting + ciber$Dendritic.cells.activated),
  Macrophage=(ciber$Monocytes + ciber$Macrophages.M0 + ciber$Macrophages.M1 + ciber$Macrophages.M2)
  )

til_ciber <- inner_join(newCiber, mannot_unique)
#dim(til_ciber)
#[1] 10241    19

#################################################################################

# Aggregate #2

newCiber <- data.frame(
  patient_barcode=ciber$pan.samplesID,
  Neutrophils=ciber$Neutrophils,
  Eosinophils=ciber$Eosinophils,
  Mast.cells=(ciber$Mast.cells.resting + ciber$Mast.cells.activated),
  Dendritic.cells=(ciber$Dendritic.cells.resting + ciber$Dendritic.cells.activated),
  Macrophage=(ciber$Macrophages.M0 + ciber$Macrophages.M1 + ciber$Macrophages.M2),
  NK.cells=(ciber$NK.cells.resting+ciber$NK.cells.activated),
  B.cells=(ciber$B.cells.naive + ciber$B.cells.memory),
  T.cells.CD8=ciber$T.cells.CD8,
  T.cells.CD4=(ciber$T.cells.CD4.naive+ciber$T.cells.CD4.memory.resting+ciber$T.cells.CD4.memory.activated),
  )

net_ciber <- inner_join(newCiber, mannot_unique)
#dim(til_ciber)
#[1] 10241    19


###################################################################################

# Aggregate #3

newCiber <- data.frame(
  patient_barcode=ciber$pan.samplesID,
  B.cells=(ciber$B.cells.naive + ciber$B.cells.memory),
  Plasma.cells=ciber$Plasma.cells,
  T.cells.CD8=ciber$T.cells.CD8,
  T.cells.CD4=(ciber$T.cells.CD4.naive+ciber$T.cells.CD4.memory.resting+ciber$T.cells.CD4.memory.activated+ciber$T.cells.follicular.helper+ciber$T.cells.regulatory..Tregs),
  T.cells.gamma.delta=ciber$T.cells.gamma.delta,
  NK.cells=(ciber$NK.cells.resting+ciber$NK.cells.activated),
  Macrophage=(ciber$Monocytes + ciber$Macrophages.M0 + ciber$Macrophages.M1 + ciber$Macrophages.M2),
  Dendritic.cells=(ciber$Dendritic.cells.resting + ciber$Dendritic.cells.activated),
  Mast.cells=(ciber$Mast.cells.resting + ciber$Mast.cells.activated),
  Neutrophils=ciber$Neutrophils,
  Eosinophils=ciber$Eosinophils
  )

gentles_ciber <- inner_join(newCiber, mannot_unique)
#dim(til_ciber)
#[1] 10241    19

#################################################################################

save(net_ciber, gentles_ciber, til_ciber, file="aggregate_cibersorts.rda")
write.table(til_ciber, file="Aggregate_Cibersort_TIL_1.tsv", sep="\t", row.names=F, quote=F)
write.table(net_ciber, file="Aggregate_Cibersort_Net_2.tsv", sep="\t", row.names=F, quote=F)
write.table(gentles_ciber, file="Aggregate_Cibersort_Gentles_3.tsv", sep="\t", row.names=F, quote=F)


#################################################################################
# comparing clusters 4 and 6
# (M2+Monocyte) and (M0+M1+CD8+Tregs)

aggCiber <- data.frame(
  patient_barcode=ciber$pan.samplesID,
  set1 =(ciber$Monocytes+ciber$Macrophages.M2),
  set2 =(ciber$Macrophages.M0 + ciber$Macrophages.M1 + ciber$T.cells.CD8 + ciber$T.cells.regulatory..Tregs)
  )

cluster_comparison <- inner_join(aggCiber, mannot_unique)

save(cluster_comparison, file="cluster_ciber_agg_comparison.rda")

################################################################
##  Subsetting Male multiome data to only contains DEGs ## 
################################################################
############ 1- Train on Multiome DEGs
Top_DEGs_df <- read.csv(".../Top_Degs_GC_PB_multiome.csv")
Top_DEGs <- Top_DEGs_df$genes
################################ 
# Path to GC and PB cells from Multiome dataset
GC <- read.csv(".../Male_Donor/GC.csv", row.names = 1)
PB <- read.csv(".../Male_Donor/PB.csv", row.names = 1)
GC_PB <- rbind(GC, PB)
shared <- intersect(Top_DEGs, colnames(GC_PB))
GC_PB_sub <- GC_PB[,shared]
write.csv(GC_PB_sub, ".../x_train_DEGs.csv", row.names = TRUE)
################################ 
##### finding the subset of scRNA as test
scRNA_GC <- read.csv(".../scRNA_Male_Donor/GC.csv", row.names = 1)
scRNA_PB <- read.csv(".../scRNA_Male_Donor/PB.csv", row.names = 1)
scRNA_GC_PB <- rbind(scRNA_GC, scRNA_PB)
shared <- intersect(Top_DEGs, colnames(scRNA_GC_PB))
scRNA_GC_PB_sub <- scRNA_GC_PB[,shared] 
write.csv(scRNA_GC_PB_sub, ".../x_test_scRNA_DEGs.csv", row.names = TRUE)
##############################################################################################
##############################################################################################
############ 2- Train on scRNA-seq
# finding the subset of scRNA as test
Top_DEGs_df <- read.csv(".../Top_Degs_GC_PB_scRNA.csv")
Top_DEGs <- Top_DEGs_df$genes

scRNA_GC <- read.csv(".../scRNA_Male_Donor/GC.csv", row.names = 1)
scRNA_PB <- read.csv(".../scRNA_Male_Donor/PB.csv", row.names = 1)
scRNA_GC_PB <- rbind(scRNA_GC, scRNA_PB)
shared <- intersect(Top_DEGs, colnames(scRNA_GC_PB))
scRNA_GC_PB_sub <- scRNA_GC_PB[,shared]  
write.csv(scRNA_GC_PB_sub, "/x_train_DEGs.csv", row.names = TRUE)
###############################################################
##### finding the subset of multiome data as test
GC <- read.csv(".../Male_Donor/GC.csv", row.names = 1)
PB <- read.csv(".../Male_Donor/PB.csv", row.names = 1)
GC_PB <- rbind(GC, PB)
shared <- intersect(Top_DEGs, colnames(GC_PB))
GC_PB_sub <- GC_PB[,shared]
write.csv(GC_PB_sub, ".../x_test_multiome_DEGs.csv", row.names = TRUE)


BiocManager::install("Seurat")
library(Seurat)
library(dplyr)
library(Matrix)
library(SingleCellExperiment)
library(data.table)
BiocManager::install("biomaRt")
library(biomaRt)
BiocManager::install("xml2")
library(xml2)

###Seurat pipeline for single-cell data analysis
#QC/ normalisation/ scaling/ dimensonality reduction/ find canonical markers/ automated celltyping/ find marker genes


#Set working directory 
setwd("~/uni_work/queen_mary/individual_project/lineage_dependant_HO_2020/")

#IF user needs to change from Ensembl ID's in features/genes list
# Load ensembl ids
genes_cellranger <- fread('data_copy/E-MTAB-8410.aggregated_filtered_counts.mtx_rows',header=F)

# Load the correct ensembl version
ensembl <- useEnsembl(biomart='ensembl',dataset='hsapiens_gene_ensembl',version=99)

new_cr <- genes_cellranger[, 1]

# Generate mapping to ensembl ids
geneMap <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position','gene_biotype'),mart=ensembl,filters='ensembl_gene_id',values=new_cr$V1)

write.csv(geneMap, file = "genemap.csv")

genemap <- read.csv("genemap.csv", header = FALSE)

genemap <- genemap[, 2:7]
colnames(genemap) <- c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position','gene_biotype')

# Remove duplicate ensembl ids (the mapping to these two genes was outdated/incorrect)
# Checked by googling

geneMap <- genemap[-which(genemap$hgnc_symbol=='C12orf74'),]
geneMap <- geneMap[-which(geneMap$hgnc_symbol=='CCL3L1'),]
geneMap$chromosome_name <- paste0('chr',geneMap$chromosome_name)

# Merge with original ensembl ids (a bit overkill to make sure no duplication)
genes_cellranger <- merge(new_cr,geneMap,by.x='V1',by.y='ensembl_gene_id')
genes_cellranger_1 <- genes_cellranger[,c(1,3:ncol(genes_cellranger))]
colnames(genes_cellranger) <- c('ensembl_gene_id','hgnc_symbol','chromosome','start','end','biotype')

gene_symbols <- genes_cellranger$hgnc_symbol

str(genes_cellranger)
write.csv(genes_cellranger, "gene_data.csv")

# Output to rds (can read back in with readRDS)
saveRDS(genes_cellranger,file='working/compiled_gene_info.rds')

#Can also use online resources to do conversions, change the feature/gene names in spreadsheet software for simplicity
#Biomart method works, but check correct HGNC/ ENSG version used. 3K of 

#Remove extra objects 
rm(counts)
rm(counts.ref)
rm(input.lineage)
rm(input)
rm(pred.hesc)
rm(ref)
rm(reference)
gc()

#Get experimental design/ metadata to add to raw counts
exp_design <- read.csv("~/uni_work/queen_mary/individual_project/lineage_dependant_HO_2020/working/ExpDesign.csv", sep = "\t")

#Read in as a data frame
exp_design <- data.frame(exp_design, row.names = 1)

#Check data structure
exp_design


#Read in 10x data
lineage.data.full <- Read10X(data.dir = "~/uni_work/queen_mary/individual_project/lineage_dependant_HO_2020/raw_data/", gene.column = 2)

gc()

#Check data
lineage.data.full
dim(lineage.data.full)
#24414 genes, 60383 cells
str(lineage.data.full)

#Make seurat object, min cells 10 (Only include features/genes found in at least n many cells ), min features 200 (Only include cells with at least n features)
lineage.data.full <- CreateSeuratObject(counts = lineage.data.full, project = "lineage", meta.data = exp_design, min.cells = 10, min.features = 200)


lineage.data.full <- subset(lineage.data.full, subset = nFeature_RNA > 200 & nFeature_RNA < 7000)

lineage.data.full[["percent.mt"]] <- PercentageFeatureSet(lineage.data.full, pattern = "^MT-")

VlnPlot(lineage.data.full, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")

gc()

lineage.data.full
dim(lineage.data.full)
#21933 genes 60383 cell
str(lineage.data.full)

#Add metadata
lineage.data.full <- AddMetaData(lineage.data.full, exp_design)

#Check metadata
lineage.data.full@meta.data


#Can split data by disease/ state/ location/ relevant metadata. Here data has been split by sample location, tumour core/ border/ non-tumour.
data.tumor <- SplitObject(lineage.data.full, split.by = "Factor.Value.sampling.site.")

#QC
#Find mitochondrial %
lineage.data.full[["percent.mt"]] <- PercentageFeatureSet(lineage.data.full, pattern = "^MT-")

#Check data now includes MT
head(lineage.data.full@meta.data)

#Violin plot for QC statistics
VlnPlot(lineage.data.full, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")

#Do MT calc for spit data
data.tumor[["percent.mt"]] <- PercentageFeatureSet(lineage.data.full, pattern = "^MT-")

#For tumour core samples
data.tumor$`tumour core`[["percent.mt"]] <- PercentageFeatureSet(lineage.data.full, pattern = "^MT-")
#Check
data.tumor$`tumour core`
#19780 cells, 21933 genes
#Plot QC stats
VlnPlot(data.tumor$`tumour core`, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")

#Same for just tumour border
data.tumor$`tumour border`[["percent.mt"]] <- PercentageFeatureSet(lineage.data.full, pattern = "^MT-")
data.tumor$`tumour border`
#20786 cells, 21933 genes
VlnPlot(data.tumor$`tumour border`, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")

#Same for just normal tissue
data.tumor$`normal tissue adjacent to neoplasm`[["percent.mt"]] <- PercentageFeatureSet(lineage.data.full, pattern = "^MT-")
data.tumor$`normal tissue adjacent to neoplasm`
#19817 cells, 21933 genes
VlnPlot(data.tumor$`normal tissue adjacent to neoplasm`, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")


#Compare QC plots for clearer image
plot1 <- FeatureScatter(lineage.data.full, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
plot2 <- FeatureScatter(lineage.data.full, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
plot1 + plot2


#Remove extra objects out of memory
rm(plot1)
rm(plot2)
gc()

install.packages("ggplot2")
library(ggplot2)

ggplot(lineage.data.full@meta.data, aes(x = lineage.data.full$percent.mt)) + geom_histogram(binwidth = 0.5, fill="orange", colour="black") + 
  ggtitle("Distribution of Percentage Mitochondrion") + geom_vline(xintercept = 10)


#Subset the data removing low-quality cells. Keep only cells with between X and Y features, < % MT ratio, RNA counts > Z 
lineage.data.full <- subset(lineage.data.full, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20 & nCount_RNA > 1000)

lineage.data.full
dim(lineage.data.full)
#21933 genes 33921 cells
str(lineage.data.full)


#Same for split data
data.tumor <- subset(lineage.data.full, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 10 & nCount_RNA > 1000)


#Same cells remaining
dim(data.tumor)
gc()

hist(colSums(lineage.data.full@assays$RNA@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")

hist(colSums(data.tumor@assays$RNA@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")


#Normalise data
lineage.data.full <- NormalizeData(lineage.data.full)

data.tumor <- NormalizeData(data.tumor)

lineage.data.full

#Find 2000 most variable genes
lineage.data.full <- FindVariableFeatures(lineage.data.full, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(lineage.data.full), 10)

write.csv(top10, "top10.csv")
saveRDS(top10, "top10.rds")


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(lineage.data.full)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2
rm(top10)
rm(plot1)
rm(plot2)
gc()

#Centre and Scale data
lineage.data.full <- ScaleData(lineage.data.full)


hist(colSums(lineage.data.full@assays$RNA@data),
     breaks = 100,
     main = "Total expression after normalisation",
     xlab = "Sum of expression")

gc()
# Run PCA
lineage.data.full <- RunPCA(lineage.data.full, features = VariableFeatures(object = lineage.data.full))

#Plot PCA
DimPlot(lineage.data.full, reduction = "pca")

#List of cell cycle genes (From Seurat)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes


# Score cell-cycle
lineage.data.full <- CellCycleScoring(lineage.data.full, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell-cycle scores and phase assignments
head(lineage.data.full[[]])

#Can plot expression of all cell cycle genes
RidgePlot(lineage.data.full, features = s.genes[31:43], ncol = 4)
#All expressed but mostly equally
#CDK1, MKI67, KIF11, CDCA3, CDC20, RANGAP1, DLGAP5, CDCA2, CKAP5, CENPE, CBX5  
#FEN1, NASP, GMNN, MSH2, E2F8

#Plot PCA, colour by various metadata
DimPlot(lineage.data.full, reduction = "pca", group.by = "Factor.Value.sampling.site.")

DimPlot(lineage.data.full, reduction = "pca", group.by = "Sample.Characteristic.organism.part.")

DimPlot(lineage.data.full, reduction = "pca", group.by = "Phase")


# Run PCA with cell-cycle data
cell.cycle.pca <- RunPCA(lineage.data.full, features = c(s.genes, g2m.genes))
DimPlot(cell.cycle.pca)

gc()
rm(cell.cycle.pca)

#Set variables to regress out as counts, mitochondrial and cell-cycle
vars_to_regress <- c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score", "Sample.Characteristic.organism.part.", "Sample.Characteristic.sex.")

#Scale data regressing out uninteresting variation
lineage.data.full <- ScaleData(lineage.data.full, vars.to.regress = vars_to_regress)

#Re-run PCA
lineage.data.full <- RunPCA(lineage.data.full, features = VariableFeatures(object = lineage.data.full))

#Plot PCA
DimPlot(lineage.data.full)
DimPlot(lineage.data.full, reduction = "pca", group.by = "Factor.Value.sampling.site.")



print(lineage.data.full[["pca"]], dims = 1:41, nfeatures = 5)

VizDimLoadings(lineage.data.full, dims = 1:5, reduction = "pca")

DimHeatmap(lineage.data.full, dims = 1:10, balanced = TRUE)

ElbowPlot(lineage.data.full, ndims = 50)

gc()
rm(exp_design)


# Find neighbours for clustering
lineage.data.full <- FindNeighbors(lineage.data.full, dims = 1:27)

#Do clustering
lineage.data.full <- FindClusters(lineage.data.full, resolution = 0.6)

install.packages("Rtsne")
library(Rtsne)

# Find tSNE 
lineage.data.full <- RunTSNE(object = lineage.data.full, dims = 1:27, check_duplicates = F)

#Plot tSNE
TSNEPlot(object = lineage.data.full)

TSNEPlot(object = lineage.data.full, group.by = "seurat_clusters", label = TRUE)

#Write out clusters and celltypes to check compatibility
write.csv(lineage.data.full$seurat_clusters, "clusters.csv")
#Used later for comparison after celltyping
#write.csv(lineage.data.full$celltypes, "celltypes.csv")


saveRDS(lineage.data.full, "lineage.data.full.rds")

lineage.data.full <- readRDS("~/uni_work/queen_mary/individual_project/lineage_dependant_HO_2020/lineage.data.full.rds")


library(reticulate)
#Old method, still works but unnecessary 
#reticulate::py_install(packages ='umap-learn')

#Run UMAP, use same dimensions as before
lineage.data.full <- RunUMAP(lineage.data.full, dims = 1:27)

#Plot UMAP, colour by clusters
DimPlot(lineage.data.full, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

#Used later after celltyping
#DimPlot(lineage.data.full, reduction = "umap", group.by = "celltypes", label = TRUE)


#saveRDS(lineage.data, "lineage.data.rds")


lineage.data.full <- readRDS("~/uni_work/queen_mary/individual_project/lineage_dependant_HO_2020/lineage.data.full.rds")

BiocManager::install('limma')
library(limma)

#Find markers for each cluster, compared to other clusters. Useful if using canonical markers, can take 1-3 mins per cluster
#cluster1.markers <- FindMarkers(lineage.data.full, ident.1 = 1, min.pct = 0.25)
#head(cluster1.markers, n = 5)

#Can choose to include/ exclude custers from comparison
#cluster1.not.20.or.0 <- FindMarkers(lineage.data, ident.1 = 1, ident.2 = c(0, 20), min.pct = 0.25)
#head(cluster1.not.20.or.0, n = 5)


#RE-RUN, Finds all markers for all custers, takes ~3min per cluster
#pbmc.markers <- FindAllMarkers(lineage.data.full, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

saveRDS(pbmc.markers, "pbmc.markers.rds")
saveRDS(top5, "top5.rds")

write.csv(top5, file = "top5.csv")
write.csv(pbmc.markers, file = "cellDEmarkers.csv")


pbmc.markers <- readRDS("~/uni_work/queen_mary/individual_project/lineage_dependant_HO_2020/pbmc.markers.rds")

str(pbmc.markers)

#Can change test/ parameters as fit 
cluster1.markers <- FindMarkers(lineage.data, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#Plot expression for canonical celltype marker genes                                  
VlnPlot(lineage.data.full, features = c("IL7R", "CCR7"))

#Plot canonical celltype marker genes using UMAP                                                
FeaturePlot(lineage.data.full, features = c(features.plot = c("IL7R", "CD14", "LYZ", "MS4A1", "CD8A", "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP")))

#Plot cannonical celltype marker genes using tsne
FeaturePlot(lineage.data.full, reduction = "tsne", features = c("MKI67", "CD14", "LYZ", "MS4A1", "CD8A", "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP"))

#Find top10 markers per cluster based on log2 fold change
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(lineage.data, features = top10$gene) + NoLegend()
         
write.csv(pbmc.markers, "pbmc.markers.full.csv")
write.csv(top10, "top10.csv")


###Automated Celltyping###

BiocManager::install("celldex")  
library(celldex)

#Get Human Cell Atlas data from celldex
hpca.se <- HumanPrimaryCellAtlasData()
#Celltypes avaliable as reference data set
hpca.se$label.main

#Remove unnecessary variables/ objects, next step is memory intensive 

gc()

###SingleR### 
BiocManager::install("SingleR")  
library(SingleR)

#Predicts celltype based on gene expression, using HPCA data as reference, takes ~10 mins for 30k sample cells on 16gb RAM, 8th gen I7 
pred.hesc <- SingleR(test = lineage.data.full@assays$RNA@counts, ref = hpca.se, assay.type.test=1, labels = hpca.se$label.main)

saveRDS(pred.hesc, "pred.hesc.rds")

pred.hesc <- readRDS("~/uni_work/queen_mary/individual_project/lineage_dependant_HO_2020/pred.hesc.rds")

#Check data
pred.hesc

#Sort by predicted celltype frequencey
sort(table(pred.hesc$pruned.labels), decreasing = TRUE)

#See how many are unpredicted
summary(pruneScores(pred.hesc))

#See scores
pruneScores(pred.hesc, get.thresholds = TRUE)

pruneScores(
  pred.hesc,
  nmads = 3,
  min.diff.med = -Inf,
  min.diff.next = 0,
  get.thresholds = TRUE
)

write.csv(pred.hesc, "pred.hesc.csv")

#Plot a score heatmap
plotScoreHeatmap(pred.hesc)

#Plot distribution of delta (Difference)
plotDeltaDistribution(pred.hesc, ncol = 3)

###CHETAH###

BiocManager::install("CHETAH")
library(CHETAH)
library(SingleCellExperiment)

#load reference data set, can be found online (https://figshare.com/s/aaf026376912366f81b6)
load("~/uni_work/queen_mary/individual_project/lineage_dependant_HO_2020/CHETAH_TME_reference.Rdata")

#Name counts variable from user dataset
counts <- lineage.data.full@assays$RNA@counts

#Name user data input, requires counts and dimensionality reduction
input <- SingleCellExperiment(assays = list(counts = counts), reducedDims = SimpleList(UMAP = lineage.data.full@reductions$umap@cell.embeddings))

#Need reference celltypes
celltypes.ref <- reference$celltypes

#Need reference counts
counts.ref <- assay(reference)

#Combine reference data
ref <- SingleCellExperiment(assays = list(counts = counts.ref),
                                     colData = DataFrame(celltypes = celltypes.ref))

#Classify. Takes ~5 mins for 30k user cells with 16gb RAM and 8th gen I7
input.lineage <- CHETAHclassifier(input = input,
                              ref_cells = ref)

#Plot CHETAH graph
PlotCHETAH(input.lineage)

PlotCHETAH(input = input.lineage, interm = TRUE)

#Keep celltype predictions
celltypes <- input.lineage$celltype_CHETAH

#Check data
celltypes 

#Sort by prediction frequency, lots of unassigned cells
sort(table(celltypes), decreasing = TRUE)

#Can change parameters to allow different confidence, lower = less confident  
input.lineage <- Classify(input = input.lineage, 0.1)

PlotCHETAH(input = input.lineage, tree = FALSE)

#Add predicted celltypes to Seurat dataframe metadata
lineage.data.full$SRcelltypes <- pred.hesc$pruned.labels
lineage.data.full$CHcelltypes <- celltypes

#Check data
lineage.data.full$SRcelltypes
lineage.data.full$CHcelltypes

saveRDS(lineage.data.full, "lineage.data.full.rds")

#Write as CSV to compare celltype predictions in spreadsheet software
write.csv(lineage.data.full$CHcelltypes, "CHcelltypes.csv")
write.csv(lineage.data.full$SRcelltypes, "SRcelltypes.csv")

#Take original predictions from paper authors to compare in spreadsheet
write.csv(lineage.data.full$Factor.Value.inferred.cell.type...authors.labels., "authorscelltypes.csv")


lineage.data.full <- readRDS("~/uni_work/queen_mary/individual_project/lineage_dependant_HO_2020/lineage.data.full.rds")

#Display predicted celltypes via umap, grouped by SingleR, CHETAH, Seurat cluster
DimPlot(lineage.data.full, reduction = "umap", group.by = "SRcelltypes", label = TRUE, repel = TRUE)
DimPlot(lineage.data.full, reduction = "umap", group.by = "CHcelltypes", label = TRUE, repel = TRUE)
DimPlot(lineage.data.full, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)


#Display predicted celltypes via tsne, grouped by SingleR, CHETAH, Seurat cluster
DimPlot(lineage.data.full, reduction = "tsne", group.by = "SRcelltypes", label = TRUE, repel = TRUE)
DimPlot(lineage.data.full, reduction = "tsne", group.by = "CHcelltypes", label = TRUE, repel = TRUE)
DimPlot(lineage.data.full, reduction = "tsne", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

#Save metadata (Includes SingleR/ CHETAH/ Seurat clustering)
write.csv(lineage.data.full@meta.data, "metadata.csv")

#Split data by sample site (Name will change for different metadata, but allows spitting of data based on where sample was taken from)
data.split <- SplitObject(lineage.data.full, split.by = "Factor.Value.sampling.site.")

#Can plot just specific samples based on location
DimPlot(data.split$`tumour core`, reduction = "umap", group.by = "SRcelltypes", label = T, repel = T)
DimPlot(data.split$`tumour border`, reduction = "umap", group.by = "SRcelltypes", label = T, repel = T)
DimPlot(data.split$`normal tissue adjacent to neoplasm`, reduction = "umap", group.by = "SRcelltypes", label = T, repel = T)

gc()

#Can assign celltypes to seurat object 
SRcelltypes <- lineage.data.full$SRcelltypes
SRcelltypes <- as.data.frame(SRcelltypes)
Idents(lineage.data.full) <- "SRcelltypes"
Idents(lineage.data.full)

#Can now use SingleR celltype predictions to find marker genes for each celltype
T.markers <- FindMarkers(lineage.data.full, ident.1 = "T_cells", min.pct = 0.25)
T.markers

# CAN NOW SUBSET INTO CELL-TYPES 
# Can recombine the various subsets to immune, stromal, tumour

#Subset into cell-type only, runs as each celltype vs all others combined!
pbmc.markers <- FindAllMarkers(lineage.data.full, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(pbmc.markers, "cell_type_markers.csv")

#Check set idents are the "active" idents
Idents(lineage.data.full)
lineage.data.full@active.ident

#Used for presentation graphics based on manual approximation of most common celltypes in each seurat cluster
#lineage.data.full <- RenameIdents(object = lineage.data.full, "0" = "T Cells/ NK Cells", "1" = "B Cells/ Pre-B Cell CD34-/ Pro-B Cell CD34+",
 #                                 "2" = "Monocyte/ Macrophage", "3" = "B Cells", "4" = "Smooth Muscle/ Fibroblast/ Chondrocytes", 
  #                                "5" = "Epithelial Cells", "6" = "Tissue Stem Cells/ Smooth Muscle Cells/ Fibroblasts",
   #                               "7" = "Epithelial Cells", "8" = "Endothelial Cells", "9" = "Epithelial Cells", "10" = "Epithelial Cells",
    #                              "11" = "Tissue Stem Cells/ Fibroblasts", "12" = "T Cells/ B Cells", "13" = "CMP/ NK Cells", 
     #                             "14" = "Neurons/ Dendritic")

##Used for presentation graphics based on manual approximation of most common celltypes in each seurat cluster
#DimPlot(lineage.data.full, reduction = "umap", label = TRUE, repel = TRUE)
#DimPlot(lineage.data.full, reduction = "tsne", label = TRUE, repel = TRUE)

#Check what celltypes are predicted as present
levels(lineage.data.full)

#Rename epithelial cells from tumour core as tumour cells (Treating all epithelial cells from tumour as being "tumour cells")
data.tumor <- SplitObject(lineage.data.full, split.by = "Factor.Value.sampling.site.")
Idents(data.tumor$`tumour core`)
#Change name
lineage.data.celltyped <- RenameIdents(object = data.tumor$`tumour core`, `Epithelial_cells` = "Tumour_Cells")
#Check name change
Idents(lineage.data.celltyped)
#Check data
sort(table(lineage.data.celltyped@active.ident), decreasing = TRUE)
sort(table(data.tumor$`tumour core`@active.ident), decreasing = TRUE)

#Merge data from tumour core back into full data
full.data <- merge(x = lineage.data.celltyped, y = list(data.tumor$`tumour border`, data.tumor$`normal tissue adjacent to neoplasm`))
#Check data
sort(table(full.data@active.ident), decreasing = TRUE)

rm(lineage.data.celltyped)
rm(data.tumor)
rm(lineage.data.full)
gc()

#Find cell markers for all data, including name changes
tumour_cell.markers <- FindAllMarkers(full.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.75)

#Write as csv for manual gene identification via spreadsheet
write.csv(tumour_cell.markers, "tumour_cell_markers2.csv")

                                    

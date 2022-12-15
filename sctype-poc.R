# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

# load processed covid dataset
scdata <- readRDS("./covid_data.rds")

# cluster and visualize
scdata <- FindNeighbors(scdata, dims = 1:10)
scdata <- FindClusters(scdata, resolution = 0.4)
scdata <- RunUMAP(scdata, dims = 1:10)
DimPlot(scdata, reduction = "umap")

# convert ENS ID to gene symbol and set gene symbols as rownames of the count matrix
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(scdata)
symbols <- NA
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                          "entrezgene_description", "hgnc_symbol"),values=genes,mart= mart)

df <- scdata[["RNA"]]@scale.data
df.conv <- merge(df,G_list,by.x="row.names",by.y="ensembl_gene_id")
df.conv.nd=df.conv[which(!duplicated(df.conv$hgnc_symbol)),]
rownames(df.conv.nd) <- df.conv.nd$hgnc_symbol
df.conv.nd <- subset(df.conv.nd, select=-c(hgnc_symbol,Row.names))
df.conv.nd <- as.matrix(sapply(df.conv.nd, as.numeric))
temp <- df.conv[which(!duplicated(df.conv$hgnc_symbol)),]
rownames(df.conv.nd) <- temp$hgnc_symbol

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = df.conv.nd, scaled = TRUE,
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix.
# In case Seurat is used, it is either scdata[["RNA"]]@scale.data (default), scdata[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or scdata[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_results = do.call("rbind", lapply(unique(scdata@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(scdata@meta.data[scdata@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scdata@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

scdata@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  scdata@meta.data$customclassif[scdata@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(scdata, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif',label.size=4)


# load auto-detection function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";

# guess a tissue type
tissue_guess = auto_detect_tissue_type(path_to_db_file = db_, seuratObject = scdata, scaled = F, assay = "RNA")  # if scaled = TRUE, make sure the data is scaled, as seuratObject[[assay]]@scale.data is used. If you just created a Seurat object, without any scaling and normalization, set scaled = FALSE, seuratObject[[assay]]@counts will be used


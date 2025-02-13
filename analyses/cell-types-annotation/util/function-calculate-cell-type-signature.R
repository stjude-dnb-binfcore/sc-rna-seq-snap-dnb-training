#' Function to calculate cell type signature
#' @param seurat_obj
#' @param results_dir
#' @param gene_markers_df
#' @param genome_name
#' 
#' @return
#' @export
#'
#' @examples
#' 
calculate_cell_type_signature <- function(seurat_obj, results_dir, gene_markers_df, genome_name) {
  #calculate_cell_type_signature <- function(Seurat.obj, Project.Path, Figure.Total, Analysis.Name, Integration.Method, Cell.Type.DF, Genome) {
  
  pdf(file = paste(results_dir, "/", "all_cell_type_signatures.pdf", sep = ""), width = length(gene_markers_df)*6, height = 6)
  
  for (i in 1:ncol(gene_markers_df)) {
    cell.type <- paste(colnames(gene_markers_df)[i], ".score", sep = "")
    gene.markers <- as.character(gene_markers_df[, i][gene_markers_df[, i] != ""])
    
    if (genome_name == "GRCh38" | genome_name == "hg19"){
      gene.markers <- toupper(gene.markers)
      
    } else if (genome_name == "mm10"){
      gene.markers <- str_to_title(gene.markers)
      
    } else if (genome_name == "Dualhg19"){
      gene.markers <- paste("hg19-", toupper(gene.markers), sep = "")
      
    } else if (genome_name == "DualGRCh38"){
      gene.markers <- paste("GRCh38-", toupper(gene.markers), sep = "")
      
    } else if (genome_name == "Dualmm10"){
      gene.markers <- paste("mm10---", str_to_title(gene.markers), sep = "") 
      
    } else if (genome_name == "GRCm39") {
      gene.markers <- paste("GRCm39---", str_to_title(gene.markers), sep = "") 
      
    } else if (genome_name == "DualGRCm39") {
      gene.markers <- paste("DualGRCm39---", str_to_title(gene.markers), sep = "") } 

    
    # sce <- as.SingleCellExperiment(seurat_obj)
    # gene.markers <- gene.markers[gene.markers %in% getGenes(sce)]
    #gene.markers <- gene.markers[gene.markers %in% rownames(seurat_obj)]
    gene.markers <- gene.markers[gene.markers %in% getGenes(seurat_obj)]
    head(gene.markers)
    seurat_obj <- AddModuleScore(seurat_obj, features = list(gene.markers), name = cell.type)
    
    # Print plots  
    print(FeaturePlot(seurat_obj, features = paste(cell.type, "1", sep=""), label = TRUE) + theme(aspect.ratio = 1))
    print(VlnPlot(seurat_obj, features = paste(cell.type, "1", sep = ""), pt.size = 0))
    print(VlnPlot(seurat_obj, features = paste(cell.type, "1", sep = ""), pt.size = 0, group.by = "ID")) }
  
  cell_cluster_scores <- seurat_obj@meta.data[, grepl(".score1", colnames(seurat_obj@meta.data))]
  
  # Predict.Cell.Score <- function(Data.Frame)
  # {
  #   Data.Frame <- sort(Data.Frame, decreasing = TRUE)
  #   if (foldchange(Data.Frame[1], Data.Frame[2]) > 1)
  #   {
  #     celltype <- sub(".score1", "", names(Data.Frame[1]))
  #   }
  #   else
  #   {
  #     celltype <- "UNDETERMINED"
  #   }
  #   
  #   return(celltype)
  # }
  # seurat_obj$predicted.cell.signature.ident <- apply(cell_cluster_scores, MARGIN = 1, FUN = Predict.Cell.Score)
  
  seurat_obj$predicted.cell.signature.ident <- sub(".score1", "", colnames(cell_cluster_scores)[max.col(cell_cluster_scores, ties.method = "first")])
  
  # Print plots
  print(DimPlot(seurat_obj, reduction = "umap", group.by = "predicted.cell.signature.ident") + theme(aspect.ratio = 1))
  print(DimPlot(seurat_obj, reduction = "umap", group.by = "Phase") + theme(aspect.ratio = 1))
  print(DimPlot(seurat_obj, reduction = "umap", group.by = "ID") + theme(aspect.ratio = 1))
  print(DimPlot(seurat_obj, reduction = "umap", split.by = "predicted.cell.signature.ident") + theme(aspect.ratio = 1))
  print(ggplot(data = data.frame(table(seurat_obj$predicted.cell.signature.ident)), aes(x = Var1, y = Freq)) + geom_bar(stat = "identity") + geom_text(aes(label = Freq)))
  
  dev.off()
  
  # Save tables
  write.table(table(seurat_obj$predicted.cell.signature.ident, seurat_obj$ID), file = paste(results_dir, "/", "all_cell_type_signatures_count_data.tsv", sep = ""), sep="\t", quote = FALSE)
  write.table(prop.table(table(seurat_obj$predicted.cell.signature.ident, seurat_obj$ID), margin = 2) *100, file = paste(results_dir, "/", "all_cell_type_signatures_percentage_data.tsv", sep = ""), sep="\t", quote=FALSE)
  
  return(seurat_obj)
}


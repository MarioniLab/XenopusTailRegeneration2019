#Load reticulate and specify python directory
library("reticulate")
use_python("/Users/hiscoc01/miniconda3/bin/python3", required = T)
py_config()

#Load other packages
library(ggplot2)
library(Matrix)
library(igraph)
library(AUCell)
library(edgeR)
library(ggrepel)
library(reshape2)
library(Seurat)

#Define python variables and functions
umap <- import("umap")
np <- import("numpy", convert = FALSE)
scp <- import("scipy")
py_run_string("
import umap
import scipy as scp
import sklearn as sk
from sklearn.utils import check_random_state
from sklearn.utils import check_array
from scipy.sparse import csr_matrix, find
import numpy as np
              ")

py_run_string("
import umap
def get_umap(y, n, minDist,seed):
    return umap.UMAP(n_neighbors=n,min_dist=minDist,metric='cosine', n_components = 2, random_state = seed).fit_transform(y)
              ")

py_run_string("
def get_umap_G(y,n, seed):
  random_state = check_random_state(seed)
  Y = check_array(y).astype(np.float64)
  G = umap.umap_.fuzzy_simplicial_set(X=Y,n_neighbors=n,random_state=random_state,metric='cosine')
  return find(G)
              ")

#Normalize by total counts
normalize <- function(counts, pseudocount = 1e4){
  librarySize <- Matrix::colSums(counts)
  return(pseudocount * t(t(counts) / librarySize ))
}

#Select HVGs based on mean expression and fano factor
compute_hvg <- function(counts.normalized, lowmean = 0.05, highmean = 0.8, fano = 0.65){
  mean.expression <- Matrix::rowMeans(counts.normalized)
  meansq.expression <- Matrix::rowMeans(counts.normalized^2)
  fano.factor <- (meansq.expression - (mean.expression^2))/mean.expression
  hvg <- counts.normalized[which((mean.expression >= quantile(mean.expression, lowmean)) & (fano.factor >= quantile(na.omit(fano.factor),fano)) & (mean.expression <= quantile(mean.expression,highmean))),]
  return(hvg)
}

#Project data using UMAP
umap_project <- function(meta, hvg, neighbours = 10, dist = 0.3, seed = 42){
  y <- np$log2(as.matrix(hvg) + 1)$T
  out <- py$get_umap(y,as.integer(neighbours),dist, as.integer(seed))
  meta$x <- out[,1]
  meta$y <- out[,2]
  return(meta)
}

#Cluster data using UMAP
umap_cluster <- function(meta, hvg, neighbours = 10, steps = 4, seed1 = 42, seed2 = 42){
  y <- np$log2(as.matrix(hvg) + 1)$T
  x <- py$get_umap_G(y,as.integer(neighbours),as.integer(seed1))
  sp <- sparseMatrix(i = (x[[1]]+1), j = (x[[2]]+1), x = as.vector(x[[3]]), dims = c(dim(hvg)[2],dim(hvg)[2]))
  g <- graph_from_adjacency_matrix(sp,weighted = TRUE,mode = "undirected")
  set.seed(seed2)
  meta$cluster <- cluster_walktrap(g,steps = steps)$membership
  return(meta)
}

#Plot UMAP projection in 2D - with or without clusters colored differently
plotMeta <- function(meta, mode = "none"){
  if(mode == "none"){
    p <- ggplot( data = meta, mapping= aes(x = x, y = y)) +
      theme_bw() +
      geom_point(size=0.5, col = "black") +
      coord_equal() +
      theme(aspect.ratio = 1)+
      theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1) 
  }
  
  else if(mode == "cluster"){
    p <- ggplot( data = meta, mapping= aes(x = x, y = y, col = factor(meta$cluster))) +
      theme_bw() +
      geom_point(size=0.5) +
      coord_equal() +
      theme(aspect.ratio = 1)+
      theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1)  +
      theme(legend.position = "none")
  }
  return(p)
}

#Plot UMAP projection in 2D, with cluster names and colours specified by "clusters"
plotAnnotation <- function(meta, clusters){
   
  p <-  ggplot( data = meta, mapping= aes(x = X, y = Y, col = factor(meta$cluster, levels = clusters$names))) +
      theme_bw() +
      geom_point(size=0.5) +
      scale_colour_manual(values = sapply(clusters$cols, as.character), labels= sapply(clusters$names, as.character), name = "") +
      coord_equal() +
      theme(aspect.ratio = 1)+
      theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1)  +
      guides(colour = guide_legend(override.aes = list(size=3,
                                                       alpha = 1))) 
  
  print(p)
  
  p <-  ggplot( data = meta, mapping= aes(x = X, y = Y, col = factor(meta$cluster, levels = clusters$names))) +
    theme_bw() +
    geom_point(size=0.5) +
    scale_colour_manual(values = sapply(clusters$cols, as.character), labels= sapply(clusters$names, as.character), name = "") +
    coord_equal() +
    theme(aspect.ratio = 1)+
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1)  +
    guides(colour = guide_legend(override.aes = list(size=3,
                                                     alpha = 1))) + theme(legend.position = "none")
  return(p)
}


#Plot a boxplot of values (e.g. gene expression levels) separated by groups (e.g. cell types)
scBox <- function(values, groups, clusters, title = "title"){
  
  df <- data.frame(score = as.vector(values),cluster = factor(groups, sapply(clusters$names, as.character)))
  p <- ggplot(df,aes(x = cluster, y=score, fill = factor(groups, sapply(clusters$names, as.character))))+
    geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=0.4, notch=FALSE)+
    theme_bw() +
    scale_fill_manual(values = sapply(clusters$cols, as.character), labels= sapply(clusters$names, as.character), name = "") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.title = element_blank()) + ggtitle(title)
  print(p)
}

#Plot cell cycle fractions for different groups (e.g. clusters)
plotCellCycle <- function(meta, clusters, title = "title"){
  output <- c()
  for (category in sort(unique(meta$CellCyclePhase))){
    output <- rbind(output,aggregate(meta$CellCyclePhase == category,list(meta$cluster),sum)$x)
  }
  rownames(output) <- sort(unique(meta$CellCyclePhase))
  colnames(output) <- aggregate((meta$CellCyclePhase == category),list(meta$cluster),sum)$Group.1
  output <- output[,match(clusters$names,colnames(output))]
  
  output <- output[c(1,3,2),]
  
  par(mfrow=c(1, 1), mar=c(10, 5, 4, 8)); barplot(prop.table(output,2),las=2,legend.text = TRUE,args.legend = list(x = "topright", bty = "n", inset=c(-0.15, 0)), cex.names = 0.65, col = c("grey","#377eb8","#e41a1c")); title(title)
}

#Calculate the fraction of cycling cells for different conditions, for each cluster. We ignore clusters with low abundances for which this estimate will be highly noisy.
getCycling <- function(meta, clusters){
  output <- c()
  for (category in sort(unique(meta$CellCyclePhase))){
    output <- rbind(output,aggregate(meta$CellCyclePhase == category,list(meta$cluster),sum)$x)
  }
  rownames(output) <- sort(unique(meta$CellCyclePhase))
  colnames(output) <- aggregate((meta$CellCyclePhase == category),list(meta$cluster),sum)$Group.1
  output <- output[,match(clusters$names,colnames(output))]
  output <- output[c(1,3,2),]
  numbers <- apply(output,2,sum)
  proliferating <- 1-prop.table(output,2)[1,]
  proliferating[numbers < 10] <- NA
  return(proliferating)
}

#Overlay different samples on the UMAP projection
plotSample <- function(m1, title, mode = "none"){
  umap <- as.data.frame(cbind(meta$X,meta$Y))
  umap1 <- as.data.frame(cbind(meta$X[m1],meta$Y[m1], meta$sample[m1]))
  if(mode == "none"){ cols <- "red"}
  else if(mode == "batch"){ cols <- factor(umap1$V3)}
  p <- ggplot( data = umap, mapping= aes(x = V1, y = V2)) +
    theme_bw() +
    geom_point(size=0.5, col = "grey") +
    scale_shape(solid = F) +
    geom_point(data = umap1, mapping = aes(x = V1, y = V2, col = cols),size=0.5, fill = NA) +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1)  +
    guides(colour = guide_legend(override.aes = list(size=3,
                                                     alpha = 1)))  + ggtitle(title) + theme(legend.position = "none")
  print(p)
}

#Plot gene expression level for a subset of cells, overlaid on the entire map
plotGeneSubset <- function(meta, countn, geneName,condition, title){
  XA <- meta$X
  YA <- meta$Y
  X <- meta$X[condition]
  Y <- meta$Y[condition]
  atlas <- as.data.frame(cbind(XA, YA))
  all <- as.data.frame(cbind(X, Y))
  gene <- countn[which(rownames(countn) == geneName), condition]
  if (length(gene) > 0) {
    log_counts = log10(gene + 1)
    Xsub <- X[log_counts > 0]
    Ysub <- Y[log_counts > 0]
    log_counts = log_counts[log_counts > 0]
    m <- order(log_counts)
    Xsub <- Xsub[m]
    Ysub <- Ysub[m]
    log_counts <- log_counts[m]
    if (length(Xsub) > 0) {
      subset <- as.data.frame(cbind(Xsub, Ysub))
      p <- ggplot(data = atlas, mapping = aes(x = XA, y = YA)) +
        geom_point(size = 0.1, col = alpha("grey", 1.0)) +
        geom_point(
          data = all,
          mapping = aes(x = X, y = Y),
          size = 1,
          shape = 21,
          colour = "black",
          stroke = 0.1,
          fill = "grey"
        ) +
        geom_point(
          data = subset,
          mapping = aes(x = Xsub, y = Ysub, fill = log_counts),
          size = 1,
          shape = 21,
          colour = "black",
          stroke = 0.1
        ) +
        scale_fill_gradient(low = "yellow",
                            high = "red",
                            limits = c(0, 1)) +
        theme_bw() +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = 1
        ) +
        ggtitle(paste0(geneName,":   ",title))
      print(p)
      
    }
  }
}




#Compute differential abundance statistics between conditions
differentialAbundance <- function(meta,labels,sample1,sample2,clusters, plot = T,minN = 5, pthresh = 0.05){
  clusterColours <- sapply(clusters$cols,as.character)
  clusterNames <- sapply(clusters$names,as.character)
  counts <- table(meta$cluster, meta$sample)
  counts <- counts[which(rowMeans(counts) > minN),]
  group <- factor(labels$Condition[match(colnames(counts),labels$Sample)])
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)

  y <- DGEList(counts=counts,group=group)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)

  m1 <- match(sample1, colnames(design))
  m2 <- match(sample2, colnames(design))
  
  contrast <- rep(0,dim(design)[2])
  contrast[m2] <- 1
  contrast[m1] <- -1
  lrt <- glmQLFTest(fit, contrast = contrast)$table
  
  subset <- lrt[which(lrt$PValue < pthresh  ),]
  if (plot){
    p <- ggplot(lrt, aes(x=logFC, y=-log10(PValue), col = factor(rownames(lrt), levels = rownames(lrt)))) + 
      theme_bw() +
      xlim(-2.5,2.5) +
      ylim(0,1.5) +
      geom_point(size=2) +
      geom_text_repel(data = subset, mapping= aes(x = logFC, y = -log10(PValue), label=rownames(subset),col = factor(rownames(subset), levels = rownames(subset))), size=3) +
      
      scale_colour_manual(values = clusterColours[match(rownames(lrt),(clusterNames))], labels= rownames(lrt), name = "") +
      theme(  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1)  +
      theme(aspect.ratio = 1) +
      guides(colour = guide_legend(override.aes = list(size=3,
                                                       alpha = 1))) +
      labs(x = "log fold change in abundance", y = "-log P-value") 
    print(p)
  }
  return(lrt)
}

#Plot gene expression values on the UMAP projection
plotGene <- function(meta, countn, geneName){
  
  X <- meta$X
  Y <- meta$Y
  atlas <- as.data.frame(cbind(X,Y))
  if(!(geneName %in% rownames(countn))){
    return()
  }
  gene <- countn[which(rownames(countn) == geneName),]
  if(max(gene) == 0){ return() }
  log_counts = log10(gene+1)
  Xsub <- X[log_counts>0]
  Ysub <- Y[log_counts>0]
  log_counts = log_counts[log_counts>0]
  m <- order(log_counts)
  Xsub <- Xsub[m]
  Ysub <- Ysub[m]
  log_counts <- log_counts[m]
  subset <- as.data.frame(cbind(Xsub, Ysub))
  p <- ggplot( data = atlas, mapping=aes(x=X,y=Y)) +
    geom_point(size=0.2,col=alpha("grey",1)) +
    geom_point(data = subset,mapping = aes(x=Xsub, y=Ysub, col = log_counts), size = 0.2) +
    theme_bw() +
    # scale_colour_gradient2(low="#ff6666", mid = "#e60000", high = "#990000") + 
    scale_colour_gradient(low="yellow", high = "red") + 
    
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", aspect.ratio = 1) +
  ggtitle(geneName)

  return(p)
  
}

#Plot UMI on the UMAP projection
plotUMI <- function(meta){
  
  X <- meta$X
  Y <- meta$Y
  atlas <- as.data.frame(cbind(X,Y))
  log_counts = log10(meta$umi+1)
  Xsub <- X[log_counts>0]
  Ysub <- Y[log_counts>0]
  log_counts = log_counts[log_counts>0]
  m <- order(log_counts)
  Xsub <- Xsub[m]
  Ysub <- Ysub[m]
  log_counts <- log_counts[m]
  subset <- as.data.frame(cbind(Xsub, Ysub))
  p <- ggplot( data = atlas, mapping=aes(x=X,y=Y)) +
    geom_point(size=0.2,col=alpha("grey",1)) +
    geom_point(data = subset,mapping = aes(x=Xsub, y=Ysub, col = log_counts), size = 0.2) +
    theme_bw() +
    scale_colour_gradient(low="yellow", high = "red") + 
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", aspect.ratio = 1) +
    ggtitle("UMI")
  print(p)
  return()
  
}

#Given a set of gene names, create .L and .S versions (long and short chromosomes, since X. Laevis is pseudotetraploid)
xenifyGenes <- function(geneList){
  
  ligands.L<-sapply(tolower((geneList)),function(x) paste0(x,".L"))
  ligands.S<-sapply(tolower(geneList),function(x) paste0(x,".S"))
  return((c(ligands.L,ligands.S)))

}

#Construct a heatmap for the average expression of different genes across different conditions/clusters
heatmap <- function(countn, genes = rownames(countn), cells = 1:length(colnames(countn)), condition, norm = "max", conditionList = "none", aspect.ratio = 1, xAngle = 45, yAngle = 0, xSize = 10, ySize = 10, xAdj = 1, plot = TRUE){
  
  
  subset <- log10(countn[na.omit(match(genes,rownames(countn))),cells]+1)
  avgd <- (t(aggregate(t(as.matrix(subset)),list(condition[cells]),mean)))
  clusters <- avgd[1,]
  avgd <- avgd[-1,]
  avgd <- matrix(as.numeric(unlist(avgd)),nrow=nrow(avgd))
  rownames(avgd) <- rownames(subset)
  colnames(avgd) <- clusters
  
  
  if (norm == "max"){
    avgd_n <- (avgd / apply(avgd,1,function(x) max(x)))
  } else if (norm == "sum"){
    avgd_n <- avgd / rowSums(avgd)
  } else {
    avgd_n <- avgd
  }
  
  if(conditionList != "none"){
    avgd_n <- avgd_n[,match(conditionList,colnames(avgd_n))]
  }

  if(length(which(is.na(rowSums(avgd_n)))) > 0){
    avgd_n <- avgd_n[-which(is.na(rowSums(avgd_n))),]
  }
  
  melted_ <- melt(t(avgd_n))
  
  p <- ggplot(data = melted_, aes(x=Var1, y=Var2, fill=value)) + 
    theme_bw() +
    # geom_tile() + 
    geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low = "white",high = "steelblue")+
    theme(axis.text.x = element_text(angle = xAngle, hjust = xAdj, size = xSize), axis.text.y = element_text(angle = yAngle, size = ySize), aspect.ratio = aspect.ratio, axis.title.x = element_blank(), axis.title.y = element_blank())
  
  if(plot == TRUE){
    print(p)
  }
  return(avgd_n)
}


#Same as heatmap above, but now selecting only one of each .L and .S pair, whichever has the highest expression
heatmap_xenify <- function(countn, genes = rownames(countn), cells = 1:length(colnames(countn)), condition, norm = "max", clusters, aspect.ratio = 1, xAngle = 45, yAngle = 0, xSize = 10, ySize = 10, xAdj = 1, plot = TRUE){
  
  
  conditionList <- clusters$names[which(clusters$names %in% unique(meta$cluster[cells]))]
  
  markers <- unique((genes))
  markers.L<-sapply(tolower((markers)),function(x) paste0(x,".L"))
  markers.S<-sapply(tolower(markers),function(x) paste0(x,".S"))
  
  
  M.L <- heatmap(countn,genes = markers.L,cells = cells,condition = condition, conditionList = conditionList, norm = "none", plot = FALSE)
  M.S <- heatmap(countn,genes = markers.S,cells = cells,condition = condition, conditionList = conditionList, norm = "none", plot = FALSE)
  
  
  max.L <- apply(M.L,1,max)
  max.S <- apply(M.S,1,max)
  
  newMarkers <- c()
  for (i in 1:length(markers)){
    
    mL <- match(markers.L[i],names(max.L))
    mS <- match(markers.S[i],names(max.S))
    
    if(!is.na(mL) & is.na(mS)){ 
      newMarkers <- c(newMarkers,markers.L[i])
    } else if(is.na(mL) & !is.na(mS)) {
      newMarkers <- c(newMarkers,markers.S[i])
    } else if(!is.na(mL) & !is.na(mS)){
      if(max.S[mS] > max.L[mL]){
        newMarkers <- c(newMarkers,markers.S[i])
      }
      else{
        newMarkers <- c(newMarkers,markers.L[i])
      }
    }
  }
  
  
  M <- heatmap(countn,genes = newMarkers,cells = cells,condition = condition, conditionList = conditionList, norm = "max", aspect.ratio = aspect.ratio, xSize = xSize, ySize = ySize, xAngle = xAngle, yAngle = yAngle)
  return(M)
}

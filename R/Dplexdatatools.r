#' id_clusters
#' provides several methods to identify UMAP clusters
#' 
#' @param bar.UMAP
#' @param method  "gmm", "CLARA_medoids", or "kmeans"
#' @param nClusters
#'
#' @return
#' @export
#'
#' @examples
id_clusters <- function(bar.UMAP, method = "CLARA_medoids" , nClusters = 10){
  require(ClusterR)
  require(ggplot2)
  #require(Rcpp) #maybe required if not included with ClusterR

  #separate counts from bar.UMAP and normalize them
  bars <- bar.UMAP[,3:ncol(bar.UMAP)]
  n.bars <- CLR(bars)
  UMAP <- bar.UMAP[,1:2]

  ##TODO! needs error handling if invalid string is entered
  clusters <- switch(method,
                     "gmm" ={ gmm <- GMM(UMAP,
                                         nClusters,
                                         dist_mode = "maha_dist",
                                         seed_mode = "random_subset",
                                         km_iter = 10,
                                         em_iter = 10)
                     pr <- predict_GMM(UMAP, gmm$centroids, gmm$covariance_matrices, gmm$weights)
                     pr$cluster_labels
                     },
                     "CLARA_medoids" ={ clm <- Clara_Medoids(UMAP,
                                                             clusters = nClusters,
                                                             distance_metric = "mahalanobis",
                                                             samples = 10,
                                                             sample_size = 0.05,
                                                             swap_phase = T)
                     clm$clusters
                     },
                     "kmeans" = { km <- KMeans_rcpp(UMAP,
                                                    clusters = nClusters,
                                                    num_init = 5,
                                                    max_iters = 100,
                                                    initializer = "kmeans++")
                     km$clusters
                     })


  #plot UMAP showing centroids
  cluster.UMAP <- as.data.frame(UMAP)
  cluster.UMAP[,"clusters"] <- factor(clusters, levels=1:nClusters)
  g <- ggplot(cluster.UMAP, aes(x=UMAP1, y=UMAP2, col=clusters), environment = environment())
  print(g + geom_point(size=1) + theme_classic() + scale_color_manual(values=rainbow(nClusters)))


  #populate with the UMI's of each centroid
  centroid.list <- lapply(1:nClusters, FUN = function(x){
    n.bars[clusters==x,]
  })

  #populate with multiplicity of each centroid
  Multiplicity = NULL

  for(i in 1:nClusters){
    L = length(centroid.list[[i]])
    Z = sum(centroid.list[[i]] < 0)
    Multiplicity[i] = (L-Z)/L
  }

  #Compare relative abundance of barcodes among clusters
  p.bar.table <- bar.UMAP[,3:ncol(bar.UMAP)]
  tot_antibody_abun <- colSums(p.bar.table) # total barcode among sequenced cells
  p.bar.table <- t(t(p.bar.table)/tot_antibody_abun) #proportion of total barcode in each cell

  #proportion of barcodes per cluster
  cluster_proportion <- Matrix::Matrix(nrow = nClusters, ncol = ncol(p.bar.table))
  for(i in 1:nClusters){
    for(j in 1:ncol(p.bar.table)){
      cluster_proportion[i,j] <- sum(p.bar.table[clusters==i,][,j])
    }
  }

  ClusterSNR <- signalNoise(cluster_proportion)
  #poor cluster SNR can mean two things. The cells took up many barcodes from solution indiscriminately i.e. dead cells and "soupy" cells.
  #Or it can indicate clusters which aggregate because of similar levels of background.
  #cells with Low levels of even background are callable and can cluster together.

  cluster <- c(1:nClusters)
  id.res <- cbind(cluster, Multiplicity, ClusterSNR, deparse.level = 1)

  #recombine UMAP data with Barcode data
  cluster.UMAP <- cbind(cluster.UMAP, bars)
  #compile results as list
  cluster.res <- list(cluster.UMAP, as.data.frame(id.res))
  return(cluster.res)
}

#' cluster.Remove
#' accepts the output of id_clusters() and returns a bar.UMAP removing all points 
#' in the specified cluster
#' @param cluster.res output of id_clusters()
#' @param bad.cluster parameter defining the cluster to be removed
#'
#' @return
#' @export
#'
#' @examples
cluster.Remove <- function(cluster.res, bad.cluster){
  cluster.UMAP <- cluster.res[[1]]
  bar.UMAP <- cluster.UMAP[(cluster.UMAP[,3] != bad.cluster),]
  bar.UMAP <- bar.UMAP[,-3]
  return(bar.UMAP)
}

#' estimateSNR
#' TODO method needs to be improved! This is a prototype function
#'
#' @param bar.table
#'
#' @return
#' @export
#'
#' @examples
estimateSNR <- function(bar.table){

  #attach signal noise ratio as column to the bar table
  snr.barTable <- cbind(bar.table, snr = signalNoise(bar.table))

  snr.min <- min(snr.barTable[,"snr"])
  snr.max <- max(snr.barTable[,"snr"])

  freq.report  <- as.data.frame(snr.barTable)
  freq.report <- density(freq.report[,"snr"])
  freq.report <- data.frame(x=freq.report$x, y = freq.report$y * freq.report$n)

  mode <- as.numeric(freq.report[freq.report[,"y"]==max(freq.report),][1])

  bar.table <- snr.barTable[snr.barTable[,"snr"] > mode,]
  bar.table <- bar.table[,1:ncol(bar.table)-1]
  message("SNR threshold: ", mode)
  return(bar.table)
}

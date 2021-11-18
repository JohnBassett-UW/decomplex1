#' Format sample barcode count matrix
#'
#'A convenience function equipped to handle the outputs of CITE-seq count or the MULTIseq Alignment Suite
#'
#'data.format() performs the following pipeline manipulations:
#'casts data input to a matrix using as.matrix()
#'transposes the input data if necessary
#'removes any uninformative columns attached during the alignment
#' @param UMI_count_matrix raw count output
#'
#' @return formatted count matrix where each column corresponds to a barcode and each row to a cell UMI
#' @export
#'
#' @examples
#' f.bar.table <- data.format(bar.table)
data.format <- function(UMI_count_matrix){
  #Check if input is type matrix array
  if(!inherits(UMI_count_matrix, "matrix")){
    message("Attempting to cast object type ", class(UMI_count_matrix), " to matrix")
    #Attempt to cast to sparse matrix
    UMI_count_matrix <- tryCatch(as.matrix(UMI_count_matrix),

                                 error = function(e){
                                   message("Object ", class(UMI_count_matrix)," could not be cast to matrix: \n", e)
                                   stop()
                                 },

                                 warning = function(w){
                                   message("Warnings detected during casting to matrix\n", w)
                                 },

                                 finally = {
                                   message("Success.")
                                 })
  }
  else{message("input class ", class(UMI_count_matrix))}
  #
  #assumes there are fewer hashes than UMIs and transposes matrix so that UMIs are expressed as rows
  if( ncol(UMI_count_matrix) > nrow(UMI_count_matrix) ){
    UMI_count_matrix <- t(UMI_count_matrix)
    message("input matrix transposed")
  }
  #
  #removes uninformative columns
  columns.remove <- c("nUMI", "nUMI_total", "unmapped") #list of known uninformative column names
  columns.retain <- setdiff(colnames(UMI_count_matrix), columns.remove)
  columns.remove <- intersect(colnames(UMI_count_matrix), columns.remove) #columns.remove now reflects columns which will be removed
  UMI_count_matrix <- UMI_count_matrix[,columns.retain]

  if(length(columns.remove) == 0 ){
    message("No bad columns detected")}
  else if(length(columns.remove) == 1){
    message("one uninformative column removed: ", paste(columns.remove, collapse = " "))}
  else{
    message(length(columns.remove), " uninformative columns removed: ", paste(columns.remove, collapse = " "))
  }

  message("rows: ", nrow(UMI_count_matrix), ", columns: ", ncol(UMI_count_matrix))
  return(UMI_count_matrix)
}

#'Perform center logratio transform on the input data
#'
#'The CLR function provided is nearly identical to the one used in the Seurat package
#'except log1p() has been replaced with log(). The rational for this change is
#'that valid barcode counts are always positive and substantially greater than 0.
#'
#' @param bar.table formatted sample barcode table
#'
#' @return center logratio transformed bar.table as class Matrixarray
#' @export
#'
#' @examples
#' n.bar.table <- CLR(bar.table)
CLR <- function(bar.table) {
  bar.table <- apply(bar.table, MARGIN = 2, FUN = function(x){
    log(x=x/exp(x=sum(log(x=x[x>0]), na.rm = TRUE)/length(x=x)))
  })
  return(bar.table)
}

#' Calculate the signal to noise ratio for each cell in a formatted sample barcode table
#'
#' @param bar.table formatted sample barcode table
#'
#' @return a vector of signal to noise ratios
#' @export
#'
#' @examples
#' snr <- signalNoise(bar.table)
signalNoise <- function(bar.table){
  snr <- NULL
  n = ncol(bar.table)
  snr <- apply(bar.table, MARGIN = 1, FUN = function(x){
    rowMax = max(x)
    noiseMax = sort(x, partial = n-1)[n-1]
    rowMax/noiseMax
  })
  return(snr)
}

#' Performs dimension reduction on a formatted sample barcode table
#'
#' @param bar.table formatted sample barcode table
#' @param normalize logical, should the return table contain CLR normalized values. defaults to False.
#' @param n_neighbors n neighbors parameter for the umap function from package uwot
#' @param min_dist minimum distance parameter for the umap function from package uwot
#'
#' @return sample barcode table with UMAP dimensions UMAP1 & UMAP2 inserted as first two columns.
#' @export
#'
#' @examples
#' bar.UMAP <- barUMAP(bar.table, normalize =F, n_neighbors =50, min_dist = 0.1)
#' bar.UMAP <- barUMAP(bar.table)
barUMAP <- function(bar.table, normalize = F, n_neighbors = 50, min_dist = 0.1){

  # Center logratio transform
  n.bar.table <- CLR(bar.table)
  #remove infinite values
  n.bar.table[is.infinite(n.bar.table[])==TRUE] = 0
  #Generate UMAP via uwot package using default parameters
  message("Performing dimension reduction...")
  UMAP.res <- uwot::umap(n.bar.table,
                    n_neighbors = n_neighbors, #Size of Local neighborhood to constrain manifold learning
                    n_components = 2, #dimensions to embed to
                    n_epochs = NULL, #default to 500 for datasets containing >= 10k vertices. 200 otherwise.
                    min_dist = min_dist, #minimum distance apart points are allowed to be
                    n_trees = 50) #number of trees to build when constructing nearest neighbors. sensible values are 10-100. larger -> better

  message("Done.")
  colnames(UMAP.res) <- c("UMAP1", "UMAP2")

  if(normalize == T){return(cbind(UMAP.res, n.bar.table))}

  return(cbind(UMAP.res, bar.table))
}

#' Covenience function for visualizing output of barUMAP()
#'
#' @param bar.UMAP output matrix from barUMAP() function
#' @param plt logical, should plots be printed to std out. defaults to True.
#' Otherwise ggplot gobs are returned as a list.
#'
#' @return If plt is FALSE returns list of ggplot grobs otherwise no return value
#' @export
#'
#' @examples
#' plotHashes(bar.UMAP)
#' plot_list <- plotHashes(bar.UMAP, plt = FALSE)
plotHashes <- function(bar.UMAP, plt = TRUE){
  require(gridExtra)
  require(ggplot2)

  #ggplot requires data.frame
  bar.UMAP <- tryCatch(as.data.frame(bar.UMAP),
                       error = function(e){
                         stop("ggplot requires type data.frame")
                       })
  #separate barcode counts from dimension reduction
  bars <- bar.UMAP[,3:ncol(bar.UMAP)]
  #plotHashes function requires normalized UMAP
  bars <- tryCatch(as.data.frame(CLR(bars)),
                   error = function(e){
                     stop("Encountered error evaluating Normalization status")
                   },
                   warning = function(w){
                     message("Normalized data detected")
                     bars
                   }#,
                   # finally = function(f){
                   #   message("Good status")
                   # }
  )

  #Remove all values below the column geometric means
  bars[bars<0] <- NaN

  #bars <- as.data.frame(bars)
  message(ncol(bars)," barcodes detected:")

  ##populate plist with barcode columns from bar.UMAP
  plist <- as.list(bars)
  #iterate over barcode columns to create list of UMAPs for each barcode with ggplot
  counter <- 0
  plist <- lapply(plist, function(x){
    counter <<- counter + 1
    ggplot(bar.UMAP, aes(x=UMAP1, y=UMAP2, col=x)) +
      geom_point(size=0.5) +
      theme_classic() +
      scale_color_gradient(low="blue", high="red") +
      theme(legend.position = "none") +
      ggtitle(names(bars[counter]))
  })
  if(plt == FALSE){
    message("Done. Output as ggplot grobs.")
    return(plist)
  }
  message("generating, plots...")
  rowsArranged = max(1, ncol(bars)%/%8)
  grid.arrange(grobs = plist , nrow=rowsArranged , ncol=min(c(ncol(bars),8))) #arrange UMAPs side by side
  message("Plotting...")
}
#' Format sample barcode count matrix
#'
#'A convenience function equipped to handle the outputs of CITE-seq count or the MULTIseq Alignment Suite
#'
#'data.format() performs the following pipeline manipulations:
#'casts data input to a matrix using as.matrix()
#'transposes the input data if necessary
#'removes any uninformative columns attached during the alignment
#' @param UMI_count_matrix raw count output
#'
#' @return formatted count matrix where each column corresponds to a barcode and each row to a cell UMI
#' @export
#'
#' @examples
#' f.bar.table <- data.format(bar.table)
data.format <- function(UMI_count_matrix){
  #Check if input is type matrix array
  if(!inherits(UMI_count_matrix, "matrix")){
    message("Attempting to cast object type ", class(UMI_count_matrix), " to matrix")
    #Attempt to cast to sparse matrix
    UMI_count_matrix <- tryCatch(as.matrix(UMI_count_matrix),

                                 error = function(e){
                                   message("Object ", class(UMI_count_matrix)," could not be cast to matrix: \n", e)
                                   stop()
                                 },

                                 warning = function(w){
                                   message("Warnings detected during casting to matrix\n", w)
                                 },

                                 finally = {
                                   message("Success.")
                                 })
  }
  else{message("input class ", class(UMI_count_matrix))}
  #
  #assumes there are fewer hashes than UMIs and transposes matrix so that UMIs are expressed as rows
  if( ncol(UMI_count_matrix) > nrow(UMI_count_matrix) ){
    UMI_count_matrix <- t(UMI_count_matrix)
    message("input matrix transposed")
  }
  #
  #removes uninformative columns
  columns.remove <- c("nUMI", "nUMI_total", "unmapped") #list of known uninformative column names
  columns.retain <- setdiff(colnames(UMI_count_matrix), columns.remove)
  columns.remove <- intersect(colnames(UMI_count_matrix), columns.remove) #columns.remove now reflects columns which will be removed
  UMI_count_matrix <- UMI_count_matrix[,columns.retain]

  if(length(columns.remove) == 0 ){
    message("No bad columns detected")}
  else if(length(columns.remove) == 1){
    message("one uninformative column removed: ", paste(columns.remove, collapse = " "))}
  else{
    message(length(columns.remove), " uninformative columns removed: ", paste(columns.remove, collapse = " "))
  }

  message("rows: ", nrow(UMI_count_matrix), ", columns: ", ncol(UMI_count_matrix))
  return(UMI_count_matrix)
}

#'Perform center logratio transform on the input data
#'
#'The CLR function provided is nearly identical to the one used in the Seurat package
#'except log1p() has been replaced with log(). The rational for this change is
#'that valid barcode counts are always positive and substantially greater than 0.
#'
#' @param bar.table formatted sample barcode table
#'
#' @return center logratio transformed bar.table as class Matrixarray
#' @export
#'
#' @examples
#' n.bar.table <- CLR(bar.table)
CLR <- function(bar.table) {
  bar.table <- apply(bar.table, MARGIN = 2, FUN = function(x){
    log(x=x/exp(x=sum(log(x=x[x>0]), na.rm = TRUE)/length(x=x)))
  })
  return(bar.table)
}

#' Calculate the signal to noise ratio for each cell in a formatted sample barcode table
#'
#' @param bar.table formatted sample barcode table
#'
#' @return a vector of signal to noise ratios
#' @export
#'
#' @examples
#' snr <- signalNoise(bar.table)
signalNoise <- function(bar.table){
  snr <- NULL
  n = ncol(bar.table)
  snr <- apply(bar.table, MARGIN = 1, FUN = function(x){
    rowMax = max(x)
    noiseMax = sort(x, partial = n-1)[n-1]
    rowMax/noiseMax
  })
  return(snr)
}

#' Performs dimension reduction on a formatted sample barcode table
#'
#' @param bar.table formatted sample barcode table
#' @param normalize logical, should the return table contain CLR normalized values. defaults to False.
#' @param n_neighbors n neighbors parameter for the umap function from package uwot
#' @param min_dist minimum distance parameter for the umap function from package uwot
#'
#' @return sample barcode table with UMAP dimensions UMAP1 & UMAP2 inserted as first two columns.
#' @export
#'
#' @examples
#' bar.UMAP <- barUMAP(bar.table, normalize =F, n_neighbors =50, min_dist = 0.1)
#' bar.UMAP <- barUMAP(bar.table)
barUMAP <- function(bar.table, normalize = F, n_neighbors = 50, min_dist = 0.1){

  # Center logratio transform
  n.bar.table <- CLR(bar.table)
  #remove infinite values
  n.bar.table[is.infinite(n.bar.table[])==TRUE] = 0
  #Generate UMAP via uwot package using default parameters
  message("Performing dimension reduction...")
  UMAP.res <- uwot::umap(n.bar.table,
                    n_neighbors = n_neighbors, #Size of Local neighborhood to constrain manifold learning
                    n_components = 2, #dimensions to embed to
                    n_epochs = NULL, #default to 500 for datasets containing >= 10k vertices. 200 otherwise.
                    min_dist = min_dist, #minimum distance apart points are allowed to be
                    n_trees = 50) #number of trees to build when constructing nearest neighbors. sensible values are 10-100. larger -> better

  message("Done.")
  colnames(UMAP.res) <- c("UMAP1", "UMAP2")

  if(normalize == T){return(cbind(UMAP.res, n.bar.table))}

  return(cbind(UMAP.res, bar.table))
}

#' Covenience function for visualizing output of barUMAP()
#'
#' @param bar.UMAP output matrix from barUMAP() function
#' @param plt logical, should plots be printed to std out. defaults to True.
#' Otherwise ggplot gobs are returned as a list.
#'
#' @return If plt is FALSE returns list of ggplot grobs otherwise no return value
#' @export
#'
#' @examples
#' plotHashes(bar.UMAP)
#' plot_list <- plotHashes(bar.UMAP, plt = FALSE)
plotHashes <- function(bar.UMAP, plt = TRUE){
  require(gridExtra)
  require(ggplot2)

  #ggplot requires data.frame
  bar.UMAP <- tryCatch(as.data.frame(bar.UMAP),
                       error = function(e){
                         stop("ggplot requires type data.frame")
                       })
  #separate barcode counts from dimension reduction
  bars <- bar.UMAP[,3:ncol(bar.UMAP)]
  #plotHashes function requires normalized UMAP
  bars <- tryCatch(as.data.frame(CLR(bars)),
                   error = function(e){
                     stop("Encountered error evaluating Normalization status")
                   },
                   warning = function(w){
                     message("Normalized data detected")
                     bars
                   }#,
                   # finally = function(f){
                   #   message("Good status")
                   # }
  )

  #Remove all values below the column geometric means
  bars[bars<0] <- NaN

  #bars <- as.data.frame(bars)
  message(ncol(bars)," barcodes detected:")

  ##populate plist with barcode columns from bar.UMAP
  plist <- as.list(bars)
  #iterate over barcode columns to create list of UMAPs for each barcode with ggplot
  counter <- 0
  plist <- lapply(plist, function(x){
    counter <<- counter + 1
    ggplot(bar.UMAP, aes(x=UMAP1, y=UMAP2, col=x)) +
      geom_point(size=0.5) +
      theme_classic() +
      scale_color_gradient(low="blue", high="red") +
      theme(legend.position = "none") +
      ggtitle(names(bars[counter]))
  })
  if(plt == FALSE){
    message("Done. Output as ggplot grobs.")
    return(plist)
  }
  message("generating, plots...")
  rowsArranged = max(1, ncol(bars)%/%8)
  grid.arrange(grobs = plist , nrow=rowsArranged , ncol=min(c(ncol(bars),8))) #arrange UMAPs side by side
  message("Plotting...")
}

#' quant_first_thresh
#' uses polynomial aproximation to get calculate derivative
#' TODO implement diff()
#' @param model.x
#'
#' @return 1st derivative zeroes
#' @export
#'
#' @examples
quant_first_thresh <- function(model.x){
  #polynomial regression model of first derivative function
  drv <- KernSmooth::locpoly(x = 1:100, y = model.x, drv = 1L, bandwidth = 1)
  drv <- data.frame(x = drv$x, y = drv$y)

  #classify maxima as right first derivitive zeroes
  extrema <- NULL
  for(k in 1:(nrow(drv)-1)){
    if( drv[,"y"][k] >= 0 && drv[,"y"][k+1] <= 0 ){
      extrema <- c(extrema,drv[,"x"][k])
    }
  }
  return(extrema)
}


#' quant_second_thresh
#'Uses polynomial approximation to calculate 2nd derivative
#'TODO diff(diff())
#' @param model.x
#'
#' @return 2nd derivative zeroes
#' @export
#'
#' @examples
quant_second_thresh <- function(model.x){
    drv_2 <- KernSmooth::locpoly(x = 1:100, y = model.x, drv = 2L, bandwidth = 1)
    drv_2 <- data.frame(x = drv_2$x, y = drv_2$y)

    #left side approach
    extreme.left <- NULL
    for(k in 1:(nrow(drv_2)-1)){
      if( drv_2[,"y"][k] <= 0 && drv_2[,"y"][k+1] >= 0 ){
        extreme.left <- c(extreme.left,drv_2[,"x"][k])
      }
    }

    #right side approach
    extreme.right <- NULL
    for(k in 1:(nrow(drv_2)-1)){
      if( drv_2[,"y"][k] >= 0 && drv_2[,"y"][k+1] <= 0 ){
        extreme.right <- c(extreme.right,drv_2[,"x"][k])
      }
    }
    ####FOR quality check#####
    #plot2 <- ggplot(drv, aes(x=x, y=y)) + geom_point() +
    #                                       geom_vline(xintercept = extreme.left, color = "green") +
    #                                       geom_vline(xintercept = extreme.right, color = "orange") +
    #                                       ggtitle(paste("1st derivitive, ",colnames(barTable.n[i])))
    #print(plot2)
    # plot3 <- ggplot(QA_out, aes(y =x, x = 1:100)) + geom_point() +
    #                                                 geom_vline(xintercept = extreme.left, color = "green") +
    #                                                 geom_vline(xintercept = extreme.right, color = "orange") +
    #                                                 ggtitle(paste("inflections, ",colnames(barTable.n[i])))
    # print(plot3)

  extrema <- c(extreme.left, extreme.right)
  return(extrema)
}


#' force_classifycells()
#'
#'Assumptions made by deMULTIplex break down when the number of barcodes is
#'small and if there is unequal representation of barcodes.
#'
#'
#'force_classifycells() provides an additional control branch to handle barcodes
#'which do not behave as assumed by deMULTIplex
#'
#'forcing cell classification will increase the number of doublets falsely called as singlets
#' @param barTable
#' @param q
#'
#' @return bc_calls, a list of singlets, doublets, and negatives called for quantile q
#' @export
#'
#' @examples
force_classifycells <- function(barTable, q, QA = 1, est.doublet) {

  ## Normalize Data: Log2 Transform, mean-center
  barTable.n <- as.data.frame(log2(barTable))
  for (i in 1:ncol(barTable)) {
    ind <- which(is.finite(barTable.n[,i]) == FALSE)
    barTable.n[ind,i] <- 0
    barTable.n[,i] <- barTable.n[,i]-mean(barTable.n[,i])
  }

  ## Pre-allocate memory to save processing time/memory
  n_BC <- ncol(barTable.n) # Number of barcodes
  n_cells <- nrow(barTable.n) # Numer of cells
  bc_calls <- rep("Negative",n_cells) # Barcode classification for each cell
  names(bc_calls) <- rownames(barTable.n)

  ## Generate thresholds for each barcode:
  ## 1. Gaussian KDE with bad barcode detection, outlier trimming
  ## 2. Define local maxima for GKDE
  ## 3. Split maxima into low and high subsets, adjust low if necessary
  ## 4. Threshold and classify cells according to user-specified inter-maxima quantile

  for (i in 1:n_BC) {
    ## Step 1: GKDE
    model <- tryCatch( { stats::approxfun(KernSmooth::bkde(barTable.n[,i], kernel="normal")) },
                       error=function(e) { print(paste0("No threshold found for ", colnames(barTable.n)[i],"...")) } )
    if (class(model) == "character") { next }
    x <-  seq(from=quantile(barTable.n[,i],0.001), to=quantile(barTable.n[,i],0.999), length.out=100)

    extrema <- quant_first_thresh(model(x))

    #set the low extreme from the first derivative
    low.extreme <- extrema[which.max(model(x)[extrema])]

    ####FOR QA####
      QA_out <- data.frame(x = model(x))
      # plot1 <- ggplot(QA_out, aes(y =x, x = 1:100)) + geom_point() +
      #   geom_vline(xintercept = extrema ) +
      #   geom_vline(xintercept = low.extreme, color = "blue") +
      #   ggtitle(colnames(barTable.n[i]))
      # print(plot1)

    ## Step 2.2 find suitable thresh for barcodes that fail Quality testing

    if (length(extrema) < 1) {
      print(paste0("No threshold found for ", colnames(barTable.n)[i],"..."))
      next
    }

    ## Step 3: Select higher maximum
    ## Assumes negative cells are largest mode
    ## Assumes positive cells diffuse throughout --will mistake some doublets as singlets

    ## compute higher threshold maxima from possible inflection points
    extrema <- quant_second_thresh(model(x))

    #deprecated control flow... This only works if the graphs have assumed deMULTIplex behavior
    #if( signif(max(extrema), digits = 2) == match(max(model(x)), model(x)) ){}

    extrema <- extrema[extrema > low.extreme]
    AUC <- numeric(length(extrema))
    AUC[1] <- integrate(f = model, lower = model(x)[extrema[1]], upper = model(x)[100] )[[1]]
      for(z in 1:length(extrema)){
        AUC[z+1] <- integrate(f = model, lower = model(x)[extrema[z]], upper = model(x)[100] )[[1]]
      }

    if(length(unique(AUC[-1])) == 1 ){
      high.extreme <- sort(extrema)[2]
      }else{
      high.extreme <- extrema[which(abs(AUC[-1]-est.doublet*AUC[1])==min(abs(AUC[-1]-est.doublet*AUC[1])))]
      }
    # plot4 <- ggplot(QA_out, aes(y = x, x = 1:100)) +
    #     geom_point() +
    #     geom_vline(xintercept = extrema ) +
    #     geom_vline(xintercept = high.extreme, color = "red") +
    #     ggtitle(paste("thresholds"), colnames(barTable.n[i]))
    # print(plot4)



    if (low.extreme == high.extreme) {print(paste0("No threshold found for ", colnames(barTable.n)[i],"...")) }

    ### FOR QA ###
    if(i == QA && q == .01){
    plot5 <- ggplot(QA_out, aes(y =x, x = 1:100)) + geom_point() + geom_vline(xintercept = c(low.extreme,high.extreme), color = c("blue", "red")) +
                                                    ggtitle(paste(colnames(barTable.n[i]), " FINAL"))
    print(plot5)
    }

    ## Step 4: Threshold and classify cells
    thresh <- quantile(c(x[high.extreme], x[low.extreme]), q)
    cell_i <- which(barTable.n[,i] >= thresh)
    n <- length(cell_i)
    if (n == 0) { next } # Skip to next barcode if no cells classified
    bc_calls[cell_i] <- sapply(bc_calls[cell_i],
                               FUN = function(x) {
                                 if (x == "Negative") {
                                   return(colnames(barTable.n)[i])
                                 } else {
                                   return("Doublet")
                                 } })
  }

  return(bc_calls)

}

#' findThresh()
#' modified acquisition of extrema
#' @param call.list
#'
#' @return
#' @export
#'
#' @examples
#'
findThresh <- function(call.list) {

  res <- as.data.frame(matrix(0L, nrow=length(call.list), ncol=4))
  colnames(res) <- c("q","pDoublet","pNegative","pSinglet")
  q.range <- unlist(strsplit(names(call.list), split="q="))
  res$q <- as.numeric(q.range[grep("0", q.range)])
  nCell <- length(call.list[[1]])

  for (i in 1:nrow(res)) {
    temp <- table(call.list[[i]])
    if ( "Doublet" %in% names(temp) == TRUE ) { res$pDoublet[i] <- temp[which(names(temp) == "Doublet")] }
    if ( "Negative" %in% names(temp) == TRUE ) { res$pNegative[i] <- temp[which(names(temp) == "Negative")] }
    res$pSinglet[i] <- sum(temp[which(names(temp) != "Doublet" & names(temp) != "Negative")])
  }

  res <- reshape2::melt(res, id.vars="q")
  res[,4] <- res$value/nCell
  colnames(res)[2:4] <- c("Subset","nCells","Proportion")

  extrema <- res[res[,"Subset"] == "pSinglet", "Proportion"]
  extrema <- max(extrema)
  extrema <- res[res[,"Proportion"] == extrema, "q"]
  return(list("extrema" = extrema, "res" = res))
}

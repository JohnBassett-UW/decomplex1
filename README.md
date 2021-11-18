# decomplex1

<h3> by John Bassett  

A fork of deMULTIplex by Chris McGinnis  
</h3>



The decomplex1 software is currently a protoype. It was developed for a specific use-case of the deMULTIplex algorithm, i.e. small numbers (<10) of antibody barcoded samples originating from frozen tumor slice culture. The necessity of this version arises from the dramatic distribution changes that occur in hashing libraries derived from antibody-barcoded frozen tumor samples as compared with PBMCs barcoded with lipid-modified oligos. 


For more information on deMULTIplex see: 

https://github.com/chris-mcginnis-ucsf/MULTI-seq1 

## Setup
### installation
    devtools::install_github('JohnBassett-UW/decomplex')

### load sample data
    load("path-to-file/bar.table1.robj")

## Formatting the data
decomplex1 will format barcode count matrices that are generated using either _CITE-seq count_ or _MULTIseq alignment suite_. 

    f.bar.table <- data.format(bar.table1)
This sample dataset was generated with _MULTI_seq alignment suite_. data.format() simply casts the data to class matrix and removes all the uninformative columns. 
``` 
Attempting to cast object type data.frame to matrix   
Success.     
2 uninformative columns removed: nUMI nUMI_total     
rows: 12145, columns: 6 
```
For data generated using other software, decomplex1 expects UMI's as rows and barcodes as columns.  

It is good practice to view the data at this stage anyways and make sure that it is formatted correctly.

    View(f.bar.table)

![sample matrix](/vignettes/Capture_matrix.PNG)

## Generate barcode space
Generate 2 dimensional embedding of the n dimensional barcode space.

    bar.UMAP <- barUMAP(f.bar.table)
    
Visualize the barcode space.

    plotHashes(bar.UMAP)

This provides a heatmap to identify barcode clusters and visualize the cross contamination profile of barcodes.

Note that color corresponds with barcode counts that are above the barcode geometric mean. Counts below the geometric mean appear gray.

![sample barcode space](/vignettes/plotHashes.png)

## Classify Nuclei
This next part is somewhat opaque, but performs barcode classifcation for each UMI. At this stage an estimation of the proportion of doublets is used to threshold barcode distributions. The sample classification workflow here is identical to deMULTIplex. The key difference is the method by which the distribution thresholds are selected. 

```
bar.table.full <- f.bar.table #preserve full data
bar.table_sweep.list <- list()
counter <- 0
est.doublet = 0.10
while (counter >= 0) {
  counter <- counter + 1
  print(paste("Round ", counter, "...",sep=""))
  bar.table_sweep.list <- list()
  n <- 0
  for (q in seq(0.01, 0.99, by=0.02)) {
    print(q)
    n <- n + 1
    bar.table_sweep.list[[n]] <- force_classifycells(f.bar.table, q=q, QA = 1, est.doublet = est.doublet)
    names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
  }
  threshold.results <- findThresh(call.list=bar.table_sweep.list)
  ####QC###
  print(ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) +
          geom_line() +
          ggtitle(paste("round", counter)) +
          theme(legend.position = "none") +
          geom_vline(xintercept=threshold.results$extrema, lty=2) +
          scale_color_manual(values=c("red","black","blue")))
  ########
  calls <- force_classifycells(f.bar.table, q=threshold.results$extrema, est.doublet = est.doublet )
  negs <- names(calls)[which(calls == "Negative")]

  if (length(negs) == 0) { break }
  if (counter == 1) { neg.cells <- negs }
  if (counter > 1) { neg.cells <- c(neg.cells, negs) }
  print(paste(length(which(calls %in% c("Negative","Doublet"))), " Singlets...", sep=""))
  print(paste(length(which(calls == "Doublet")), " Negatives...", sep=""))
  print(paste(length(which(calls == "Negative")), " Negatives...", sep=""))
  f.bar.table <- f.bar.table[-which(rownames(f.bar.table) %in% neg.cells), ]
} 
```

In this example a single distinct mode is visibile before classification. The majority of the classifiable nuclei are distributed broadly. The minimum threshold is now selected by the mode and the max threshold is selected to adequately classify doublets. 

![before thresholding](/vignettes/Sample_barcode_original.png)



After classification is finished and negatives have been removed, two modes are clearly visible. The larger mode (blue line) contains true negatives for the specified barcode and the smaller mode contains the stable barcode calls.

![Sample thresholding](/vignettes/SampleBarcodeThresholding.png)



This step ends with the final optimization of Singlets to doublets to negatives.

![optimization](/vignettes/Optimization.png)



Already if classifications are visualized, barcodes for most nuclei can be accurately called.



![Classifications](/vignettes/Classification.png)



## Reclassify Nuclei

However, this can still be improved by performing negative reclassification as in deMULTIplex, but now using appropriately selected extrema. 

```
final.calls <- c(calls, rep("Negative",length(neg.cells)))
names(final.calls)<- c(names(calls), neg.cells)

reclass.cells <- findReclassCells(bar.table.full, names(final.calls)[which(final.calls=="Negative")])
reclass.res <- rescueCells(bar.table.full, final.calls, reclass.cells)
```

Visualizing this shows a substantial recovery of negatives.

```
bar.UMAP.final <- as.data.frame(bar.UMAP)
bar.UMAP.final[,"Stable"] <- rep("unknown")
bar.UMAP.final[names(final.calls.rescued), "Stable"] <- final.calls.rescued

ggplot(bar.UMAP.final, aes(x=UMAP1, y=UMAP2, color=Stable)) + geom_point(size=0.5) + theme_classic() + scale_color_manual(values=c("dodgerblue","goldenrod","darkred","seagreen","darkorchid3","pink","black","grey"))
```

![Negative Reclassification](/vignettes/Reclassification.png)


## Save calls and generate QC statistics

Finally statistics can be generated and classifications are saved.

```
#generate hash calling summary table
barcode <- c(paste("Bar", 1:6, sep =""), "Doublet", "Negative")
cells.recovered <- data.frame(identity = barcode, cells = rep("empty"), fraction = rep("empty"))
for(i in 1:length(barcode)){
  cells.recovered[i,2] <- nrow(bar.UMAP.final[which(bar.UMAP.final[,"Stable"] == barcode[i]),])
}
cells.recovered[,2] <- as.numeric(cells.recovered[,2])
cells.recovered[,3] <- cells.recovered[,2]/sum(cells.recovered[,2])
cells.recovered

#generate sample call table
sampleCalls <- cbind(rownames(bar.UMAP.final), bar.UMAP.final[,"Stable"])
colnames(sampleCalls) <- c("Cell UMI", "group")

sampleGroups <- c("Day 0", "Day 2", "Day 5", "Day 7", "Day 10", "Day 14", "Doublet", "Negative")
sampleCalls <- as.data.frame(sampleCalls)
sampleCalls[,"sample"] <- rep("empty")
for(i in 1:nrow(sampleCalls)){
  for(k in 1: length(barcode)){
    if(sampleCalls[i,2] == barcode[k]){
      sampleCalls[i,3] <- sampleGroups[k]
    }
  }
}
head(sampleCalls)
save(sampleCalls, file = "sampleCalls.robj")
```

This example shows reasonable levels of doublets and negatives. A total of 16% of the nuclei are counted as doublets and 17% as negatives.

      identity cells   fraction
    1     Bar1  2075 0.17085220
    2     Bar2  1589 0.13083573
    3     Bar3  1170 0.09633594
    4     Bar4  1424 0.11724990
    5     Bar5   668 0.05500206
    6     Bar6  1112 0.09156031
    7   Doublet  1953 0.16080692
    8 Negative  2154 0.17735694


#Load sample data
load("bar.table1.robj")

###generate barcode space embedding###
f.bar.table <- data.format(bar.table1)
View(f.bar.table) #good practice: double check formatting of raw data
bar.UMAP <- barUMAP(f.bar.table)
plotHashes(bar.UMAP) #good practice: visualize barcode space

###Classify Cells###
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

###Reclassify Negatives###
final.calls <- c(calls, rep("Negative",length(neg.cells)))
names(final.calls)<- c(names(calls), neg.cells)

reclass.cells <- findReclassCells(bar.table.full, names(final.calls)[which(final.calls=="Negative")])
reclass.res <- rescueCells(bar.table.full, final.calls, reclass.cells)

###Visualize Class Stability###
ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) +
  geom_point() + xlim(c(nrow(reclass.res)-1,1)) +
  ylim(c(0,1.05)) +
  geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
  geom_vline(xintercept = 10, color = "blue") #This value will be dataset-specific

###assemble call lists based on class stability###
final.calls.rescued <- final.calls
rescue.ind <- which(reclass.cells$ClassStability >= 4) ## Note: Value will be dataset-specific
final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]

final.calls.rescued_all <- final.calls
rescue.ind <- which(reclass.cells$ClassStability >= 0) ## Note: Value will be dataset-specific
final.calls.rescued_all[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]

bar.UMAP.final <- as.data.frame(bar.UMAP)
bar.UMAP.final[,"Classified"] <- rep("unknown")
bar.UMAP.final[,"Stable"] <- rep("unknown")
bar.UMAP.final[,"All_Reclassed_Cells"] <- rep("unknown")
bar.UMAP.final[names(final.calls), "Classified"] <- final.calls
bar.UMAP.final[names(final.calls.rescued), "Stable"] <- final.calls.rescued
bar.UMAP.final[names(final.calls.rescued_all), "All_Reclassed_Cells"] <- final.calls.rescued_all

###Visualize calls###
ggplot(bar.UMAP.final, aes(x=UMAP1, y=UMAP2, color=Classified)) + geom_point(size=0.5) + theme_classic() + scale_color_manual(values=c("dodgerblue","goldenrod","darkred","seagreen","darkorchid3","pink","black","grey"))
ggplot(bar.UMAP.final, aes(x=UMAP1, y=UMAP2, color=Stable)) + geom_point(size=0.5) + theme_classic() + scale_color_manual(values=c("dodgerblue","goldenrod","darkred","seagreen","darkorchid3","pink","black","grey"))
ggplot(bar.UMAP.final, aes(x=UMAP1, y=UMAP2, color=All_Reclassed_Cells)) + geom_point(size=0.5) + theme_classic() + scale_color_manual(values=c("dodgerblue","goldenrod","darkred","seagreen","darkorchid3","pink","black","grey"))

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

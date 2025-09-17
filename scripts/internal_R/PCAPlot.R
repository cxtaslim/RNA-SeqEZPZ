# modified 2021/02/02 by CT to determine nsub and avoid error in vst
PCAPlot <- function (object, group=NULL,counts.trans,varInt,typeTrans, ntop = min(500, nrow(counts.trans)), 
                     col, batch=NULL,outfile = TRUE,batchRem=FALSE) 
{
  if (typeTrans == "VST") {
	# calculate the number of rows with counts >5 to determine nsub and avoid error
	nsub=sum(rowMeans(assays(object)$counts)>5)
	if(nsub>=1000){
	 deseq.trans <- vst(object,blind=FALSE,nsub=1000)
	}else{
	 deseq.trans <- varianceStabilizingTransformation(object,blind=FALSE)
	}
  }else  {
    deseq.trans <- rlog(object,blind=FALSE)
  } 

  if (outfile & !batchRem){ 

    png(filename = "figures/PCA.png", width = cairoSizeWrapper(1800 * 
         2), height = cairoSizeWrapper(1800), res = 300)
  }  
  if (outfile & batchRem){ 
    png(filename = "figures/PCA_batchRem.png", width = cairoSizeWrapper(1800 * 
         2), height = cairoSizeWrapper(1800), res = 300)
  }
  if(batchRem){
    mat <- assay(deseq.trans)
    mat <- limma::removeBatchEffect(mat,colData(deseq.trans)[,batch])
    assay(deseq.trans) <- mat
  }
  
  pcaData <- DESeq2::plotPCA(deseq.trans, intgroup = varInt ,ntop=ntop,returnData=TRUE)
  pcaData$rep <- colData(deseq.trans)[,"rep"]
  # as factor so rep1 is not followed by rep10
  pcaData$rep <- factor(pcaData$rep,levels=gtools::mixedsort(unique(pcaData$rep)))	
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ngroup=length(unique(group))
  # shape values for plotting
  shape_val <- c(19,15,17,18,setdiff(0:100,c(19,15,17,18)))
  if(!is.null(batch)) {
    pcaData$batch = colData(deseq.trans)[,batch]
    p<-ggplot(pcaData, aes(x = PC1, y = PC2, color = group, shape = batch)) +
      geom_point(size =4,stroke=2) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      # coord_fixed() +
      # Adding plotting space to fit all points
      scale_y_continuous(expand=expansion(mult=c(0.1,0.1)))+
	  scale_color_manual(values=col)+
      # need to specify shapes manually otherwise there is no shape after 6	
      scale_shape_manual(values = shape_val) +
      # make group legend on top
      guides(
       color =	guide_legend(order = 1),
       shape =	guide_legend(order =2,
		# shrink shape on legend so they don't overlap
                override.aes=list(size=3,stroke=2))
      )

  } else{
    p<-ggplot(pcaData, aes(x = PC1, y = PC2, color = group,shape=rep)) +
      geom_point(size =4,stroke=2) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
     # coord_fixed() +
     # Adding plotting space to fit all points
     scale_y_continuous(expand=expansion(mult=c(0.1,0.1)))+
	  scale_color_manual(values=col)+
      # need to specify shapes manually otherwise there is no shape after 6	
      scale_shape_manual(values = shape_val) +
      # make group legend on top
      guides(
	color = guide_legend(order = 1),
        shape = guide_legend(order =2,
		# shrink shape on legend so they don't overlap
		override.aes=list(size=3,stroke=2))
      )     
  }
  print(p)
  rv = apply(counts.trans, 1, var, na.rm = TRUE)
  pca = prcomp(t(counts.trans[order(rv, decreasing = TRUE), 
                              ][1:ntop, ]))
  if (outfile) 
    dev.off()
  return(invisible(pcaData$x))
}


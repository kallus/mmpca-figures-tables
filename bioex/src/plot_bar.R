plot_bar <- function(z,annotAssoc,...){
    
  zOrderD <- order(z,decreasing=T)
  zPlot <- z[zOrderD]
  
  zMethGenes <- sapply(names(zPlot),function(n){
    ind <- which(annotAssoc$Name==n)
    if(length(ind)>0){
      geneNames <- unique(annotAssoc$UCSC_RefGene_Name[ind])
      if(length(geneNames)>0){
        return(paste(geneNames,collapse=";"))  
      }else{
        return("")
      }
    }else{
      return("")
    }
  })
  
  names(zPlot)[zMethGenes!=""] <- zMethGenes[zMethGenes!=""] #paste(names(y),":",yAppend,sep="")
  colArray <- array("#cccccc",length(zPlot))
  colArray[zMethGenes!=""] <- "purple"
  
  h <- barplot(zPlot, horiz=TRUE,names.arg=NA,cex.names=1,col=colArray,xlab="Comp",...)
  
  zPlotNegInd <- zPlot < 0
  zPlotPosInd <- zPlot >= 0
  
  text(0,h[zPlotNegInd],names(zPlot)[zPlotNegInd],cex=0.6,pos=4,xpd=T)
  text(0,h[zPlotPosInd],names(zPlot)[zPlotPosInd],cex=0.6,pos=2,xpd=T)

}
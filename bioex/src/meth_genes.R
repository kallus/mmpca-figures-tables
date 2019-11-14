meth_genes <- function(nameArray,annotAssoc){
  geneNames <- sapply(nameArray,function(n){
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
  
  nameArray[geneNames!=""] <- geneNames[geneNames!=""]
  nameArray
  
}


export_gsea <- function(thisRun,outputPathModule,origMat){
  
  outDat <- list()
  sumFiles <- list()
  
  pThresh <- 1.1 # Include everything
  
  library(XML)
  for(i in 1:ncol(origMat)){
    
    print(i)
    
    wdF <- gsub('/','\\\\',paste(outputPathModule,"\\output_",thisRun,"\\output_",i,sep=""))
    wdFS <- gsub('/','\\\\',paste(thisRun,"\\output_",i,sep=""))
    
    file.names <- dir(wdF)
    
    for(j in 1:length(file.names)){
      
      thisId <- str_split(file.names[j],'\\.')[[1]][3]
      
      sumFiles[i] <- paste(wdFS,'\\',file.names[j],'\\index.html',sep="")
      
      fNeg <- paste(wdF,'\\',file.names[j],'\\gsea_report_for_cluster_',i,'_repos_neg_',thisId,'.html',sep="")
      fPos <- paste(wdF,'\\',file.names[j],'\\gsea_report_for_cluster_',i,'_repos_pos_',thisId,'.html',sep="")
      
      fNeg2 <- paste(wdF,'\\',file.names[j],'\\gsea_report_for_cluster_',i,'_neg_',thisId,'.html',sep="")
      fPos2 <- paste(wdF,'\\',file.names[j],'\\gsea_report_for_cluster_',i,'_pos_',thisId,'.html',sep="")
      
      if(!file.exists(fNeg)){
        fNeg <- fNeg2
      }
      
      if(!file.exists(fPos)){
        fPos <- fPos2
      }
      
      if(file.exists(fNeg)){
        dfNeg <- readHTMLTable(fNeg,stringsAsFactors=F,header=F)[[1]]
        dfNegSub <- dfNeg[as.numeric(dfNeg[,8])<pThresh,c(2,5,8)]
      }else{
        dfNegSub <- array(0,c(0,3))
      }
      
      if(file.exists(fPos)){
        dfPos <- readHTMLTable(fPos,stringsAsFactors=F,header=F)[[1]]
        dfPosSub <- dfPos[as.numeric(dfPos[,8])<pThresh,c(2,5,8)]
      }else{
        dfPosSub <- array(0,c(0,3))
      }
      
      tabOut <- rbind(dfNegSub,dfPosSub)
      
      outDat[[length(outDat)+1]] <- tabOut
      
    }
    
  }
  
  names(outDat) <- colnames(gseaMat)
  
  allPath <- unique(unlist(lapply(outDat,function(x){ x[x[,3]<pThresh,1] })))
  pathMat <- array(NA,c(length(outDat),length(allPath)))
  pathMat2 <- array(NA,c(length(outDat),length(allPath)))
  
  
  colnames(pathMat) <- allPath
  rownames(pathMat) <- names(outDat)
  
  for(i in 1:length(outDat)){
    
    thisDat <- outDat[[i]]
    thisDat <- thisDat[thisDat[,3]<pThresh,,drop=F]
    
    if(nrow(thisDat)>0){
      
      pathMat[i,match(thisDat[,1],allPath)] <- as.numeric(thisDat[,3])
      pathMat2[i,match(thisDat[,1],allPath)] <- as.numeric(thisDat[,2])
      
      
    }
    
  }
  
  distMethod <- "euclidean"
  
  thisDist <- as.dist(dist_na(t(pathMat2),method=distMethod))
  thisDist[is.na(thisDist)] <- mean(thisDist,na.rm=T)
  hc1 <- hclust(thisDist)
  
  thisDist <- as.dist(dist_na(pathMat2,method=distMethod))
  thisDist[is.na(thisDist)] <- mean(thisDist,na.rm=T)
  hc2 <- hclust(thisDist)
  
  pathMatO <- pathMat[order.dendrogram(as.dendrogram(hc1)),order.dendrogram(as.dendrogram(hc2))]
  pathMat2O <- pathMat2[order.dendrogram(as.dendrogram(hc1)),order.dendrogram(as.dendrogram(hc2))]
  sumFilesO <- sumFiles[order.dendrogram(as.dendrogram(hc1))]
  
  #----- 
  
  printer = file(paste(outputPathModule,'/GSEA_',thisRun,'.html',sep=""),"w") 
  
  pathMatO[is.na(pathMatO)] <- ''
  
  pathMatCol <- array("#ffffff",dim(pathMatO))
  pathMatCol[pathMat2O>0] <- "#fddbc7"
  pathMatCol[pathMat2O<0] <- "#d1e5f0"
  
  writeLines(paste('<html><head><link rel="stylesheet" href="style.css"></head><body>',sep=""),con=printer)
  writeLines(paste('<table style="width:100%">',sep=""),con=printer)
  
  writeLines(paste("<tr><th></th><th class='rotate'><div><span>",paste(colnames(pathMatO),collapse="</span></div></th><th class='rotate'><div><span>"),"</span></div></tr>",sep=""),con=printer)
  
  
  makeTd <- function(v1,v2){
    
    paste("<td style='background:",v1,"'>",v2,"</td>",sep="")
    
  }
  
  for(i in 1:ncol(gseaMat)){
    
    writeLines(paste('<tr><td><a target="_blank" href="',sumFilesO[i],'">',rownames(pathMatO)[i],"</a></td>",paste(lapply(1:ncol(pathMatO),function(j){makeTd(pathMatCol[i,j],pathMatO[i,j])}),collapse=""),'</tr>',sep=""),con=printer)
    
  }
  
  writeLines(paste('</table></body></html>',sep=""),con=printer)
  
  close(printer) 
  
  save(pathMat,file=paste(outputPathModule,"/pathMat_",thisRun,".RData",sep=""))
  
}







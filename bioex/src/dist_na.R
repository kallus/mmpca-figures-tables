dist_na <- function(inMat,method,doAbs=F){
  
  distOut <- matrix(NA,nrow=dim(inMat)[2],ncol=dim(inMat)[2])
  rownames(distOut) <- colnames(inMat)
  colnames(distOut) <- colnames(inMat)
  
  for(i in 1:dim(inMat)[2]){
    for(j in 1:dim(inMat)[2]){
      
      if(i!=j){
        x<- inMat[,i]
        y<- inMat[,j]
        isNa <- is.na(x)|is.na(y)
        
        thisDist <- NA
        
        if(sum(!isNa)>0){
          
          if(method=="pearson"){
            if(doAbs){
              thisDist<- 1-abs(cor(x[!isNa],y[!isNa], method = "pearson"))
            }else{
              thisDist<- 1-cor(x[!isNa],y[!isNa], method = "pearson")
            }
          }else if(method=="spearman"){
            
            if(doAbs){
              thisDist<- 1-abs(cor(x[!isNa],y[!isNa], method = "spearman"))
            }else{
              thisDist<- 1-cor(x[!isNa],y[!isNa], method = "spearman")
            }
            
          }else if(method=="euclidean"){
            thisDist <- sqrt(sum((x[!isNa]-y[!isNa])^2))
          }
          
          if(!is.na(thisDist)){
            distOut[i,j] <- thisDist
          }
        }
      }
    }
  }
  
  distOut
}
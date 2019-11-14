
plot_scatter <- function(x,y,xlab='',ylab='',mtext.cex=1,pointColors=NA){
  
  
  if(is.na(pointColors)){
    colorArray <- '#444444'
  }else{
    colorArray <- sapply(pointColors,color_function)
  }
  
  plot(x,y,xlab=xlab,ylab=ylab,pch=16,col=colorArray)
  
  corP <- cor.test(x,y)$p.value
  
  mtext(paste(c("p: ",format.pval(corP)),sep="",collapse=""),cex=mtext.cex)

  corP
  
}

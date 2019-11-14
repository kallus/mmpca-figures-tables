
plot_boxpoints <- function(xBox,yBox,ylab='',xlab='',mtext.cex=1,pointColors=NA){
  
  boxplot(yBox~xBox,ylab=ylab,xlab=xlab,xpd=T, bty="n", outline=F)
  
  if(is.na(pointColors)){
    colorArray <- sapply(xBox,color_function)
  }else{
    colorArray <- sapply(pointColors,color_function)
  }
  points(as.numeric(factor(xBox))+rnorm(0,0.1,n=length(xBox)),yBox,pch=16,cex=1,col=colorArray)
  
  lmSubtype <- lm(yBox~1+as.factor(xBox))
  summary(lmSubtype)
  anovaSubtype <- anova(lmSubtype)
  
  anovaP <- anovaSubtype[[5]][1]
  kruskalP <- kruskal.test(yBox~as.factor(xBox))$p.value
  
  mtext(paste(c("p: ",format.pval(kruskalP)),sep="",collapse=""),cex=mtext.cex)
  #mtext(paste(c("ANOVA p: ",format.pval(anovaP),", Kruskal-Wallis p: ",format.pval(kruskalP)),sep="",collapse=""),cex=mtext.cex)

  kruskalP
  
}

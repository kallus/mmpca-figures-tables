plot_survival <- function(y,clinical,z,i,annotAssoc){
  
  yNA <- is.na(y)	
  
  ySub <- y[!yNA]
  survivalSub <- clinical$survival[!yNA]
  vitalstatSub <- clinical$vitalstat[!yNA]
  ageSub <- clinical$age[!yNA]
  
  cM2 <- summary(coxph(Surv(survivalSub, vitalstatSub) ~ ySub))
  cM3 <- summary(coxph(Surv(survivalSub, vitalstatSub) ~ ySub + ageSub))
  
  sexSub <- clinical$sex[!yNA]
  sexSub[sexSub==""] <- NA 
  
  dataS <- data.frame(surv=survivalSub,vital=vitalstatSub,y=ySub,age=ageSub,sex=sexSub)
  dataS <- dataS[complete.cases(dataS),]
  
  dataS$sex <- as.numeric(factor(dataS$sex,levels=c("MALE","FEMALE")))-1
  
  cM4 <- summary(coxph(with(dataS,Surv(surv, vital)) ~ y + age + sex,data=dataS))
  
  
  stratPat <- array("Intermediate",length(ySub))
  stratPat[order(ySub)[1:floor(length(stratPat)/4)]] <- "Low"
  stratPat[order(ySub,decreasing=T)[1:floor(length(stratPat)/4)]] <- "High"
  
  include <- stratPat!="Intermediate"
  stratVect <- relevel(as.factor(stratPat[include]),ref="Low")
  
  cM <- coxph(Surv(survivalSub[include], vitalstatSub[include])~stratVect)
  modelP <- summary(cM)$logtest[3]
  
  dataS <- data.frame(surv=survivalSub[include],vital=vitalstatSub[include],strat=stratVect)
  dataS$SurvObj <- with(dataS, Surv(surv,vital))
  kmStrat<- npsurv(SurvObj ~ strat, data = dataS)
  
  survplot(fit  = kmStrat,
           conf = c("none","bands","bars")[2],
           xlab = "Days", ylab = "Survival",
           ## xlim(0,100),
           col=c("skyblue","tomato"),
           lwd=2,
           lty=1,
           label.curves = list(keys = "lines"),  # legend instead of direct label
           levels.only  = T,                    # show only levels, no label
           abbrev.label = FALSE,                    # if label used, abbreviate
           ## fun = function(x) {1 - x},            # Cumulative probability plot         
           loglog   = FALSE,                        # log(-log Survival) plot
           logt     = FALSE,                        # log time
           time.inc = 300,                           # time increment
           dots     = F,                        # dot grid
           n.risk   = F,                         # show number at risk
           ## srt.n.risk = 0, sep.n.risk = 0.056, adj.n.risk = 1,
           y.n.risk = -0.2, cex.n.risk = 0.6
  )
  
  mtext(paste("p:",format.pval(cM2$coefficients[,5])),line=0, cex=0.8)
  
  zOrderD <- order(z,decreasing=T)
  zOrderI <- order(z)
  
  #legend("topright",inset=c(0,0.1),horiz=F,legend=names(z[zOrderD[1:10]]),fill=c('#ffffff'), border=FALSE, bty="n", cex=0.7)
  #legend("topleft",inset=c(0,0.1),horiz=F,legend=names(z[zOrderI[1:10]]),fill=c('#ffffff'), border=FALSE, bty="n", cex=0.7)
  
  nTop <- 20
  
  zPlot <- c( z[zOrderD[1:nTop]], rev(z[zOrderI[1:nTop]]))
  
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
  
  #h <- barplot(zPlot, horiz=TRUE,names.arg=NA,cex.names=1,col=colArray,xlab="Comp")
  
  zPlotNegInd <- zPlot < 0
  zPlotPosInd <- zPlot >= 0
  
  #text(0,h[zPlotNegInd],names(zPlot)[zPlotNegInd],cex=0.6,pos=4,xpd=T)
  #text(0,h[zPlotPosInd],names(zPlot)[zPlotPosInd],cex=0.6,pos=2,xpd=T)
  
  return(list(cM2$coefficients[,5],cM3$coefficients[,5],cM4$coefficients[,5]))
  
}

plot_bar_signature <- function(y,clinical,z,i,annotAssoc, ...){
  
  yNA <- is.na(y)	
  
  ySub <- y[!yNA]
  survivalSub <- clinical$survival[!yNA]
  vitalstatSub <- clinical$vitalstat[!yNA]
  ageSub <- clinical$age[!yNA]
  
  cM2 <- summary(coxph(Surv(survivalSub, vitalstatSub) ~ ySub))
  cM3 <- summary(coxph(Surv(survivalSub, vitalstatSub) ~ ySub + ageSub))
  
  sexSub <- clinical$sex[!yNA]
  sexSub[sexSub==""] <- NA 
  
  dataS <- data.frame(surv=survivalSub,vital=vitalstatSub,y=ySub,age=ageSub,sex=sexSub)
  dataS <- dataS[complete.cases(dataS),]
  
  dataS$sex <- as.numeric(factor(dataS$sex,levels=c("MALE","FEMALE")))-1
  
  cM4 <- summary(coxph(with(dataS,Surv(surv, vital)) ~ y + age + sex,data=dataS))
  
  
  stratPat <- array("Intermediate",length(ySub))
  stratPat[order(ySub)[1:floor(length(stratPat)/4)]] <- "Low"
  stratPat[order(ySub,decreasing=T)[1:floor(length(stratPat)/4)]] <- "High"
  
  include <- stratPat!="Intermediate"
  stratVect <- relevel(as.factor(stratPat[include]),ref="Low")
  
  cM <- coxph(Surv(survivalSub[include], vitalstatSub[include])~stratVect)
  modelP <- summary(cM)$logtest[3]
  
  dataS <- data.frame(surv=survivalSub[include],vital=vitalstatSub[include],strat=stratVect)
  dataS$SurvObj <- with(dataS, Surv(surv,vital))
  kmStrat<- npsurv(SurvObj ~ strat, data = dataS)
  
#  survplot(fit  = kmStrat,
#           conf = c("none","bands","bars")[2],
#           xlab = "Days", ylab = "Survival",
#           ## xlim(0,100),
#           col=c("skyblue","tomato"),
#           lwd=2,
#           lty=1,
#           label.curves = list(keys = "lines"),  # legend instead of direct label
#           levels.only  = T,                    # show only levels, no label
#           abbrev.label = FALSE,                    # if label used, abbreviate
#           ## fun = function(x) {1 - x},            # Cumulative probability plot         
#           loglog   = FALSE,                        # log(-log Survival) plot
#           logt     = FALSE,                        # log time
#           time.inc = 300,                           # time increment
#           dots     = TRUE,                        # dot grid
#           n.risk   = F,                         # show number at risk
#           ## srt.n.risk = 0, sep.n.risk = 0.056, adj.n.risk = 1,
#           y.n.risk = -0.2, cex.n.risk = 0.6
#  )
  
  #mtext(paste("Comp.", i, ":",format.pval(cM2$coefficients[,5]),paste(format.pval(cM3$coefficients[,5]),collapse=" "),sep=", "),line=1,cex=0.8)
  
  zOrderD <- order(z,decreasing=T)
  zOrderI <- order(z)
  
  #legend("topright",inset=c(0,0.1),horiz=F,legend=names(z[zOrderD[1:10]]),fill=c('#ffffff'), border=FALSE, bty="n", cex=0.7)
  #legend("topleft",inset=c(0,0.1),horiz=F,legend=names(z[zOrderI[1:10]]),fill=c('#ffffff'), border=FALSE, bty="n", cex=0.7)
  
  nTop <- 20
  
  zPlot <- c( z[zOrderD[1:nTop]], rev(z[zOrderI[1:nTop]]))
  
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
  
  h <- barplot(zPlot, horiz=TRUE,names.arg=NA,cex.names=1,col=colArray, ...)
  par(xpd=T)
  legend('bottomleft', inset=c(-0.3, 0), lwd=4, lty=1, col=c('purple', '#cccccc'), legend=c('Meth.', 'GE'))
  par(xpd=F)
  mtext('Max. magnitude variables',line=0, cex=0.8)
  
  zPlotNegInd <- zPlot < 0
  zPlotPosInd <- zPlot >= 0
  
  text(0,h[zPlotNegInd],names(zPlot)[zPlotNegInd],cex=0.6,pos=4,xpd=T)
  text(0,h[zPlotPosInd],names(zPlot)[zPlotPosInd],cex=0.6,pos=2,xpd=T)
  
  return(list(cM2$coefficients[,5],cM3$coefficients[,5],cM4$coefficients[,5]))
  
}

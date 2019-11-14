
plot_survival_clusters <- function(y,clinical,clusterColors){
  
  yNA <- is.na(y)	
  
  ySub <- y[!yNA]
  survivalSub <- clinical$survival[!yNA]
  vitalstatSub <- clinical$vitalstat[!yNA]
  
  # Individual comparisons:
  uClust <- unique(ySub)
  clustP <- sapply(uClust,function(c){
    strat <- array(0,length(ySub))
    strat[ySub==c] <- 1
    
    dataS <- data.frame(surv=survivalSub,vital=vitalstatSub,strat=as.factor(strat))
    dataS$SurvObj <- with(dataS, Surv(surv,vital))
    
    summary(coxph(SurvObj ~ strat, data = dataS))$coefficients[,5]
    
  }
  )
  
  #Full model comparison:
  dataS <- data.frame(surv=survivalSub,vital=vitalstatSub,strat=as.factor(ySub))
  dataS$SurvObj <- with(dataS, Surv(surv,vital))
  anovaP <- anova(coxph(SurvObj ~ strat, data = dataS))$`Pr(>|Chi|)`[2]
  
  kmStrat<- npsurv(SurvObj ~ strat, data = dataS)
  
  kmStrataNames <- names(kmStrat$strata)
  survplot(fit  = kmStrat,
           conf = c("none","bands","bars")[1],
           xlab = "Days", ylab = "Survival",
           ## xlim(0,100),
           col=clusterColors,
           lwd=2,
           lty=1,
           label.curves = F,#list(keys = "lines"),  # legend instead of direct label
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
  
  mtext(paste("p:",format.pval(min(clustP))),line=0, cex=0.8)
  
}

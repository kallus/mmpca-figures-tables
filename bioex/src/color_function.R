color_function <- function(inD){
  
  require(RColorBrewer)
  thisPal <- brewer.pal(9,"Set1")
  bPal2 <- brewer.pal(12,"Paired")
  
  if(is.na(inD)){
    
    return('#ffffff')
    
  }else{
    
    if(inD==''){
      return("#cccccc")
    }
    
    if(inD=='HGCC'){
      return("#d85218")
    }
    if(inD=='TCGA'){
      return('#0072bc')
    }
    
    if(inD=='x'){
      return('#ffffff') 
    }else if(inD==''){
      return('#ffffff')
    }else if(inD=='NA'){
      return('#ffffff')
    }
    
    if(inD=='+'){
      return('#ff0000')
    }else if(inD=='-'){
      return('#00ff00')
    }
    
    if(inD=='G-CIMP'){
      return('#984EA3')
    }else if(inD=='NON G-CIMP'){
      return('#444444')
    }
    
    if(inD=='PN'){
      return('#984EA3') # purple
    }else if(inD=='CL'){
      return('#377EB8')  # blue
    }else if(inD=='NL'){
      return('#4DAF4A')  # green
    }else if(inD=="MS"){
      return('#E41A1C') # red
    }
    
    
    if(inD=='RTK I'){
      return('#984EA3') # purple
    }else if(inD=='RTK II'){
      return('#377EB8')  # blue
    }else if(inD=='MES'){
      return('#4DAF4A')  # green
    }
    
    if(inD=='I'){
      return('#ff0000')
    }else if(inD=='H'){
      return('#00ff00')
    }else if(inD=='L'){
      return('#0000ff')
    }
    
    if(inD=="F"){
      return('#E41A1C')
    }else if(inD=="M"){
      return('#377EB8')
    }
    
    if(is.numeric(inD)){
      return(thisPal[inD])
    }else{
      
      inD <- as.numeric(inD)
      if(is.numeric(inD)){
        return(thisPal[inD])
      }
    }
    
  }
  
  return('#ffffff')
}
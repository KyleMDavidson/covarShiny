
#sets things up to have 10 discrete breaks in the synthethized color gradient on heatmap
#max is, by default, max of all 
covImageui=function(inputcov,filename,maxColor){
  if (missing(maxColor)){
    maxColor = max(inputcov)
  }
  
  
  
  colorStep=maxColor / 10
  colors = c(seq(0,colorStep,length=100),seq(colorStep+.001,colorStep*2,length=100),seq(colorStep*2+.001,colorStep*3,length=100),
             seq(colorStep*3+.001,colorStep*4,length=100), seq(colorStep*4+.001,colorStep*5,length=100),
             seq(colorStep*5.001,colorStep*6,length=100), seq(colorStep*6+.001,colorStep*7,length=100),
             seq(colorStep*7.001,colorStep*8,length=100), seq(colorStep*8+.001,colorStep*9,length=100),
             seq(colorStep*9.001,colorStep*10,length=100))
  my_palette <- colorRampPalette(c("black", "blue"))(n = 999)
  
  library('gplots')
  
  
  condition=missing(filename)
  x=getwd()
  if (!condition){
    if (dir.exists('./covarImages')){
      setwd(file.path('.','./covarImages'))
    }else{
      dir.create('./covarImages',showWarnings=FALSE)
      setwd(file.path('.','./covarImages'))
    }
    jpeg(file=filename)
  }
  #cat (file=stderr(), nrow(inputcov))
  heatmap.2(as.matrix(inputcov[1:nrow(inputcov),1:ncol(inputcov)]), col=my_palette, dendrogram="none", 
            Rowv=FALSE, Colv=FALSE,breaks=colors, density.info="none", trace="none", 
            symm=F,symkey=F,symbreaks=T, scale="none")
  
  if (!condition){
    dev.off()
    setwd('./..')
    
  }
  
  
  inputcov
}

seqsAll=function(seqs, pos){
  for (i in 1:length(seqs)){
    print(seqs[[i]][pos])
  } 
}
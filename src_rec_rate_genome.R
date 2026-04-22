#### MapMarey

#######################
### Obtain recombination rate per windows using Loess (or lowess for LOcally WEighted Scatterplot Smoothing)
### function from MareyMap (Clement Rezvoy, Gabriel Marais)

### map = file with physical and genetical position
### span = the size of the window
### degree = the degree of the fitted curve
### data = a dataset of markers where recombination rates need to be estimated


recLoess<-function(map, data, span, degree){
  
  model <- loess(map$geneticalDistances ~ map$physicalPositions, span=span, degree=degree)
  
  #physical positions +1
  pp1 <- data$physicalPositions+1
  #physical positions -1
  pm1 <- data$physicalPositions-1
  #predicted genetical positions @ phys+1
  gp1 <- predict(model,newdata=pp1)
  #predicted genetical positions @ phys-1
  gm1 <- predict(model,newdata=pm1)
  Prates <- mapply(function(xa,xb,ya,yb){
    round((yb-ya)/((xb-xa)/1000000),2)
  },pm1,pp1,gm1,gp1)
  Prates
}

setwd("C:/Users/mgabrielli/Documents/Doctorat/full_genomes/rec_rate/")
data1=read.table("positions.recomb.txt",header=TRUE)

pdf("rec_chr.pdf",14,7)
par(mfrow=c(3,6))
m=matrix("", nrow = 0, ncol = 3)
# do loess fit for each chromosome
for (chr in levels(data1$map)){
  
  setwd("C:/Users/mgabrielli/Documents/Doctorat/rec_rate/")
  map=read.table("Backstrom_Recombination_Map_forMaryMap.txt",header=TRUE)
  colnames(map)=c("set","mkr","geneticalDistances","map","physicalPositions")
  map=map[which(map$map==chr),]
  map=map[order(map$physicalPositions),]
  
  if (nrow(map)>10){
    data=read.table("positions.recomb.txt",header=TRUE)
    colnames(data)=c("map","physicalPositions")
    data=data[which(data$map==chr),]
    data=data[order(data[,2]),]
    
    # span = 0.4
    # degree = 2 (polynomiale degree 2)
    recomb=recLoess(map,data,0.4,2)
    recomb[recomb<0]=0
    m=rbind(m,cbind(chr,data$physicalPositions,recomb))
    plot(recomb~data$physicalPositions,main=chr,xlab="pos",ylab="rec",cex=0.1)
    } else {
    if (nrow(map)>2){
      recomb=max(map$geneticalDistances)/(max(map$physicalPositions)/10^6)
      m=rbind(m,cbind(chr,data$physicalPositions,recomb))
    }
  }
}
dev.off()
colnames(m)=c("chr","pos","rec")
write.table(m,"rec_chr.txt",col.names=TRUE,row.names=FALSE,quote=FALSE)

############## figure 5

#### Z. Borbonicus and highlight high Fst with high and low recombination rate

### high Fst = top 0.05%
### high rec rata = greater than 1cM

lighten <- function(color, factor=0.5){
  col <- col2rgb(color)
  col <- col+(255-col)*factor
  col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
  col
}
library(qqman)

pdf("/Users/mgabrielli/Documents/Doctorat/full_genomes/scans_diversity/filt/noF/scans.Fst.chr.all.zobo.col.highlight.rec.highlow.pdf",width=12,height=10)
par(mfrow=c(6,1), mai = c(0.3, 0.3, 0.2, 0.1))
pops=c("W_S","W_NE","S_NE","NE_PN","W_PN","S_PF")
title=c("LBHB vs. BNB","GHB vs. LBHB","GHB vs. BNB","HIGH-North vs. GHB","HIGH-North vs. LBHB","HIGH-South vs. BNB")
cols=c("orchid","orchid","orchid","#7570b3","#7570b3","#7570b3")
for (i in 1:length(pops)){
  setwd("/Users/mgabrielli/Documents/Doctorat/full_genomes/scans_diversity/filt")
  scan=read.csv("scans_pi_dxy_fst_SNP.tab.chr.header.pseudoK",sep="\t")
  scan=scan[,c(ncol(scan)-1,ncol(scan),grep(paste0("Fst_",pops[i]),colnames(scan)))]
  
  if (length(which(scan[,3]=="NaN"))>0){
    scan=scan[-which(scan[,3]=="NaN"),] # remove sites with na when calculating scan
  }
  scan=as.data.frame(scan)
  scan$chr=gsub("chromosome_","chr",scan$chr)
  scan$code=paste(scan$chr,scan$pos,sep="_")
  scan$code=gsub(" ","",scan$code)
  
  # rec
  setwd("/Users/mgabrielli/Documents/Doctorat/full_genomes/rec_rate")
  rec=read.table("rec_chr.txt",header=TRUE)
  rec=rec[!is.na(rec[,3]),]
  rec$code=paste(rec$chr,rec$pos,sep="_")
  tab=merge(scan,rec,by="code")

  ## remove Z chromosome
  tab.auto=tab[!tab$chr.x%in%"chrZ",]
  out=tab.auto[tab.auto[,4]>quantile(tab.auto[,4],0.995),]
  out1=out[out[,7]>1,]
  out2=out[!out[,7]>1,]
  
  ## for Z chromosome
  tab.sex=tab[tab$chr.x%in%"chrZ",]
  outb=tab.sex[tab.sex[,4]>quantile(tab.sex[,4],0.995),]
  out1b=outb[outb[,7]>1,]
  out2b=outb[!outb[,7]>1,]
  
  out1=rbind(out1,out1b)
  out2=rbind(out2,out2b)
  
  scan[,1]=gsub(",", "", scan[,1])
  scan[,1]=gsub("chr","",scan[,1])
  scan[grep("1A",scan[,1]),1]="1.1"
  scan[grep("4A",scan[,1]),1]="4.4"
  scan[grep("z",scan[,1]),1]="30"
  scan[grep("Z",scan[,1]),1]="30"
  scan[grep("LGE22",scan[,1]),1]="29"
  
  scan[grep("-",scan[,3]),3]="0"
  
  
  scansubset<-scan[complete.cases(scan),]
  SNP<-c(1:(nrow(scansubset)))
  colnames(scansubset)=c("CHR","BP","P","code")
  mydf<-data.frame(SNP,scansubset)
  
  mydf$CHR=as.numeric(as.character(mydf$CHR))
  mydf$BP=as.numeric(as.character(mydf$BP))
  mydf$P=as.numeric(as.character(mydf$P))
  
  mydf=mydf[order(mydf$CHR, mydf$BP), ]
  
  out1=merge(out1,mydf,by="code")
  out1=out1$SNP
  
  out2=merge(out2,mydf,by="code")
  out2=out2$SNP
  
  if (i<length(pops)){
    manhattan2(mydf,chr="CHR",bp="BP",p="P",snp="SNP",logp=FALSE,chrlabs=c(1,"1A",c(2:4),"4A",c(5:27), "LGE22","Z"),main=title[i],ylab="Fst",cex=0.6,yaxt="n",xaxt="n",ylim=c(0,1.05),col=c(cols[i],lighten(cols[i])),font.main=3,highlight1 = out1, highlight2=out2)
    axis(2,at=c(0,0.5,1))
  } else {
    manhattan2(mydf,chr="CHR",bp="BP",p="P",snp="SNP",logp=FALSE,chrlabs=c(1,"1A",c(2:4),"4A",c(5:27), "LGE22","Z"),main=title[i],ylab="Fst",cex=0.6,yaxt="n",ylim=c(0,1.05),col=c(cols[i],lighten(cols[i])),font.main=3,highlight1 = out1, highlight2=out2)
    axis(2,at=c(0,0.5,1))
  }
}
dev.off()

##### function to highlight two sets in two colours
manhattan2<-function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10",
                                                                               "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05),
                      genomewideline = -log10(5e-08), highlight1 = NULL, highlight2 = NULL, logp = TRUE,
                      ...)
{
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x)))
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x)))
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x)))
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x)))
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]]))
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]]))
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]]))
    stop(paste(p, "column should be numeric."))
  d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]]))
    d = transform(d, SNP = x[[snp]])
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    options(scipen = 999)
    d$pos = d$BP/1e+06
    ticks = floor(length(d$pos))/2 + 1
    xlabel = paste("Chromosome", unique(d$CHR), "position(Mb)")
    labs = ticks
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in sort(unique(d$index))) {
      if (i == 1) {
        d[d$index == i, "pos"] <- d[d$index == i, "BP"]
      } else {
        lastbase <- lastbase + tail(subset(d, index == i - 1)$BP, 1)
        d[d$index == i, "pos"] <- d[d$index == i, "BP"] + lastbase
      }
      ticks <- c(ticks, mean(range(d[d$index == i, "pos"], na.rm = TRUE)))
    }
    xlabel <- "Chromosome"
    labs <- sort(unique(d$CHR))
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i",
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0,
                                                                     ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log10))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in%
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep(col, max(d$CHR))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      with(d[d$index == unique(d$index)[i], ], points(pos,
                                                      logp, col = col[icol], pch = 20, ...))
      icol = icol + 1
    }
  }
  if (suggestiveline)
    abline(h = suggestiveline, col = "blue")
  if (genomewideline)
    abline(h = genomewideline, col = "red")
  if (!is.null(highlight1)) {
    if (any(!(highlight1 %in% d$SNP)))
      warning("You're trying to highlight1 SNPs that don't exist in your results.")
    d.highlight1 = d[which(d$SNP %in% highlight1), ]
    with(d.highlight1, points(pos, logp, col = "green3", pch = 20,
                              ...))
  }
  if (!is.null(highlight2)) {
    if (any(!(highlight2 %in% d$SNP)))
      warning("You're trying to highlight2 SNPs that don't exist in your results.")
    d.highlight2 = d[which(d$SNP %in% highlight2), ]
    with(d.highlight2, points(pos, logp, col = "red", pch = 20,
                              ...))
  }
}

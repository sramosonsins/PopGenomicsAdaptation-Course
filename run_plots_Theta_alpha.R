library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0){
  setwd("./")
  show(sprintf("No filename assigned. Exiting"))
  quit()
}else{
  filename=args[1]
}
show(sprintf("filename: %s", filename))
filename.theta.neutral <- sprintf("%s.out.ms.neutral_variability.txt",filename)
filename.theta.functional <- sprintf("%s.out.ms.functional_variability.txt",filename)
filename.div.neutral <- sprintf("%s.out.ms.neutral_divergence.txt",filename)
filename.div.functional <- sprintf("%s.out.ms.functional_divergence.txt",filename)

Ps <- read.table(header=T,file=filename.theta.neutral)
Pn <- read.table(header=T,file=filename.theta.functional)
Ds <- as.numeric(read.table(header=F,file=filename.div.neutral))
Dn <- as.numeric(read.table(header=F,file=filename.div.functional))

Ps <-  Ps[1,c(3,1,2,4)]
Pn <-  Pn[1,c(3,1,2,4)]

alpha <- 1 - Pn/Ps * Ds/Dn

pdf(sprintf("%s.alpha.pdf",filename))
plot(x=c(1,2,3,4),y=as.matrix(alpha),type="b",ylim=c(-2,1),
     main=sprintf("%s. \nalpha from Theta",filename),
     xlab="Theta (FuLi,Watt,Taj,FayWu",ylab="alpha",
     col="red")

dev.off()

#a <- read.table(header=T,file="table_Theta_alpha.txt")
#am <- as.matrix(a[,-1])
#plot(am[1,],type="b",ylim=c(-2,1),main="alpha from Theta",xlab="Theta (FuLi,Watt,Taj,FayWu",ylab="alpha")
#for(i in 2:15) {lines(am[i,],type="b",col=rainbow(20)[i])}
#legend(x="bottomright", col = rainbow(20)[1:15],legend=as.array(a[,1]),lty=1)
#a

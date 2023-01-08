setwd("~/Desktop/Papers/Gene-editing/Submission - Figure/")
## BMI, height and IQ
load("summary_stats_GWS_Yengo_Savage_2018_v2.Rdata") # BMI, height and IQ

png("Figure_2.png",width=2500,height=2800,res=300)
plotQT <- function(dt,desired_direction){
  ## Trait Increasing Allele Frequency (TIAF)
  TIAF <- ifelse(dt[,"BETA"]>0,dt[,"EAF"],1-dt[,"EAF"])
  BETA <- abs(dt[,"BETA"])
  if(desired_direction==1){
    contribution_to_mean <- 2 * BETA * (1 - TIAF)
    Y <- cumsum(sort(contribution_to_mean,decreasing = TRUE))
  }else{
    contribution_to_mean <- 2 * BETA * TIAF
    Y <- cumsum(sort(contribution_to_mean,decreasing = TRUE))
  }
  return(Y)
}
pQT_ht    <- plotQT(dt=height,desired_direction=1)
pQT_iq    <- plotQT(dt=iq,desired_direction=1)

nset <- c(length(pQT_ht),length(pQT_iq))

qts   <- c("ht","iq")
qtTab <- matrix(NA,nrow=max(nset),ncol=length(nset))
for(i in 1:length(qts)){
  x <- get(paste0("pQT_",qts[i]))
  qtTab[1:length(x),i] <- x
}
colnames(qtTab) <- casefold(qts,upper = TRUE)
write.table(qtTab,"qtTab_HT_and_IQ.txt",quote=F,row.names = F,sep="\t")

Cols <- c(ht="dodgerblue",iq="coral1")
par(mar=c(5,5,3,1))
MaxAlleles <- 100
plot(c(1,MaxAlleles),c(0.2,10),type="n",log="xy",cex.lab=1.3,
     xlab="Number of Edited Allele(s)\n(Ranked from largest to smallest expected effect size)",
     ylab="Expected Increase in Population Mean (in trait SD)",
     axes=FALSE)
axis(1,at=c(1,2,5,10,20,50,MaxAlleles),line = -0.5)
axis(2,at=c(0.2,0.5,1,2,5,10),
     labels=c("0.2","0.5","1","2","5","10"))
for(i in 1:length(qts)){
  x <- qtTab[1:MaxAlleles,i]
  points(1:MaxAlleles,x,pch=14+i,col=Cols[qts[i]],type="b")
  lines(1:MaxAlleles,x,col=Cols[qts[i]],lwd=2)
  #abline(h=max(x,na.rm = T),col=Cols[qts[i]],lty=2)
  segments(x0=1,y0=max(x,na.rm = T),
           x1=MaxAlleles,y1=max(x,na.rm = T),
           col=Cols[qts[i]],lty=2)
}

legend(45,0.5,legend=c("Height","IQ"),
       col=Cols,pch=14+(1:10),
       box.lty=0,cex=1.1,lwd=2,lty=1)
legend(20,0.3,legend="Highest Reduction",lty=2,box.lty=0)
dev.off()

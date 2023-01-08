setwd("~/Desktop/Papers/Gene-editing/Submission - Figure/")
## FG -- ref. Lagou et al. NC (2021)
## WHR -- ref. Pullit et al. HMG (2018)
## Lipids  -- ref. Graham et al. Nature (2021)
## SBP -- Warren et al. (2021)

fg    <- read.table("FG_36snps.txt",h=T,stringsAsFactors = F,sep="\t")
Vp_fg <- mean( fg$N * 2 * fg$EAF * (1-fg$EAF) * (fg$SE^2) )
sd_fg <- sqrt(Vp_fg)
FG    <- fg[,c("SNP","EAF","BETA")]
FG[,"BETA"] <-  FG[,"BETA"] / sd_fg

whr  <- read.table("WHRadjBMI_Pullit_et_al2018.txt",h=T,stringsAsFactors = F)
WHR  <- whr[,c("SNP","frqA1.whr","beta.whr")]; colnames(WHR) <- c("SNP","EAF","BETA")
WHR  <- WHR[which(WHR$EAF>0.01 & WHR$EAF<0.99),]

glgc <- read.table("GLGC_Graham2021.txt",h=T,stringsAsFactors = F,sep="\t")
glgc <- glgc[which(glgc$Ancestry=="EUR"),]
ldl  <- glgc[which(glgc$Lipid=="LDL"),]
tg   <- glgc[which(glgc$Lipid=="logTG"),]

sd_ldl <- sqrt( mean( ldl$N * 2 * ldl$AF * (1-ldl$AF) * (ldl$Effectsize_SD^2) ) )
sd_tg  <- sqrt( mean( tg$N * 2 * tg$AF * (1-tg$AF) * (tg$Effectsize_SD^2) ) )

LDL <- ldl[,c("rs_dbSNP150","AF","Effectsize")]; colnames(LDL) <- c("SNP","EAF","BETA")
TG  <- tg[,c("rs_dbSNP150","AF","Effectsize")]; colnames(TG) <- c("SNP","EAF","BETA")

LDL <- LDL[which(LDL$EAF>0.01 & LDL$EAF<0.99),]
TG  <- TG[which(TG$EAF>0.01 & TG$EAF<0.99),]

icbp <- read.table("ICBP_Warren2022.txt",h=T,stringsAsFactors = F,sep="\t",comment.char = "!")
icbp <- icbp[which(icbp$CHR%in%(1:22)),]

sbp <- na.omit(icbp[which(icbp$Trait=="SBP"),])
dbp <- na.omit(icbp[which(icbp$Trait=="DBP"),])

sd_sbp <- sqrt( mean( sbp$Neff * 2 * sbp$EAF * (1-sbp$EAF) * (sbp$SE^2) ) )
sd_dbp <- sqrt( mean( dbp$Neff * 2 * dbp$EAF * (1-dbp$EAF) * (dbp$SE^2) ) )

SBP <- sbp[,c("rsID","EAF","Effect")]; colnames(SBP) <- c("SNP","EAF","BETA"); SBP[,"BETA"] <- SBP[,"BETA"] / sd_sbp
DBP <- dbp[,c("rsID","EAF","Effect")]; colnames(DBP) <- c("SNP","EAF","BETA"); DBP[,"BETA"] <- DBP[,"BETA"] / sd_dbp

## BMI, height and IQ -- focus on BMI
load("summary_stats_GWS_Yengo_Savage_2018_v2.Rdata") # BMI, height and IQ

png("Figure_1.png",width=5000,height=2800,res=300)
op <- par(mfrow=c(1,2))
## Get data
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
pQT_bmi    <- plotQT(dt=bmi,desired_direction=-1)
pQT_whr    <- plotQT(dt=WHR,desired_direction=-1)

pQT_fg     <- plotQT(dt=FG,desired_direction=-1)

pQT_ldl    <- plotQT(dt=LDL,desired_direction=-1)
pQT_tg     <- plotQT(dt=TG,desired_direction=-1)

pQT_sbp    <- plotQT(dt=SBP,desired_direction=-1)
pQT_dbp    <- plotQT(dt=DBP,desired_direction=-1)


nset <- c(length(pQT_bmi),length(pQT_whr),
          length(pQT_fg) ,length(pQT_ldl),
          length(pQT_tg),length(pQT_sbp),
          length(pQT_dbp)
)

qts <- c("bmi","whr","fg","ldl","tg","sbp","dbp")
qtTab <- matrix(NA,nrow=max(nset),ncol=length(nset))
for(i in 1:length(qts)){
  x <- get(paste0("pQT_",qts[i]))
  qtTab[1:length(x),i] <- x
}
#matplot(1:100,-qtTab[1:100,],type="b",pch=19)
colnames(qtTab) <- casefold(qts,upper = TRUE)

write.table(qtTab,"qtTab_09Jun2022.txt",quote=F,row.names = F,sep="\t")

Cols <- c(height="dodgerblue",bmi="coral1",iq="goldenrod")
Cols <- c("dodgerblue","coral1","goldenrod","darkred","lightgreen",
          "violet","pink")#,"darkgreen","grey","darkblue")
names(Cols) <- qts
par(mar=c(5,5,3,1))
MaxAlleles <- 100
plot(c(1,MaxAlleles),c(10,0.1),type="n",log="xy",cex.lab=1.3,
     xlab="Number of Edited Allele(s)\n(Ranked from largest to smallest expected effect size)",
     ylab="Expected Reduction Population Mean (in trait SD)",
     axes=FALSE)
axis(1,at=c(1,2,5,10,20,50,MaxAlleles),line = -0.5)
#axis(2)
axis(2,at=c(0.1,0.2,0.5,1,2,5,10),
     labels=c("-0.1","-0.2","-0.5","-1","-2","-5","-10"))
for(i in 1:length(qts)){
  #x <- get(paste0("pQT_",qts[i]))
  x <- qtTab[1:MaxAlleles,i]
  points(1:MaxAlleles,x,pch=14+i,col=Cols[qts[i]],type="b")
  lines(1:MaxAlleles,x,col=Cols[qts[i]],lwd=2)
  #abline(h=max(x,na.rm = T),col=Cols[qts[i]],lty=2)
  Ymax <- max(x,na.rm = T) 
  Xmax <- which.max(x)
  segments(x0=1,y0=Ymax,
           x1=Xmax,y1=Ymax,
           col=Cols[qts[i]],lty=2)
}

legend(50,0.5,legend=casefold(qts,upper = TRUE),col=Cols,pch=14+(1:10),
       box.lty=0,cex=1.1,lwd=2,lty=1)
legend(25,0.15,legend="Highest Reduction",lty=2,box.lty=0)

## Diseases
ad  <- read.table(file= "ad_marioni2018_ukbf.txt", header=TRUE)
mdd <- read.table(file= "mdd_howard2019_ukbf.txt", header=TRUE)
scz <- read.table(file= "scz_ng2014_ukbf.txt", header=TRUE) # need to add T2D and CAD later
t2d <- read.table(file="t2d_Vujkovic2020.txt",h=T,stringsAsFactors = F) # https://www.nature.com/articles/s41588-020-0637-y#Sec45
t2d <- t2d[which(t2d$P<5e-16),]
t2d <- t2d[,c("rsid","EAF","Beta")]
colnames(t2d) <- colnames(ad)

## CAD: https://www.medrxiv.org/content/10.1101/2021.05.24.21257377v1.full.pdf
cad <- read.table("cad_aragam2021.txt",h=T,stringsAsFactors = F)
cad <- cad[which(cad[,"P"]<5e-14),]

plotDisease <- function(dt,K){
  t    <- qnorm(1-K)
  z    <- dnorm(t)
  ## Risk Lowering Allelic Frequency
  RLAF <- ifelse(dt[,"Beta"]>0,dt[,"freq_A1"],1-dt[,"freq_A1"])
  BETA <- abs(dt[,"Beta"]) * K * (1 - K) / z
  contribution_to_mean <- 2 * RLAF * BETA
  Y    <- cumsum(sort(contribution_to_mean,decreasing = TRUE))
  K_edited <- pnorm(t,mean=-Y,sd=1,lower.tail = FALSE) # same variance?
  return(K_edited)
}

K_ad     <- 0.05
K_mdd    <- 0.15
K_scz    <- 0.01
K_t2d    <- 0.10
K_cad    <- 0.06

p_ad     <- plotDisease(dt=ad,K=K_ad)
p_mdd    <- plotDisease(dt=mdd,K=K_mdd)
p_scz    <- plotDisease(dt=scz,K=K_scz)
p_t2d    <- plotDisease(dt=t2d,K=K_t2d)
p_cad    <- plotDisease(dt=t2d,K=K_cad)

nsnpmax  <- max(c(length(p_ad),
                  length(p_mdd),
                  length(p_scz),
                  length(p_t2d),
                  length(p_cad)))
diseaseTab <- cbind.data.frame(T2D=c(p_t2d,rep(NA,nsnpmax-length(p_t2d))),
                               CAD=c(p_cad,rep(NA,nsnpmax-length(p_cad))),
                               SCZ=c(p_scz,rep(NA,nsnpmax-length(p_scz))),
                               MDD=c(p_mdd,rep(NA,nsnpmax-length(p_mdd))),
                               AD=c(p_ad,rep(NA,nsnpmax-length(p_ad))))

write.table(diseaseTab,"diseaseTab_v2.txt",quote=F,row.names = F,sep="\t")

ymin  <- min(c(p_ad,p_mdd,p_scz,p_t2d,p_cad))
ymax  <- 0.2
ColsD <- c(ad="dodgerblue",mdd="coral1",scz="goldenrod",t2d="darkred",cad="lightgreen")
par(mar=c(5,5,3,1))
plot(c(1,100),c(ymin,ymax),type="n",log="xy",cex.lab=1.2,
     xlab="Number of Edited Allele(s)\n(Ranked from largest to smallest expected effect size)",
     ylab="Expected Prevalence among Edited Genomes",
     axes=FALSE)
axis(1,line=-0.5)
axis(2,at=10**seq(-1,-8,by=-1),
     labels=c("10%","1%",".1%",".01%",".001%",".0001%",
              ".00001%",".000001%"))

for(i in 1:ncol(diseaseTab)){
  x <- diseaseTab[1:MaxAlleles,i]
  points(1:MaxAlleles,x,pch=14+i,col=ColsD[casefold(colnames(diseaseTab)[i],upper = FALSE)],type="b")
  lines(1:MaxAlleles,x,col=ColsD[casefold(colnames(diseaseTab)[i],upper = FALSE)],lwd=2)
  #abline(h=min(x,na.rm = T),col=ColsD[casefold(colnames(diseaseTab)[i],upper = FALSE)],lty=2)
  Ymin <- min(x,na.rm = T)
  Xmin <- which.min(x)
  segments(x0=1,y0=Ymin,
           x1=Xmin,y1=Ymin,
           col=ColsD[casefold(colnames(diseaseTab)[i],upper = FALSE)],lty=2)
}

# points(1:length(p_ad),c(p_ad),pch=15,col=ColsD["ad"],type="b")
# points(1:length(p_mdd),c(p_mdd),pch=16,col=ColsD["mdd"],type="b")
# points(1:length(p_scz),c(p_scz),pch=17,col=ColsD["scz"],type="b")
# points(1:length(p_t2d),c(p_t2d),pch=18,col=ColsD["t2d"],type="b")
# points(1:length(p_cad),c(p_cad),pch=19,col=ColsD["cad"],type="b")

# abline(h=min(p_ad),col=ColsD["ad"],lty=2)
# abline(h=min(p_mdd),col=ColsD["mdd"],lty=2)
# abline(h=min(p_scz),col=ColsD["scz"],lty=2)
# abline(h=min(p_t2d),col=ColsD["t2d"],lty=2)
# abline(h=min(p_cad),col=ColsD["cad"],lty=2)

legend(1,1e-5,legend=c("Alzheimer's (5%)","MDD (15%)","Schizophrenia (1%)",
                       "Type 2 Diabetes (10%)","Coronary Artery Disease (6%)"),col=ColsD,pch=15:19,
       box.lty=0,cex=1.1,title="Disease (Prevalence)",title.adj = .1,title.col = 2)
legend(1,1e-7,legend="Lowest value",lty=2,box.lty=0)
par(op)
dev.off()

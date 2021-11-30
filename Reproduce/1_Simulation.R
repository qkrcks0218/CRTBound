library("CRTBound")

###################################
# Data cleaning
###################################

dir <- "http://web.hku.hk/~bcowling/influenza/data/HongKongNPIstudyV3/"

LabPCR <- read.csv(paste(dir, "home_pcr.csv", sep=""))
HHchar <- read.csv(paste(dir, "hchar_h.csv", sep=""))
HHadhere <- read.csv(paste(dir, "adherence_m.csv", sep=""))

mark <- data.frame(hhID = unique(HHadhere$hhID))                                            # unique hhID
LabPCR <- merge(LabPCR,mark,by="hhID",all.y=TRUE)                                           # merge all hhID in HHadhere
LabPCR <- LabPCR[order(LabPCR$hhID,LabPCR$member,LabPCR$visit),]                            # sort
AData <- data.frame(hhID = HHadhere$hhID, member = HHadhere$member)                         # Analyzing data

hc.temp <- reshape(LabPCR, timevar="visit", idvar=c("hhID","member"), direction="wide", v.names="PCR")
AData <- merge(AData,hc.temp, by=c("hhID","member"), all.x=TRUE)
names(AData) <- c("hhID","member","V1","V2","V3")      # V1, V2, V3 : Outcome at visiting 1,2,3, 1337 contacts, 322 households

## exd_index: none of V0/V1 culture is A/B; Contact exclusion: V1 culture is A/B

for (i in 1:nrow(AData)){
    if(AData$member[i]==0 & ( (!is.na(AData$V1[i]) & AData$V1[i]==0)|is.na(AData$V1[i]) )) {
        AData$exd_index[i] <- 1
    } else {
        AData$exd_index[i] <- 0
    }
    
    if(AData$member[i]!=0 &  !is.na(AData$V1[i]) & AData$V1[i]!=0) {
        AData$d_contact[i] <- 1
    } else {
        AData$d_contact[i] <- 0
    }
}

# exd_index = 1 : Exclude index member because they were not infected or NA
# d_contact = 1 : Exclude other members because they were infected

d_contactid <- unique(AData$hhID[AData$d_contact==1|AData$exd_index==1]) # for excluding co-index households
d_contact <- data.frame(hhID=d_contactid)
d_contact$d_contact <- 1

AData <- merge(AData[,1:5],d_contact,all.x=TRUE)
AData.valid <- AData[is.na(AData$d_contact),-6]         # 1053 contacts, 259 households

# Following Cowling et al. (2009), replace NA values of V1,V2,V3 with 0
AData.valid$V1[is.na(AData.valid$V1)] <- 0
AData.valid$V2[is.na(AData.valid$V2)] <- 0
AData.valid$V3[is.na(AData.valid$V3)] <- 0

# Inserting intervention
intervention <- read.csv(paste(dir, "randomarm_407.csv", sep=""))
intervention$hhID <- as.numeric(substr(intervention$hhID,5,7))
AData.valid <- merge(AData.valid,intervention[,1:2],by="hhID")
AData.valid <- merge(AData.valid,HHadhere,by=c("hhID","member"),all.x=TRUE)

sum(unique(cbind(AData.valid$hhID,AData.valid$intervention))[,2]==1)    # 90 clusters under control
sum(unique(cbind(AData.valid$hhID,AData.valid$intervention))[,2]==3)    # 86 clusters under hand sanitizer
sum(unique(cbind(AData.valid$hhID,AData.valid$intervention))[,2]==4)    # 83 clusters under Facemask + hand sanitizer

# Remove households that include contacts whose self-reported mask and/or hand sanitizer usage are NA
errorHH <- data.frame(hhID=unique(AData.valid[is.na(AData.valid$mask+AData.valid$handrub),1]),error=1)
AData.valid.noerror <- merge(AData.valid,errorHH,by="hhID",all.x=TRUE)
AData.valid.noerror <- AData.valid.noerror[is.na(AData.valid.noerror$error),]
AData.valid <- AData.valid.noerror[,-dim(AData.valid.noerror)[2]]       # 937 contacts, 232 households

# Remove households that include contacts whose covariates are missing
errorHH <- data.frame(hhID=unique(AData.valid[is.na(AData.valid$age+AData.valid$male+AData.valid$vaccine08),1]),error=1)
AData.valid.noerror <- merge(AData.valid,errorHH,by="hhID",all.x=TRUE)
AData.valid.noerror <- AData.valid.noerror[is.na(AData.valid.noerror$error),]
AData.valid <- AData.valid.noerror[,-dim(AData.valid.noerror)[2]]       # 919 contacts, 227 households

# Focusing on A=0 (control) vs A=1 (Facemask + Hand sanitizer) by removing households that were recommended to use facemasks only
AData.valid <- AData.valid[AData.valid$intervention!=3,]                # 615 contacts, 151 households

# Remove redundant variables
AData.valid <- AData.valid[,c("hhID","member","V2","V3","intervention","mask","handrub","male","age","vaccine08")]

###################################
# Generate Simulated Data
###################################

probCS <- function(X){
    age <- X[,3]
    gender <- X[,2]
    probNT <- probAT <- probCO <- rep(0,dim(X)[1])
    for(jj in 1:dim(X)[1]){
        if(gender[jj]==1){
            PCO <- exp( -(age[jj]-20)*(age[jj]-50)/150 )
            PAT <- exp( -(age[jj]-40)*(age[jj]-65)/150 )
            PNT <- 1
            probCO[jj] <- PCO/(PNT+PAT+PCO)
            probNT[jj] <- PNT/(PNT+PAT+PCO)
            probAT[jj] <- PAT/(PNT+PAT+PCO)
        } else {
            PCO <- exp( -(age[jj]-20)*(age[jj]-50)/20 )
            PAT <- exp( -(age[jj]-40)*(age[jj]-65)/20 )
            PNT <- 1
            probCO[jj] <- PCO/(PNT+PAT+PCO)
            probNT[jj] <- PNT/(PNT+PAT+PCO)
            probAT[jj] <- PAT/(PNT+PAT+PCO)
        }
    }
    
    return(cbind(probNT,probAT,probCO))
    
}

###################################
# Compliance Probability
###################################

# png("Compliance.png",width=7,height=2,unit="in",res=500)

X.grid <- expand.grid(1,c(0,1),seq(0,100,1))
PROB <- probCS(X.grid)

layout(matrix(c(1,2,3,4,5,6),1,3,byrow=T),widths=c(3,3,1),heights=c(1))
par(mar=c(3,3,2,1))
plot(X.grid[X.grid[,2]==1,3], PROB[X.grid[,2]==1,1],xlim=c(0,100),ylim=c(0,1),type='l',xlab="",ylab="")
par(new=T)
plot(X.grid[X.grid[,2]==1,3], PROB[X.grid[,2]==1,2],xlim=c(0,100),ylim=c(0,1),type='l',col=2,xlab="",ylab="")
par(new=T)
plot(X.grid[X.grid[,2]==1,3], PROB[X.grid[,2]==1,3],ylim=c(0,1),type='l',col=4,xlab="",ylab="")
title(ylab="Prob",line=2)
title(main="Male",line=1)

plot(X.grid[X.grid[,2]==0,3], PROB[X.grid[,2]==0,1],xlim=c(0,100),ylim=c(0,1),type='l',xlab="",ylab="")
par(new=T)
plot(X.grid[X.grid[,2]==0,3], PROB[X.grid[,2]==0,2],xlim=c(0,100),ylim=c(0,1),type='l',col=2,xlab="",ylab="")
par(new=T)
plot(X.grid[X.grid[,2]==0,3], PROB[X.grid[,2]==0,3],xlim=c(0,100),ylim=c(0,1),type='l',col=4,xlab="",ylab="")
title(main="Female",line=1)

plot.new()
segments(0,0.7,0.3,0.7,col=1)
text(0.4,0.7,"NT",pos=4)
segments(0,0.5,0.3,0.5,col=2)
text(0.4,0.5,"AT",pos=4)
segments(0,0.3,0.3,0.3,col=4)
text(0.4,0.3,"CO",pos=4)

# dev.off()

###################################
# Reform Data
###################################

HHID <- AData.valid$hhID

XX <- AData.valid[,c("male","age","vaccine08")]
ff <- formula("~male+age+I(age^2)+vaccine08")

Pre.Cov <- model.matrix(ff,XX)
colnames(Pre.Cov)[1] <- "Const"

PROB <- probCS(Pre.Cov)
probNT <- PROB[,1]
probAT <- PROB[,2]
probCO <- PROB[,3]

set.seed(8888)

CompSt <- rep("NT",dim(Pre.Cov)[1])
for(jj in 1:dim(Pre.Cov)[1]){
    CompSt[jj] <- sample(c("NT","AT","CO"),1,prob=c(probNT[jj],probAT[jj],probCO[jj]))
}

A0 <- as.numeric(CompSt=="AT")
A1 <- as.numeric(CompSt!="NT")

prop.of.trt1 <- prop.of.trt0 <- rep(0,dim(Pre.Cov)[1])
for(jj in 1:dim(Pre.Cov)[1]){
    groupind <- which(HHID==HHID[jj])
    groupind <- setdiff(groupind,jj)
    if(length(groupind)>0){
        prop.of.trt1[jj] <- sum(CompSt[groupind]!="NT")/length(groupind)
        prop.of.trt0[jj] <- sum(CompSt[groupind]=="AT")/length(groupind)
    } else {
        prop.of.trt1[jj] <- 0
        prop.of.trt0[jj] <- 0
    }
    
}

para <- c(-3,2,4,2)
probA0 <- rep(para[1],dim(Pre.Cov)[1]) +                                   # Base (-3)
    para[2]*as.numeric(CompSt=="AT") +                               # Treatment effect (2)
    para[3]*AData.valid[,"vaccine08"]*as.numeric(CompSt=="AT") +          # Additional treatment effect on individuals that have been vaccinated (4)
    para[4]*prop.of.trt0*as.numeric(CompSt!="AT")                    # spillover effect from AT to nonAT (2)

probA1 <- rep(para[1],dim(Pre.Cov)[1]) +                                   # Base (-3)
    para[2]*as.numeric(CompSt!="NT") +                               # Treatment effect (2)
    para[3]*AData.valid[,"vaccine08"]*as.numeric(CompSt!="NT") +          # Additional treatment effect on individuals that have been vaccinated (4)
    para[4]*prop.of.trt1*as.numeric(CompSt=="NT")                    # spillover effect from nonNT to NT (2)


Y0 <- rbinom(dim(Pre.Cov)[1],1,1/(1+exp(-probA0)))
Y1pot <- rbinom(dim(Pre.Cov)[1],1,1/(1+exp(-probA1)))
Y1 <- apply(cbind(Y1pot,Y0),1,max)

set.seed(1234)

Zc <- rep(0,length(unique(HHID)))
Zc[sample(151,72)] <- 1
Z <- rep(0,length(HHID))
for(ii in 1:151){
    Z[which(HHID==sort(unique(HHID))[ii])] <- Zc[ii]
}

Simulated.Data <- Simulation.Reform(Y0,Y1,Z,A0,A1,HHID,Pre.Cov,seed=8888)

###################################
# To use the same simulated data set in Park and Kang (2021), run
###################################

# load("Simulation.RData")

###################################
# True Effects
###################################

TrueITT <- mean( Simulated.Data$Y1-Simulated.Data$Y0 )
TrueNT  <- sum( (Simulated.Data$Y1-Simulated.Data$Y0)*as.numeric( Simulated.Data$CompSt=="NT" ) )/ sum(as.numeric( Simulated.Data$CompSt=="NT" ))
TrueAT  <- sum( (Simulated.Data$Y1-Simulated.Data$Y0)*as.numeric( Simulated.Data$CompSt=="AT" ) )/ sum(as.numeric( Simulated.Data$CompSt=="AT" ))
TrueCO  <- sum( (Simulated.Data$Y1-Simulated.Data$Y0)*as.numeric( Simulated.Data$CompSt=="CO" ) )/ sum(as.numeric( Simulated.Data$CompSt=="CO" ))

###################################
# Simulation
###################################

TotalIter <- 10 # Number of total iteration
SSsize <- 1000    # Number of bootstrap CIs for bounds

RESULT11 <- matrix(0,TotalIter,4)
RESULT21 <- matrix(0,TotalIter,20)
RESULT3B1 <- matrix(0,TotalIter,15)
RESULT3B2 <- matrix(0,TotalIter,15)
RESULT3B3 <- matrix(0,TotalIter,15)
RESULT3B4 <- matrix(0,TotalIter,15)

colnames(RESULT11) <- c("Estimate","True","Pvalue","Coverage")
colnames(RESULT21) <- c(sprintf("Estimate_%0.1d",0:4),
                        sprintf("True_%0.1d",0:4),
                        sprintf("Pvalue_%0.1d",0:4),
                        sprintf("Coverage_%0.1d",0:4))
colnames(RESULT3B1) <-
    colnames(RESULT3B2) <-
    colnames(RESULT3B3) <-
    colnames(RESULT3B4) <- c(sprintf("Estimate_%s",c("NTLB","NTUB","ATLB","ATUB","COLB","COUB")),
                             sprintf("True_%s",c("NTLB","NTUB","ATLB","ATUB","COLB","COUB")),
                             sprintf("Coverage_%s",c("NT","AT","CO")))

COVER1 <- function(p,s,t){
    as.numeric(p-s*qnorm(0.975)<=t & p+s*qnorm(0.975)>=t )
}
COVER2 <- function(p,t){
    as.numeric(p[1]<=t[1] & p[2]>=t[2] )
}

for(iter in 1:TotalIter){
    
    set.seed(iter)
    
    Zc <- rep(0,length(unique(HHID)))
    Zc[sample(151,72)] <- 1
    Z <- rep(0,length(HHID))
    for(ii in 1:151){
        Z[which(HHID==sort(unique(HHID))[ii])] <- Zc[ii]
    }
    
    Simulated.Data <- Simulation.Reform(Y0,Y1,Z,A0,A1,HHID,Pre.Cov,seed=8888)
    
    
    TrueITT <- mean( Simulated.Data$Y1-Simulated.Data$Y0 )
    TrueNT  <- sum( (Simulated.Data$Y1-Simulated.Data$Y0)*as.numeric( Simulated.Data$CompSt=="NT" ) )/ sum(as.numeric( Simulated.Data$CompSt=="NT" )) # 41/248
    TrueAT  <- sum( (Simulated.Data$Y1-Simulated.Data$Y0)*as.numeric( Simulated.Data$CompSt=="AT" ) )/ sum(as.numeric( Simulated.Data$CompSt=="AT" )) # 14/90
    TrueCO  <- sum( (Simulated.Data$Y1-Simulated.Data$Y0)*as.numeric( Simulated.Data$CompSt=="CO" ) )/ sum(as.numeric( Simulated.Data$CompSt=="CO" )) # 92/277
    
    ################################################
    # ITT effect analysis
    ################################################
    
    ITTest <- ITT(Simulated.Data, Input.Type="Sim")
    
    RESULT11[iter,] <- c(ITTest[c(1,4)], 2-2*pnorm(abs(ITTest[3])), COVER1(ITTest[1],ITTest[2],ITTest[4]))
    
    HTEest <- HTE(Simulated.Data, Xvar=1:5, constant=TRUE , Input.Type="Sim")
    
    RESULT21[iter,] <- c(HTEest$Estimate[1,],HTEest$Estimate[4,],
                         2-2*pnorm(abs(HTEest$Estimate[3,])),
                         COVER1(HTEest$Estimate[1,],HTEest$Estimate[2,],HTEest$Estimate[4,]))
    
    ################################################
    # Compliance group effect analysis
    ################################################
    
    BOUND1 <- SharpBound(Simulated.Data , paraC=1:5, method="Linear", CIcalc=TRUE, SSsize=SSsize, level=0.95, Input.Type="Sim")
    BOUND2 <- SharpBound(Simulated.Data , paraC=1:5, method="Linear", CIcalc=TRUE, SSsize=SSsize, level=0.95, Input.Type="Sim")
    BOUND3 <- LongHudgens(Simulated.Data , paraC=c(2,5), CIcalc=TRUE, SSsize=SSsize, level=0.95, Input.Type="Sim")
    BOUND4 <- Bound.Intersect(BOUND2, BOUND3, level=0.95)
    
    RESULT3B1[iter,] <- as.numeric(c(BOUND1$Bound[1,],BOUND1$Bound[2,],
                                     COVER2(BOUND1$BootCIBound[1:2],BOUND1$Bound[2,1:2]),
                                     COVER2(BOUND1$BootCIBound[3:4],BOUND1$Bound[2,3:4]),
                                     COVER2(BOUND1$BootCIBound[5:6],BOUND1$Bound[2,5:6])))
    RESULT3B2[iter,] <- as.numeric(c(BOUND2$Bound[1,],BOUND2$Bound[2,],
                                     COVER2(BOUND2$BootCIBound[1:2],BOUND2$Bound[2,1:2]),
                                     COVER2(BOUND2$BootCIBound[3:4],BOUND2$Bound[2,3:4]),
                                     COVER2(BOUND2$BootCIBound[5:6],BOUND2$Bound[2,5:6])))
    RESULT3B3[iter,] <- as.numeric(c(BOUND3$Bound.Adj[1,],BOUND3$Bound.Adj[2,],
                                     COVER2(BOUND3$BootCIBound.Adj[1:2],BOUND3$Bound.Adj[2,1:2]),
                                     COVER2(BOUND3$BootCIBound.Adj[3:4],BOUND3$Bound.Adj[2,3:4]),
                                     COVER2(BOUND3$BootCIBound.Adj[5:6],BOUND3$Bound.Adj[2,5:6])))
    RESULT3B4[iter,] <- as.numeric(c(BOUND4$Bound[1,],BOUND4$Bound[2,],
                                     COVER2(BOUND4$BootCIBound[1:2],BOUND4$Bound[2,1:2]),
                                     COVER2(BOUND4$BootCIBound[3:4],BOUND4$Bound[2,3:4]),
                                     COVER2(BOUND4$BootCIBound[5:6],BOUND4$Bound[2,5:6])))
    
    print(iter)
    
}

###################################
# Summary
###################################

# ITT
hist(RESULT11[,1]) ; abline(v=RESULT11[1,2],col=2) # ITT estimator is consistent
mean(RESULT11[,3]) # Average p-value
mean(RESULT11[,4]) # Coverage >= 0.95

# HTT
par(mfrow=c(2,3))
for(jj in 1:5){
    hist(RESULT21[,jj]) # Heterogeneous parameters are consistently estimated
    abline(v=RESULT21[,jj+5],col=2)
}
apply(RESULT21[,11:15],2,FUN="mean") # Average p-value
apply(RESULT21[,16:20],2,FUN="mean") # Coverage >= 0.95

# Bound
rbind(RESULT3B1[1,7:12],
      RESULT3B2[1,7:12],
      RESULT3B3[1,7:12],
      RESULT3B4[1,7:12]) # compare bounds

BB <- RESULT3B2 # RESULT3B2, RESULT3B3, RESULT3B4
for(jj in 1:6){
    hist(BB[,jj],breaks=seq(0,1,length=101)) # Bound Estimators are consistent
    abline(v=BB[,jj+6],col=2)
}
apply(BB[,13:15],2,FUN="mean") # Coverage >= 0.95


###################################
# Summary of the original simulation in Park and Kang (2021)
###################################

# ITT

RESULT11 <- read.csv("1_SimResult/RESULT11.csv")
RESULT21 <- read.csv("1_SimResult/RESULT21.csv")
RESULT3B1 <- read.csv("1_SimResult/RESULT3B1.csv")
RESULT3B2 <- read.csv("1_SimResult/RESULT3B2.csv")
RESULT3B3 <- read.csv("1_SimResult/RESULT3B3.csv")
RESULT3B4 <- read.csv("1_SimResult/RESULT3B4.csv")



RR <- rbind( as.numeric(c(RESULT11$True[1],RESULT21[1,6:10])),
             as.numeric(c(mean(RESULT11$Estimate),apply(RESULT21[,1:5],2,mean))),
             as.numeric(c(mean(RESULT11$Estimate),apply(RESULT21[,1:5],2,mean)))-as.numeric(c(RESULT11$True[1],RESULT21[1,6:10])),
             as.numeric(c(sd(RESULT11$Estimate),apply(RESULT21[,1:5],2,sd))),
             as.numeric(c(mean(RESULT11$Pvalue),apply(RESULT21[,11:15],2,mean))),
             as.numeric(c(mean(RESULT11$Coverage),apply(RESULT21[,16:20],2,mean))))

colnames(RR) <- c("Overall ITT","Const","Male","Age","Age^2","Vaccine08")
rownames(RR) <- c("True","Estimate","Bias","Empirical SE","Average p-value","Coverage")

RR

# True Bounds

RR <- rbind( as.numeric(c(TrueNT, RESULT3B1[1, 7: 8], RESULT3B2[1, 7: 8], RESULT3B3[1, 7: 8], RESULT3B4[1, 7: 8])),
             as.numeric(c(TrueAT, RESULT3B1[1, 9:10], RESULT3B2[1, 9:10], RESULT3B3[1, 9:10], RESULT3B4[1, 9:10])),
             as.numeric(c(TrueCO, RESULT3B1[1,11:12], RESULT3B2[1,11:12], RESULT3B3[1,11:12], RESULT3B4[1,11:12])) )

colnames(RR) <- c("True","(C1)LB","(C1)UB","(C2)LB","(C2)UB","(C3)LB","(C3)UB","(C4)LB","(C4)UB")
rownames(RR) <- c("NT","AT","CO")

RR

# Estimated Bounds

RR <- rbind( as.numeric( c( apply(RESULT3B1[,1:2],2,mean)-RESULT3B1[1, 7: 8],
                            apply(RESULT3B2[,1:2],2,mean)-RESULT3B2[1, 7: 8],
                            apply(RESULT3B3[,1:2],2,mean)-RESULT3B3[1, 7: 8],
                            apply(RESULT3B4[,1:2],2,mean)-RESULT3B4[1, 7: 8])),
             as.numeric( c( apply(RESULT3B1[,1:2],2,sd),
                            apply(RESULT3B2[,1:2],2,sd),
                            apply(RESULT3B3[,1:2],2,sd),
                            apply(RESULT3B4[,1:2],2,sd)) ),
             c(mean(RESULT3B1$Coverage_NT),mean(RESULT3B1$Coverage_NT),
               mean(RESULT3B2$Coverage_NT),mean(RESULT3B2$Coverage_NT),
               mean(RESULT3B3$Coverage_NT),mean(RESULT3B3$Coverage_NT),
               mean(RESULT3B4$Coverage_NT),mean(RESULT3B4$Coverage_NT)),
             as.numeric( c( apply(RESULT3B1[,3:4],2,mean)-RESULT3B1[1, 9:10],
                            apply(RESULT3B2[,3:4],2,mean)-RESULT3B2[1, 9:10],
                            apply(RESULT3B3[,3:4],2,mean)-RESULT3B3[1, 9:10],
                            apply(RESULT3B4[,3:4],2,mean)-RESULT3B4[1, 9:10])),
             as.numeric( c( apply(RESULT3B1[,3:4],2,sd),
                            apply(RESULT3B2[,3:4],2,sd),
                            apply(RESULT3B3[,3:4],2,sd),
                            apply(RESULT3B4[,3:4],2,sd)) ),
             c(mean(RESULT3B1$Coverage_AT),mean(RESULT3B1$Coverage_AT),
               mean(RESULT3B2$Coverage_AT),mean(RESULT3B2$Coverage_AT),
               mean(RESULT3B3$Coverage_AT),mean(RESULT3B3$Coverage_AT),
               mean(RESULT3B4$Coverage_AT),mean(RESULT3B4$Coverage_AT)),
             as.numeric( c( apply(RESULT3B1[,5:6],2,mean)-RESULT3B1[1,11:12],
                            apply(RESULT3B2[,5:6],2,mean)-RESULT3B2[1,11:12],
                            apply(RESULT3B3[,5:6],2,mean)-RESULT3B3[1,11:12],
                            apply(RESULT3B4[,5:6],2,mean)-RESULT3B4[1,11:12])),
             as.numeric( c( apply(RESULT3B1[,5:6],2,sd),
                            apply(RESULT3B2[,5:6],2,sd),
                            apply(RESULT3B3[,5:6],2,sd),
                            apply(RESULT3B4[,5:6],2,sd)) ),
             c(mean(RESULT3B1$Coverage_CO),mean(RESULT3B1$Coverage_CO),
               mean(RESULT3B2$Coverage_CO),mean(RESULT3B2$Coverage_CO),
               mean(RESULT3B3$Coverage_CO),mean(RESULT3B3$Coverage_CO),
               mean(RESULT3B4$Coverage_CO),mean(RESULT3B4$Coverage_CO)))

colnames(RR) <- c("(C1)LB","(C1)UB","(C2)LB","(C2)UB","(C3)LB","(C3)UB","(C4)LB","(C4)UB")
rownames(RR) <- c("NT_Bias","NT_SE","NT_Coverage",
                  "AT_Bias","AT_SE","AT_Coverage",
                  "CO_Bias","CO_SE","CO_Coverage")
RR

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

# Focusing on households that (intervention-onset)<=36hr
# Removing index contact
# Using index contact's covariates as cluster-level covariates
threshold <- 36
AData.subgroup <- merge(AData.valid,HHchar[,c(1,4,9)],by="hhID")
AData.subgroup.index <- AData.subgroup[AData.subgroup$member==0,c(1,8,9,10)]
colnames(AData.subgroup.index) <- c("hhID","i.male","i.age","i.vaccine08")
AData.subgroup <- merge(AData.subgroup,AData.subgroup.index,by="hhID")
AData.subgroup <- AData.subgroup[abs(AData.subgroup$onset_v1_delay)<=threshold,]
HHsize <- AData.subgroup[AData.subgroup$member==0,]$house_size
AData.subgroup <- AData.subgroup[AData.subgroup$member>0,]              # 290 non-index contacts, 96 households

####################################################
# Balance Check
####################################################

Y <- 1-apply(cbind(AData.subgroup$V2,AData.subgroup$V3),1,max)
Z <- as.numeric(AData.subgroup$intervention==4)
C <- AData.subgroup$hhID
A <- as.numeric(AData.subgroup$mask<4&AData.subgroup$handrub<4)
Class <- unique(C)
N <- length(Y)
n <- rep(0,N)
for(jj in 1:N){
  n[jj] <- sum(C==C[jj])
}
AData.subgroup$n <- n

XX <- AData.subgroup[,c("male","age","vaccine08","i.male","i.age","i.vaccine08","house_size","n")]
ff <- formula("~male+age+I(age^2)+vaccine08+
              i.male+i.age+I(i.age^2)+i.vaccine08+
              I(n==3)+I(n==4)+I(n>=5)+house_size")
X <- model.matrix(ff,XX)[,1:13]

Reformed.Data <- Data.Reform(Y,Z,A,C,X,seed=1234)

c1 <- c("","Number of clusters","Number of individuals","Number of individuals who actually took treatment","Average of outcomes",
        "Average cluster size","Proportion of male","Average of age","Proportion of vaccinated individuals",
        "Proportion of male index individuals","Average of age of index individuals","Proportion of vaccinated index individuals",
        "House size")
c2 <- c("Treated",sum(Reformed.Data$Zc),sum(Reformed.Data$Z),sum(Reformed.Data$A[Reformed.Data$Z==1]),round(mean(Reformed.Data$Y[Reformed.Data$Z==1]),4),
        round(mean(Reformed.Data$nc[Reformed.Data$Zc==1]),4),round(mean(Reformed.Data$X[Reformed.Data$Zc==1,2]),4),round(mean(Reformed.Data$X[Reformed.Data$Zc==1,3]),4),round(mean(Reformed.Data$X[Reformed.Data$Zc==1,5]),4),
        round(mean(Reformed.Data$X[Reformed.Data$Zc==1,6]),4),round(mean(Reformed.Data$X[Reformed.Data$Zc==1,7]),4),round(mean(Reformed.Data$X[Reformed.Data$Zc==1,9]),4),round(mean(Reformed.Data$X[Reformed.Data$Zc==1,13]),4))
c3 <- c("Control",sum(1-Reformed.Data$Zc),sum(1-Reformed.Data$Z),sum(Reformed.Data$A[Reformed.Data$Z==0]),round(mean(Reformed.Data$Y[Reformed.Data$Z==0]),4),
        round(mean(Reformed.Data$nc[Reformed.Data$Zc==0]),4),round(mean(Reformed.Data$X[Reformed.Data$Zc==0,2]),4),round(mean(Reformed.Data$X[Reformed.Data$Zc==0,3]),4),round(mean(Reformed.Data$X[Reformed.Data$Zc==0,5]),4),
        round(mean(Reformed.Data$X[Reformed.Data$Zc==0,6]),4),round(mean(Reformed.Data$X[Reformed.Data$Zc==0,7]),4),round(mean(Reformed.Data$X[Reformed.Data$Zc==0,9]),4),round(mean(Reformed.Data$X[Reformed.Data$Zc==0,13]),4))
c4 <- c("t-statistic comparing treated and control groups","-","-","-","-",
        round(t.test(Reformed.Data$nc[Reformed.Data$Zc==0],Reformed.Data$nc[Reformed.Data$Zc==1])$statistic,4),
        round(t.test(Reformed.Data$X[Reformed.Data$Zc==0,2],Reformed.Data$X[Reformed.Data$Zc==1,2])$statistic,4),
        round(t.test(Reformed.Data$X[Reformed.Data$Zc==0,3],Reformed.Data$X[Reformed.Data$Zc==1,3])$statistic,4),
        round(t.test(Reformed.Data$X[Reformed.Data$Zc==0,5],Reformed.Data$X[Reformed.Data$Zc==1,5])$statistic,4),
        round(t.test(Reformed.Data$X[Reformed.Data$Zc==0,6],Reformed.Data$X[Reformed.Data$Zc==1,6])$statistic,4),
        round(t.test(Reformed.Data$X[Reformed.Data$Zc==0,7],Reformed.Data$X[Reformed.Data$Zc==1,7])$statistic,4),
        round(t.test(Reformed.Data$X[Reformed.Data$Zc==0,9],Reformed.Data$X[Reformed.Data$Zc==1,9])$statistic,4),
        round(t.test(Reformed.Data$X[Reformed.Data$Zc==0,13],Reformed.Data$X[Reformed.Data$Zc==1,13])$statistic,4))

CC <- cbind(c1," & ",c2," & ",c3," & ",c4," \\\\ \\hline")

print(data.frame(CC),row.names=FALSE)

####################################################
# Estimation
####################################################
############################
# ITT effect analysis
############################

ITTest <- ITT(Reformed.Data, Input.Type="Data")

############################
# Heterogeneous ITT effect analysis
############################

HTEest <- HTE(Reformed.Data, Xvar=1:13, constant = TRUE, Input.Type="Data")

############################
# Compliance grooup effect analysis
############################

set.seed(12321)

# Sharp bound under (C1)
Bound1 <- SharpBound(Reformed.Data, paraC=1:13, method="Linear", CIcalc=TRUE, SSsize=1000, level=0.95, Input.Type="Data")

# Sharp bound under (C2)
Bound2 <- SharpBound(Reformed.Data, paraC=1:13, method="Logistic", CIcalc=TRUE, SSsize=1000, level=0.95, Input.Type="Data")

# Bound using (Unadjusted), [male]x[cluster size] (Adjusted); approach from Long and Hudgens (2013)
Bound3 <- LongHudgens(Reformed.Data, paraC=c(2,10,11,12), CIcalc=TRUE, SSsize=1000, level=0.95, Input.Type="Data")

# Intersecting (C2) and (Adjusted)
Bound4 <- Bound.Intersect(Bound2,Bound3, level=0.95)

####################################################
# Report
####################################################

# Overall/Heterogeneous ITT effects

CI <- function(p,t){
  cbind(p-qnorm(0.975)*t,p+qnorm(0.975)*t)
}

RR <- cbind(c(ITTest[1],HTEest$Estimate[1,1:13]),
            c(ITTest[2],HTEest$Estimate[2,1:13]),
            CI(c(ITTest[1],HTEest$Estimate[1,1:13]),
               c(ITTest[2],HTEest$Estimate[2,1:13])),
            c(ITTest[3],HTEest$Estimate[3,1:13])^2,
            1-pchisq(c(ITTest[3],HTEest$Estimate[3,1:13])^2,1))
colnames(RR) <- c("Estimate","SE","LB","UB","chi1-statistic","p-value")
rownames(RR) <- c("Overall ITT",colnames(HTEest$Estimate))

RR

# Bounds

RR <- rbind(as.numeric(c( Bound1$Bound[1:2],
                          Bound2$Bound[1:2],
                          Bound3$Bound.Adj[1:2],
                          Bound4$Bound[1:2])),
            as.numeric(c( Bound1$BootCIBound[1:2],
                          Bound2$BootCIBound[1:2],
                          Bound3$BootCIBound.Adj[1:2],
                          Bound4$BootCIBound[1:2])),
            as.numeric(c( Bound1$Bound[3:4],
                          Bound2$Bound[3:4],
                          Bound3$Bound.Adj[3:4],
                          Bound4$Bound[3:4])),
            as.numeric(c( Bound1$BootCIBound[3:4],
                          Bound2$BootCIBound[3:4],
                          Bound3$BootCIBound.Adj[3:4],
                          Bound4$BootCIBound[3:4])),
            as.numeric(c( Bound1$Bound[5:6],
                          Bound2$Bound[5:6],
                          Bound3$Bound.Adj[5:6],
                          Bound4$Bound[5:6])),
            as.numeric(c( Bound1$BootCIBound[5:6],
                          Bound2$BootCIBound[5:6],
                          Bound3$BootCIBound.Adj[5:6],
                          Bound4$BootCIBound[5:6])))
colnames(RR) <- c("(C1)LB","(C1)UB","(C2)LB","(C2)UB","(C3)LB","(C3)UB","(C4)LB","(C4)UB")
rownames(RR) <- c("NT_Est","NT_CI","AT_Est","AT_CI","CO_Est","CO_CI")

RR

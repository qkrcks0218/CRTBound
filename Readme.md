# CRTBound

This Github repository contains CRTBound R package that implements the methodologies in [Park and Kang (2021)](https://www.tandfonline.com/doi/full/10.1080/01621459.2021.1983437 "IVNet"). We consider a cluster randomized trial (CRT) when both partial interference and noncompliance are present. In particular, the package focuses on the three types of the causal estimands:
* Overall intent-to-Treat (ITT) Effect 
* Heteogeneous ITT effect
* Compliance group ITT effects (i.e. network effects) arising from CRTs under the presence of interference and noncompliance

To infer the overall, heterogeneous ITT effects, we propose a modest extension of a nonparametric, regression-esque method based on Ding et al. (2019). In CRTBound package, the inference of the overall, heterogeneous ITT effects can be done by running `ITT(Data)` and `HTE(Data)` functions, respecitvely.

To infer the network effects, we use the classification-based bounding methods via linear programming. In CRTBound package, this method is implemented in `SharpBound(Data)` function. Additionally, we propose extensions of bounds by Grilli and Mealli (2008) and Long and Hudgens (2013) to CRT settings with interference. In CRTBound package, these methods are implemented in `LongHudgens` function. Lastly, by intersecting the bound estimates from the above two approaches, narrower bounds can be obtained. In CRTBound package, these methods are implemented in `Bound.Intersect` function. 


This package is currently in beta.

## Installation

To install CRTBound package in R, run the commands below:

```{r}
library(devtools)
install_github("qkrcks0218/CRTBound")
```

## Example Usage

Here are some examples of analyzing a simulated CRT:

```{r}
library(CRTBound)
#################################
# Generate Population
#################################

J <- 100 ; m <- 60
nc <- sample(c(3,4,5),J,replace=T)
N <- sum(nc)
C <- rep(1:J,nc)
n <- rep(nc,nc)
X1 <- rnorm(N)
X2 <- rbinom(N,1,0.3)
X3 <- rep(rnorm(J),nc)
A0 <- rbinom(N,1,logistic(-2+2*X3))
A1 <- apply(cbind(A0,rbinom(N,1,logistic(-2+3*X1+3*X2+2*X3))),1,max)
OtherA1 <- (rep(aggregate(A1~C,FUN="sum")[,2],nc)-A1)/(n-1)
NT <- which(A1==0 & A0==0)
AT <- which(A1==1 & A0==1)
CO <- which(A1==1 & A0==0)
Y0 <- Y1 <- rep(0,N)
for(jj in NT){
    Y0[jj] <- rbinom(1,1,logistic(-2+2*OtherA1[jj]))
    Y1[jj] <- max(Y0[jj],rbinom(1,1,logistic(-2+2*OtherA1[jj])))
}
for(jj in AT){
    Y0[jj] <- rbinom(1,1,logistic(2+X1[jj]+X2[jj]))
    Y1[jj] <- max(Y0[jj],rbinom(1,1,logistic(2+X1[jj]+X2[jj])))
}
for(jj in CO){
    Y0[jj] <- rbinom(1,1,logistic(-2+2*OtherA1[jj]))
    Y1[jj] <- max(Y0[jj],rbinom(1,1,logistic(2+X1[jj]+X2[jj])))
}
X <- cbind(1,X1,X2,X3,n)

Zc <- rep(0,J)
Zc[sort(sample(J,m))] <- 1
Z <- rep(Zc,nc)
A <- Z*A1 + (1-Z)*A0
Y <- Z*Y1 + (1-Z)*Y0

#################################
# Reform the Data
#################################

Reformed.Data <- Data.Reform(Y,Z,A,C,X,seed=1)

#################################
# Bounds
#################################

Bound1 <- SharpBound(Reformed.Data,paraC=c(1,2,3,4),method="Logistic",CIcalc=TRUE,SSsize=100,level=0.95,seed=1)
Bound2 <- LongHudgens(Reformed.Data,paraC=c(3),CIcalc=TRUE,SSsize=100,level=0.95,seed=1)
Bound3 <- Bound.Intersect(Bound1,Bound2,level=0.95)
```







## References

Grilli Leonardo & Mealli Fabrizia (2008) **Nonparametric Bounds on the Causal Effect of University Studies on Job Opportunities Using Principal Stratification**, _Journal of Educational and Behavioral Statistics_, 33(1), 111-130 [[link](https://journals.sagepub.com/doi/abs/10.3102/1076998607302627 "GM")]

Dustin M. Long & Michael G. Hudgens (2013) **Sharpening Bounds on Principal Effects with Covariates**, _Biometrics_, 69(4), 812-819 [[link](https://onlinelibrary.wiley.com/doi/10.1111/biom.12103 "LH")]

Peng Ding, Avi Feller & Luke Miratrix (2019) **Decomposing Treatment Effect Variation**, _Journal of the American Statistical Association_, 114(525), 304-317 [[link](https://www.tandfonline.com/doi/full/10.1080/01621459.2017.1407322 "Ding")]

Chan Park & Hyunseung Kang (2021+) **Assumption-Lean Analysis of Cluster Randomized Trials in Infectious Diseases for Intent-to-Treat Effects and Network Effects**, _Journal of the American Statistical Association_, In press [[link](https://www.tandfonline.com/doi/full/10.1080/01621459.2021.1983437 "IVNet")]
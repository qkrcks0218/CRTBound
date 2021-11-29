\name{SharpBound}
\alias{SharpBound}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Inference of the Compliance Group ITT Effects Using Classifier-based Approach}
\description{
Analyze the data (obtained from \code{\link{Data.Reform}} function) to infer the compliance group ITT effect using classifier-based approach; see Section 4 of Park and Kang (2021) for details.
}
\usage{
SharpBound(Reformed.Data, paraC=NULL, method="Linear", CIcalc=FALSE, SSsize=1000, level=0.95, seed=1){
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Reformed.Data}{A list obtained from \code{\link{Data.Reform}} function.}
  \item{paraC}{\emph{(Optional)} Indices of the pre-treatment covariate matrix Reformed.Data$X that are used in inference of the compliance group ITT effect. The default option is to use all pre-treatment covariates.}
  \item{method}{\emph{(Optional)} The type of the transformation function. "Linear" and "Logistic" are available where "Linear" and "Logistic" use the linear and logistic transformations. The default option is "Linear"; see Section 4 of Park and Kang (2021) for details.}
  \item{CIcalc}{\emph{(Optional)} A logical value whether the bootstrap-based confidence intervals are obtained. The dafault option is TRUE.}
  \item{SSsize}{\emph{(Optional)} The total number of the bootstrap resamples used for the construction of the confidence intervals. The dafault option is 1000.}
  \item{level}{\emph{(Optional)} The confidence level of the confidence interval. The default option is 0.95.}
  \item{seed}{\emph{(Optional)} A random seed to generate the bootstrap resamples used for the construction of the confidence intervals.}
}
}
\value{
\item{Bound}{Bound estimates of the compliance group ITT effects.}
\item{Violation}{Summary whether the linear program was feasible.}
\item{Resample.Bound}{Bootstrap bound estimates of the compliance group ITT effects.}
\item{BootCIBound}{Confidence intervals of the bounds of the compliance group ITT effects.}
}
\references{
Chan Park & Hyunseung Kang (2021+) \strong{Assumption-Lean Analysis of Cluster Randomized Trials in Infectious Diseases for Intent-to-Treat Effects and Network Effects}, \emph{Journal of the American Statistical Association} \href{https://www.tandfonline.com/doi/full/10.1080/01621459.2021.1983437}{[Link]}
}
\examples{
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
# Sharp Bounds
#################################

Bound1 <- SharpBound(Reformed.Data,paraC=c(1,2,3,4),method="Logistic",CIcalc=TRUE,SSsize=100,level=0.95,seed=1)
Bound2 <- LongHudgens(Reformed.Data,paraC=c(3),CIcalc=TRUE,SSsize=100,level=0.95,seed=1)
Bound3 <- Bound.Intersect(Bound1,Bound2,level=0.95)
}

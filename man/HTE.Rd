\name{HTE}
\alias{HTE}
\title{Inference of the Heterogeneous ITT Effect
}
\description{
Analyze the data (obtained from \code{\link{Data.Reform}} or \code{\link{Simulation.Reform}} functions) to infer the heterogeneous ITT effect; see Section 3 of Park and Kang (2021) for details.
}
\usage{HTE(Data, Xvar=NULL, constant=TRUE, Input.Type="Data")}
\arguments{
  \item{Data}{A list obtained from \code{\link{Data.Reform}} and \code{\link{Simulation.Reform}} functions.}
  \item{Xvar}{\emph{(Optional)} Indices of the pre-treatment covariate matrix Data$X that are used in inference of the heterogeneous ITT effect. The default option is to use all pre-treatment covariates.}
  \item{constant}{\emph{(Optional)} A logical value whether the constant intercept is included in the analysis.}
  \item{Input.Type}{\emph{(Optional)} A variable whether Data is from \code{\link{Data.Reform}} (Input.Type="Data") or \code{\link{Simulation.Reform}} (Input.Type="Sim") function. The default The default option is "Data".}
}
\value{
\item{Estimate}{A matrix containing the point estimates of the heterogeneity parameters, the standard errors, and the z-statistics. When Input.Type="Sim", the true heterogeneity parameters and the true standard errors are available.}
\item{NonConst.Chi.Sq.Statistic}{A chi-square statistic obtained from testing the overall non-intercept null hypothesis.}
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
Simulated.Data <- Simulation.Reform(Y0,Y1,Z,A0,A1,C,X,seed=1)

#################################
# HT Effect
#################################

HT.Data <- HTE(Reformed.Data,Xvar=1:5,Input.Type="Data")
HT.Sim.Data <- HTE(Simulated.Data,Xvar=1:5,Input.Type="Sim")
}

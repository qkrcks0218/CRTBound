\name{Bound.Intersect}
\alias{Bound.Intersect}
\title{Sharpening Bounds of the Compliance Group ITT Effects Using Multiple Approaches}
\description{
Sharpening bounds of the compliance group ITT effects using the approaches in Grilli and Mealli (2008), Long and Hudgens (2013), and Park and Kang (2021).}
\usage{Bound.Intersect(SharpBound.Object, LongHudgens.Object, level=0.95)}
\arguments{
  \item{SharpBound.Object}{An object obtained from \code{\link{SharpBound}} function.}
  \item{LongHudgens.Object}{An object obtained from \code{\link{LongHudgens}} function.}
  \item{level}{\emph{(Optional)} The confidence level of the confidence interval. The default option is 0.95.}
}
\value{
\item{Bound}{Bound estimates of the compliance group ITT effects. When SharpBound.Object and LongHudgegns.Object have the true bounds, the true bounds are availab.e}
\item{BootCIBound}{Confidence intervals of the bounds of the compliance group ITT effects.}
}
\references{
Grilli Leonardo & Mealli Fabrizia (2008) \strong{Nonparametric Bounds on the Causal Effect of University Studies on Job Opportunities Using Principal Stratification}. \emph{Journal of Educational and Behavioral Statistics} \href{https://journals.sagepub.com/doi/abs/10.3102/1076998607302627}{[Link]}

Dustin M. Long & Michael G. Hudgens (2013) \strong{Sharpening Bounds on Principal Effects with Covariates}. \emph{Biometrics} \href{https://onlinelibrary.wiley.com/doi/10.1111/biom.12103}{[Link]}

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
# Bounds
#################################

Bound1 <- SharpBound(Reformed.Data,paraC=1:5,method="Logistic",CIcalc=TRUE,SSsize=100,level=0.95,seed=1,Input.Type="Data")
Bound2 <- LongHudgens(Reformed.Data,paraC=c(3),CIcalc=TRUE,SSsize=100,level=0.95,seed=1,Input.Type="Data")
Bound3 <- Bound.Intersect(Bound1,Bound2,level=0.95,Input.Type="Data")

Bound1.Sim <- SharpBound(Simulated.Data,paraC=1:5,method="Logistic",CIcalc=TRUE,SSsize=100,level=0.95,seed=1,Input.Type="Sim")
Bound2.Sim <- LongHudgens(Simulated.Data,paraC=c(3),CIcalc=TRUE,SSsize=100,level=0.95,seed=1,Input.Type="Sim")
Bound3.Sim <- Bound.Intersect(Bound1.Sim,Bound2.Sim,level=0.95,Input.Type="Sim")
}

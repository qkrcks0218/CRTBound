\name{Data.Reform}
\alias{Data.Reform}
\title{Reform the data
}
\description{
Reform the data to a compatible form for \code{\link{ITT}}, \code{\link{HTE}}, \code{\link{SharpBound}}, and \code{\link{LongHudgens}} functions.
}
\usage{Data.Reform(Y,Z,A,C,X,seed=1)}
\arguments{
  \item{Y}{\eqn{N-by-1} vector containing the outcome of individuals where \eqn{N} is the number of individuals in the population.}
  \item{Z}{\eqn{N-by-1} 0-1 vector containing the treatment assignment where Z=1 indicates an individual was assigned to treatment.}
  \item{A}{\eqn{N-by-1} 0-1 vector containing the actual treatment usage of individuals where A=1 indicates an individual actually used treatment.}
  \item{C}{\eqn{N-by-1} numeric vector containing cluster ID where individuals from the same cluster have the same C value.}
  \item{X}{\eqn{N-by-p} numeric matrix or data frame containing the pre-treatment covariates.}
  \item{seed}{\emph{(Optional)} A random seed to generate the random noise used in \code{\link{SharpBound}}; see Section 4.2 of Park and Kang (2021) for details.}
}
\value{
\item{Y}{\eqn{N-by-1} vector; same as input Y.}
\item{Z}{\eqn{N-by-1} vector; same as input Z.}
\item{A}{\eqn{N-by-1} vector; same as input A.}
\item{C}{\eqn{N-by-1} vector; same as input C.}
\item{X}{\eqn{N-by-p} matrix or data frame; same as input X.}
\item{N}{The number of individuals in the population.}
\item{J}{The number of clusters in the population.}
\item{Class}{\eqn{J-by-1} vector of cluster ID.}
\item{nc}{\eqn{J-by-1} vector of cluster sizes.}
\item{n}{\eqn{N-by-1} vector of cluster sizes and the same as \code{rep(nc,nc)}.}
\item{Zc}{\eqn{J-by-1} vector of treatment assignment and the simplified form of input Z.}
\item{m}{The number of clusters that were assigned to treatment.}
\item{error}{\eqn{N-by-3} vector of the random noises and are used in \code{\link{SharpBound}}; see Section 4.2 of Park and Kang (2021) for details.}
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

}

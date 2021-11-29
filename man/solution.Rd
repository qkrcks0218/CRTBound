\name{solution}
\alias{solution}
\title{Find the Calibration Value
}
\description{
Find the calibration value \emph{q} satisfying \eqn{T = \sum_{ji} I_{c,h} (L_{ji} - q)} from a grid-search.}
\usage{solution(L,T,h,c,thr)}
\arguments{
  \item{L}{A vector in the surrogate indicator function \eqn{\mathcal{I}_{c,h}}.}
  \item{T}{A number in the left hand side.}
  \item{h,c}{Parameters of the surrogate indicator function \eqn{\mathcal{I}_{c,h}}; see Section 4.2 of Park and Kang (2021) for details.}
  \item{thr}{A number used as the threshold in the grid-search. The true solution \emph{q}* satisfies \eqn{|q*-q|<thr}. }
}
\value{
  A solution \emph{q} to \eqn{T = \sum_{ji} \mathcal{I}_{c,h} (L_{ji} - q)}.
}
\references{
Chan Park & Hyunseung Kang (2021+) \strong{Assumption-Lean Analysis of Cluster Randomized Trials in Infectious Diseases for Intent-to-Treat Effects and Network Effects}, \emph{Journal of the American Statistical Association} \href{https://www.tandfonline.com/doi/full/10.1080/01621459.2021.1983437}{[Link]}
}
\examples{
L <- c(1.1, 1.8, 2.7, 4.1, 5.2,
       6.1, 7.4, 7.9, 8.8, 9.9)
T <- 5
c <- h <- 10^(-4)
thr <- 10^(-4)
solution(L,T,h,c,thr) # Returns 5.65
}

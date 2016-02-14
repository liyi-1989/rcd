\name{rcd-package}
\docType{package}
\title{
\packageTitle{rcd}
}
\description{
\packageDescription{rcd}
}
\details{

%The DESCRIPTION file:
%\packageDESCRIPTION{rcd}
%\packageIndices{rcd}

This package includes the implementation of the robust copula dependence (RCD). The main function rcd is used with dependence measure bwtween two or more random variables. 
}
\author{
\packageAuthor{rcd}

Maintainer: \packageMaintainer{rcd}
}
\references{
Yale Chang, Yi Li, Adam Ding, Jennifer Dy, A Robust-Equitable Copula Dependence Measure for Feature Selection , AISTATS 2016.
}

\keyword{ Copula, Dependence measure, Equitability }

\examples{
  n <- 1000
  x <- runif(n)
  y <- runif(n)
  rcdkde(x,y)
}

\name{multiplsrog_loading}
\alias{multiplsrog_loading}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiset partial least squares rank order of groups loading}
\description{
This function computes multiset partial least squares rank order of groups (PLS-ROG) loading from the result of the "multipls_rog" function.
}
\usage{
multiplsrog_loading(multiplsrog)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multiplsrog}{

  The following variables (P,T,Q and U) are included in the multiplsrog object.

  P : A list of matrix with Multiset PLS-ROG coefficients for the explanatory variables in each column for each dataset

  T : A list of matrix with Multiset PLS-ROG scores for the explanatory variables in each column for each dataset

  Q : A matrix with Multiset PLS-ROG coefficients for the response variable in each column

  U : A matrix with Multiset PLS-ROG scores for the response variable in each column

  tau : Matrix for strength parameter of the connection between omics datasets or between omics dataset and group information.

  }
}
\details{
The multiset PLS-ROG loading and the p-value by statistical hypothesis testing is added to the multiplsrog object returned by the "multipls_rog" function in this package.

}
\value{
The return value is a list object that contains the following elements:

Score : A list of matrix with Multiset PLS score for auxiliary variable in each column for each dataset

R : A list of matrix with Multiset PLS-ROG loadings in each column for each dataset

p.value : A list of matrix with p-value of multiset PLS-ROG loading in each column for each dataset

}
\references{
Yamamoto H. (2022) Multiset partial least squares with rank order of groups for integrating multi-omics data, bioRxiv.
}
\author{
Hiroyuki Yamamoto
}
\examples{
data(whhl)
X <- whhl$X
Y <- whhl$Y
D <- whhl$D
tau <- whhl$tau

multiplsrog <- multipls_rog(X,Y,tau,D)
# multiplsrog <- multipls_rog(X,Y,tau,D, kappa=0.999)

multiplsrog <- multiplsrog_loading(multiplsrog)
multiplsrog$loading$R
multiplsrog$loading$p.value
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }


\name{multipls_loading}
\alias{multipls_loading}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiset partial least squares loading}
\description{
This function computes multiset partial least squares (PLS) loading from the result of the "multipls_geigen" function.
}
\usage{
multipls_loading(multipls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multipls}{

  The following variables (P,T,Q and U) are included in the multipls object.

  P : A list of matrix with Multiset PLS coefficients for the explanatory variables in each column for each dataset

  T : A list of matrix with Multiset PLS scores for the explanatory variables in each column for each dataset

  Q : A matrix with Multiset PLS coefficients for the response variable in each column

  U : A matrix with Multiset PLS scores for the response variable in each column

  }
}
\details{
The multiset PLS loading and the p-value by statistical hypothesis testing is added to the multipls object returned by the "multipls_geigen" function in this package.

}
\value{
The return value is a list object that contains the following elements:

R : A list of matrix with Multiset PLS loadings in each column for each dataset

p.value : A list of matrix with p-value of multiset PLS loading in each column for each dataset

}
\references{
Yamamoto H. (2022) Multiset partial least squares with rank order of groups for integrating multi-omics data, bioRxiv.
}
\author{
Hiroyuki Yamamoto
}
\examples{
# COVID19 multiomics data
data(covid19)

X <- covid19$X
Y <- covid19$Y
tau <- covid19$tau

multipls <- multipls_geigen(X,Y,tau)

multipls <- multipls_loading(multipls)
multipls$loading$R
multipls$loading$p.value
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }


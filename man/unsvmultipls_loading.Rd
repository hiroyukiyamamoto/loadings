\name{unsvmultipls_loading}
\alias{unsvmultipls_loading}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Unsupervised multiset partial least squares loading}
\description{
This function computes unsupervised multiset partial least squares loading from the result of the "unsv_multipls" function.
}
\usage{
unsvmultipls_loading(unsvmultipls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{unsvmultipls}{

  The following variables (P,T) are included in the unsvmultipls object.

  P : A list of matrix with Multiset PLS coefficients in each column for each dataset

  T : A list of matrix with Multiset PLS scores in each column for each dataset

  }
}
\details{
The unsupervised multiset PLS loading and the p-value by statistical hypothesis testing is added to the unsvmultipls object returned by the "unsv_multipls" function in this package.
}
\value{
The return value is a list object that contains the following elements:

R : A list of matrix with Unsupervised multiset PLS loadings in each column for each dataset

p.value : A list of matrix with p-value of Unsupervised multiset PLS loading in each column for each dataset

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
tau <- rbind(c(0,1/2),c(1/2,0))

unsvmpls <- unsv_multipls(X,tau)

unsvmultipls <- unsvmultipls_loading(unsvmpls)
unsvmultipls$loading$R
unsvmultipls$loading$p.value

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }


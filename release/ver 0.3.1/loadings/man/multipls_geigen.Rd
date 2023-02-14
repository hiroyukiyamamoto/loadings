\name{multipls_geigen}
\alias{multipls_geigen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiset PLS : Multiset partial least squares}
\description{
This function performs Multiset partial least squares (Multiset PLS).
In this function, data matrix is automatically scaled to zero mean and unit variance (i.e. autoscaling) for each variables.
}
\usage{
multipls_geigen(X,Y,tau)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
	List of data matrix that include variables in each columns.
  }
  \item{Y}{
	Dummy matrix that include group information 0,1 in each columns.
  }
  \item{tau}{
	Matrix for strength parameter of the connection between omics datasets or between omics dataset and group information.
  }
}
\details{
Diagonal elements of matrix tau must be 0.
}
\value{
The return value is a list object that contains the following elements:

P : A list of matrix with Multiset PLS coefficients for the explanatory variables in each column for each dataset

T : A list of matrix with Multiset PLS scores for the explanatory variables in each column for each dataset

Q : A matrix with Multiset PLS coefficients for the response variable in each column

U : A matrix with Multiset PLS scores for the response variable in each column

tau : Matrix for strength parameter of the connection between omics datasets or between omics dataset and group information (same as input value).

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
tau <- whhl$tau

multipls <- multipls_geigen(X,Y,tau)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }


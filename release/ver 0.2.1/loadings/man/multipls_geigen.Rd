\name{multipls_geigen}
\alias{multipls_geigen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiset Partial Least Squares}
\description{
This function performs Multiset Partial Least Squares (Multiset PLS).
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
	Parameter matrix for strength of the connection between omics datasets or between omics dataset and group information.
  }
}
\details{
diagonal elements are set to 0. The kappa represents the degree of smoothing. The value of kappa increases, the strength of the smoothing increases.
}
\value{
The return value is a list object that contains the following elements:

P : A matrix with OS-PC loading in each column

T : A matrix with OS-PC score in each column

MT : A matrix with averaging OS-PC score for repeated data in each column
(If not for repeated data, the matrix MT equals to the matrix T)

Q : A matrix with OS-PC loading for auxiliary variable in each column

U : A matrix with OS-PC score for auxiliary variable in each column

}
\references{
Yamamoto H., Nakayama Y., Tsugawa H. (2021) OS-PCA: Orthogonal Smoothed Principal Component Analysis Applied to Metabolome Data, Metabolites, 11(3):149.
}
\author{
Hiroyuki Yamamoto
}
\examples{
# metabolic turnover data
data(covid19)

X <- covid19$X
Y <- covid19$Y
tau <- covid19$tau

multipls <- multipls_geigen(X,Y,tau)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }


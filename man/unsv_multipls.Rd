\name{unsv_multipls}
\alias{unsv_multipls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Unsupervised multiset partial least squares}
\description{
This function performs Unsupervised multiset partial least squares.
In this function, data matrix is automatically scaled to zero mean and unit variance (i.e. autoscaling) for each variables.

}
\usage{
unsv_multipls(X,tau)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
	List of data matrix that include variables in each columns.
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

P : A list of matrix with Multiset PLS coefficients in each column for each dataset

T : A list of matrix with Multiset PLS scores in each column for each dataset

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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }


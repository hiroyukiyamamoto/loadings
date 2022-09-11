\name{multipls_rog}
\alias{multipls_rog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiset PLS-ROG : Multiset partial least squares with rank order of groups}
\description{
This function performs multiset partial least squares with rank order of groups (Multiset PLS-ROG).
In this function, data matrix is automatically scaled to zero mean and unit variance (i.e. autoscaling) for each variables.
}
\usage{
multipls_rog(X,Y,tau,D,kappa)
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
  \item{D}{
	Differential matrix.
  }
  \item{kappa}{
	The smoothing parameter (default : 0.999).
  }
}
\details{
Diagonal elements of matrix tau must be 0.
}
\value{
The return value is a list object that contains the following elements:

P : A list of matrix with Multiset PLS-ROG coefficients for the explanatory variables in each column for each dataset

T : A list of matrix with Multiset PLS-ROG scores for the explanatory variables in each column for each dataset

Q : A matrix with Multiset PLS-ROG coefficients for the response variable in each column

U : A matrix with Multiset PLS-ROG scores for the response variable in each column

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
D <- covid19$D

multiplsrog <- multipls_rog(X,Y,tau,D)
# multiplsrog <- multipls_rog(X,Y,tau,D, kappa=0.999)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }


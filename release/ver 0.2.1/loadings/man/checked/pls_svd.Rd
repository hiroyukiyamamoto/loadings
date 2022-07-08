\name{pls_svd}
\alias{pls_svd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Partial least squares}
\description{
This function  performs partial least squares. In this function, data matrix for explanatory variable is automatically scaled to zero mean and unit variance (i.e. autoscaling) for each variables.
}
\usage{
pls_svd(X,Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
	Data matrix of explanatory variables that include variables in each columns.
  }
  \item{Y}{
	Dummy matrix that include group information 0,1 in each columns.
  }
}
\details{
This function is wrapper function of "pls_rog" fucntion that the smoothing parameter kappa=0.
}
\value{
The return value is a list object that contains the following elements:

P : A matrix with PLS loading for explanatory variable in each column

T : A matrix with PLS score for explanatory variable in each column

Q : A matrix with PLS loading for response variable in each column

U : A matrix with PLS score for response variable in each column
}
\references{
Barker, M. and Rayens, W. (2003) Partial Least Squares for Discrimination. Journal of Chemometrics, 17, 166-173.
}
\author{
Hiroyuki Yamamoto
}
\examples{
data(whhl)
X <- whhl$X
Y <- whhl$Y

pls <- pls_svd(X,Y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }


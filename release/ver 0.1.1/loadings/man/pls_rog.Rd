\name{pls_rog}
\alias{pls_rog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Partial least squares rank order of groups}
\description{
This function performs partial least squares rank order of groups (PLS-ROG).
In this function, data matrix is automatically scaled to zero mean and unit variance (i.e. autoscaling) for each variables.
}
\usage{
pls_rog(X,Y,kappa)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
	Data matrix of explanatory variables that include variables in each columns.
  }
  \item{Y}{
	Dummy matrix that include group information 0,1 in each culumns.
  }
  \item{kappa}{
	The smoothing parameter (default : kappa = 0.999).
  }
}
\details{
The kappa represents the degree of smoothing. The value of kappa increases, the strength of the smoothing increases.
}
\value{
P : A matrix with PLS loading for explanatory variable in each column

T : A matrix with PLS score for explanatory variable in each column

Q : A matrix with PLS loading for response variable in each column

U : A matrix with PLS score for response variable in each column
}
\references{
Yamamoto, H. (2017) PLS-ROG: Partial least squares with rank order of groups., Journal of Chemometrics, 31(3) (2017) e2883.
}
\author{
Hiroyuki Yamamoto
}
\examples{
data(whhl)

pls <- pls_rog(X,Y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pca/msea }


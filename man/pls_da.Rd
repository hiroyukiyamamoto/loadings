\name{pls_da}
\alias{pls_da}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Partial least squares discriminant analysis (PLS-DA)}
\description{
This function performs partial least squares discriminant analysis (PLS-DA). In this function, data matrix for explanatory variable is automatically scaled to zero mean and unit variance (i.e. autoscaling) for each variables.
}
\usage{
pls_da(X,Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
	Data matrix of explanatory variables that include variables in each columns.
  }
  \item{Y}{
	Dummy matrix that include group information 0,1 in each columns.
  }
  \item{k}{
	Number of components.
  }
}
\details{
This function calculates PLS-DA. For PLS, use the 'pls_svd' function for PLS.
}
\value{
The return value is a list object that contains the following elements:

Wx : A matrix with PLS-DA weight for explanatory variable in each column

P: A matrix containing the PLS-DA loadings for each explanatory variable in the columns, before transformation.

T : A matrix with PLS-DA score for explanatory variable in each column

Wy : A matrix with PLS-DA weight for explanatory variable in each column

S : A matrix with PLS score for response variable in each column

}
\references{
Yamamoto, H. et al., Dimensionality reduction for metabolome data using PCA, PLS, OPLS, and RFDA with differential penalties to latent variables", Chemom. Intell. Lab. Syst., 98 (2009)
}
\author{
Hiroyuki Yamamoto
}
\examples{
data(whhl)
X <- whhl$X$liver
Y <- whhl$Y

plsda <- pls_da(X,Y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

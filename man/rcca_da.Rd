\name{rcca_da}
\alias{rcca_da}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Regularized canonical correlation analysis for discriminant analysis (RCCA-DA)}
\description{
This function performs regularized canonical correlation analysis for discriminant analysis (RCCA-DA). In this function, data matrix for explanatory variable is automatically scaled to zero mean and unit variance (i.e. autoscaling) for each variables.
}
\usage{
rcca_da(X,Y,lambda,k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
	Data matrix of explanatory variables that include variables in each columns.
  }
  \item{Y}{
	Dummy matrix that include group information 0,1 in each columns.
  }
  \item{lambda}{
	The regularized parameter has a value in the range [0, 1), meaning it can be 0 but is less than 1."
  }
  \item{k}{
	Number of components.
  }
}
\details{
RCCA-DA is equivalent to Regularized Fisher discriminant analysis, theoretically.
}
\value{
The return value is a list object that contains the following elements:

Wx : A matrix with RCCA-DA weight for explanatory variable in each column

P: A matrix containing the RCCA-DA loadings for each explanatory variable in the columns, before transformation.

T : A matrix with RCCA-DA score for explanatory variable in each column

Wy : A matrix with RCCA-DA weight for explanatory variable in each column

S : A matrix with RCCA-DA score for response variable in each column

}
\references{
Yamamoto, H. et al., Canonical correlation analysis for multivariate regression and its application to metabolic fingerprinting", Biochem. Eng. Journal, 40 (2008) 199-204.
}
\author{
Hiroyuki Yamamoto
}
\examples{
data(whhl)
X <- whhl$X$liver
Y <- whhl$Y

rcca <- rcca_da(X,Y,0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

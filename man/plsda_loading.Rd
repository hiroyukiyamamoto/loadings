\name{plsda_loading}
\alias{plsda_loading}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Partial least squares discriminant analysis (PLS-DA) loading}
\description{
This function computes partial least squares discriminant analysis (PLS-DA) loading from the result of the 'pls_da' function.
}
\usage{
plsda_loading(plsda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plsrog}{

  The following variables (Wx,P,T,Wy and S) are included in the plsda object.

  Wx : A matrix with PLS-DA weight for explanatory variable in each column

P: A matrix containing the PLS-DA loadings for each explanatory variable in the columns, before transformation.

T : A matrix with PLS-DA score for explanatory variable in each column

Wy : A matrix with PLS-DA weight for explanatory variable in each column

S : A matrix with PLS score for response variable in each column

  }
}
\details{
The PLS-DA loading and the p-value by statistical hypothesis testing is added to the pls object returned by the 'pls_da' function in this package.
}
\value{

  The return value is a list object that contains the following elements:

  R : PLS-DA loading (Correlation coefficient between PLS-DA score of explanatory variable and each variables in data matrix.)

  p.value : p-value of PLS-DA loading

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

plsda <- plsda_loading(plsda)

plsda$loading$R
plsda$loading$p.value

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }


\name{plsrog_loading}
\alias{plsrog_loading}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Partial least squares rank order of groups loading}
\description{
This function computes partial least squares rank rder of groups (PLS-ROG) loading from the result of the "pls_rog" function.
PLS-ROG loading can also be computed from the result of the "pls_rog" function.
}
\usage{
pls_loading(pls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pls}{

  The following variables (P,T,Q and U) are included in the pls object.

  P : A matrix with PLS-ROG loading for explanatory variable in each column

  T : A matrix with PLS-ROG score for explanatory variable in each column

  Q : A matrix with PLS-ROG loading for response variable in each column

  U : A matrix with PLS-ROG score for response variable in each column

  }
}
\details{
The PLS-ROG loading and the p-value by statistical hypothesis testing is added to the pls object returned by the "pls_rog" function in this package.
}
\value{

  The return value is a list object that contains the following elements:

  R : PLS-ROG loading (Correlation coefficient between PL-ROG score of response variable and each variables in data matrix.)

  p.value : p-value of PLS-ROG loading

}
\references{
Yamamoto, H. (2017) PLS-ROG: Partial least squares with rank order of groups., Journal of Chemometrics, 31(3) (2017) e2883.
}
\author{
Hiroyuki Yamamoto
}
\examples{
data(whhl)
X <- whhl$X$liver
Y <- whhl$Y
D <- whhl$D

plsrog <- pls_rog(X,Y,D)

plsrog <- pls_loading(plsrog)
plrogs$loading$R
plsrog$loading$p.value

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }


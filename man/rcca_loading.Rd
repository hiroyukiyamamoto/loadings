\name{rcca_loading}
\alias{rcca_loading}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Regularized canonical correlation analysis (RCCA-DA) loading}
\description{
This function computes regularized canonical correlation analysis for discriminant analysis (RCCA-DA) loading from the result of the 'rcca_da' function.
}
\usage{
rcca_loading(rcca)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rcca}{

  The following variables (Wx,P,T,Wy and S) are included in the rcca object.

  Wx : A matrix with RCCA-DA weight for explanatory variable in each column

P: A matrix containing the RCCA-DA loadings for each explanatory variable in the columns, before transformation.

T : A matrix with RCCA-DA score for explanatory variable in each column

Wy : A matrix with RCCA-DA weight for explanatory variable in each column

S : A matrix with RCCA-DA score for response variable in each column

  }
}
\details{
The RCCA-DA loading and the p-value by statistical hypothesis testing is added to the pls object returned by the 'rcca_da' function in this package.
}
\value{

  The return value is a list object that contains the following elements:

  R : RCCA-DA loading (Correlation coefficient between PLS-DA score of explanatory variable and each variables in data matrix.)

  p.value : p-value of RCCA-DA loading

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
rcca <- rcca_da(X,Y,0.5,2)

rcca <- rcca_loading(rcca)

rcca$loading$R
rcca$loading$p.value

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }


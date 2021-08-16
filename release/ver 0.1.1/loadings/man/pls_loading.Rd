\name{pls_loading}
\alias{pls_loading}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Partial least squares loading}
\description{
This function computes partial least squares (PLS) loading from the result of the "pls_svd" function.
PLS loading can also be computed from the result of the "pls_eigen" function in chemometrics package.
}
\usage{
pls_loading(pls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pls}{
P : A matrix with PLS loading for explanatory variable in each column

T : A matrix with PLS score for explanatory variable in each column

Q : A matrix with PLS loading for response variable in each column

U : A matrix with PLS score for response variable in each column

}
}
\details{
The PLS loading and the p-value by statistical hypothesis testing is added to the pls object returned by the "pls_svd" function in this pacjage or "pls_eigen" function in chemometrics package.
}
\value{
loadings$R : PLS loading for explanatory variable

loadings$p.value : p-value of PLS loading for explanatory variable
}
\references{
Yamamoto, H. (2017) PLS-ROG: Partial least squares with rank order of groups., Journal of Chemometrics, 31(3) (2017) e2883.
}
\author{
Hiroyuki Yamamoto
}
\examples{
data(whhl)

pls <- pls_svd(X,Y)

pls <- pls_loading(pls)
pls$loading$R
pls$loading$p.value

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pca/msea }

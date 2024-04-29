\name{onekpca_loading}
\alias{onekpca_loading}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{One-sided kernel principal component loading}
\description{
This function computes one-sided kernel principal component (one-sided KPC) loading from the result of the "one_kpca" function.
}
\usage{
onekpca_loading(onekpca)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{onekpca}{

  The following variables (P,T and U) are included in the onekpca object.


  P : A matrix with one-sided KPCA loading in each column

  T : A matrix with one-sided KPCA score for linear side in each column

  U : A matrix with one-sied KPCA score for nonlinear side in each column

  }
}
\details{
The one-sided KPC loading and the p-value by statistical hypothesis testing is added to the onekpca object returned by the "one_kpca" function.
}
\value{

  The return value is a list object that contains the following elements:

  R : one-sided KPC loading (Correlation coefficient between one-sided KPC score of nonlinear variable (kernel matrix) and each variables in data matrix.)

  p.value : p-value of one-sided KPC loading

}
\references{
Yamamoto H. (2023) One-sided Kernel Principal Component Analysis, Jxiv, <doi:10.51094/jxiv.262>.
}
\author{
Hiroyuki Yamamoto
}
\examples{
# COVID-19 metabolome data
data(covid19)

X <- covid19$X$metabolomics
K <- covid19$K$metabolomics

onekpca <- one_kpca(X,K)

onekpca <- onekpca_loading(onekpca)

onekpca$loading$R
onekpca$loading$p.value

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }


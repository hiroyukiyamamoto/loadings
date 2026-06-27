\name{one_kpca}
\alias{one_kpca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{One-sided kernel principal component analysis}
\description{
This function performs one-sided kernel principal component analysis (one-sided KPCA).
In this function, data matrix is automatically scaled to zero mean and unit variance (i.e. autoscaling) for each variables.
}
\usage{
one_kpca(X,K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
	Data matrix that include variables in each columns.
  }
  \item{K}{
	Kernel matrix computed from the data matrix X.
  }
}
\details{
The kernel matrix K, which is the argument of the one_kpca function, must be centered.
}
\value{
The return value is a list object that contains the following elements:

P : A matrix with one-sided KPCA loading in each column

T : A matrix with one-sided KPCA score for linear side in each column

U : A matrix with one-sied KPCA score for nonlinear side in each column

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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }


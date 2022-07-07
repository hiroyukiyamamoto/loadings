\name{pca_loading}
\alias{pca_loading}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Principal component loading}
\description{
This function computes principal component (PC) loading from the result of the "prcomp" function.
(The "princomp" function is not supported. For "princomp" function, the "loadings" function in stats package should be used.)
In this function, data matrix is should be scaled to zero mean and unit variance (i.e. autoscaling) for each variables.
}
\usage{
pca_loading(pca)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pca}{
	 The return object of the "prcomp" function.
  }
}
\details{
The PC loading and the p-value by statistical hypothesis testing is added to the pca object returned by the "prcomp" function.
}
\value{

  The return value is a list object that contains the following elements:

  R : PC loading (Correlation coefficient between PC score and each variables in data matrix.)

  p.value : p-value of PC loading

}
\references{
Yamamoto H., Fujimori T., Sato H., Ishikawa G., Kami K., Ohashi Y. (2014). "Statistical hypothesis testing of factor loading in principal component analysis and its application to metabolite set enrichment analysis". BMC Bioinformatics, (2014) 15(1):51.
}
\author{
Hiroyuki Yamamoto
}
\examples{
data(fasting)
X <- fasting$X

pca <- prcomp(X, scale=TRUE)

pca <- pca_loading(pca)
pca$loading$R # PC loading
pca$loading$p.value # p-value
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ principal component analysis / loading }


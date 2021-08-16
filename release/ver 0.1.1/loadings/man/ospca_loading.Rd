\name{ospca_loading}
\alias{ospca_loading}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Orthogonal smoothed principal component loading}
\description{
This function computes orthogonals smoothed principal component (OS-PC) loading from the result of the "os_pca" function.
In this function, data matrix is should be scaled to zero mean and unit variance (i.e. autoscaling) for each variables.
}
\usage{
ospca_loading(ospca)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ospca}{
	 The return object of the "os_pca" function.
}
}
\details{
The OS-PC loading and the p-value by statistical hypothesis testing is added to the ospca object returned by the "os-pca" function.
}
\value{
list of PC scores, factor loadings, p-value and q-value (by Benjamini and Hochberg),
 and contribution ratios
}
\references{
Yamamoto H., Nakayama Y., Tsugawa H. (2021) OS-PCA: Orthogonal Smoothed Principal Component Analysis Applied to Metabolome Data, Metabolites, 11(3):149.
}
\author{
Hiroyuki Yamamoto
}
\examples{
data(GreenTea)

ospca <- os_pca(X)

ospca <- ospca_loading(ospca)
ospca$loading$R
ospca$loading$p.value
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pca/msea }


\name{ospca_loading}
\alias{ospca_loading}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Orthogonal smoothed principal component loading}
\description{
This function computes orthogonals smoothed principal component (OS-PC) loading from the result of the "os_pca" function.
In this function, data matrix is should be scaled to zero mean and unit variance (i.e. autoscaling) for each variables.
This is also wrapper function of pls_loading function.
}
\usage{
ospca_loading(ospca)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ospca}{

  The following variables (P,T,Q and U) are included in the ospca object.

	P : A matrix with OS-PC loading in each column.

  T : A matrix with OS-PC score in each column

  Q : A matrix with OS-PC loading for auxiliary variable in each column.

  U : A matrix with OS-PC score for auxiliary variable in each column

  }
}
\details{
The OS-PC loading and the p-value by statistical hypothesis testing is added to the ospca object returned by the "os-pca" function.
}
\value{
  \item{loadings}{

  The following variables (R and p.value) are included in the loadings object.

  R : OS-PCA loading (Correlation coefficient between OS-PC score of auxiliary variable and each variables in data matrix.)

  p.value : p-value of OS-PCA loading

  }
}
\references{
Yamamoto H., Nakayama Y., Tsugawa H. (2021) OS-PCA: Orthogonal Smoothed Principal Component Analysis Applied to Metabolome Data, Metabolites, 11(3):149.
}
\author{
Hiroyuki Yamamoto
}
\examples{
# metabolic turnover data
data(turnover)
X <- turnover$X
D <- turnover$D

ospca <- os_pca(X,D,0.999)

ospca <- ospca_loading(ospca)
ospca$loading$R
ospca$loading$p.value

# metabolome data

data(greentea)
X <- greentea$X
D <- greentea$D
M <- greentea$M

ospca <- os_pca(X,D,0.999,M)

ospca <- ospca_loading(ospca)
ospca$loading$R
ospca$loading$p.value
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ orthogonal smoothed principal component analysis / loading }


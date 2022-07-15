\name{os_pca}
\alias{os_pca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Orthogonal smoothed principal component analysis}
\description{
This function performs orthogonal smoothed principal component analysis (OS-PCA).
In this function, data matrix is automatically scaled to zero mean and unit variance (i.e. autoscaling) for each variables.
}
\usage{
os_pca(X,D,kappa,M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
	Data matrix that include variables in each columns.
  }
  \item{D}{
	Differential matrix.
  }
  \item{kappa}{
	The smoothing parameter (default : 0.999).
  }
  \item{M}{
	Averaging matrix for repeated data (default : Identity matrix).
  }
}
\details{
The kappa represents the degree of smoothing. The value of kappa increases, the strength of the smoothing increases.
}
\value{
The return value is a list object that contains the following elements:

P : A matrix with OS-PC loading in each column

T : A matrix with OS-PC score in each column

MT : A matrix with averaging OS-PC score for repeated data in each column
(If not for repeated data, the matrix MT equals to the matrix T)

Q : A matrix with OS-PC loading for auxiliary variable in each column

U : A matrix with OS-PC score for auxiliary variable in each column

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

# metabolome data
data(greentea)

X <- greentea$X
D <- greentea$D
M <- greentea$M

ospca <- os_pca(X,D,0.1,M)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }


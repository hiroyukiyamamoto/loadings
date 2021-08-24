\name{turnover}
\alias{turnover}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Metabolomic Turnover analysis of the yeast culture mediums}
\description{
This includes differential isotopomer ratio from metabolomic turnover analysis of yeast cultures by gas chromatgraphy mass spectrometry. The data matrix contains 11 samping time for three groups of different yeast strains and culture condition (33 samples total) and 60 peaks.
}
\usage{
data(turnover)
}
%- maybe also 'usage' for other objects documented here.
\format{
The list object "turnover" contains the following elements:

	X : Data matrix of explanatory variables that include variables in each columns.

	Y : Dummy matrix that include group information 0,1 in each culumns.

  D : Differential matrix.

}

\references{
Nakayama Y., Tamada Y., Tsugawa H., Bamba T., Fukusaki E. Novel Strategy for Non-Targeted Isotope-Assisted Metabolomics by Means of Metabolic Turnover and Multivariate Analysis. Metabolites. 2014;4:722-739.
}
\examples{
data(turnover)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ metabolome data }


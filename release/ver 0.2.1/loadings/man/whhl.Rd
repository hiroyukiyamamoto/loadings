\name{whhl}
\alias{whhl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Metabolome analysis of plasma, liver, heart and brain samples from rabbits with hyperlipidemia.}
\description{
This is metabolome data from plasma, liver, heart and brain samples of normal, hyperlipidemic, and statin-treated rabbits measured by capillary electrophoresis time-of-flight mass spectrometry. The data matrix contains 9 samples (3 samples in each groups).
}
\usage{
data(whhl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
The list object "whhl" contains the following elements:

  X : List of data matrix that include samples in each rows and metabolites in each columns for plasma, liver, heart and brain samples.

  Y : Dummy matrix that include group information 0,1 in each columns.

  D : Example of differential matrix used in PLS-ROG or multiset PLS-ROG.

  tau : Example of strength parameter matrix used in mutliset PLS or multiset PLS-ROG. All weights were set as the same value of 0.1.

}
\references{
Ooga, T., Sato, H., Nagashima, A., Sasaki, K., Tomita, M., Soga T. and Ohashi, Y. (2011) Metabolomic Anatomy of Animal Model Revealing Homeostatic Imbalances in Dyslipidemia. Y. Mol. BioSyst.,7, 1217-23.
}
\examples{
data(whhl)

liver <- whhl$X$liver


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datasets }


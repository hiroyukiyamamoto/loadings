\name{covid19}
\alias{covid19}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Proteome and metabolome data of COVID-19 serum samples.}
\description{
This is proteome and metabolome data of COVID-19 serum samples [B. Shen et al, Cell. 182 (2020) 59-72.e15.]. This dataset includes three groups of healthy control, mild and severe COVID-19.
}
\usage{
data(covid19)
}
%- maybe also 'usage' for other objects documented here.
\format{
The list object "covid19" contains the following elements:

  X : List of data matrix that include samples in each rows and protein or metabolite in each columns for proteomics and metabolomics.

  Y : Dummy matrix that include group information 0,1 in each columns.

  D : Example of differential matrix used in PLS-ROG or multiset PLS-ROG.

  tau : Example of strength parameter matrix used in mutliset PLS or multiset PLS-ROG. All weights were set as the same value of 0.1.

  K : List of kernel matrix by hyperbolic tangent kernel (scale=1, offset=1).

}
\source{
B. Shen, et al, Proteomic and Metabolomic Characterization of COVID-19 Patient Sera, Cell. 182 (2020) 59-72.e15.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datasets }


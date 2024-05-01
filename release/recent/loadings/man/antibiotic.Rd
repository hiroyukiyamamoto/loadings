\name{antibiotic}
\alias{antibiotic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Metagenome and metabolome data from antibiotic-treated mouse feces.}
\description{
This is metagenome and metabolome data of antibiotic-treated mouse fecal samples [Y. Wakita, et al, BMC Microbiology, 18:188 (2018)]. This dataset includes six groups, comprising three types (Clea, CR, and Slc) of control mice and their antibiotic-treated counterparts.
}
\usage{
data(antibiotic)
}
%- maybe also 'usage' for other objects documented here.
\format{
The list object "antibiotic" contains the following elements:

  data$X$microbiome$dist: Weighted UniFrac distance matrix for the metagenome.

  data$X$microbiome$comp: Composition data for the metagenome.

  data$X$metabolome : Metabolome data.

  data$Y: A dummy matrix that includes group information with 0 or 1 in each column.

  annotation$metabolome : Compound annotation for metabolome data.

}
\source{
Y. Wakita, et al, Taxonomic classification for microbiome analysis, which correlates well with the metabolite milieu of the gut, BMC microbiology, 18:188 (2018)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datasets }

\name{MBmca-package}
\alias{MBmca-package}
\alias{MBmca}
\docType{package}
\title{
Microbead Surface Nucleic Acid Melting Curve Analysis
}
\description{ Nucleic acid Melting Curve Analysis is a powerful method to
investigate the interaction of double stranded nucleic acids. The MBmca package provides data sets 
and lightweight utilities for nucleic acid melting curve analysis and presentation on microbead
surfaces. Methods include melting curve data preprocessing (smooth, normalize, rotate, background
subtraction), data inspection (comparison of multiplex melting curves) with location parameters
(mean, median), deviation parameters (standard deviation, median absolute deviation, Coefficient of
Variation) and calculation of the melting peaks including the second derivative. The
second derivative melting peaks is implemented as parameter to further
characterize the melting behavior. Plot functions to illustrate data quality,
smoothed curves and derivatives are available too.
}
\details{
\tabular{ll}{
Package: \tab MBmca\cr
Type: \tab Package\cr
Version: \tab 0.0.2-1\cr
Date: \tab 2013-05-30\cr
License: GPL (>= 2)
}

}
\author{
Stefan Roediger <stefan_roediger@gmx.de>
}
\references{ 
A Highly Versatile Microscope Imaging Technology Platform for the Multiplex Real-Time Detection of 
Biomolecules and Autoimmune Antibodies. S. Roediger, P. Schierack, A. Boehm, J. Nitschke, I. Berger, 
U. Froemmel, C. Schmidt, M. Ruhland, I. Schimke, D. Roggenbuck, W. Lehmann and C. Schroeder.
\emph{Advances in Biochemical Bioengineering/Biotechnology}. N(N):1--40, 2012. 
\url{http://www.ncbi.nlm.nih.gov/pubmed/22437246}

Surface Melting Curve Analysis with R. S. Roediger, A. Boehm and I. Schimke. \emph{The R Journal}. 
N(N):N--NN, 2013. \url{http://journal.r-project.org/}
}

\keyword{ package }
\examples{
require(MBmca)
}

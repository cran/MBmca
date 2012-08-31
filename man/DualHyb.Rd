\name{DualHyb}
\alias{DualHyb}
\docType{data}
\title{ Surface melting curve data from direct hybridization experiment of short oligonucleotide 
probes with bimodal melting curve pattern.
}
\description{
	     A melting curve experiment with four microbead populations and short oligonucleotide probes 
(direct hybridization) and longer probes (dual-hybridization probes) capture probe. Detection probes 
for human VIM (vimentin), MLC-2v (myosin regulatory light chain 2, ventricular/cardiac muscle isoform) 
and SERCA2 (Atp2a2 - ATPase, Calcium-transporting ATPase sarcoplasmic reticulum type, slow twitch skeletal muscle isoform). 
One sequence of VIM contained a mutation at position 41.
}
\usage{data(DualHyb)}
\format{
  A data frame with the melting curves of three different capture and detection probe pairs for HRPT1 and MLC-2v. First column contains 
  the temperature (in degree Celsius, 0.5 degree Celsius per step) followed by melting curves of HRPT1 on 12 microbead populations and  melting curves of MLC-2v 
  on 12 microbead populations.
  \describe{
    \item{\code{T}}{a numeric vector, Temperature in degree Celsius.}
    \item{\code{MLC2v}}{a numeric vector, MLC-2v with quencher-labeled detection probe}
    \item{\code{SERCA2}}{a numeric vector, SERCA2 without quencher-labeled detection probe}
    \item{\code{VIM.w.Mutation}}{a numeric vector, mutated VIM with quencher-labeled detection probe}
    \item{\code{VIM.wo.Mutation}}{a numeric vector, native VIM with quencher-labeled detection probe}
  }
}

\details{ The melting curve was conducted with short oligonucleotide probes
(direct hybridization) and longer probes (dual-hybridization probes) on the
surface of microbeads (sequences and materials according to Roediger et al.
(2012)) using the VideoScan platform by Roediger et al. (2012). The dyes and
quencher used were Atto 647N and BHQ2.
}

\source{ Original data were measured with the VideoScan platform by Roediger et al. (2012).
A Highly Versatile Microscope Imaging Technology Platform for the Multiplex Real-Time Detection of Biomolecules and Autoimmune Antibodies.
S. Roediger, P. Schierack, A. Boehm, J. Nitschke, I. Berger, U. Froemmel, C. Schmidt, M. Ruhland, I. Schimke, D. Roggenbuck, W. Lehmann and C. Schroeder.
\emph{Advances in Biochemical Bioengineering/Biotechnology}. N(N):1--40, 2012. \url{http://www.ncbi.nlm.nih.gov/pubmed/22437246}
}
\seealso{ \code{\link{MFIerror}}, \code{\link{mcaSmoother}}, \code{\link{diffQ}}, \code{\link{diffQ2}}, \code{\link{DMP}}, \code{\link{MultiMelt}}}
\examples{
data(DualHyb)
}
\keyword{datasets}

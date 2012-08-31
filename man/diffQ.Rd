\name{diffQ}
\alias{diffQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculation of the melting temperature (Tm) from the first derivative
}
\description{ \code{diffQ} is used to calculate the melting temperature (Tm)
but also for elementary graphical operations (e.g., show the Tm or the
derivative). It does not require smoothed data for the MCA. First the
approximate Tm is determined as the \code{min()} and/or \code{max()} from the
first derivative. This approximate peak value is the starting-point for an
accurate calculation. The function takes a defined number n (maximum 8) of the
left and the right neighbor values and fits a quadratic polynomial. The
quadratic regression of the X (temperature) against the Y (fluorescence)
range gives the coefficients. The optimal quadratic polynomial is chosen based
on the highest adjusted R-squared value. \code{diffQ} returns an
objects of the class \code{list}. To accessing components of lists is done as
described elsewhere either be name or by number.
}
\usage{
diffQ(xy, FCT = min, fws = 8, col = 2, PEAK = FALSE, negDeriv = TRUE, 
	  deriv = FALSE, derivlimits = FALSE, derivlimitsline = FALSE,
	  vertiline = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xy}{ is a \code{data.frame} containing in the first column the temperature and in the second column the fluorescence values. Preferably the output from \code{mcaSmoother} is used. }
  \item{FCT}{ accepts \code{min} or \code{max} as option and is used to define whether to find a local minimum (``negative peak'') or local maximum (``positive peak''). }
  \item{fws}{ defines the number (n) of left an right neighbors to use for the calculation of the quadratic polynomial. }
  \item{col}{ is a graphical parameter used to define the length of the line used in the plot. }
  \item{PEAK}{ shows the peak in the plot. }
  \item{negDeriv}{ uses the positive first derivative instead of the negative. }
  \item{deriv}{ shows the first derivative with the color assigned to \code{col}. }
  \item{derivlimits}{ shows the neighbors (fws) used to calculate the Tm as points in the plot. }
  \item{derivlimitsline}{ shows the neighbors (fws) used to calculate the Tm as line in the plot. }
  \item{vertiline}{ draws a vertical line at the Tms.}
}

\value{
    \item{diffQ() }{ returns a comprehensive list with results from the first derivative. The list includes a \code{data.frame} 
		      of the derivative ("xy"). The temperature range ("limits_xQ") and fluorescence range ("limits_diffQ") to 
		      calculate the peak value. "fluo_x" is the approximate fluorescence at the approximate melting temperature. 
		      The calculated melting temperature ("Tm") with the corresponding fluorescence intensity ("fluoTm"). 
		      The number of neighbors ("fws") and the adjusted R-squared ("adj.R.squ") to fit.
	  }

}

\author{
Stefan Roediger
}

\seealso{ \code{\link{diffQ2}}, \code{\link{mcaSmoother}}}
\examples{
data(MultiMelt)
HPRT1 <- matrix(NA,12,4,
	dimnames = list(colnames(MultiMelt[,c(2:13)]),
	    c("Fluo", "Tm", "Tm1D2", "Tm2D2")))
for (i in 2:13) {
       tmp <- mcaSmoother(MultiMelt[,1], 
                     MultiMelt[,i])
       tmpTM <- diffQ2(tmp)
       HPRT1[i-1,1] <- max(tmp[["y.sp"]])
       HPRT1[i-1,2] <- tmpTM[["TmD1"]][["Tm"]]
       HPRT1[i-1,3] <- tmpTM[["Tm1D2"]][["Tm"]]
       HPRT1[i-1,4] <- tmpTM[["Tm2D2"]][["Tm"]]
}
plot(HPRT1[,1], HPRT1[,2], 
       xlab="refMFI", ylab="T", 
       main="HPRT1", xlim=c(2.1,2.55),
       ylim=c(72,82), pch=19,
       col=1:12, cex=1.8)
points(HPRT1[,1], HPRT1[,3], pch=15)
points(HPRT1[,1], HPRT1[,4], pch=15)
abline(lm(HPRT1[,2] ~ HPRT1[,1]))
abline(lm(HPRT1[,3] ~ HPRT1[,1]))
abline(lm(HPRT1[,4] ~ HPRT1[,1]))
}

\keyword{ melting }
\keyword{ Tm }


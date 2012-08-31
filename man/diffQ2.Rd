\name{diffQ2}
\alias{diffQ2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Calculation of the melting temperatures (Tm, Tm1D2 and Tm2D2) from the first and the second derivative
}
\description{ \code{diffQ2()} calls instances of \code{diffQ()} to calculate the
Tm1D2 and Tm2D2. The options are similar to \code{diffQ()}. Both
\code{diffQ()} and \code{diffQ2()} return objects of the class \code{list}. To
accessing components of lists is done as described elsewhere either be name or
by number.
}
\usage{
diffQ2(xy, FCT = max, fws = 8, col = 2, PEAK = FALSE, deriv = FALSE, 
	   negDeriv = TRUE, derivlimits = FALSE, 
	   derivlimitsline = FALSE, vertiline = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xy}{ is a \code{data.frame} containing in the first column the temperature and in the second column the fluorescence values. Preferably the output from \code{mcaSmoother} is used. }
  \item{FCT}{ accepts \code{min} or \code{max} as option and is used to define whether to find a local minimum (``negative peak'') or local maximum (``positive peak''). }
  \item{fws}{ defines the number (n) of left an right neighbors to use for the calculation of the quadratic polynomial. }
  \item{col}{ is a graphical parameter used to define the length of the line used in the plot. }
  \item{PEAK}{ shows the peak in the plot. }
  \item{negDeriv}{ calculates the negative derivative (default). If \code{FALSE} the positive first negative is calculated. }
  \item{deriv}{ shows the first derivative with the color assigned to \code{col}. }
  \item{derivlimits}{ shows the number (n) used to calculate the Tm as points in the plot. }
  \item{derivlimitsline}{ shows the number (n) used to calculate the Tm as line in the plot. }
  \item{vertiline}{ draws a vertical line at the Tms. }
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{TmD1 }{\code{TmD1} returns a comprehensive list with results from the first derivative. The list includes a \code{data.frame} of the derivative ("xy"). 
	      The temperature range ("limits_xQ") and fluorescence range ("limits_diffQ") to calculate the peak value. "fluo_x" is the approximate fluorescence at the approximate melting temperature.
	      The calculated melting temperature ("Tm") with the corresponding fluorescence intensity ("fluoTm"). 
	      The number of points ("fws") and the adjusted R-squared ("adj.R.squ") to fit.
		}
  \item{Tm1D2 }{\code{Tm1D2} returns comprehensive a list (start of melting process) from the second derivative similar to TmD1.}
  \item{Tm2D2 }{\code{Tm2D2} returns comprehensive a list (end of melting process) from the second derivative similar to TmD1.}
  \item{xTm1.2.D2 }{\code{xTm1.2.D2} returns the combined temperature of Tm1D2 and Tm2D2 (start and end of melting process) from the second derivative.}
  \item{yTm1.2.D2 }{\code{yTm1.2.D2} returns the combined intensity of Tm1D2 and Tm2D2 (start and end of melting process) from the second derivative.}
}

\author{ Stefan Roediger }
\seealso{ \code{\link{diffQ}}, \code{\link{mcaSmoother}}}

\examples{
data(MultiMelt)
HPRT1 <- matrix(NA,12,4,
	dimnames = list(colnames(MultiMelt[,c(2:13)]),
	    c("Fluo", "Tm", "Tm1D2", "Tm2D2")))
for (i in 2:13) {
       tmp <- mcaSmoother(MultiMelt[,1], 
                     MultiMelt[,i])
       tmpTM <- diffQ2(tmp)
       HPRT1[i-1,1] <- max(tmp$y)
       HPRT1[i-1,2] <- tmpTM$TmD1$Tm
       HPRT1[i-1,3] <- tmpTM$Tm1D2$Tm
       HPRT1[i-1,4] <- tmpTM$Tm2D2$Tm
}
HPRT1
}

\keyword{ Tm }
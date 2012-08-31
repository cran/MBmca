\name{MFIerror}
\alias{MFIerror}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiple comparison of the temperature dependent variance of the refMFI
}
\description{
 MFIerror is used for a fast multiple comparison of the
temperature dependent variance of the refMFI. MFIerror returns an
object of the class data.frame with columns ``Temperature'', ``Location''
(Mean, Median), ``Deviation'' (Standard Deviation, Median Absolute Deviation)
and ``Coefficient of Variation''.
}
\usage{
MFIerror(x, y, CV = FALSE, RSD = FALSE, 
	 rob = FALSE, PLOT = TRUE, type = "p", 
	 pch = 19, length = 0.05, col = "black")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	    is the column of a data frame for the temperature.
}
  \item{y}{
	    are multiple columns of fluorescence values from a
\code{data.frame} (e.g., [,c(1:n)]).
}
  \item{CV}{ If \code{CV} is true the coefficient of variation (CV) is plotted.
}
  \item{RSD}{ Setting the option \code{RSD=TRUE} shows the relative standard
deviation (RSD) in percent.
}
  \item{rob}{ Using the option \code{rob} as TRUE the median and the median
absolute deviation (MAD) is plotted instead of the mean and standard deviation.
}
  \item{PLOT}{sets \code{MFIerror()} to plot
the results (default). In the default setting (\code{CV=FALSE}) the mean with the standard
deviations is plotted.
}
  \item{type}{ is a graphical parameter setting the plot use lines,
points or both (see \code{\link{plot}}).
}
  \item{pch}{ is a graphical parameter used to definde the symbol
used in the plot.
}
  \item{length}{ \code{length} is a graphical parameter used to definde the
length of the error bar used in the plot.
}
  \item{col}{ \code{col} is a graphical parameter used to define the length of
the error bar used in the plot.
}
}

\value{
  \item{res }{
      returns a \code{data.frame} containing the "Temperature", "Location" (mean, median), "Deviation" (standard deviation, median absolute deviation), "Coefficient of Variance" (CV, RSD) sequential in the columns.}
}

\author{
Stefan Roediger
}
\seealso{ \code{\link{mcaSmoother}}}
\examples{
data(MultiMelt)
MFIerror(MultiMelt[,1], MultiMelt[,c(2:13)])
}

\keyword{ deviation }
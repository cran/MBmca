\name{mcaSmoother}
\alias{mcaSmoother}
\title{
Function to preprocess melting curve data.
}
\description{
The \code{mcaSmoother()} is used for data preprocessing. Measurements from
experimental systems may occasionally include missing values
(NA). \code{mcaSmoother()} uses \code{approx()} to fill up NAs under the
assumption that all measurements were equidistant. The original data remain
unchanged and only the NAs are substituted. Following it calls
\code{smooth.spline()} to smooth the curve. Different strengths can be set using
the option \code{df\_fact} (f default~0.95). Internally it takes the degree of
freedom value from the spline and multiplies it with a factor between 0.6 and 1.1.
Values lower than 1 result in stronger smoothed curves. A temperature range for
a simple linear background correction. The linear trend is estimated by a robust
linear regression using \code{lmrob()}. In case criteria for a robust linear
regression are violated \code{lm()} is automatically used. Optionally a Min-Max
normalization between 0 and 1 can be used by setting the option \code{MinMax} to
\code{TRUE}. This is useful if the fluorescence values between samples vary
considerably, for example due to high background, particularly in solution.
}

\usage{ mcaSmoother(x, y, bgadj = FALSE, bg = NULL, MinMax = FALSE, df_fact = 0.95) }

\arguments{
  \item{x}{ \code{x} is the column of a data frame for the temperature. }
  \item{y}{ \code{y} is the column of a data frame for the fluorescence values. }
  \item{bgadj}{ \code{bgadj} is used to adjust the background signal. This causes \code{mcaSmoother} to use the data subset defined by \code{bg} for the linear regression and background correction. }
  \item{bg}{ \code{bg} is used to define the range for the background reduction (e.g., \code{bg=c(1:21)}). }
  \item{MinMax}{ \code{MinMax} is used to scale the fluorescence a Min-Max normalization between 0 and 1 can be used by setting the option \code{MinMax} to \code{TRUE}. }
  \item{df_fact}{ \code{df_fact} is a factor to smooth the curve. Different strengths can be set using the option \code{df\_fact} (f default~0.95). Internally it takes the degree of freedom value from the spline and multiplies it with a factor between 0.6 and 1.1. Values lower than 1 result in stronger smoothed curves. }
}
\value{ \item{xy }{returns a \code{data.frame} with the temperature ("x") in the first and the preprocessed fluorescence values ("y.sp") in the second column.} }

\author{ Stefan Roediger }

\seealso{ \code{\link{MFIerror}}, \code{\link{lmrob}}, \code{\link{smooth.spline}}, \code{\link{lm}}}

\examples{
data(DMP)

plot(DMP[,1], DMP[,6], 
      xlim=c(20,95), xlab="T",
      ylab="refMFI", pch=19, col=8)
f <- c(0.6, 0.8, 1.0)
for (i in c(1:3)) { 
 	lines(mcaSmoother(DMP[,1],
         DMP[,6],df_fact=f[i]),
         col=i, lwd=2)
}
legend(20, 1.5, paste("f", f, sep=": "),
      cex=1.2, col= 1:3, bty="n",
      lty=1, lwd=4)
}

\keyword{ smooth }
\keyword{ background }
\keyword{ normalization }

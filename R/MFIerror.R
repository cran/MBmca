"MFIerror" <-
function (x, y, CV=FALSE, RSD=FALSE, rob=FALSE, PLOT=TRUE, type="p", pch=19, length=0.05, col="black")
{	#Define if "robust" or standard function should be used as measures
	options(warn = -1)
	#Test if x and y exist.
	if (is.null(x)) 
	   stop("Enter temperature")
	if (is.null(y)) 
	   stop("Enter fluorescence data")

	if (rob) {
		      loc.fct <- median
		      dev.fct <- mad
			} else {
				loc.fct <- mean
				dev.fct <- sd
				}
	y.m	<-	apply(y, 1, loc.fct)
	y.sd	<-	apply(y, 1, dev.fct)
	if (RSD) {
		  y.cv	<-	(y.sd / y.m) * 100
	   
		  } else {
			  y.cv	<-	y.sd / y.m
			    }
	res	<-	data.frame(x, y.m, y.sd, y.cv)
	names(res)	<- c("Temperature", "Location", "Deviation", "Coefficient of Variance")
	#Plot the Coefficient of Variance
	if (PLOT) {
		   if (CV) {
			    plot(res[,1], res[,4], xlab="T", ylab="CV", col=col, pch=pch)
		      #Plot the location with error bars.
		      } else {
			      plot(res[,1], res[,2], ylim=c(min(res[,2] - res[,3]), max(res[,2] + res[,3])), xlab="T", ylab="MFI", type=type)
			      arrows(res[,1], res[,2] + res[,3], res[,1], res[,2] - res[,3], angle=90, code=3, length=length, col=col)
			      }
		  }
	#res is the an object of the type data.frame containing the temperature, location, deviation and coefficient of variance.
	return(res = res)
}

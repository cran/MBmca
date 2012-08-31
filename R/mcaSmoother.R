"mcaSmoother" <-
function(x, y, bgadj=FALSE, bg=NULL, MinMax=FALSE, df_fact=0.95)
{	
	options(warn = -1)
	#Test if df_fact is within a meaningful range.
	if (df_fact < 0.6 || df_fact > 1.1) 
	    stop("df_fact size must be within 0.6 and 1.1.")
	#Test if x and y exist and have identical lengths.
	if (is.null(x)) 
	   stop("Enter temperature")
	if (is.null(y)) 
	   stop("Enter fluorescence data")
	if (length(x) != length(y)) 
	   stop("Use temperature and fluorescence data with same number of elements")
	#Test if background adjustment was set and background vector is not empty.
	if ((is.null(bg)) && (bgadj==TRUE))
	   stop("Enter temperature background range.")

	#Test if y contains missing values. In case of missing values a regression is 
	#used to estimate the missing value.
	if (length(which(is.na(y) == TRUE)) > 0) { 
	y[which(is.na(y))] <- approx(x, y, n = length(x))$y[c(which(is.na(y)))]}

	#Smooth the curve with a cubic spline. Takes first the degree of freedom from the cubic spline.
	#The degree of freedom is than used to smooth the curve by a user defined factor.
	df_tmp <- data.frame(smooth.spline(x,y)$df) 
	y.sp <- smooth.spline(x, y, df = (df_tmp * df_fact))$y

	# If the argument bgadj is set TRUE, bg must be  used to define a temperature range for a linear 
	#background correction. The linear trend is estimated by a robust linear regression using lmrob().
	#In case criteria for a robust linear regression are violated lm() is automatically used.
	if (bgadj) {
	  	if (class(try(lmrob(y.sp[bg] ~ x[bg]), silent = T)) == "try-error") { 
		coefficients <- data.frame(lm(y.sp[bg] ~ x[bg])[1]) 
		  	} else { 
		    	lmrob_control <- suppressWarnings(lmrob(y.sp[bg] ~ x[bg])) 
		    	if ((class(lmrob_control) != "try-error") && (lmrob_control$converged == TRUE)) { 
		      		coefficients <- data.frame(lmrob(y.sp[bg] ~ x[bg])[1]) 
		    	} else { 
		      		coefficients <- data.frame(lm(y.sp[bg] ~ x[bg])[1]) 
		    		} 
		  	} 
		  y_norm <- y.sp - (coefficients[2,1] * x + coefficients[1,1]) #Subtracts the linear trend from the smoothed values.
	} else {y_norm <- data.frame(y.sp)}
		
	if (MinMax) {y_norm <- (y_norm - min(y_norm)) / (max(y_norm) - min(y_norm))} #Performs a "Min-Max Normalization" between 0 and 1.
	
	#Returns an object of the type data.frame containing the temperature in the first column 
	#and the preprocessed fluorescence data in the second column.
	return(xy = data.frame(x,y_norm))
}
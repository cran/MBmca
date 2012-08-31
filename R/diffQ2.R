"diffQ2" <-
function (xy,FCT=max,fws=8,col=2,PEAK=FALSE,deriv=FALSE,negDeriv=TRUE,derivlimits=FALSE,derivlimitsline=FALSE,vertiline=FALSE)
{
	    #Test if fws (number of neighbors) is within a meaningful range.
	    options(warn = -1)
	      fws <- round(fws)
	      if (fws < 1 || fws > 8) 
	      stop("Fit window size must be within 1 and 8.")
	    
	    #Calls instances of diffQ to calculate the Tm of the first and the second derivatives.
	    #The arguments are similar to diffQ
	    #Calculates the first derivative and provides the starting information for the second derivative.
	    TmD1	<-	diffQ(xy,fws=8,col=2,PEAK=PEAK,deriv=deriv,negDeriv=negDeriv,derivlimits=derivlimits,derivlimitsline=derivlimitsline,vertiline=vertiline)
	    #Calculates the second derivative and from the first derivative for the minimum of the melting curve.
	    Tm1D2	<-	diffQ(TmD1$xy,FCT=min,fws=fws,col=col,PEAK=PEAK,deriv=deriv,negDeriv=FALSE,derivlimits=derivlimits,derivlimitsline=derivlimitsline,vertiline=vertiline)
	    #Calculates the second derivative and from the first derivative for the maximum of the melting curve.
	    Tm2D2	<-	diffQ(TmD1$xy,FCT=max,fws=fws,col=col,PEAK=PEAK,deriv=deriv,negDeriv=FALSE,derivlimits=derivlimits,derivlimitsline=derivlimitsline,vertiline=vertiline)

	    DATA	<-	TmD1[["xy"]]
	    LIMIT	<-	Tm1D2[["limits_xQ"]]
	    REG1	<- 	subset(DATA, DATA[["xQ"]] >= LIMIT[1] 
					& DATA[["xQ"]] <= LIMIT[length(LIMIT)], select=c(xQ,diffQ))

	    LIMIT	<-	Tm2D2[["limits_xQ"]]
	    REG2	<-	subset(DATA, DATA[["xQ"]] >= LIMIT[1] 
					& DATA[["xQ"]] <= LIMIT[length(LIMIT)], select=c(xQ,diffQ))

	    xQ 		<- NULL; rm(xQ); 

	    #Calculates the melting temperatures based on a local regression with a quadratic polynomial for the Tm1D2 peak value.
	    lm2.1	<-	lm(REG1[,2] ~ REG1[,1] + I(REG1[,1]^2))
	    coeflm2.1	<-	data.frame(lm2.1[1])[c(1:3),]
	    coeflm2_y.1 <-	coeflm2.1[1] + coeflm2.1[2] * Tm1D2[["Tm"]] + coeflm2.1[3] * Tm1D2[["Tm"]]^2

	    #Calculates the melting temperatures based on a local regression with a quadratic polynomial for the Tm2D2 peak value.
	    lm2.2	<-	lm(REG2[,2] ~ REG2[,1] + I(REG2[,1]^2))
	    coeflm2.2	<-	data.frame(lm2.2[1])[c(1:3),]
	    coeflm2_y.2	<-	coeflm2.2[1] + coeflm2.2[2] * Tm2D2[["Tm"]] + coeflm2.2[3] * Tm2D2[["Tm"]]^2
	    
	    #Vectors of the two melting temperatures of the second derivative.
	    x <- c(Tm1D2[["Tm"]], Tm2D2[["Tm"]])
	    #Vectors of the two intensities at the melting temperatures of the second derivative.
	    y <- c(coeflm2_y.1, coeflm2_y.2)

	    #Returns an object of the type list containing the data and data.frames from above including the approximate 
	    #difference quotient values, melting temperatures of the first derivative and the second derivative, intensities and used neighbors.
	    return(list(TmD1=TmD1, Tm1D2=Tm1D2, Tm2D2=Tm2D2, xTm1.2.D2=x, yTm1.2.D2=y))
}
"diffQ" <-
function (xy,FCT=min,fws=8,col=2,PEAK=FALSE,negDeriv=TRUE,deriv=FALSE,derivlimits=FALSE,derivlimitsline=FALSE,vertiline=FALSE)
{ 
	  #Test if fws (number of neighbors) is within a meaningful range.
	  options(warn = -1)
	    fws <- round(fws)
	    if (fws < 1 || fws > 8) 
	    stop("Fit window size must be within 1 and 8.")
	    
	  list_res <- list()

	  #Take the x and y values from the object of type data.frame.
	  x <- xy[,1]
	  y <- xy[,2]

	  options(warn = -1)
	  #Test if x and y exist.
	  if (is.null(x)) 
	      stop("Enter temperature")
	  if (is.null(y)) 
	      stop("Enter fluorescence data")

	  #Function to calculate the approximate derivative.
	  if (negDeriv) {y <- y * -1}
	  N <- length(x); low <- 1:2; up <- (N-1):N
	  xQ <- x[3:length(x) - 1] 
	  diffQ <- (y[-low] - y[-up]) / (x[-low] - x[-up])

	  #Optional draw line of calculated approximate derivative.
	  if (deriv) {lines(xQ, diffQ, col = col)} 

	  #First: Find approximate range of neighbors for minimum or maximum of temperature peak.
	  #Second: Calculate quadratic polynomial for ranges of minimum or maximum temperature peaks.
	  #Return the results in an object of the type list.
	  for (i in 2:9) {
	    suppressMessages(RANGE <- which(diffQ == FCT(diffQ[(i + 1):(N-i)]))) 
	    if (class(try(na.omit(xQ[(RANGE - i):(RANGE + i)]), silent = T)) == "try-error") {
	      list_res[[i-1]] <- NA
	    } else {
	      limits_xQ <- na.omit(xQ[(RANGE-i):(RANGE + i)])
	      limits_diffQ <- na.omit(diffQ[(RANGE-i):(RANGE + i)])
	      fluo_x <- limits_diffQ[length(limits_diffQ) - i]
	      lm2 <- lm(limits_diffQ ~ limits_xQ + I(limits_xQ^2))
	      coeflm2 <- data.frame(lm2[1])[c(1:3),]
	      coeflm2_y <- coeflm2[1] + coeflm2[2] * limits_xQ + coeflm2[3] * limits_xQ^2
	      lm2sum <- summary(lm(coeflm2_y ~ limits_diffQ))
	      list_res[[i-1]] <- list(i,limits_xQ,limits_diffQ,fluo_x,lm2,coeflm2,coeflm2_y,lm2sum)
	    }
	  }

	  #Determine the optimal fitted quadratic polynomial for ranges of minimum or maximum temperature peaks
	  #bases on the adjusted  R squared.
	  Rsq <- matrix(NA, nrow = fws, ncol = 2)
	  colnames(Rsq) <- c("fw", "Rsqr")
	      for (i in 1:fws) {
		Rsq[i,1] <- list_res[[i]][[1]]
		Rsq[i,2] <- list_res[[i]][[8]]$adj.r.squared
	      }
	      list_res <- list_res[[which(Rsq[,2] == max(na.omit(Rsq[,2])))]]
	      names(list_res) <- list("fw","limits_xQ","limits_diffQ","fluo_x","lm2","coeflm2","coeflm2_y","lm2sum")

	      limits_xQ <- list_res$limits_xQ
	      limits_diffQ <- list_res$limits_diffQ
	      fluo_x <- list_res$fluo_x
	      lm2 <- list_res$lm2
	      coeflm2 <- list_res$coeflm2
	      coeflm2_y <- list_res$coeflm2_y
	      lm2sum <- list_res$lm2sum
	      fw <- list_res$fw

	      if (derivlimits) {points(limits_xQ, limits_diffQ, cex = 1, pch = 19, col = col)}
	      if (derivlimitsline) {lines(spline(limits_xQ[1:length(limits_xQ)], lm2[["fitted.values"]][1:length(lm2[["fitted.values"]])]), col = "orange", lwd = 2)}
	      abl <- -lm2$coefficients[2] / (2 * lm2$coefficients[3])
	      y <-	coeflm2[1] + coeflm2[2] * abl + coeflm2[3] * abl^2
	      if (vertiline) {abline(v = abl, col = "grey")}
	      if (PEAK) {points(abl,y, cex = 1, pch = 19, col = col)}

	      #Returns an object of the type list containing the data and data.frames from above including the approximate 
	      #difference quotient values, melting temperatures, intensities and used neighbors.
	      return(list(xy = data.frame(xQ,diffQ), fluo_x = fluo_x, limits_xQ = limits_xQ, limits_diffQ = limits_diffQ,
		   	Tm = abl, fluoTm = y, adj.R.squ = lm2sum$adj.r.squared, fws = fw))
}
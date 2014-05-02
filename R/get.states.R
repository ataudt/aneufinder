get.states <- function(model, threshold=0.5, separate.zeroinflation=FALSE) {
	numbins <- nrow(model$posteriors)
	states <- rep(NA,numbins)
	if (is(model,"chromStar.univariate.model")) {
		if (separate.zeroinflation) {
			states[ model$posteriors[,3]<=threshold & model$posteriors[,2]<=model$posteriors[,1] ] <- 1
			states[ model$posteriors[,3]<=threshold & model$posteriors[,2]>=model$posteriors[,1] ] <- 2
			states[ model$posteriors[,3]>threshold ] <- 3
			states <- as.factor(c("zeroinflation","unmodified","modified"))[states]
		} else {
			states <- ifelse(model$posteriors[,3]>threshold, 2, 1)
			states <- as.factor(c("unmodified","modified"))[states]
		}
		return(states)
	} else if (is(model,"chromStar.multivariate.model")) {
		states <- model$stateorder[apply(model$posteriors, 1, which.max)]
		return(states)
	} else {
		stop("Supply either a univariate or multivariate chromStar model")
	}
	
}


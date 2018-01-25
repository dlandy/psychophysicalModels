# classicPsychophysics models



#' stevensPowerLaw
#' @param stimuli a vector of stimuli, between 0 and inf
#' @param beta the power to raise the stimulus to, equivalent to the ratio of taus in the reconstucted version 
#' @param delta proportionality constant--the constant multiple on the stimulus value
#' @param responses an optional vector of responses.
#' @param mode a mode. can validly be "prediction", or "logLikelihoodOfResponses". If the latter, responses must be given, and tau must be given
#' @param tau if log-likelihood is to be calculated, we must decide the predicted variability of responses. This is captured in tau, the precision
#' @return A vector the transformed stimuli 
#' @seealso gonzalezWu, spence
#' @export
#' @examples
#' stevensPowerLaw(1:100)
#' stevensPowerLaw(1:100, beta=0.5, delta=2)
#' stevensPowerLaw(1:100, beta=0.9, delta=2, responses=2*(1:100)^.9, mode="logLikelihoodOfResponses)
stevensPowerLaw <- function(stimuli
                                   , beta = 1
                                   , delta = 1
                                   , responses = NULL
                                   , tau = 1
                                   , mode="prediction"){
  predictions <- delta * stimuli^beta
  if(mode=="prediction"){
    return(predictions)
  } else{ 
    return(0-sum(log(dnorm(responses-predictions, sd=1/sqrt(tau)))))
  }
}







#' spence
#' @param stimuli a vector of stimuli, between 0 and 1
#' @param beta the power to raise the stimulus to, equivalent to the ratio of taus in the reconsturcted version 
#' @param responses an optional vector of responses.
#' @param mode a mode. can validly be "prediction", or "logLikelihoodOfResponses". If the latter, responses must be given, and tau must be given
#' @param tau if log-likelihood is to be calculated, we must decide the predicted variability of responses. This is captured in tau, the precision
#' @return A vector the transformed stimuli 
#' @seealso classicGonzalezWu, spence
#' @export
#' @examples
#' spence(seq(0.00, 1, 0.01))
#' spence(seq(0.00, 1, 0.01), beta=0.5)
#' spence(seq(0.00, 1, 0.01), beta=0.9, responses=spence(seq(0.00, 1, 0.01), 0.8), mode="logLikelihoodOfResponses)
spence <- function(stimuli
                                   , beta = 1
                                   , responses = NULL
                                   , tau = 1
                                   , mode="prediction"){
  predictions <- stimuli^beta/(stimuli^beta + (1-stimuli)^beta)
  if(mode=="prediction"){
    return(predictions)
  } else{ 
    return(0-sum(log(dnorm(responses-predictions, sd=1/sqrt(tau)))))
  }
}


#' gonzalezWu
#' @param stimuli a vector of stimuli, between 0 and 1
#' @param beta the power to raise the stimulus to, equivalent to the ratio of taus in the reconsturcted version 
#' @param delta the power to raise the stimulus to, equivalent to the ratio of taus in the reconsturcted version 
#' @param responses an optional vector of responses.
#' @param mode a mode. can validly be "prediction", or "logLikelihoodOfResponses". If the latter, responses must be given, and tau must be given
#' @param tau if log-likelihood is to be calculated, we must decide the predicted variability of responses. This is captured in tau, the precision
#' @return A vector the transformed stimuli 
#' @seealso classicGonzalezWu, spence
#' @export
#' @examples
#' gonzalezWu(seq(0.00, 1, 0.01))
#' gonzalezWu(seq(0.00, 1, 0.01), beta=0.5)
#' gonzalezWu(seq(0.00, 1, 0.01), beta=0.9, responses=gonzalezWu(seq(0.00, 1, 0.01), 0.8), mode="logLikelihoodOfResponses)
gonzalezWu <- function(stimuli
                   , beta = 1
                   , delta = 1
                   , responses = NULL
                   , tau = 1
                   , mode="prediction"){
  predictions <- delta*stimuli^beta/(delta*stimuli^beta + (1-stimuli)^beta)
  if(mode=="prediction"){
    return(predictions)
  } else{ 
    return(0-sum(log(dnorm(responses-predictions, sd=1/sqrt(tau)))))
  }
}


#' landyBrowerBiasedSpence
#' @param stimuli a vector of stimuli, between 0 and 1
#' @param beta the power to raise the stimulus to, equivalent to the ratio of taus in the reconsturcted version 
#' @param bias a fixed shift (in log odds space) of the stimuli
#' @param responses an optional vector of responses.
#' @param mode a mode. can validly be "prediction", or "logLikelihoodOfResponses". If the latter, responses must be given, and tau must be given
#' @param tau if log-likelihood is to be calculated, we must decide the predicted variability of responses. This is captured in tau, the precision
#' @return A vector the transformed stimuli 
#' @seealso classicGonzalezWu, spence
#' @export
#' @examples
#' landyBrowerBiasedSpence(seq(0.00, 1, 0.01))
#' landyBrowerBiasedSpence(seq(0.00, 1, 0.01), beta=0.5)
#' landyBrowerBiasedSpence(seq(0.00, 1, 0.01), beta=0.9, bias=1, responses=gonzalezWu(seq(0.00, 1, 0.01), 0.8), mode="logLikelihoodOfResponses)
landyBrowerBiasedSpence <- function(stimuli
                       , beta = 1
                       , bias = 0
                       , responses = NULL
                       , tau = 1
                       , mode="prediction"){
  predictions <- bias^beta*stimuli^beta/(bias^beta*stimuli^beta + (1-stimuli)^beta)
  if(mode=="prediction"){
    return(predictions)
  } else{ 
    return(0-sum(log(dnorm(responses-predictions, sd=1/sqrt(tau)))))
  }
}







#' Create Predicted data patterns from a 'middle out' asymmetric gonzalez & Wu model
#'
#' This model calculates, given a set of parameters, the predicted results of an 'asymmetric' HD GW model. This 
#' model uses a central (but maybe off-center) reference point, and applies the G&W model from that central point ot each end
#' Each end may be separately scaled.
#' @param stimuli a Vector of stimuli
#' @param beta The curvature of the space, corresponding to a stevens beta. <1 yields overestimation of small values.
#' @param delta The location of the prior / the elevation
#' @param scaling The overall distance to the priors, used when the scaling should be symmetric. The value should go from 0 to 1. 0 means that the reference point is at infinity.
#' @param leftScaling The scaling of the lefthand side, or lower values. 
#' @param rightScaling The scaling of the righthand side, or higher values
#' @keywords proportionJudgments
#' @export
#' @examples
#' HollandsDyreGonzalezWuAsymmeric(c(1,2))
HollandsDyreGonzalezWuAsymmeric <- function(stimuli
                                            , beta=1
                                            , delta=1
                                            , scaling=1
                                            , leftScaling=scaling
                                            , rightScaling=scaling
                                            , center=0.5){
  scaling <- leftScaling*(stimuli<center) + rightScaling*(stimuli>=center)
  x <- center-(center-stimuli)*scaling
  dist <- abs(x-center)/((x>center)*(1-center)+(x<=center)*(center)) # Proportional distance from center to closest edge, ranges from 0 to 1.
  perception <- delta*dist^beta/((delta*dist^beta)+(1-dist)^beta) # scaled version of x
  result <- center+
    sign(x-center)*((x>center)*perception*(1-center) + (x<center)*perception*(center))/scaling
  result
}



#' Create Predicted data patterns from a vanilla gonzalez & Wu model, scaled only on the right.  Suitable for one-sided bounds like line length
#'
#' This model calculates, given a set of parameters, the predicted results of an 'asymmetric' HD GW model. This 
#' model uses a central (but maybe off-center) reference point, and applies the G&W model from that central point ot each end
#' Each end may be separately scaled.
#' @param stimuli a Vector of stimuli
#' @param gamma The curvature of the space, corresponding to a stevens beta. <1 yields overestimation of small values.
#' @param delta The location of the prior / the elevation
#' @keywords proportionJudgments
#' @export
#' @examples
GonzalezWuRightScaled <- function(stimuli,beta, delta, scaling=1){
  
  x <- stimuli*scaling
  perception <- (delta*x^beta) / 
    ((delta*x^beta + (1-x)^beta))
  result <- perception/scaling
}



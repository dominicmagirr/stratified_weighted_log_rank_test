#' Weighted log-rank test
#'
#' @param df A data frame. Assume standard structure for time-to-event data.
#' @param trt_colname A character string. The name of the treatment column in \code{df}.
#' @param time_colname A character string. The name of the column in \code{df} with the survival times.
#' @param event_colname A character string. The name of the column in \code{df} with the event status. Assumes 1 means event; 0 means censored.
#' @param wlr The type of weighted log-rank test. Either the default "lr" for a standard log-rank test, "mw" for a modestly-weighted log-rank test, or "fh" for the Fleming-Harrington rho-gamma family.
#' @param s_star This is a parameter for the "mw"-test. Either \code{s_star} or \code{t_star} must be specified. The weights are defined as w(t) = 1 / max(S(t-), \code{s_star}), where S is the Kaplan-Meier estimate in the pooled data.
#' @param t_star This is a parameter for the "mw"-test. Either \code{s_star} or \code{t_star} must be specified. The weights are defined as w(t) = 1 / max(S(t-), S(\code{t_star}), where S is the Kaplan-Meier estimate in the pooled data.
#' @param rho Rho parameter in the "fh"-test. The weights are defined as w(t) = S(t-)^\code{rho} (1 - S(t-))^\code{gamma}, where S is the Kaplan-Meier estimate in the pooled data.
#' @param gamma Gamma parameter in the "fh"-test. The weights are defined as w(t) = S(t-)^\code{rho} (1 - S(t-))^\code{gamma}, where S is the Kaplan-Meier estimate in the pooled data.
#' @return A data frame. The outcome of the weighted log-rank test. There is a column \code{o_minus_e_trt} to indicate which treatment the "obs - exp" refers to.
#' @export

wlrt_fast <- function(df,
                      trt_colname,
                      time_colname,
                      event_colname,
                      wlr = "lr",
                      s_star = NULL,
                      t_star = NULL,
                      rho = NULL,
                      gamma = NULL){
  
  if(!all(is.character(c(trt_colname,
                         time_colname,
                         event_colname)))) stop("trt_colname, time_colname and event_colname must all be character strings")
  
  if (any(is.na(df[,c(trt_colname,
                      time_colname,
                      event_colname)]))) stop("NA's in data set. wlrt doesn't have a default for missing data.")
  
  if (!(wlr %in% c("lr", "fh", "mw"))) stop("wlr must be one of: 'lr', 'fh', 'mw'")
  
  #### fit pooled data
  fit_pool <- survival::survfit(survival::Surv(eval(as.name(time_colname)),
                                               eval(as.name(event_colname))) ~ 1,
                                data= df,
                                timefix = FALSE)
  
  
  ### get survival probabilities for the pooled data
  fail <- fit_pool$time[fit_pool$n.event > 0]
  km_pool <- fit_pool$surv[fit_pool$n.event > 0]
  km_pool_minus <- c(1, km_pool[1:(length(km_pool) - 1)])
  
  
  if (wlr == "lr"){
    ### standard log-rank test weights
    w <- rep(1, length(km_pool_minus))
  }
  else if (wlr == "fh"){
    if (is.null(rho) || is.null(gamma)) stop("must specify rho and gamma")
    if (rho < 0 || gamma < 0) stop("rho and gamma must be non-negative")
    
    w <- km_pool_minus ^ rho * (1 - km_pool_minus) ^ gamma
    
  }
  else{
    if (is.null(t_star) && is.null(s_star)) stop("must specify either t_star or s_star")
    if (!is.null(t_star) && !is.null(s_star)) stop("must specify either t_star or s_star (not both)")
    ### modest weights
    if (!is.null(t_star)){
      
      if(any(fail < t_star)){
        
        w <- pmin(1 / km_pool_minus,
                  1 / km_pool[max(which(fail < t_star))])
      }
      else {
        w <- rep(1, length(km_pool_minus))
      }
      
    }
    else {
      w <- pmin(1 / km_pool_minus,
                1 / s_star)
    }
  }
  
  ### produce risk table
  trt <- df[[trt_colname]]
  u_trt <- sort(unique(trt))
  k <- length(u_trt) ### always equal to 2 in this package
  times <- df[[time_colname]]
  status <- df[[event_colname]]
  
  neventg <- table(trt[status > 0], times[status > 0])
  nevent <- colSums(neventg)
  
  nriskg <- matrix(1, length(fail), k)
  for (i in 1:k) nriskg[, i] <- colSums(matrix(rep(fail,each = sum(trt == u_trt[i])),,length(fail)) <= times[trt == u_trt[i]])
  nrisk <- rowSums(nriskg)
  
  ## test statistics
  
  observed <- w %*% t(neventg)
  expected <- w %*% (nriskg * (nevent/nrisk))
  
  u <- (observed - expected)[, 2]
  
  v_u <- w^2 * (nevent * (nrisk - nevent)/(nrisk - 1))
  v_u[nrisk == 1] <- 0
  v_u <- (diag(c(v_u %*% (nriskg/nrisk))) - t(nriskg/nrisk) %*% ((nriskg/nrisk) * v_u))[2,2]
  
  z = u / sqrt(v_u)
  
  data.frame(u = (observed - expected)[, 2],
             v_u = v_u,
             z = z,
             o_minus_e_trt = names((observed - expected)[, 2]))
  
}
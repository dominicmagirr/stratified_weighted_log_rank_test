#' Stratified weighted log-rank test
#'
#' @param df A data frame. Assume standard structure for time-to-event data.
#' @param trt_colname A character string. The name of the treatment column in \code{df}.
#' @param time_colname A character string. The name of the column in \code{df} with the survival times.
#' @param event_colname A character string. The name of the column in \code{df} with the event status. Assumes 1 means event; 0 means censored.
#' @param strat_colname A character string. The name of the stratifying-variable column in \code{df}. Will be converted into a factor. Care is needed here -- the function doesn't do any checks that there are sufficient numbers in each category.
#' @param wlr The type of weighted log-rank test. Either the default "lr" for a standard log-rank test, "mw" for a modestly-weighted log-rank test, or "fh" for the Fleming-Harrington rho-gamma family.
#' @param s_star This is a parameter for the "mw"-test. Either \code{s_star} or \code{t_star} must be specified. The weights are defined as w(t) = 1 / max(S(t-), \code{s_star}), where S is the Kaplan-Meier estimate in the pooled data.
#' @param t_star This is a parameter for the "mw"-test. Either \code{s_star} or \code{t_star} must be specified. The weights are defined as w(t) = 1 / max(S(t-), S(\code{t_star}), where S is the Kaplan-Meier estimate in the pooled data.
#' @param rho Rho parameter in the "fh"-test. The weights are defined as w(t) = S(t-)^\code{rho} (1 - S(t-))^\code{gamma}, where S is the Kaplan-Meier estimate in the pooled data.
#' @param gamma Gamma parameter in the "fh"-test. The weights are defined as w(t) = S(t-)^\code{rho} (1 - S(t-))^\code{gamma}, where S is the Kaplan-Meier estimate in the pooled data.
#' @return A list containing the overall test result and a data frame with the outcome of the weighted log-rank test in each strata. There is a column \code{o_minus_e_trt} to indicate which treatment the "obs - exp" refers to.
#' @export

swlrt_fast <- function(df,
                       trt_colname,
                       time_colname,
                       event_colname,
                       strat_colname = NULL,
                       wlr = "lr",
                       s_star = NULL,
                       t_star = NULL,
                       rho = NULL,
                       gamma = NULL){
  
  if (is.null(strat_colname)) return(wlrt_fast(df,
                                               trt_colname,
                                               time_colname,
                                               event_colname,
                                               wlr,
                                               s_star,
                                               t_star,
                                               rho,
                                               gamma))
  
  if (!is.character(strat_colname)) stop("strat_colname must be a character string")
  
  if(!is.factor(df[[strat_colname]])){
    df[[strat_colname]] <- as.factor(df[[strat_colname]])
  }
  
  dfs <- purrr::map(levels(df[[strat_colname]]),
                    function(x) df[df[[strat_colname]] == x,])
  
  uvs <- purrr::map_df(dfs,
                       wlrt_fast,
                       trt_colname = trt_colname,
                       time_colname = time_colname,
                       event_colname = event_colname,
                       wlr = wlr,
                       s_star = s_star,
                       t_star = t_star,
                       rho = rho,
                       gamma = gamma)
  
  rownames(uvs) <- paste0(strat_colname, levels(df[[strat_colname]]))
  
  list(by_strata = uvs,
       z = sum(uvs$u) / sqrt(sum(uvs$v_u)))
  
  
}
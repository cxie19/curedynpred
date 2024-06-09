#' Predict conditional survival probability using a joint model
#' @description Predicts an individual's conditional probability of not experiencing the event of
#' interest in an additional time \eqn{t_{hor}} given that the patient remains risk-free at least until
#' time \eqn{L} using a joint model with a flexible-hazards cure model for survival data (JMFHC) or
#' a joint model with a proportional hazards cure model for survival data (JMPHC)
#'
#' @param landmark_surv_vec a vector of landmark times for prediction of survival probabilities in the same unit as the observed survival times and the original form of measurement times
#' in the data set.
#' @param thor_list a list of vectors of the time horizons of prediction in the same unit as the parameter landmark_surv_vec in the order of the corresponding landmark times.
#' @param object an object of jmfhc_point_est() function estimated by JMFHC or JMPHC for the training data.
#' @param test_dat a data frame containing observations for prediction
#' @param no_cores the number of cores used during the prediction. The default is 7.
#'
#' @return a list containing sublists for each assigned landmark time from landmark_surv_vec and the object of jmfhc_point_est().
#'         Within each sublist at each assigned landmark time, there are data frames for time horizons of predictions containing id and
#'         corresponding predicted conditional survival probabilities.
#' @export
#'
#' @examples
#' jmfhc_est <- jmfhc::jmfhc_point_est(data=train_dat, event_time="event.time", event_status="event",
#'                                     id="id", beta_variable=c("x1","x2"), gamma_variable=c("x1","x2"),
#'                                     fu_measure_original="measure",fu_measure="measure",
#'                                     fu_time_original="measure.time",fu_time_fixed_variable="measure.time",
#'                                     fu_time_random_variable="measure.time",
#'                                     baseline_var_lmm=c("x1","x2"),no_cores=7)
#' # jmfhc_est <- readRDS("jmfhc_estresult.rds")
#' predict_surv_jmfhc  <- pred_surv_joint_model(landmark_surv_vec=c(10,20),
#'                                              thor_list=list(c(5,10),c(5,10)),
#'                                              object=jmfhc_est,
#'                                              test_dat=test_dat,
#'                                              no_cores=7)
#' predict_surv_jmfhc$`landmark=10`$`thor=10`
#' predict_surv_jmfhc$`landmark=20`$`thor=10`
#'
#' @import parallel
#' @import foreach
#' @import doSNOW
#' @import cubature
pred_surv_joint_model <- function(landmark_surv_vec,
                                  thor_list,
                                  object,
                                  test_dat,
                                  no_cores=7){

  # variables in the data set
  event_time <- object$setting$event_time
  measure_time <- object$setting$fu_time_original
  id <- object$setting$id
  fu_measure <- object$setting$fu_measure
  fu_time_fixed_variable <- object$setting$fu_time_fixed_variable
  fu_time_random_variable <- object$setting$fu_time_random_variable
  baseline_var_lmm <- object$setting$baseline_var_lmm
  beta_variable <- object$setting$beta_variable
  gamma_variable <- object$setting$gamma_variable
  # the covariance matrix of the random effects
  Sigma <- object$re_cov
  # number of random effects in the analysis
  length_random_var <- nrow(Sigma)

  base_dat <- object$dat_baseline
  base_dat <- base_dat[order(base_dat[,event_time]),]
  est <- object$coef

  compute_surv <- function(){

    # numerator of the computation for conditional survival probability
    ind_num_surv <- function(x) {

      # density function of random effects
      X <- matrix(x,ncol=1)
      Q <- (-1/2)*t(X)%*%solve(Sigma)%*%(X)
      f <- exp(Q)

      # G(b_i(L))
      fixed_effects <- matrix(est[,which(grepl("fixed_",colnames(est)))],ncol=1)
      if (is.null(baseline_var_lmm)){
        ind_dat$true_Y <- matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_fixed_variable]),ncol=nrow(fixed_effects))%*%fixed_effects+matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_random_variable]),ncol=length_random_var)%*%X
      } else{ind_dat$true_Y <- matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_fixed_variable],unlist(ind_dat[,baseline_var_lmm])),ncol=nrow(fixed_effects))%*%fixed_effects+matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_random_variable]),ncol=length_random_var)%*%X}
      ind_dat$diff_Y <- ind_dat[,fu_measure]-ind_dat$true_Y
      G_b <- prod(exp(-ind_dat$diff_Y^2/(2*est[1,"error_sd"])))

      # Pr(T>=L)
      if (is.null(gamma_variable)){
        survival_thorL <- exp(-exp(unlist(c(1,ind_dat[1,beta_variable],X))%*%(est[,which(grepl("beta_",colnames(est)))]))*
                                (F0L_hor))
      }else{
        survival_thorL <- exp(-exp(unlist(c(1,ind_dat[1,beta_variable],X))%*%(est[,which(grepl("beta_",colnames(est)))]))*
                                (F0L_hor)^(exp(unlist(ind_dat[1,gamma_variable])%*%(est[,which(grepl("gamma_",colnames(est)))]))))
      }

      return(G_b*survival_thorL*f)
    }

    # denominator of the computation for the estimated survival rate
    ind_deno_surv <- function(x) {

      # density function of random effects
      X <- matrix(x,ncol=1)
      Q <- (-1/2)*t(X)%*%solve(Sigma)%*%(X)
      f <- exp(Q)

      # G(b_i(L))
      fixed_effects <- matrix(est[,which(grepl("fixed_",colnames(est)))],ncol=1)
      if (is.null(baseline_var_lmm)){
        ind_dat$true_Y <- matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_fixed_variable]),ncol=nrow(fixed_effects))%*%fixed_effects+matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_random_variable]),ncol=length_random_var)%*%X
      } else{ind_dat$true_Y <- matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_fixed_variable],unlist(ind_dat[,baseline_var_lmm])),ncol=nrow(fixed_effects))%*%fixed_effects+matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_random_variable]),ncol=length_random_var)%*%X}
      ind_dat$diff_Y <- ind_dat[,fu_measure]-ind_dat$true_Y
      G_b <- prod(exp(-ind_dat$diff_Y^2/(2*est[1,"error_sd"])))

      # Pr(T>=L)
      if (is.null(gamma_variable)){
        survival_L <- exp(-exp(unlist(c(1,ind_dat[1,beta_variable],X))%*%(est[,which(grepl("beta_",colnames(est)))]))*
                            (F0L))
      }else{
        survival_L <- exp(-exp(unlist(c(1,ind_dat[1,beta_variable],X))%*%(est[,which(grepl("beta_",colnames(est)))]))*
                            (F0L)^(exp(unlist(ind_dat[1,gamma_variable])%*%(est[,which(grepl("gamma_",colnames(est)))]))))
      }

      return(G_b*survival_L*f)
    }

    evaluate <- function(){
      num <- tryCatch({
        num <- withTimeout({
          arg=ind_num_surv
          bound <- 0.1
          result_1 <- adaptIntegrate(arg, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
          bound <- 0.2
          result_2 <- adaptIntegrate(arg, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
          while ((abs(result_2-result_1)/result_1 > 1e-2|result_1==0) & bound < 6){
            bound <- bound+0.1
            result_1 <- result_2
            result_2 <- adaptIntegrate(arg, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
          }
          value <- result_2
          return(value)
        }, timeout = 60*1)
      }, TimeoutException = function(ex) {
        return(result_1)
      })

      deno <- tryCatch({
        deno <- withTimeout({
          arg=ind_deno_surv
          bound <- 0.1
          result_1 <- adaptIntegrate(arg, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
          bound <- 0.2
          result_2 <- adaptIntegrate(arg, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
          while ((abs(result_2-result_1)/result_1 > 1e-2|result_1==0) & bound < 6){
            bound <- bound+0.1
            result_1 <- result_2
            result_2 <- adaptIntegrate(arg, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
          }
          value <- result_2
          return(value)
        }, timeout = 60*1)
      }, TimeoutException = function(ex) {
        return(result_1)
      })

      return(c(num,deno))
    }

    num <- pcubature(ind_num_surv, lowerLimit = c(-Inf, -Inf), upperLimit = c(+Inf,+Inf))$integral
    deno <-pcubature(ind_deno_surv, lowerLimit = c(-Inf, -Inf), upperLimit = c(+Inf,+Inf))$integral
    comment <- "pcubature"

    if (is.na(num)|is.na(deno)){
      comment <- "NA exists. hcubature"
      num <- hcubature(ind_num_surv, lowerLimit = c(-Inf, -Inf), upperLimit = c(+Inf,+Inf))$integral
      deno <-hcubature(ind_deno_surv, lowerLimit = c(-Inf, -Inf), upperLimit = c(+Inf,+Inf))$integral

      if (num/deno > 1){

        eval_result <- evaluate()
        num <- eval_result[1]
        deno <- eval_result[2]

        bound <- 3
        while ((num/deno>1 & (num/deno-1)>1e-5) & bound >0.1){
          num <- adaptIntegrate(ind_num_surv, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
          deno <- adaptIntegrate(ind_deno_surv, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
          bound <- bound - 0.1
        }
      }
    }else if (num==0 | deno==0){
      comment <- "0 exists."
      eval_result <- evaluate()
      num <- eval_result[1]
      deno <- eval_result[2]

      bound <- 3
      while ((num/deno>1 & (num/deno-1)>1e-5) & bound >0.1){
        num <- adaptIntegrate(ind_num_surv, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
        deno <- adaptIntegrate(ind_deno_surv, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
        bound <- bound - 0.1
      }
    } else if (num/deno>1 & (num/deno-1)>1e-5){
      comment <- "greater than 1."
      eval_result <- evaluate()
      num <- eval_result[1]
      deno <- eval_result[2]

      bound <- 3
      while ((num/deno>1 & (num/deno-1)>1e-5 ) & bound >0.1){
        num <- adaptIntegrate(ind_num_surv, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
        deno <- adaptIntegrate(ind_deno_surv, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
        bound <- bound - 0.1
      }
    }

    return(num/deno)
    #return(c(num/deno,comment))
  }

  result <- foreach(L=landmark_surv_vec)%do%{

    # the data set containing patients who are still at risk at the landmark time
    landmark_dat <- test_dat[test_dat[,event_time]>= L & test_dat[,measure_time] <= L,]
    # obtain F_0(L) and F_0(L+t_hor)
    F0L <- ifelse(sum(base_dat[,event_time]<=L)!=0,base_dat$base_cdf[sum(base_dat[,event_time]<=L)],0)

    pred_surv_L_thor <- foreach (t_hor=thor_list[[which(landmark_surv_vec==L)]])%do%{
      F0L_hor <- ifelse(sum(base_dat[,event_time]<=(L+t_hor))!=0,base_dat$base_cdf[sum(base_dat[,event_time]<=(L+t_hor))],0)

      cl <- makeCluster(no_cores)
      registerDoSNOW(cl)
      ind_survival <- foreach(i=unique(landmark_dat[,id]),.combine="rbind",.packages=c("cubature","R.utils"))%dopar%{
        ind_dat <- landmark_dat[landmark_dat[,id]==i,]
        compute <- compute_surv()
        return(c(i,compute))
      }
      stopCluster(cl)
      if (unique(landmark_dat[,id])==1){
        ind_survival <- matrix(ind_survival,ncol=2,nrow=1)
      }
      ind_survival <- data.frame(id=ind_survival[,1],con_surv=ind_survival[,2])
      return(ind_survival)
    }
    pred_surv_L_thor <- append(pred_surv_L_thor,list(landmark_dat))
    names(pred_surv_L_thor) <- c(paste0("thor=",thor_list[[which(landmark_surv_vec==L)]]),"test_data_landmark")
    return(pred_surv_L_thor)
  }

  result <- append(result,list(object))
  names(result) <- c(paste0("landmark=",landmark_surv_vec),"train_est")
  return(result)

}

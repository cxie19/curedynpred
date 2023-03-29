#' Predict conditional probability of being cured
#' @description Predicts an individual's the conditional probability of being cured at landmark time \eqn{L}
#'
#' @param L a value of the landmark time in the same unit as the observed survival times and the original form of measurement times
#' in the data set.
#' @param predict.id a character for the type of prediction. It could be "all","one", or "new".
#' predict.id="all" is to do prediction for all patients who are still at risk at time \eqn{L} in the data set.
#' predict.id="one" is to do prediction for a specific patient who is still at risk at landmark time \eqn{L} in the data set, where the
#' patient id is specified in the parameter predict.id.one.
#' predict.id="new" is to do prediction for a new patient who is not in the data set and at risk at landmark time \eqn{L}. This new
#' patient's new values are specified in the parameters new.fu_measure, new.fu_time_variable, new.baseline_value_lmm, new.z_value and new.x_value.
#' @param predict.id.one a patient's id number from the data set when predict.id="one",
#' which could be numeric and character depending on the form of patient id column in the data set. By default predict.id.one = NULL.
#' @param new.fu_measure a vector of the new patient's biomarker measurements up to time \eqn{L} when predict.id="new". All measurements are
#' in the same form specified for the parameter fu_measure in the fitted longitudinal submodel of jmfhc_point_est().
#' By default new.fu_measure = NULL.
#' @param new.fu_time_variable a vector of the new patient's measurement time points in the forms specified for the parameter
#' fu_time_variable in the fitted longitudinal submodel of jmfhc_point_est() when predict.id="new". When fu_time_variable is more than one variable,
#' the elements in new.fu_time_variable is placed in the order of variables in the fu_time_variable.
#' By default new.fu_time_variable = NULL.
#' @param new.baseline_value_lmm value(s) of the baseline covariate(s) in the order specified for the parameter baseline_var_lmm in the fitted
#' longitudinal submodel of jmfhc_point_est() when predict.id="new". When the baseline_var_lmm=NULL in jmfhc_point_est(), new.baseline_value_lmm
#' can be set as NULL. By default new.baseline_value_lmm = NULL.
#' @param new.z_value value(s) of the long-term baseline covariate(s) in the order specified for the parameter beta_variable in the fitted
#' cure submodel of jmfhc_point_est() when predict.id="new". By default new.z_value = NULL.
#' @param new.x_value value(s) of the short-term baseline covariate(s) in the order specified for the parameter gamma_variable in the fitted
#' cure submodel of jmfhc_point_est() when predict.id="new". When the gamma_variable=NULL in jmfhc_point_est(), new.x_value can be set as NULL.
#' By default new.x_value = NULL.
#' @param object an object of jmfhc_point_est() function.
#' @param no_cores the number of cores used during the estimation. The default is 7.
#'
#' @return a data.frame with id for the predicted patients who are still at risk at time \eqn{L} and predicted conditional probabilities of being
#' cured at time \eqn{L}.
#' @export
#'
#' @examples
#' result_jmfhc <- jmfhc::jmfhc_point_est(data=longdat,
#'                                event_time="event.time", event_status="event",
#'                                id="patient.id",beta_variable="trt", gamma_variable="trt",
#'                                fu_measure="measure", fu_time_original="mes.times",
#'                                fu_time_variable="mes.times")
#' predict_cure <- est_cure_L(L=10,predict.id="all",object=result_jmfhc)
#' predict_cure <- est_cure_L(L=10,predict.id="one",predict.id.one=3,object=result_jmfhc)
#' predict_cure <- est_cure_L(L=10,predict.id="new",
#'                       new.fu_measure=c(5.2,5.1,1.6,0.9,-1.1,-2.6,-5.3,-8.0,-7.5,-11,-12),
#'                       new.fu_time_variable= 0:10,
#'                       new.baseline_value_lmm=NULL,
#'                       new.z_value=0,
#'                       new.x_value=0,
#'                       object=result_jmfhc)
#' predict_cure
#'
#' @import jmfhc
#' @import dplyr
#' @import parallel
#' @import foreach
#' @import doSNOW

est_cure_L <- function(L,
                       predict.id,
                       predict.id.one=NULL,
                       new.fu_measure=NULL,
                       new.fu_time_variable=NULL,
                       new.baseline_value_lmm=NULL,
                       new.z_value=NULL,
                       new.x_value=NULL,
                       object,
                       no_cores=7){

  # variables in the data set
  event_time <- object$setting$event_time
  measure_time <- object$setting$fu_time_original
  id <- object$setting$id
  fu_measure <- object$setting$fu_measure
  fu_time_variable <- object$setting$fu_time_variable
  baseline_var_lmm <- object$setting$baseline_var_lmm
  beta_variable <- object$setting$beta_variable
  gamma_variable <- object$setting$gamma_variable

  updated_dat <- object$dat_long
  base_dat <- object$dat_baseline
  est <- object$coef

  # the data set containing patients who are still at risk at the landmark time
  landmark_dat <- updated_dat[updated_dat[,event_time]>= L & updated_dat[,measure_time] <= L,]

  if(predict.id=="one"){
    if(!predict.id.one %in% unique(landmark_dat[,id])){
      stop("The entered patient id for prediction is not in the set of patients who are still at risk at the prespecified landmark time")
    }
  }

  # the covariance matrix of the random effects
  Sigma <- object$re_cov

  # number of random effects in the analysis
  length_random_var <- nrow(Sigma)

  # obtain F_0(L)
  base_dat <- base_dat[order(base_dat[,event_time]),]
  F0L <- ifelse(sum(base_dat[,event_time]<=L)!=0,base_dat$base_cdf[sum(base_dat[,event_time]<=L)],0)

  compute_cure <- function(){

    # numerator of the computation for the estimated cure rate for a subject at a pre-specified time point
    ind_num_cure <- function(x) {

      # density function of random effects
      X <- matrix(x,ncol=1)
      Q <- (-1/2)*t(X)%*%solve(Sigma)%*%(X)
      f <- exp(Q)
      # f <- (1/(2*pi))^(length_random_var/2)*((det(Simga))^(-1/2))*exp(Q)

      # G(b_i(L))
      fixed_effects <- matrix(est[,which(grepl("fixed_",colnames(est)))],ncol=1)
      if (is.null(baseline_var_lmm)){
        ind_dat$true_Y <- matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_variable]),ncol=nrow(fixed_effects))%*%fixed_effects+matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_variable]),ncol=length_random_var)%*%X
      } else{ind_dat$true_Y <- matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_variable],ind_dat[,baseline_var_lmm]),ncol=nrow(fixed_effects))%*%fixed_effects+matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_variable]),ncol=length_random_var)%*%X}
      ind_dat$diff_Y <- ind_dat[,fu_measure]-ind_dat$true_Y
      G_b <- prod(exp(-ind_dat$diff_Y^2/(2*est[1,"error_sd"])))

      # Pr(cure, T>=L)
      cure_base <- exp(-exp(matrix(c(1,ind_dat[1,beta_variable],X),nrow=1)%*%matrix(est[,which(grepl("beta_",colnames(est)))],ncol=1)))

      return(G_b*cure_base*f)
    }

    # denominator of the computation for the estimated cure rate for a subject at a pre-specified time point
    ind_deno_cure <- function(x) {

      # density function of random effects
      X <- matrix(x,ncol=1)
      Q <- (-1/2)*t(X)%*%solve(Sigma)%*%(X)
      f <- exp(Q)

      # G(b_i(L))
      fixed_effects <- matrix(est[,which(grepl("fixed_",colnames(est)))],ncol=1)
      if (is.null(baseline_var_lmm)){
        ind_dat$true_Y <- matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_variable]),ncol=nrow(fixed_effects))%*%fixed_effects+matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_variable]),ncol=length_random_var)%*%X
      } else{ind_dat$true_Y <- matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_variable],ind_dat[,baseline_var_lmm]),ncol=nrow(fixed_effects))%*%fixed_effects+matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_variable]),ncol=length_random_var)%*%X}
      ind_dat$diff_Y <- ind_dat[,fu_measure]-ind_dat$true_Y
      G_b <- prod(exp(-ind_dat$diff_Y^2/(2*est[1,"error_sd"])))

      # Pr(T>=L)
      if (is.null(gamma_variable)){
        survival <- exp(-exp(matrix(c(1,ind_dat[1,beta_variable],X),nrow=1)%*%matrix(est[,which(grepl("beta_",colnames(est)))],ncol=1))*
                          (F0L))
      }else{
        survival <- exp(-exp(matrix(c(1,ind_dat[1,beta_variable],X),nrow=1)%*%matrix(est[,which(grepl("beta_",colnames(est)))],ncol=1))*
                          (F0L)^(exp(matrix(ind_dat[1,gamma_variable],nrow=1)%*%matrix(est[,which(grepl("gamma_",colnames(est)))],ncol=1))))
      }

      return(G_b*survival*f)
    }

    evaluate <- function(){
      num <- tryCatch({
        num <- withTimeout({
          arg=ind_num_cure
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
          arg=ind_deno_cure
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

    num <- pcubature(ind_num_cure, lowerLimit = c(-Inf, -Inf), upperLimit = c(+Inf,+Inf))$integral
    deno <-pcubature(ind_deno_cure, lowerLimit = c(-Inf, -Inf), upperLimit = c(+Inf,+Inf))$integral
    comment <- "pcubature"

    if (is.na(num)|is.na(deno)){
      comment <- "NA exists. hcubature"
      num <- hcubature(ind_num_cure, lowerLimit = c(-Inf, -Inf), upperLimit = c(+Inf,+Inf))$integral
      deno <-hcubature(ind_deno_cure, lowerLimit = c(-Inf, -Inf), upperLimit = c(+Inf,+Inf))$integral

      if (num/deno > 1){

        eval_result <- evaluate()
        num <- eval_result[1]
        deno <- eval_result[2]

        bound <- 3
        while ((num/deno>1 & (num/deno-1)>1e-5) & bound >0.1){
          num <- adaptIntegrate(ind_num_cure, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
          deno <- adaptIntegrate(ind_deno_cure, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
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
        num <- adaptIntegrate(ind_num_cure, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
        deno <- adaptIntegrate(ind_deno_cure, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
        bound <- bound - 0.1
      }
    } else if (num/deno>1 & (num/deno-1)>1e-5){
      comment <- "greater than 1."
      eval_result <- evaluate()
      num <- eval_result[1]
      deno <- eval_result[2]

      bound <- 3
      while ((num/deno>1 & (num/deno-1)>1e-5 ) & bound >0.1){
        num <- adaptIntegrate(ind_num_cure, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
        deno <- adaptIntegrate(ind_deno_cure, lowerLimit= c(-bound,-bound),  upperLimit=c(bound,bound))$integral
        bound <- bound - 0.1
      }
    }

    return(num/deno)
    # return(c(num/deno,comment))

  }

  if (predict.id=="all"){
    cl <- makeCluster(no_cores)
    registerDoSNOW(cl)
    cure <- foreach(i=unique(landmark_dat[,id]),.combine="rbind",.packages=c("cubature","R.utils"))%dopar%{
      ind_dat <- landmark_dat[landmark_dat[,id]==i,]
      compute <- compute_cure()
      return(c(i,compute))
    }
    stopCluster(cl)
    all_cure <- data.frame(id=cure[,1],con_cure=cure[,2])
    return(all_cure)
  }else if(predict.id=="one"){
    ind_dat <- landmark_dat[landmark_dat[,id]==predict.id.one,]
    compute <- compute_cure()
    one_cure <- data.frame(id=predict.id.one,con_cure=compute)
    return(one_cure)
  } else if (predict.id=="new"){
    new.fu.matrix <- matrix(c(new.fu_measure,new.fu_time_variable),ncol=1+length(fu_time_variable))
    new.baseline.matrix <- matrix(c(new.baseline_value_lmm,new.z_value,new.x_value),ncol=length(c(new.baseline_value_lmm,new.z_value,new.x_value)),nrow=nrow(new.fu.matrix))
    new.variable <- cbind(new.fu.matrix,new.baseline.matrix)
    colnames(new.variable) <- c(fu_measure,fu_time_variable,baseline_var_lmm,beta_variable,gamma_variable)
    if(sum(duplicated(colnames(new.variable)))!=0){
      ind_dat <- as.data.frame(new.variable[,-which(duplicated(colnames(new.variable)))])
    }else{
      ind_dat <-as.data.frame(new.variable)
    }
    compute <- compute_cure()
    new_cure <- data.frame(id="new",con_cure=compute)
    return(new_cure)
  }
}

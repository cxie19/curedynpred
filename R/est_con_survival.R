#' Predict conditional survival probability
#' @description Predicts an individual's conditional probability of not experiencing the event of
#' interest in an additional time \eqn{t_{hor}} given that the patient remains risk-free at least until
#' time \eqn{L} and evaluates the fitted joint cure model's predictive performance by the time-dependent
#' area under the curve (AUC) of the receiver operating characteristic (ROC) and Brier score.
#'
#' @param L a value of the landmark time in the same unit as the observed survival times and the original form of measurement times
#' in the data set.
#' @param t_hor a value of the time horizon of prediction in the same unit as the parameter L.
#' The summation of L and t_hor needs to be less than the maximum of observed survival times in the data set.
#' @param predict.id a character for the type of prediction. It could be "all","one", or "new".
#' predict.id="all" is to do prediction for all patients who are still at risk at time \eqn{L} in the data set.
#' predict.id="one" is to do prediction for a specific patient who is still at risk at landmark time \eqn{L} in the data set, where the
#' patient id is specified in the parameter predict.id.one.
#' predict.id="new" is to do prediction for a new patient who is not in the data set and at risk at landmark time \eqn{L}. This new
#' patient's new values are specified in the parameters new.fu_measure, new.fu_time_fixed_variable, new.fu_time_random_variable,
#' new.baseline_value_lmm, new.z_value and new.x_value.
#' @param predict.id.one a patient's id number from the data set when predict.id="one",
#' which could be numeric and character depending on the form of patient id column in the data set. By default predict.id.one = NULL.
#' @param new.fu_measure a vector of the new patient's biomarker measurements up to time \eqn{L} when predict.id="new". All measurements are
#' in the same form specified for the parameter fu_measure in the fitted longitudinal submodel of jmfhc_point_est().
#' By default new.fu_measure = NULL.
#' @param new.fu_time_fixed_variable a vector of the new patient's measurement time points in the forms specified for the parameter
#' fu_time_fixed_variable in the fitted longitudinal submodel of jmfhc_point_est() when predict.id="new". When fu_time_fixed_variable is more than one variable,
#' the elements in new.fu_time_fixed_variable is placed in the order of variables in the fu_time_fixed_variable.
#' By default new.fu_time_fixed_variable = NULL.
#' @param new.fu_time_random_variable a vector of the new patient's measurement time points in the forms specified for the parameter
#' fu_time_random_variable in the fitted longitudinal submodel of jmfhc_point_est() when predict.id="new". When fu_time_random_variable is more than one variable,
#' the elements in new.fu_time_random_variable is placed in the order of variables in the fu_time_random_variable.
#' By default new.fu_time_random_variable = NULL.
#' @param new.baseline_value_lmm value(s) of the baseline covariate(s) in the order specified for the parameter baseline_var_lmm in the fitted
#' longitudinal submodel of jmfhc_point_est() when predict.id="new". When the baseline_var_lmm=NULL in jmfhc_point_est(), new.baseline_value_lmm
#' can be set as NULL. By default new.baseline_value_lmm = NULL.
#' @param new.z_value value(s) of the long-term baseline covariate(s) in the order specified for the parameter beta_variable in the fitted
#' cure submodel of jmfhc_point_est() when predict.id="new". By default new.z_value = NULL.
#' @param new.x_value value(s) of the short-term baseline covariate(s) in the order specified for the parameter gamma_variable in the fitted
#' cure submodel of jmfhc_point_est() when predict.id="new". When the gamma_variable=NULL in jmfhc_point_est(), new.x_value can be set as NULL.
#' By default new.x_value = NULL.
#' @param AUC a logical value (True or False) to determine whether to evaluate the AUC. By default AUC=FALSE.
#' @param Brier a logical value (True or False) to determine whether to evaluate the Brier score. By default Brier=FALSE.
#' @param object an object of jmfhc_point_est() function.
#' @param no_cores the number of cores used during the estimation. The default is 7.
#'
#' @return a list containing the prediction results depending on the arguments for the parameters predict.id, AUC, and Brier.
#'         The following items are returned.
#' \item{con_surv} a data.frame with a column of patient id and a column of predicted conditional probabilities of not experiencing the event of
#' interest in an additional time \eqn{t_{hor}} at time \eqn{L}. When predict.id="all", the data.frame contains predicted probabilities for
#' all patients who are still at risk at time \eqn{L}. When predict.id="one", the data.frame contains the specified patient's predicted probability.
#' When predict.id="new", the data.frame contains this new patient's predicted probability.
#' \item{AUC} the computed time-independent AUC value if AUC=TRUE.
#' \item{Brier_score} the computed Brier score if Brier=TRUE.
#' \item{Landmark} the value of landmark time specified in the function argument L.
#' \item{horizon} the value of time horizon of prediction specified in the function argument t_hor.
#' @export
#'
#' @examples
#' result_jmfhc <- jmfhc::jmfhc_point_est(data=jmfhc_dat, event_time="event.time", event_status="event",
#'                                        id="patient.id", beta_variable="trt", gamma_variable="trt",
#'                                        fu_measure_original="measure",fu_measure="measure",
#'                                        fu_time_original="mes.times",fu_time_fixed_variable="mes.times",
#'                                        fu_time_random_variable="mes.times")
#' predict_surv <- est_con_survival(L=10,t_hor=5,predict.id="all",AUC=TRUE,Brier=TRUE,object=result_jmfhc)
#' predict_surv <- est_con_survival(L=10,t_hor=5,predict.id="one",predict.id.one=3,AUC=FALSE,Brier=FALSE,object=result_jmfhc)
#' predict_surv <- est_con_survival(L=10,t_hor=5,predict.id="new",
#'                       new.fu_measure=c(5.2,5.1,1.6,0.9,-1.1,-2.6,-5.3,-8.0,-7.5,-11,-12),
#'                       new.fu_time_fixed_variable= 0:10,
#'                       new.fu_time_random_variable= 0:10,
#'                       new.baseline_value_lmm=NULL,
#'                       new.z_value=0,
#'                       new.x_value=0,
#'                       AUC=FALSE,Brier=FALSE,
#'                       object=result_jmfhc)
#' predict_surv
#'
#' @import jmfhc
#' @import cubature
#' @import parallel
#' @import foreach
#' @import doSNOW
#' @import R.utils
#' @import dplyr
#' @import tdROC
#' @import survival

est_con_survival <- function(L,
                             t_hor,
                             predict.id,
                             predict.id.one,
                             new.fu_measure=NULL,
                             new.fu_time_fixed_variable=NULL,
                             new.fu_time_random_variable=NULL,
                             new.baseline_value_lmm=NULL,
                             new.z_value=NULL,
                             new.x_value=NULL,
                             AUC=FALSE,
                             Brier=FALSE,
                             object,
                             no_cores=7){

  # variables in the data set
  event_time <- object$setting$event_time
  event_status <- object$setting$event_status
  measure_time <- object$setting$fu_time_original
  id <- object$setting$id
  fu_measure <- object$setting$fu_measure
  fu_time_fixed_variable <- object$setting$fu_time_fixed_variable
  fu_time_random_variable <- object$setting$fu_time_random_variable
  baseline_var_lmm <- object$setting$baseline_var_lmm
  beta_variable <- object$setting$beta_variable
  gamma_variable <- object$setting$gamma_variable

  updated_dat <- object$dat_long
  base_dat <- object$dat_baseline
  est <- object$coef

  # the data set containing patients who are still at risk at the landmark time
  landmark_dat <- updated_dat[updated_dat[,event_time]>= L & updated_dat[,measure_time] <= L,]
  landmark_dat$residual <- landmark_dat[,event_time]-L

  if(predict.id=="one"){
    if(!predict.id.one %in% unique(landmark_dat[,id])){
      stop("The entered patient id for prediction is not in the set of patients who are still at risk at the prespecified landmark time")
    }
  }

  # the covariance matrix of the random effects
  Sigma <- object$re_cov

  # number of random effects in the analysis
  length_random_var <- nrow(Sigma)

  # obtain F_0(L) and F_0(L+t_hor)
  base_dat <- base_dat[order(base_dat[,event_time]),]
  F0L <- ifelse(sum(base_dat[,event_time]<=L)!=0,base_dat$base_cdf[sum(base_dat[,event_time]<=L)],0)
  F0L_hor <- ifelse(sum(base_dat[,event_time]<=(L+t_hor))!=0,base_dat$base_cdf[sum(base_dat[,event_time]<=(L+t_hor))],0)

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
      } else{ind_dat$true_Y <- matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_fixed_variable],ind_dat[,baseline_var_lmm]),ncol=nrow(fixed_effects))%*%fixed_effects+matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_random_variable]),ncol=length_random_var)%*%X}
      ind_dat$diff_Y <- ind_dat[,fu_measure]-ind_dat$true_Y
      G_b <- prod(exp(-ind_dat$diff_Y^2/(2*est[1,"error_sd"])))

      # Pr(T>=L)
      if (is.null(gamma_variable)){
        survival_thorL <- exp(-exp(matrix(c(1,ind_dat[1,beta_variable],X),nrow=1)%*%matrix(est[,which(grepl("beta_",colnames(est)))],ncol=1))*
                                (F0L_hor))
      }else{
        survival_thorL <- exp(-exp(matrix(c(1,ind_dat[1,beta_variable],X),nrow=1)%*%matrix(est[,which(grepl("beta_",colnames(est)))],ncol=1))*
                               (F0L_hor)^(exp(matrix(ind_dat[1,gamma_variable],nrow=1)%*%matrix(est[,which(grepl("gamma_",colnames(est)))],ncol=1))))
      }

      return(G_b*survival_thorL*f)
    }

    # denominator of the computation for the estimated cure rate for a subject at a pre-specified time point
    ind_deno_surv <- function(x) {

      # density function of random effects
      X <- matrix(x,ncol=1)
      Q <- (-1/2)*t(X)%*%solve(Sigma)%*%(X)
      f <- exp(Q)

      # G(b_i(L))
      fixed_effects <- matrix(est[,which(grepl("fixed_",colnames(est)))],ncol=1)
      if (is.null(baseline_var_lmm)){
        ind_dat$true_Y <- matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_fixed_variable]),ncol=nrow(fixed_effects))%*%fixed_effects+matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_random_variable]),ncol=length_random_var)%*%X
      } else{ind_dat$true_Y <- matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_fixed_variable],ind_dat[,baseline_var_lmm]),ncol=nrow(fixed_effects))%*%fixed_effects+matrix(c(rep(1,nrow(ind_dat)),ind_dat[,fu_time_random_variable]),ncol=length_random_var)%*%X}
      ind_dat$diff_Y <- ind_dat[,fu_measure]-ind_dat$true_Y
      G_b <- prod(exp(-ind_dat$diff_Y^2/(2*est[1,"error_sd"])))

      # Pr(T>=L)
      if (is.null(gamma_variable)){
        survival_L <- exp(-exp(matrix(c(1,ind_dat[1,beta_variable],X),nrow=1)%*%matrix(est[,which(grepl("beta_",colnames(est)))],ncol=1))*
                            (F0L))
      }else{
        survival_L <- exp(-exp(matrix(c(1,ind_dat[1,beta_variable],X),nrow=1)%*%matrix(est[,which(grepl("beta_",colnames(est)))],ncol=1))*
                            (F0L)^(exp(matrix(ind_dat[1,gamma_variable],nrow=1)%*%matrix(est[,which(grepl("gamma_",colnames(est)))],ncol=1))))
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

  if (AUC==TRUE|Brier==TRUE|predict.id=="all"){

    # compute all at-risk patients' conditional survival probabilities
    cl <- makeCluster(no_cores)
    registerDoSNOW(cl)
    ind_survival <- foreach(i=unique(landmark_dat[,id]),.combine="rbind",.packages=c("cubature","R.utils"))%dopar%{
      print(i)
      ind_dat <- landmark_dat[landmark_dat[,id]==i,]
      compute <- compute_surv()
      print(compute)
      return(c(i,compute))
    }
    stopCluster(cl)
    ind_survival <- data.frame(id=ind_survival[,1],con_surv=ind_survival[,2])

    colnames(landmark_dat)[colnames(landmark_dat)==id] <- "id"
    landmark_base <-  landmark_dat %>%
      group_by(id) %>%
      arrange(measure_time) %>%
      slice_head(n = 1)
    for_evaluate <- merge(ind_survival,landmark_base,by="id")

    if (AUC==TRUE){
      fm <- tdROC( X = 1-for_evaluate$con_surv, Y = for_evaluate$residual, delta = for_evaluate[,event_status], tau = t_hor, span = 0.1, nboot = 0)
      auc_value <- fm$AUC[1,1]
    }

    if (Brier==TRUE){
      # Kaplan-Meier of residual censoring time
      for_evaluate <- for_evaluate[order(for_evaluate[,event_time]),]
      for_evaluate$censor <- ifelse(for_evaluate[,event_status]==1,0,1)
      km_censor <- survfit(Surv(residual, censor) ~ 1, data=for_evaluate)
      S_thor <- ifelse(sum(km_censor$time<=t_hor)==0,1,km_censor$surv[sum(km_censor$time<=t_hor)])
      for_evaluate$S_res <- sapply(seq(nrow(for_evaluate)),function(x)
        km_censor$surv[sum(km_censor$time<=for_evaluate$residual[x])])
      for_evaluate$S_res[for_evaluate$S_res==0]<- Inf

      # category 1 (y_i-L <= t_hor and event=1)
      for_evaluate$cat1 <- for_evaluate$residual<=t_hor & for_evaluate[,event_status]==1
      # category 2 (y_i-L > t_hor)
      for_evaluate$cat2 <- for_evaluate$residual>t_hor

      ind_brier <-(((0-for_evaluate$con_surv)^2*for_evaluate$cat1)/for_evaluate$S_res) + (((1-for_evaluate$con_surv)^2*for_evaluate$cat2)/S_thor)
      brier_score <- mean(ind_brier)
    }

    if (predict.id=="all"){
      if (AUC==FALSE & Brier==FALSE){
        return(list(con_surv=ind_survival,Landmark=L, horizon=t_hor))
      }else if (AUC==TRUE & Brier==FALSE){
        return(list(con_surv=ind_survival,AUC=auc_value,Landmark=L, horizon=t_hor))
      }else if (AUC==FALSE & Brier==TRUE){
        return(list(con_surv=ind_survival,Brier_score=brier_score,Landmark=L, horizon=t_hor))
      }else if (AUC==TRUE & Brier==TRUE){
        return(list(con_surv=ind_survival,AUC=auc_value,Brier_score=brier_score,Landmark=L, horizon=t_hor))
      }
    }
  }

  if(predict.id=="one"){
    if (!exists("ind_survival")){
      ind_dat <- landmark_dat[landmark_dat[,id]==predict.id.one,]
      compute <- compute_surv()
    }else{
      compute <- ind_survival[ind_survival[,"id"]==predict.id.one,"con_surv"]
    }

    one_survival <- data.frame(id=predict.id.one,con_surv=compute)

    if (AUC==FALSE & Brier==FALSE){
      return(list(con_surv=one_survival,Landmark=L, horizon=t_hor))
    }else if (AUC==TRUE & Brier==FALSE){
      return(list(con_surv=one_survival,AUC=auc_value,Landmark=L, horizon=t_hor))
    }else if (AUC==FALSE & Brier==TRUE){
      return(list(con_surv=one_survival,Brier_score=brier_score,Landmark=L, horizon=t_hor))
    }else if (AUC==TRUE & Brier==TRUE){
      return(list(con_surv=one_survival,AUC=auc_value,Brier_score=brier_score,Landmark=L, horizon=t_hor))
    }

  } else if (predict.id=="new"){
    new.fu.matrix <- matrix(c(new.fu_measure,new.fu_time_fixed_variable,new.fu_time_random_variable),ncol=1+length(fu_time_fixed_variable)+length(fu_time_random_variable))
    new.baseline.matrix <- matrix(c(new.baseline_value_lmm,new.z_value,new.x_value),ncol=length(c(new.baseline_value_lmm,new.z_value,new.x_value)),nrow=nrow(new.fu.matrix))
    new.variable <- cbind(new.fu.matrix,new.baseline.matrix)
    colnames(new.variable) <- c(fu_measure,fu_time_fixed_variable,fu_time_random_variable,baseline_var_lmm,beta_variable,gamma_variable)
    if(sum(duplicated(colnames(new.variable)))!=0){
      ind_dat <- as.data.frame(new.variable[,-which(duplicated(colnames(new.variable)))])
    }else{
      ind_dat <-as.data.frame(new.variable)
    }
    compute <- compute_surv()
    new_survival <- data.frame(id="new",con_surv=compute)

    if (AUC==FALSE & Brier==FALSE){
      return(list(con_surv=new_survival,Landmark=L, horizon=t_hor))
    }else if (AUC==TRUE & Brier==FALSE){
      return(list(con_surv=new_survival,AUC=auc_value,Landmark=L, horizon=t_hor))
    }else if (AUC==FALSE & Brier==TRUE){
      return(list(con_surv=new_survival,Brier_score=brier_score,Landmark=L, horizon=t_hor))
    }else if (AUC==TRUE & Brier==TRUE){
      return(list(con_surv=new_survival,AUC=auc_value,Brier_score=brier_score,Landmark=L, horizon=t_hor))
    }
  }
}

#' Predict conditional cure or survival probability using a landmark model
#' @description Predicts an individual's conditional probability of being cured at a landmark time \eqn{L} and/or
#' conditional probability of not experiencing the event of interest in an additional time \eqn{t_{hor}}
#' given that the patient remains risk-free at least until time \eqn{L} using
#' a landmark flexible-hazards cure (LFHC) model or a landmark proportional hazards cure (LPHC) model
#' incorporating the longitudinal biomarker values by including baseline and the most current biomarker values as covariates.
#'
#' @param cure_prob a logical value (True or False) to determine whether to predict conditional cure probabilities for subject(s)
#' in the data set for prediction (test_prep_data).
#' @param landmark_cure_vec a vector of landmark times for prediction of cure probabilities in the same unit as the observed survival times and the original form of measurement times
#' in the data sets. By default landmark_cure_vec = NULL.
#' @param surv_prob a logical value (True or False) to determine whether to predict conditional survival probabilities for subject(s)
#' in the data set for prediction (test_prep_data).
#' @param landmark_surv_vec a vector of landmark times for prediction of survival probabilities in the same unit as the observed survival times and the original form of measurement times
#' in the data sets. By default landmark_surv_vec = NULL.
#' @param thor_list a list of vectors of the time horizons of prediction in the same unit as the parameter landmark_surv_vec in the order of the corresponding landmark times.
#' @param train_prep_data a data frame used as a training data for estimation.
#' @param test_prep_data a data frame containing subject(s) for prediction.
#' @param id the variable name corresponding to id in the train_prep_data and test_prep_data.
#' @param event_time the variable name corresponding to the event time in the train_prep_data and test_prep_data.
#' @param event_status the variable name corresponding to the event status in the train_prep_data and test_prep_data.
#' @param fu_time_original the variable name corresponding to the original form of measurement time points in the train_prep_data and test_prep_data.
#' Variables event_time and fu_time_original are in the same time unit.
#' @param fu_measure_original the variable name corresponding to the longitudinal biomarker measurements in the original form.
#' @param fu_measure the variable name corresponding to the longitudinal biomarker measurements in the original form or any transformation used in the model.
#' @param beta_variable the names of the long-term covariates besides the baseline and most current biomarker values.
#' @param gamma_include_measure a logical value (True or False) to determine whether to include the baseline and most current biomarker values as short-term covariates.
#' By default gamma_include_measure = FALSE.
#' @param gamma_variable the names of the short-term covariates besides the baseline and most current biomarker values. By default gamma_variable = NULL.
#' @param no_cores the number of cores used during the prediction. The default is 7.
#'
#' @return a list containing sublists for all assigned landmark times from landmark_cure_vec if cure_prob=TRUE and landmark_surv_vec if surv_prob=TRUE.
#'         Each sublist at each assigned landmark time contains the following items.
#'         \item{train_est} the estimation result of the landmark model containing coefficients (coef), the number of iterations (iter), updated data with
#'         estimated baseline hazard function and baseline cumulative hazard function (data), and a list containing all the specified parameters in the argument for the fitted landmark model (setting).
#'         \item{test_data_landmark} a data frame with risk-free subjects at landmark time.
#'         \item{pred_cure} returns if cure_prob=TRUE. A data frame contains id and its corresponding predicted conditional cure probability.
#'         \item{pred_surv} returns if surv_prob=TRUE. Data frames for time horizons of predictions contain id and
#'         corresponding predicted conditional survival probabilities.
#' @export
#'
#' @examples
#' lfhc_pred_all <- pred_landmark_model(cure_prob=T,landmark_cure_vec=c(10,20),
#'                                      surv_prob=T,landmark_surv_vec=c(10,20),
#'                                      thor_list=list(c(5,10),c(5,10)),
#'                                      train_prep_data=train_dat,test_prep_data=test_dat,
#'                                      id="id",event_time="event.time",event_status = "event",
#'                                      fu_time_original="measure.time",
#'                                      fu_measure_original="measure",fu_measure="measure",
#'                                      beta_variable = c("x1","x2"),
#'                                      gamma_include_measure=FALSE,
#'                                      gamma_variable= c("x1","x2"))
#' # predicted conditional cure probabilities for the test data at landmark 10 and 20
#'pred_cure_10 <- lfhc_pred_all$`landmark=10`$pred_cure
#'pred_cure_20 <- lfhc_pred_all$`landmark=20`$pred_cure
#'
#'# predicted conditional survival probabilities for the test data
#'# landmark 10 with thor=5 and 10
#'pred_surv_L10_thor5 <- lfhc_pred_all$`landmark=10`$pred_surv$`thor=5`
#'pred_surv_L10_thor10 <- lfhc_pred_all$`landmark=10`$pred_surv$`thor=10`
#'# landmark 20 with thor=5 and 10
#'pred_surv_L20_thor5 <- lfhc_pred_all$`landmark=20`$pred_surv$`thor=5`
#'pred_surv_L20_thor10 <- lfhc_pred_all$`landmark=20`$pred_surv$`thor=10`
#'
#' @import parallel
#' @import foreach
#' @import doSNOW
#' @import dplyr
#' @import fhc

pred_landmark_model <- function(cure_prob,landmark_cure_vec=NULL,
                                surv_prob,landmark_surv_vec=NULL,thor_list=NULL,
                                train_prep_data,test_prep_data,
                                id,event_time,event_status,
                                fu_time_original,fu_measure_original,fu_measure,
                                beta_variable,
                                gamma_include_measure=FALSE,gamma_variable=NULL,
                                no_cores=7){

  if(cure_prob & surv_prob){
    landmark_vec <- unique(c(landmark_cure_vec,landmark_surv_vec))
  }else if(cure_prob==TRUE & surv_prob==FALSE){
    landmark_vec <- landmark_cure_vec
  }else if(cure_prob==FALSE & surv_prob==TRUE){
    landmark_vec <- landmark_surv_vec
  }


  beta_variable <- c(beta_variable,"measure","current_measure")
  if(gamma_include_measure){
    gamma_variable <- c(gamma_variable,"measure","current_measure")
  }

  result <- foreach(landmark=landmark_vec)%do%{

    cat(paste("Start to estimate model at time", landmark,"\n"))
    train_dat_base <- prep_data_landmark(dat=train_prep_data,landmark=landmark,id,event_time,fu_time_original,fu_measure)

    # estimation
    result_est_L <- fhcmodel(data=train_dat_base,event_status=event_status,event_time=event_time,id=id,
                             beta_variable=beta_variable,
                             gamma_variable=gamma_variable,
                             se=F)
    setting <-  list(beta_variable=beta_variable,gamma_variable=gamma_variable,
                     gamma_include_measure=gamma_include_measure,
                     event_status=event_status,event_time=event_time,id=id,
                     fu_time_original=fu_time_original,
                     fu_measure_original=fu_measure_original,fu_measure=fu_measure)
    result_est_L <- append(result_est_L,list(setting))
    names(result_est_L)[length(result_est_L)] <- "setting"
    coef <- result_est_L$coef
    train_dat_base <- result_est_L$data
    train_dat_base <- train_dat_base[order(train_dat_base$event.time),]

    test_dat_base <- prep_data_landmark(dat=test_prep_data,landmark=landmark,id,event_time,fu_time_original,fu_measure)

    # predict conditional cure probabilities
    if(cure_prob & landmark%in%landmark_cure_vec){
      cat(paste("Predict conditional cure probabilities","at time",landmark,"\n"))
      cl <- makeCluster(no_cores)
      registerDoSNOW(cl)
      pred_cure_L <- foreach(i=unique(test_dat_base$id),.combine="rbind")%dopar%{

        pred_cure <- function(z_value){
          exp(-exp(coef[1]+coef[2:(1+length(z_value))]%*%t(z_value)))
        }
        ind_dat <- test_dat_base[test_dat_base$id==i,]
        return(c(i,pred_cure(ind_dat[,beta_variable])))
      }
      stopCluster(cl)

      pred_cure_L[,1] <- as.numeric(pred_cure_L[,1])
      pred_cure_L[,2] <- as.numeric(pred_cure_L[,2])
      colnames(pred_cure_L) <- c("id","cure")
    }

    # predict conditional survival probabilities
    if(surv_prob & landmark%in%landmark_surv_vec){
      pred_surv_L_thor <- foreach (t_hor=thor_list[[which(landmark_surv_vec==landmark)]])%do%{
        cat(paste("Predict conditional survival probabilities for another",t_hor,"at time",landmark,"\n"))
        # probability of being at risk for another t_hor at L
        F0t_horL <- ifelse(sum(train_dat_base$event.time<=t_hor)!=0,train_dat_base$base_cdf[sum(train_dat_base$event.time<=t_hor)],0)

        # one subject's predicted survival probability
        pred_surv <- function(patient.id,z_var,x_var){
          ind_dat <- test_dat_base[test_dat_base$id==patient.id,]
          z_value <- ind_dat[,z_var]
          x_value <- ind_dat[,x_var]
          cond_sur <- exp(-exp(coef[1]+coef[2:(1+length(z_value))]%*%t(z_value))*(F0t_horL^exp(coef[(2+length(z_value)):(length(coef))]%*%t(x_value))))
          return(cond_sur)
        }

        cl <- makeCluster(no_cores)
        registerDoSNOW(cl)
        pred_surv_L <- foreach(i=unique(test_dat_base$id),.combine="rbind")%dopar%{
          return(c(i,pred_surv(i,z_var=beta_variable,x_var=gamma_variable)))
        }
        stopCluster(cl)

        pred_surv_L[,1] <- as.numeric(pred_surv_L[,1])
        pred_surv_L[,2] <- as.numeric(pred_surv_L[,2])
        colnames(pred_surv_L) <- c("id","surv")

        return(pred_surv_L)

      }
      names(pred_surv_L_thor)<-paste0("thor=",thor_list[[which(landmark_surv_vec==landmark)]])
    }
    if(exists("pred_cure_L")&exists("pred_surv_L_thor")){
      return(list(train_est=result_est_L,test_data_landmark=test_dat_base,pred_cure=pred_cure_L,pred_surv=pred_surv_L_thor,test_data_long=test_prep_data))
    }else if(!exists("pred_cure_L")&exists("pred_surv_L_thor")){
      return(list(train_est=result_est_L,test_data_landmark=test_dat_base,pred_surv=pred_surv_L_thor,test_data_long=test_prep_data))
    }else if(exists("pred_cure_L")&!exists("pred_surv_L_thor")){
      return(list(train_est=result_est_L,test_data_landmark=test_dat_base,pred_cure=pred_cure_L,test_data_long=test_prep_data))
    }else if(!exists("pred_cure_L")&!exists("pred_surv_L_thor")){
      return(list(train_est=result_est_L,test_data_landmark=test_dat_base,test_data_long=test_prep_data))
    }
  }

  names(result) <- paste0("landmark=",landmark_vec)

  return(result)

}

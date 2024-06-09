#' Plot a patient's observed biomarker values and predicted conditional survival function
#' @description Plots a patient's observed biomarker values up to a landmark time \eqn{L} and predicted
#' individual conditional survival function at \eqn{L}
#'
#' @param landmark a value of the landmark time included in the object assigned below.
#' @param predict.id a subject id who is still risk free at landmark time in the test dataset.
#' @param model.type a character for the type of model used for prediction. It could be "joint model" or "landmark model".
#' @param object an object of pred_surv_joint_model() function for JMFHC model or pred_landmark_model() function for the LFHC model.
#' @param original_biomarker_form a logical value (True or False) to determine whether to use the originial form of biomarker values biomarker values
#' presented at the y-axis of the plot.
#' @param no_cores the number of cores used during the estimation. The default is 7.
#'
#' @return a plot containing a patient's observed biomarker measurments up to the time \eqn{L} and the predicted conditional survival probability at \eqn{L}
#' @export
#'
#' @examples
#' # Example of JMFHC
#' jmfhc_est <- readRDS("jmfhc_estresult.rds")
#' predict_surv_jmfhc  <- pred_surv_joint_model(landmark_surv_vec=c(10,20),
#'                                              thor_list=list(c(5,10),c(5,10)),
#'                                              object=jmfhc_est,
#'                                              test_dat=test_dat,
#'                                              no_cores=7)
#' plot_con_surv(landmark=10,predict.id=1,model.type="joint model",object=predict_surv_jmfhc,original_biomarker_form=TRUE,no_cores=7)
#'
#' # Example of the LFHC model
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
#' plot_con_surv(landmark=10,predict.id=1,model.type="landmark model",object=lfhc_pred_all,original_biomarker_form=TRUE,no_cores=7)
#'
plot_con_surv <- function(landmark,
                          predict.id,
                          model.type,
                          object, # the estimation and prediction
                          original_biomarker_form, #TRUE or FALSE
                          no_cores=7){

  twoord.stackplot <- function(lx, rx, ldata, rdata, lcol, rcol, ltype, rtype,
                               border, rylab, lylab, xlab, mar = c(6, 6, 4, 6),
                               lylimits,rylimits,main,font.main=1,
                               lwd.left=3,lwd.right=rep(3,5),cex=1.2){
    #convert to matrix
    if(is.vector(ldata)){ldata <- as.matrix(ldata)}
    if(is.vector(rdata)){rdata <- as.matrix(rdata)}
    xlimits <- range(lx, rx)
    oldmar <- par("mar")
    par(mar = mar)

    #------------------------------------------------------------------------
    #left y-axis plot
    plot(lx, ldata[, 1], xlim=xlimits, ylim=lylimits, col=lcol[1],
         type=ltype[1], axes=FALSE, ylab="", xlab="",main=main,font.main=font.main,
         lwd=lwd.left,cex.lab=cex,
         cex.axis=3,
         cex.main=2)
    xylim <- par("usr")
    mtext(lylab, 2, line=3, col = "black",cex=cex)
    title(xlab="Time", line=3, cex.lab=2)
    axis(1,cex.axis=1.5) #x axis
    axat <- axis(2, col="black", labels=FALSE) #left y axis
    abline(v=xylim[1], col="black")
    mtext(axat, 2, 1, at = axat, col = "black",cex=cex)
    box()

    #------------------------------------------------------------------------
    #right y-axis plot
    par(new=TRUE)
    plot(rx, rdata[, 1], xlim=xlimits, ylim=rylimits, col=rcol[1],
         type=rtype[1], axes=FALSE, ylab="", xlab="",lwd=lwd.right[1])
    if(ncol(rdata) > 1)	{
      for(i in 2:ncol(rdata)){
        lines(rx, rdata[, i], col=rcol[i], type=rtype[i],lwd=lwd.right[i])
      }
    }
    axat <- axis(4, col="black", labels=FALSE) #right y axis
    abline(v=xylim[1], col="black")
    mtext(axat, 4, 1, at = axat, col = "black",cex=cex)
    mtext(rylab, 4, line=3, col = "black",cex=cex)
    abline(v=landmark,lty = 2, lwd = lwd.left-0.5)
    par(mar = oldmar)
  }

  if (model.type=="landmark model"){
    object <- object[[which(names(object)==paste0("landmark=",landmark))]]
    train_base_dat <- object$train_est$data
    landmark_dat <- object$test_data_landmark
    measure_time <- object$train_est$setting$fu_time_original
    id <- object$train_est$setting$id
    if (original_biomarker_form){
      landmark_dat_long <- object$test_data_long
      biomarker_value <- object$train_est$setting$fu_measure_original
      bio_label <- "Original biomarker value"
    }else{
      biomarker_value <- object$train_est$setting$fu_measure
      bio_label <- "Transformed biomarker value"
    }
    biomarker <- landmark_dat_long[landmark_dat_long[,id]==predict.id & landmark_dat_long[,measure_time]<=landmark,c(biomarker_value,measure_time)]

  }else if(model.type=="joint model"){
    train_base_dat <- object$train_est$dat_baseline
    landmark_dat <- object[[which(names(object)==paste0("landmark=",landmark))]]$test_data_landmark
    measure_time <- object$train_est$setting$fu_time_original
    id <- object$train_est$setting$id
    if (original_biomarker_form){
      biomarker_value <- object$train_est$setting$fu_measure_original
      bio_label <- "Original biomarker value"
    }else{
      biomarker_value <- object$train_est$setting$fu_measure
      bio_label <- "Transformed biomarker value"
    }
    biomarker <- landmark_dat[landmark_dat[,id]==predict.id & landmark_dat[,measure_time]<=landmark,c(biomarker_value,measure_time)]

  }
  colnames(biomarker) <- c("bio_value","time")
  biomarker_range <- range(biomarker[,1])

  # variables in the data set
  event_time <- object$train_est$setting$event_time
  event_status <- object$train_est$setting$event_status

  if(!predict.id %in% unique(landmark_dat[,id])){
    stop("The entered patient id for prediction is not in the set of patients who are still at risk at the prespecified landmark time")
  }else{
    id_label = paste0("ID=",predict.id)
  }

  base_dat <- train_base_dat[order(train_base_dat[,event_time]),]
  # remove event times with duplicated F_0(t)
  unique.event_time <- base_dat[!duplicated(base_dat$base_cdf),event_time]
  if (landmark < max(unique.event_time)){
    select.time <- c(quantile(unique.event_time[unique.event_time>landmark], probs = seq(0, 1, by = .1)))
    if (model.type=="landmark model"){
      ind_con_surv <- foreach (t_hor=select.time-landmark,.combine = "c")%do%{
        # probability of being at risk for another t_hor at L
        F0t_horL <- ifelse(sum(base_dat$event.time<=t_hor)!=0,base_dat$base_cdf[sum(base_dat$event.time<=t_hor)],0)

        # one subject's predicted survival probability
        pred_surv <- function(patient.id,z_var,x_var){
          ind_dat <- landmark_dat[landmark_dat$id==patient.id,]
          z_value <- ind_dat[,z_var]
          x_value <- ind_dat[,x_var]
          coef <- object$train_est$coef
          cond_sur <- exp(-exp(coef[1]+coef[2:(1+length(z_value))]%*%t(z_value))*(F0t_horL^exp(coef[(2+length(z_value)):(length(coef))]%*%t(x_value))))
          return(cond_sur)
        }

        pred_surv_L <- pred_surv(patient.id=predict.id,z_var=object$train_est$setting$beta_variable,
                                 x_var=object$train_est$setting$gamma_variable)

        return(pred_surv_L)
      }
    }else if(model.type=="joint model"){
      prob <- pred_surv_joint_model(landmark_surv_vec=landmark, thor_list=list(select.time-landmark),
                                    object=object$train_est,test_dat=landmark_dat[landmark_dat[,id]==predict.id,],no_cores=no_cores)
      ind_con_surv <- sapply(seq(prob[[1]][seq(select.time)]),function(x) prob[[1]][seq(select.time)][[x]][,2])
    }
    # conditional survival probabilities starting from L
  }else{
    select.time <- landmark
    ind_con_surv <-  rep(1,length(select.time))
  }

  return(twoord.stackplot(lx=biomarker$time,rx=c(landmark,select.time),
                          ldata=biomarker$bio_value,lylimits=biomarker_range,
                          rdata=c(1,ind_con_surv),
                          lcol="red",
                          rcol="blue",
                          ltype="b", rtype="l",
                          lylab=bio_label,
                          rylab="Predicted conditional survival probability",
                          xlab="Time",
                          main=paste0(id_label," at L=",landmark),
                          font.main=1,rylimits=c(0,1),
                          lwd.right=3))
}

#' Plot a patient's observed biomarker values and predicted conditional survival function
#' @description Plots a patient's observed biomarker values up to a landmark time \eqn{L} and predicted
#' individual conditional survival function at \eqn{L}
#'
#' @param L a value of the landmark time in the same unit as the observed survival times and the original form of measurement times
#' in the data set.
#' @param predict.id a character for the type of prediction. It could be one" or "new".
#' predict.id="one" is to do prediction for a specific patient who is still at risk at landmark time \eqn{L} in the data set, where the
#' patient id is specified in the parameter predict.id.one.
#' predict.id="new" is to do prediction for a new patient who is not in the data set and at risk at landmark time \eqn{L}. This new
#' patient's new values are specified in the parameters new.fu_measure, new.fu_time_fixed_variable, new.fu_time_random_variable,
#' new.baseline_value_lmm, new.z_value and new.x_value.
#' @param predict.id.one a patient's id number from the data set when predict.id="one",
#' which could be numeric and character depending on the form of patient id column in the data set.
#' By default predict.id.one = NULL.
#' @param biomarker_form a character for the type of biomarker values presented as the y-axis of the plot.
#' It could be "original" or "transformed".
#' @param new.fu_measure_original a vector of the new patient's longitudinal measurement times in the original form.
#' When biomarker_form="transformed", new.fu_measure_original can be set as NULL.
#' By default new.fu_measure_original = NULL.
#' @param new.fu_measure a vector of the new patient's biomarker measurements up to time \eqn{L} when predict.id="new". All measurements are
#' in the same form specified for the parameter fu_measure in the fitted longitudinal submodel of jmfhc_point_est().
#' By default new.fu_measure = NULL.
#' @param new.fu_time_original a vector of the new patient's longitudinal measurement times in the original form, which is in the same
#' unit of the landmark time.
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
#' @param object an object of jmfhc_point_est() function.
#' @param no_cores the number of cores used during the estimation. The default is 7.
#'
#' @return a plot containing a patient's observed biomarker measurments up to the time \eqn{L} and the predicted conditional survival probability at \eqn{L}
#' @export
#'
#' @examples
#' result_jmfhc <- jmfhc::jmfhc_point_est(data=jmfhc_dat, event_time="event.time", event_status="event",
#'                                        id="patient.id", beta_variable="trt", gamma_variable="trt",
#'                                        fu_measure_original="measure",fu_measure="measure",
#'                                        fu_time_original="mes.times",fu_time_fixed_variable="mes.times",
#'                                        fu_time_random_variable="mes.times")
#' plot_con_surv(L=10,predict.id="one",predict.id.one=3,biomarker_form="original",object=result_jmfhc)
#' plot_con_surv(L=10,predict.id="new",biomarker="original",
#'               new.fu_measure_original=c(5.2,5.1,1.6,0.9,-1.1,-2.6,-5.3,-8.0,-7.5,-11,-12),
#'               new.fu_measure=c(5.2,5.1,1.6,0.9,-1.1,-2.6,-5.3,-8.0,-7.5,-11,-12),
#'               new.fu_time_original= 0:10,
#'               new.fu_time_fixed_variable= 0:10,
#'               new.fu_time_random_variable= 0:10,
#'               new.baseline_value_lmm=NULL,
#'               new.z_value=0,
#'               new.x_value=0,
#'               object=result_jmfhc)
#'
#' @import parallel
#' @import foreach
#' @import doSNOW
#'
plot_con_surv <- function(L,
                          predict.id,
                          predict.id.one=NULL,
                          biomarker_form,
                          new.fu_measure_original=NULL,
                          new.fu_measure=NULL,
                          new.fu_time_original=NULL,
                          new.fu_time_fixed_variable=NULL,
                          new.fu_time_random_variable=NULL,
                          new.baseline_value_lmm=NULL,
                          new.z_value=NULL,
                          new.x_value=NULL,
                          object,
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
    abline(v=L,lty = 2, lwd = lwd.left-0.5)
    par(mar = oldmar)
  }

  # variables in the data set
  event_time <- object$setting$event_time
  event_status <- object$setting$event_status
  measure_time <- object$setting$fu_time_original
  id <- object$setting$id
  updated_dat <- object$dat_long
  base_dat <- object$dat_baseline
  # the data set containing patients who are still at risk at the landmark time
  landmark_dat <- updated_dat[updated_dat[,event_time]>= L & updated_dat[,measure_time] <= L,]

  if (predict.id=="one"){
    if(!predict.id.one %in% unique(landmark_dat[,id])){
      stop("The entered patient id for prediction is not in the set of patients who are still at risk at the prespecified landmark time")
    }
    id_label = paste0("ID=",predict.id.one)
  }else{
    id_label = "A new patient"
  }

  if (biomarker_form=="original"){
    biomarker_value <- object$setting$fu_measure_original
    biomarker_range <- range(updated_dat[, biomarker_value])
    if (predict.id=="new"){
      biomarker_value <- new.fu_measure_original
    }
    bio_label <- "Original biomarker value"
  }else{
    biomarker_value <- object$setting$fu_measure
    biomarker_range <- range(updated_dat[, biomarker_value])
    if (predict.id=="new"){
      biomarker_value <- new.fu_measure
    }
    bio_label <- "Transformed biomarker value"
  }

  base_dat <- base_dat[order(base_dat[,event_time]),]
  # remove event times with duplicated F_0(t)
  unique.event_time <- base_dat[!duplicated(base_dat$base_cdf),event_time]
  if (L < max(unique.event_time)){
    select.time <- c(quantile(unique.event_time[unique.event_time>L], probs = seq(0, 1, by = .1)))

    # conditional survival probabilities starting from L
    cl <- makeCluster(no_cores)
    registerDoSNOW(cl)
    if (predict.id=="one"){
      ind_con_surv <-foreach(thor.time=select.time,.combine="c",.packages=c("cubature","R.utils"),.export= "est_con_survival")%dopar%{
        prob <- est_con_survival(L=L,t_hor=thor.time-L,predict.id="one",predict.id.one=predict.id.one,object=object)
        print(thor.time)
        print(prob)
        return(prob$con_surv[1,2])
      }
    }else{
      ind_con_surv <-foreach(thor.time=select.time,.combine="c",.packages=c("cubature","R.utils"),.export= "est_con_survival")%dopar%{
        prob <- est_con_survival(L=L,t_hor=thor.time-L,predict.id="new",
                                 new.fu_measure=new.fu_measure,
                                 new.fu_time_fixed_variable=new.fu_time_fixed_variable,
                                 new.fu_time_random_variable=new.fu_time_random_variable,
                                 new.baseline_value_lmm=new.baseline_value_lmm,
                                 new.z_value=new.z_value,
                                 new.x_value=new.x_value,
                                 object=object)
        print(thor.time)
        print(prob)
        return(prob$con_surv[1,2])
      }
    }
    stopCluster(cl)
  }else{
    select.time <- L
    ind_con_surv <- 1
  }

  if(predict.id=="one"){
    biomarker <- updated_dat[updated_dat[,id]==predict.id.one & updated_dat[,measure_time]<=L,c(biomarker_value,measure_time)]
    colnames(biomarker) <- c("bio_value","time")
  }else{
    biomarker <- data.frame(bio_value=biomarker_value[new.fu_time_original<=L],time=new.fu_time_original[new.fu_time_original<=L])
  }

  return(twoord.stackplot(lx=biomarker$time,rx=c(L,select.time,max(base_dat[,"event.time"])),
                          ldata=biomarker$bio_value,lylimits=biomarker_range,
                          rdata=c(1,ind_con_surv,ind_con_surv[length(ind_con_surv)]),
                          lcol="red",
                          rcol="blue",
                          ltype="b", rtype="l",
                          lylab=bio_label,
                          rylab="Predicted conditional survival probability",
                          xlab="Time",
                          main=paste0(id_label," at L=",L),
                          font.main=1,rylimits=c(0,1),
                          lwd.right=3))
}

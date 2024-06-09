prep_data_landmark <- function(dat,landmark,id,event_time,fu_time_original,fu_measure){

  colnames(dat)[colnames(dat)==id] <- "id"
  colnames(dat)[colnames(dat)==event_time] <- "event.time"
  colnames(dat)[colnames(dat)==fu_time_original] <- "measure.time"
  colnames(dat)[colnames(dat)==fu_measure] <- "measure"

  # extract the most current longitudinal biomarker value at the landmark time (last value carried forward method)
  dat_most_current_landmark <-
    dat %>%
    group_by(id) %>%
    arrange(measure.time) %>%
    filter(measure.time<=landmark & event.time>=landmark)

  dat_most_current_landmark <- dat_most_current_landmark[
    with(dat_most_current_landmark, order(id, measure.time)),
  ]

  # baseline data #
  dat_base <- dat_most_current_landmark[dat_most_current_landmark$measure.time==0,]

  # add the most current biomarker value
  dat_most_current_one <-
    dat_most_current_landmark %>% group_by(id) %>% filter(row_number()==n()) %>% select(id,measure)
  colnames(dat_most_current_one)[colnames(dat_most_current_one)=="measure"] <- "current_measure"
  dat_base <- merge(dat_base,dat_most_current_one,by="id")

  return(dat_base)
}

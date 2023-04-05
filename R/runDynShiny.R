#' Shiny Web App for indiviudal dynamic prediction
#'
#' @return predicted values
#' @export
#'
#' @examples runDynPred()
#' # For demenstration, jmfhc_estresult.rds is an example file for RDS.
#' # demo_new_patient.csv is a new patient's records up to the landmark time=10
#'
#' @import shiny
runDynPred <- function(){
  runApp('curedynpred_shiny')
}

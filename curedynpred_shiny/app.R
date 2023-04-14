#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(curedynpred)

# Define UI for application
ui <- fluidPage(

  titlePanel("Individual Dynamic Prediction for Cure and Survival Based on Longitudinal Data"),

    sidebarPanel(

      # import the rds file
      fileInput("file",
                "Upload the object (.RDS file) of the fitted JMFHC from jmfhc_point_est() function"),

      tags$hr(),

      # set the landmark time with a slider
      uiOutput("slider_landmark"),

      # select patients for prediction
      radioButtons("predict.id", label = "Perform prediction on:",
                   choices = list("A specifc at-risk patient from the data set" = 1,
                                  "A new at-risk patient" = 2),
                   selected = 1),

      # select the item to predict
      radioButtons("predict.item", "What to predict",
                   choices = list("Indivudal conditional cure rate" = 1,
                                  "Individual conditional survival probability" = 2),
                   selected = 1),

      checkboxInput("plot","Plot observed longitudinal measurements and individual conditional survival function?"),

    ),

    # only show this panel if predict.item==2. create a slider for time horizon of prediction
    conditionalPanel(
      condition = "input.predict.item==2",
      uiOutput("slider_t_hor"),
    ),

    # only show this panel if predict.id==1. select a patient id from a select box
    conditionalPanel(
      condition = "input.predict.id==1",
      uiOutput("predict.id.one"),
      uiOutput("predict.new")
    ),

    # Display
    mainPanel(
      textOutput("resultlabel"),
      tableOutput("predict_table"),
      plotOutput("predict_plot")
    )

)



# Define server logic
server <- function(input, output) {

  options(shiny.maxRequestSize=50*1024^2)

  object <- reactiveValues(dat_baseline=NULL,setting=NULL)
  observeEvent(input$file, {
    file_to_read <- input$file
    object$dat_baseline <- readRDS(file_to_read$datapath)$dat_baseline
    object$setting <- readRDS(file_to_read$datapath)$setting
  })

  output$slider_landmark <- renderUI({
    req(object$dat_baseline)
    min_event_time <- round(range(object$dat_baseline$event.time)[1])
    max_event_time <-  round(range(object$dat_baseline$event.time)[2])
    sliderInput("slider_landmark", "Set a landmark time for prediction:",
                min = min_event_time,
                max = max_event_time, value = min_event_time, step = 1)
  })

  output$slider_t_hor <- renderUI({
    req(object$dat_baseline,input$predict.item)
    if (input$predict.item==2){
      half_event_time <-  round(range(object$dat_baseline$event.time)[2]/2)
      sliderInput("slider_t_hor", "Set a time horizon for prediction:",
                  min = 0,
                  max = half_event_time, value = 0, step = 0.5)
    }
  })

  require.var <- reactive({
    req(input$file,input$predict.id)
    if(input$predict.id==2){
      req(input$file)
      var <- unique(c(object$setting$fu_measure_original,object$setting$fu_measure,
                      object$setting$fu_time_original,object$setting$fu_time_fixed_variable,
                      object$setting$fu_time_random_variable,object$setting$baseline_var_lmm,
                      object$setting$beta_variable,object$setting$gamma_variable))
      paste0(var,collapse = ',')
    }
  })

  output$predict.id.one <- renderUI({
    req(object$dat_baseline,input$predict.id)
    if(input$predict.id==1){
      unique.id <- unique(object$dat_baseline[object$dat_baseline[,"event.time"]> input$slider_landmark,"patient.id"])
      selectInput("predict.id.one", "Select one at-risk patient's ID from the data set",
                  choices = unique.id)
    }
  })

  output$predict.new <- renderUI({
    req(object$dat_baseline,input$predict.id)
    if (input$predict.id==2){
      fileInput("predict.new",
                paste("Upload the new patient's records up to the landmark time in the long format via CSV File including variables","[",require.var(),"]."),
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv"))
    }
  })

  output$resultlabel <- renderText({
    if (input$predict.id==1){
      if (input$predict.item==1){
        paste("Predicted conditional cure probability for patient ID =",input$predict.id.one,"at landmark time =", input$slider_landmark)
      }else if (input$predict.item==2){
        paste("Predicted conditional survival probability for patient ID =",input$predict.id.one,"at landmark time =", input$slider_landmark, "with the time horizon of prediction =", input$slider_t_hor)
      }
    }else if (input$predict.id==2){
      if (input$predict.item==1){
        paste("Predicted conditional cure probability for the new patient at landmark time =", input$slider_landmark)
      }else if (input$predict.item==2){
        paste("Predicted conditional survival probability for the new patient at landmark time =", input$slider_landmark, "with the time horizon of prediction =", input$slider_t_hor)
      }
    }
  })

  result <- reactive({
    req(input$file,input$predict.item,input$predict.id,input$slider_landmark)
    file <- readRDS(input$file$datapath)
    if(input$predict.id==1 && input$predict.item==1){
      req(input$predict.id.one)
      est_cure_L(L=input$slider_landmark,predict.id="one",predict.id.one=input$predict.id.one,
                 object=file)
    }else if (input$predict.id==1 && input$predict.item==2){
      req(input$predict.id.one,input$slider_t_hor)
      est_con_survival(L=input$slider_landmark,t_hor=input$slider_t_hor,predict.id="one",
                       predict.id.one=input$predict.id.one,
                       object=file)$con_surv
    }
    else if(input$predict.id==2){
      req(input$predict.new)
      newpatient_to_read <- input$predict.new
      newpatient <- read.table(newpatient_to_read$datapath,sep=",",header=T)
      fu_measure <- newpatient[,object$setting$fu_measure]
      fu_time_fixed_variable <- unlist(newpatient[,object$setting$fu_time_fixed_variable])
      fu_time_random_variable <- unlist(newpatient[,object$setting$fu_time_random_variable])
      if (length(newpatient[1,object$setting$baseline_var_lmm])==0){
        baseline_var_lmm <- NULL
      } else{
        baseline_var_lmm <- newpatient[1,object$setting$baseline_var_lmm]
      }
      z_value <- newpatient[1,object$setting$beta_variable]
      if (length(newpatient[1,object$setting$gamma_variable])==0){
        x_value <- NULL
      } else{
        x_value <- newpatient[1,object$setting$gamma_variable]
      }

      if(input$predict.item==1){
        est_cure_L(L=input$slider_landmark,predict.id="new",
                   new.fu_measure=fu_measure,
                   new.fu_time_fixed_variable=fu_time_fixed_variable,
                   new.fu_time_random_variable=fu_time_random_variable,
                   new.baseline_value_lmm=baseline_var_lmm,
                   new.z_value=z_value,
                   new.x_value=x_value,
                   object=file)
      }else if (input$predict.item==2){
        req(input$slider_t_hor)
        est_con_survival(L=input$slider_landmark,t_hor=input$slider_t_hor,predict.id="new",
                         new.fu_measure=fu_measure,
                         new.fu_time_fixed_variable=fu_time_fixed_variable,
                         new.fu_time_random_variable=fu_time_random_variable,
                         new.baseline_value_lmm=baseline_var_lmm,
                         new.z_value=z_value,
                         new.x_value=x_value,
                         object=file)$con_surv
      }
    }
  })

  output$predict_table <- renderTable({
    result()
  },digits=4)

  output$predict_plot <- renderPlot({
    req(input$file,input$predict.id,input$slider_landmark,input$plot)
    file <- readRDS(input$file$datapath)
    if (input$plot==TRUE){
      if (input$predict.id==1){
        req(input$predict.id.one)

          plot_con_surv(L=input$slider_landmark,predict.id="one",predict.id.one=input$predict.id.one,
                        biomarker_form="original",object=file)
      }else if(input$predict.id==2){
          req(input$predict.new)
          newpatient_to_read <- input$predict.new
          newpatient <- read.table(newpatient_to_read$datapath,sep=",",header=T)
          fu_measure_original <- newpatient[,object$setting$fu_measure_original]
          fu_measure <- newpatient[,object$setting$fu_measure]
          fu_time_original <- newpatient[,object$setting$fu_time_original]
          fu_time_fixed_variable <- unlist(newpatient[,object$setting$fu_time_fixed_variable])
          fu_time_random_variable <- unlist(newpatient[,object$setting$fu_time_random_variable])
          if (length(newpatient[1,object$setting$baseline_var_lmm])==0){
            baseline_var_lmm <- NULL
          } else{
            baseline_var_lmm <- newpatient[1,object$setting$baseline_var_lmm]
          }
          z_value <- newpatient[1,object$setting$beta_variable]
          if (length(newpatient[1,object$setting$gamma_variable])==0){
            x_value <- NULL
          } else{
            x_value <- newpatient[1,object$setting$gamma_variable]
          }

            plot_con_surv(L=input$slider_landmark,predict.id="new",
                          biomarker_form="original",
                          new.fu_measure_original=fu_measure_original,
                          new.fu_measure=fu_measure,
                          new.fu_time_original=fu_time_original,
                          new.fu_time_fixed_variable=fu_time_fixed_variable,
                          new.fu_time_random_variable=fu_time_random_variable,
                          new.baseline_value_lmm=baseline_var_lmm,
                          new.z_value=z_value,
                          new.x_value=x_value,
                          object=file)
        }
    }
  })

}


# Run the application
shinyApp(ui = ui, server = server)

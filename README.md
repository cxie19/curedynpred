# curedynpred
Perform individual dynamic prediction on cure and survival
curedynpred is an R package to predict individual conditional cure rates and conditional 
survival probabilities. It can also evaluate the fitted joint model's predictive performance 
on conditional survival probabilities via the time-dependent area under 
the receiver operating characteristic curve (AUC) and Brier score. 
The fitted joint model is our proposed joint cure model that simultaneously deals 
with longitudinal data in a linear mixed-effects submodel 
and flexible patterns of hazard ratios over time in a promotion time cure submodel.
It is named a joint model with a flexible-hazards cure submodel (JMFHC).
The model estimation for this proposed model can be found  at https://github.com/cxie19/jmfhc. <br />

## How to get started

Install the R package using the following commands on the R console:

```{r}
install.packages("devtools")
devtools::install_github("cxie19/curedynpred")
library(curedynpred)
```

An example data set called *long_dat* is provided in this package. It
is a longitudinal and cure-survival data set. Its documentation can be 
seen by using the following command.

```{r}
help(long_dat)
```

The function *est_cure_L* is called to predicts individual conditional probability of being cured at a landmark time.
Its documentation can be seen by using the following command.

```{r}
help(est_cure_L)
```

The function *est_con_survival* is called to predict individual conditional probability of not experiencing the event of
interest in an additional time given that a patient remains risk-free at least until a landmark time.
This function evaluates the fitted joint cure model's predictive performance by the time-dependent
AUC and Brier score.
Its documentation can be seen by using the following command.

```{r}
help(est_con_survival)
```

## Example
For example, we want to fit a JMFHC for the example data *long_dat*.
This joint model has repeatedly measured biomarker values as the outcome of the 
longitudinal submodel with measurement times as the 
covariate and treatment as the short- and long-term covariate in the cure 
submodel. These two submodels share individual random effects.
The time unit for the time related variables is month.

We call the function *jmfhc_point_est* for point estimation from the R package cxie219/jmfhc, and the following command is used.

```{r}
result_coef <- jmfhc_point_est(data=jmfhc_dat, 
                               event_time="event.time", event_status="event", 
                               id="patient.id", 
                               beta_variable="trt", gamma_variable="trt", 
                               fu_measure="measure", fu_time_original="mes.times",                                        fu_time_variable="mes.times")
```
The point estimation could be found as a file called jmfhc_estresult.rds.

To predict conditional cure rates at 10 months for all patients who are still at risk at 10 months, we call the function *est_cure_L*.

```{r}
predict_cure <- est_cure_L(L=10,predict.id="all",object=result_jmfhc)
```

To predict all patients' conditional probability of not experiencing the event of interest in an additional 5 months 
given that the patient remains risk-free at least 10 months, we call the function *est_con_survival*. We also compute the AUC 
and Brier score.

```{r}
predict_surv <- est_con_survival(L=10,t_hor=5,predict.id="all",AUC=TRUE,Brier=TRUE,object=result_jmfhc)
```



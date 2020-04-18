# MeasurementErrorAnalysis
The function calculates the parameters estimates of the measurement error regression model. Developed as part of an undergraduate research project at Moravian College. 

## Environment Requirements
- R (version 3.4.3)
- R libraries 
  - knitr
  - tidyverse
  - mosaic
  - MASS
  
 ## Functions
 The MEM_functionMult function requires three inputs passed as .csv files: a nm×1 matrix for the response variable, a nm×p matrix for the predictor variables, and a nm×1 matrix with unique indentification variables for the participants. Additionally, the function requires a confidence level as a numeric value. Here, n is the number of participants,  <img src="https://render.githubusercontent.com/render/math?math=m_i"> is the number of replications for each individual i, and p is the number of predictors for the model. The function returns a a px1 for the estimator of <img src="https://render.githubusercontent.com/render/math?math=\beta_1">, a single value for <img src="https://render.githubusercontent.com/render/math?math=\beta_0">, and a pxp <img src="https://render.githubusercontent.com/render/math?math=Var(\beta_1)"> matrix. In addition to the parameter estimates, the funtion returns the standard error of each <img src="https://render.githubusercontent.com/render/math?math=Var(\beta_1)"> value, the corresponding test statistic, and the confidence interval. 
  
## Running the Model
 R --slave --no-restore --file=MEM_function.R --args PredictorMatrix.csv ResponseMatrix.csv IDMatrix.csv CIlevel(numeric) 
 
 


 

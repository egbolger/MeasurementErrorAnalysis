library(knitr)
library(mosaic)
library(MASS)
library(tidyverse)




MEM_functionMult <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  x2 <-args[1]
  y2 <-args[2]
  id2<-args[3]
  CIlevel<-args[4]
  
  x2 <- read.csv(file = x2, header = TRUE)
  y2 <- read.csv(file = y2, header = TRUE)
  id2 <- read.csv(file = id2, header = TRUE)
  
  dat2<-x2
  dat2$y<-y2
  dat2$id2<-as.numeric(as.factor(unlist(id2)))
  
  #Transforming x,y to matrices to use in for loops later
  x2<-as.matrix(x2)
  y2<-as.matrix(y2)
  
  #Number participants, predictors, responses
  numpart<- length(unique(unlist(id2)))
  numpred <- ncol(x2)
  numresp <- ncol(y2)
  
  #Keeping participants with more than one day of data
  dat2_rep <- dat2 %>% dplyr::group_by(id2) %>% filter(n()>1)
  x2_repdf <- dat2_rep %>% dplyr::select(id2,colnames(x2))
  y2_repdf <- dat2_rep %>% dplyr:: select(id2,starts_with("y"))
  x2_rep<-as.matrix(x2_repdf[,c(colnames(x2))])
  y2_rep<-as.matrix(unlist(y2_repdf[,2]))
  #y2_rep<- as.matrix(y2_rep[,1])
  numpart_rep <- length(unique(dat2_rep$id2))
  
  #Calculation for the average across the total number of days for each individual
  x_bardotmatrix<-matrix(ncol = numpred, nrow = numpart)
  for(j in 1:numpred){
    for(i in 1:length(unique(unlist(id2)))){
      x_bardotmatrix[i,j] <- mean(x2[as.numeric(as.factor(unlist(id2)))==i,j])
    }
  }
  
  #Calculation for average of the averages above, ie total average per predictor
  x_dotdotmatrix<-matrix(ncol = numpred, nrow = 1)
  for (k in 1:numpred){
    x_dotdotmatrix[,k] <- mean(x_bardotmatrix[,k])
  }
  
  #Calculation for average response for each individual
  y_bardotmatrix<-matrix(ncol = numresp, nrow = numpart)
  for(j in 1:numresp){
    for(i in 1:length(unique(unlist(id2)))){
      y_bardotmatrix[i,j] <- mean(y2[as.numeric(as.factor(unlist(id2)))==i,j])
    }
  }
  
  #Calculation for total average for each response
  y_dotdotmatrix<-matrix(ncol = numresp, nrow = 1)
  for (k in 1:numresp){
    y_dotdotmatrix[,k] <- mean(y_bardotmatrix[,k])
  }
  
  #Initialization of the transpose of M_xx and M_xy
  Mxx_sdxt <- matrix(nrow=numpred, ncol = numpart)
  Mxx_sdx <- matrix(nrow=numpart, ncol = numpred)
  Mxx_sdmatrix<- matrix(nrow = numpred, ncol = numpred)
  for (j in 1:numpred){
    for (i in 1:numpart){
      #take each individual predictor value minus the total average of that predictor
      Mxx_sdx[i,j] <- x_bardotmatrix[i,j] - x_dotdotmatrix[,j]
      Mxx_sdxt <- t(Mxx_sdx)
      Mxx_sdmatrix <- Mxx_sdxt %*% Mxx_sdx
    }
  }
  M_xxmatrix <- (1/(numpart-1))*Mxx_sdmatrix
  
  #Initialization of the transpose of M_xy and M_xy
  Mxy_sdy <- matrix(nrow=numpart, ncol = numresp)
  Mxy_sdmatrix<- matrix(nrow = numpred, ncol = numpred)
  for (j in 1:ncol(Mxy_sdy)){
    for (i in 1:numpart){
      #take each individual response value minus the total average of that response
      Mxy_sdy[i,j] <- y_bardotmatrix[i,j] - y_dotdotmatrix[,j]
      Mxx_sdxt <- t(Mxx_sdx)
      Mxy_sdmatrix <- Mxx_sdxt %*% Mxy_sdy
    }
  }
  M_xymatrix <- (1/(numpart-1))*Mxy_sdmatrix
  
  Myy_sdyt <- matrix(nrow=numresp, ncol = numpart)
  Myy_sdy <- matrix(nrow=numpart, ncol = numresp)
  Myy_sdmatrix<- matrix(nrow = numresp, ncol = numresp)
  for (j in 1:numresp){
    for (i in 1:numpart){
      #take each individual response value minus the total average of that predictor
      Myy_sdy[i,j] <- y_bardotmatrix[i,j] - y_dotdotmatrix[,j]
      Myy_sdyt <- t(Myy_sdy)
      Myy_sdmatrix <- Myy_sdyt %*% Myy_sdy
    }
  }
  M_yymatrix <- (1/(numpart-1))*Myy_sdmatrix
  
  #Day to day variability
  #The variability for each person for each predictor
  #Sum of 'squares' of each individual predictor minus the average predictor for the person
  siguu_i<-array(dim = c(numpred, numpred, numpart_rep))
  for (i in 1:numpart_rep){
    #for each individual's predictor value compute the variance
    siguu_i[,,i]<-var(x2_rep[dat2_rep$id2 ==i,])
  }
  Sig_uu <- apply(siguu_i, c(1,2), mean, na.rm=T)
  
  #The variability for each person for each response
  #Sum of 'squares' of each individual response minus the average response for the person
  sigww_i<-array(dim = c(numresp, numresp, numpart_rep))
  for (i in 1:numpart_rep){
    #for each individual's response value compute the variance
    sigww_i[,,i]<-var(y2_rep[dat2_rep$id2 ==i,])
  }
  Sig_ww <- apply(sigww_i, c(1,2), mean, na.rm=T)
  
  #The covariance matrix for each predictor and response
  #Each individual predictor minus average predictor times each individual response minus average response
  sigwu_i<-array(dim=c(numpred,numresp,numpart_rep))
  for (i in 1:numpart_rep){
    #for each individual's value compute the covariance
    sigwu_i[,,i] <- cov(y2_rep[dat2_rep$id2 == i], x2_rep[dat2_rep$id2 == i,])
  }
  Sig_wu <- apply(sigwu_i, c(1,2), mean, na.rm=T)
  
  #Computing the beta matrix (1xnumpred) and the value for beta_0
  beta_matrix<- ginv(M_xxmatrix-Sig_uu)%*%(M_xymatrix-Sig_wu)
  beta_0mult<- y_dotdotmatrix - (x_dotdotmatrix %*% beta_matrix)
  
  s_vvmult<-(1/(numpart-numpred))*((t(y_bardotmatrix - rep(beta_0mult,numpart) - (x_bardotmatrix%*%beta_matrix)))%*% (y_bardotmatrix - rep(beta_0mult,numpart) - (x_bardotmatrix%*%beta_matrix)))
  
  s_rrmult <- Sig_ww - (2*(t(beta_matrix)%*%(Sig_wu))) + (t(beta_matrix)%*%(Sig_uu)%*%(beta_matrix))
  
  s_vvmult<-as.numeric(s_vvmult)
  s_rrmult<-as.numeric(s_rrmult)
  #Variance for error in the equation
  Sig_qq<-s_vvmult - s_rrmult
  
  #Error in e is independent of x and u
  Sig_eu <-array(0, dim = c(1, numpred))
  
  #Variance of error in the response
  Sig_ee <- M_yymatrix - 2*t(M_xymatrix)%*%beta_matrix + t(beta_matrix)%*%M_xxmatrix%*%beta_matrix + 2*Sig_eu%*%beta_matrix - t(beta_matrix)%*%Sig_uu%*%beta_matrix
  
  #Variance for the beta estimates
  var_betamult<-
    (1/numpart)*((ginv(M_xxmatrix-Sig_uu)*s_vvmult)
                 
                 + (ginv(M_xxmatrix-Sig_uu)%*%(Sig_uu*s_vvmult+((Sig_wu - Sig_uu%*%beta_matrix)%*%(t(Sig_wu - Sig_uu%*%beta_matrix))))%*%ginv(M_xxmatrix-Sig_uu))) +
    
    (1/(numpart-numpred))*(ginv(M_xxmatrix-Sig_uu)) %*% (Sig_uu*s_rrmult +((Sig_wu - Sig_uu%*%beta_matrix)%*%(t(Sig_wu - Sig_uu%*%beta_matrix)))) %*%(ginv(M_xxmatrix-Sig_uu))
  
  #Finding the Standard Error
  diag_variance<- diag(var_betamult)
  SE <-matrix(nrow = nrow(var_betamult), 1)
  for (i in 1:nrow(var_betamult)){
    SE[i,] <- round(sqrt(diag_variance[i]),2)
  }
  
  #Finding the test statistics
  tstat_matrix <- matrix(nrow = nrow(beta_matrix), ncol = 1)
  for (i in 1:nrow(beta_matrix)){
    tstat_matrix[i,1] <- round(beta_matrix[i,1]/ sqrt(diag_variance[i]),2)
  }
  
  #Finding the 95% Confidence Intervals
  CIlevel<-as.numeric(CIlevel)
  CI_matrix <- array(dim = c(2, 1, nrow = nrow(beta_matrix)))
  for (i in 1:nrow(beta_matrix)){
    left<-round(beta_matrix[i,1]-qt(CIlevel,df=numpart-nrow(beta_matrix))*sqrt(diag_variance[i]), 2)
    right<-round(beta_matrix[i,1]+qt(CIlevel,df=numpart-nrow(beta_matrix))*sqrt(diag_variance[i]),2)
    CI<- matrix(c(left, right))
    CI_matrix[,,i] <- CI
  }
  
  #Calculating the Residual Values and BLUP for x
  nu<-y_bardotmatrix - x_bardotmatrix%*%beta_matrix - rep(beta_0mult, numpart)
  
  x_blupmult<- t(matrix(rep(x_dotdotmatrix,numpart), ncol = numpart)) + t((ginv(M_xxmatrix)%*% (M_xxmatrix - (1/7)*Sig_uu)) %*% (t(x_bardotmatrix - t(matrix(rep(x_dotdotmatrix,numpart), ncol = numpart)))))
  
  
  res <- list(beta_matrix = beta_matrix, beta_0mult = beta_0mult, var_betamult=var_betamult, SE=SE, tstat_matrix=tstat_matrix, CI_matrix=CI_matrix, nu=nu, x_blupmult=x_blupmult)
  
  return(res)
}

MEM_functionMult()




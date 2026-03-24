library(nnet)
library(neuralnet)
library(dplyr)
#library(LatticeKrig) 
library(maps)
library(plot3D)
library(far)
library(fda)
library(autoFRK)
library(akima)
library(ggplot2) 
#library(mapview)

source("FNN_source.R")
source("functions_source.R")
load("MidwestData.RData") # request by the email
load("loc.info.RData")

year = c(1999:2020)
n = length(year)
spatial = FALSE
kriging = FALSE
### functional data grids
timepts = seq(0,1, length=365) #c(1:365)

### response and scalar variables 
precip.type=3  # use the monthly precipitation data
if(precip.type==1){
  regdat.sub = regdat.V2[ ,c("Yield.c", "avgPRCP")]
}else if (precip.type==2){
  regdat.sub = data.frame(Yield.c = regdat.V2$Yield.c, precip.quart)
}else if (precip.type==3){
  regdat.sub = data.frame(Yield.c = regdat.V2$Yield.c, precip.month)
}



############################################
### hyperparameters and tuning parameters
###
###     this combination is from the tune study
############################################
n.hidden.layer = 4
n.neuron.layer = c(64, 64, 64, 64)
link.ftn = c("sigmoid", "sigmoid","sigmoid","sigmoid")
nbasis = 7

### spatial basis for coef estimation and kriging
load("MRTS.dat.10.RData")
MRTS.dat.coef = MRTS.dat

load("MRTS.dat.150.RData")
MRTS.dat.cov = MRTS.dat


# 
MRTS.merge.coef = c();
MRTS.merge.cov = c()
for( k1 in 1:nrow(regdat.V2)){
    tmp.k1 = regdat.V2[k1, "CountyI"]
    MRTS.merge.coef = rbind(MRTS.merge.coef, MRTS.dat.coef[ which(loc.info$CountyI == tmp.k1), ])
    MRTS.merge.cov = rbind(MRTS.merge.cov, MRTS.dat.cov[ which(loc.info$CountyI == tmp.k1), ])
}


### define the fourier basis to represent functional data
f.nbasis = 21
fr_obj = create.fourier.basis(c(0,1), f.nbasis)


dim.K=2  # daily max and min temperatuer trajectories

fun_fd1 =  Data2fd(timepts, t(fundat[ ,1:365 ]), fr_obj)
fun_fd2 =  Data2fd(timepts, t(fundat[,366:730 ]), fr_obj)
func_cov_1 = fun_fd1$coefs
func_cov_2 = fun_fd2$coefs

fun_data = array(dim = c((f.nbasis), nrow(fundat), dim.K))
fun_data[,,1] = func_cov_1;  
fun_data[,,2] = func_cov_2

##############################################################
### cross-validation
##############################################################
ndivison=10;
set.seed(123)
groups <- gen.groups(nrow(regdat), ndivison)

pe.mat=c(); w.pe.mat=c()
pe=c(); w.pe=c()

spatial= TRUE   # TRUE: spatially varying network parameters/ FALSE: constant network parameter
kriging = TRUE   # TRUE: spatial random process estimation/ FAKSE: no spatial process estimation

## The SVD-funNet is the model with "spatial=TRUE" and "kriging=TRUE".

pred.list <- true.list <- list()
for(vv in 1:ndivison){
  
  cv.id <- groups[[vv]]
  
  dat_train = regdat.sub[-cv.id, ]
  dat_test = regdat.sub[cv.id,]
  
  weight_train <- regdat.V2[ -cv.id, "Area"]
  weight_test <- regdat.V2[ cv.id, "Area"]
  
  fundat_train = array(dim = c((f.nbasis), nrow(dat_train), dim.K))
  fundat_test = array(dim = c((f.nbasis), nrow(dat_test), dim.K))
  fundat_train[,,1:dim.K] = fun_data[, -cv.id, 1:dim.K]
  fundat_test[,,1:dim.K] = fun_data[, cv.id, 1:dim.K]
  
  if(spatial==TRUE){
    
    spatial_train_coef = as.matrix(MRTS.merge.coef[-cv.id, ])
    spatial_test_coef = as.matrix(MRTS.merge.coef[cv.id, ] )
    
    if(kriging==TRUE){
      spatial_train_cov = MRTS.merge.cov[-cv.id, ]
      spatial_test_cov = MRTS.merge.cov[cv.id, ]
      
      train.FNN <- FNN_sp_V2_dat(resp = dat_train$Yield.c,
                                 func_cov = fundat_train,
                                 scalar_cov = as.matrix(dat_train[,c(2:ncol(dat_train))]),
                                 spatial_basis= spatial_train_coef,
                                 spatial_basis_cov = spatial_train_cov,
                                 basis_choice = c("fourier","fourier" ),
                                 num_basis = c(nbasis,nbasis), 
                                 hidden_layers = n.hidden.layer,
                                 neurons_per_layer = n.neuron.layer,
                                 activations_in_layers = link.ftn,
                                 domain_range = list(range(timepts), range(timepts)),
                                 epochs = 500,
                                 output_size = 1,
                                 loss_choice = "mse",
                                 metric_choice = list("mean_squared_error"),
                                 val_split = 0.1,  # original code 0.15 oe 0.2
                                 learn_rate = 0.002,
                                 patience_param = 15,
                                 early_stop = T, # original code TRUE
                                 print_info = F,
                                 batch_size = 32
      )
      
      # Predicting
      pred_fnn = FNN_Predict_sp_V2_dat(train.FNN,
                                       fundat_test,
                                       scalar_cov = as.matrix(dat_test[,c(2:ncol(dat_test))]),
                                       spatial_basis= spatial_test_coef,
                                       spatial_basis_cov = spatial_test_cov,
                                       basis_choice = c("fourier","fourier"),
                                       num_basis = c(nbasis,nbasis),
                                       domain_range = list(range(timepts), range(timepts))
      )
    }else if(kriging==FALSE){
      train.FNN <- FNN_sp_V3_dat(resp = dat_train$Yield.c,
                                 func_cov = fundat_train,
                                 scalar_cov = as.matrix(dat_train[,c(2:ncol(dat_train))]),
                                 spatial_basis= spatial_train_coef,
                                 #spatial_basis_cov = spatial_train_cov,
                                 basis_choice = c("fourier","fourier" ),
                                 num_basis = c(nbasis,nbasis), 
                                 hidden_layers = n.hidden.layer,
                                 neurons_per_layer = n.neuron.layer,
                                 activations_in_layers = link.ftn,
                                 domain_range = list(range(timepts), range(timepts)),
                                 epochs = 500,
                                 output_size = 1,
                                 loss_choice = "mse",
                                 metric_choice = list("mean_squared_error"),
                                 val_split = 0.1,  # original code 0.15 oe 0.2
                                 learn_rate = 0.002,
                                 patience_param = 15,
                                 early_stop = T, # original code TRUE
                                 print_info = F,
                                 batch_size = 32
      )
      
      # Predicting
      pred_fnn = FNN_Predict_sp_V3_dat(train.FNN,
                                       fundat_test,
                                       scalar_cov = as.matrix(dat_test[,c(2:ncol(dat_test))]),
                                       spatial_basis= spatial_test_coef,
                                       #spatial_basis_cov = spatial_test_cov,
                                       basis_choice = c("fourier","fourier"),
                                       num_basis = c(nbasis,nbasis),
                                       domain_range = list(range(timepts), range(timepts))
      )
      
    }

  }else if(spatial==FALSE){
    train.FNN <- FNN(resp = dat_train$Yield.c,
                     func_cov = fundat_train,
                     scalar_cov = as.matrix(dat_train[,c(2:ncol(dat_train))]),
                     basis_choice = c("fourier","fourier"),
                     num_basis = c(nbasis,nbasis), 
                     hidden_layers = n.hidden.layer,
                     neurons_per_layer = n.neuron.layer,
                     activations_in_layers = link.ftn,
                     domain_range = list(range(timepts), range(timepts)),
                     epochs = 500,
                     output_size = 1,
                     loss_choice = "mse",
                     metric_choice = list("mean_squared_error"),
                     val_split = 0.15,
                     learn_rate = 0.002,
                     patience_param = 15,
                     early_stop = T,
                     print_info = F)    

    
    # Predicting
    pred_fnn = FNN_Predict(train.FNN,
                           fundat_test,
                           scalar_cov = as.matrix(dat_test[,c(2:ncol(dat_test))]),
                           basis_choice = c("fourier","fourier"),
                           num_basis = c(nbasis,nbasis),
                           domain_range = list(range(timepts), range(timepts))
    )  
    }
 
  pr.test = pred_fnn

  pred.list[[vv]] <- pr.test
  true.list[[vv]] <- dat_test$Yield.c
  
  (test.err = mean((dat_test$Yield.c - pr.test)^2, na.rm = T))
  (w.test.err = sum(weight_test*(dat_test$Yield.c - pr.test)^2/sum(weight_test) ))
  
  print(paste0(vv,"th: ",test.err, w.test.err))
  pe= c(pe, test.err)
  w.pe = c(w.pe, w.test.err)
}

save(pred.list, true.list, pe, w.pe, 
     file="...")



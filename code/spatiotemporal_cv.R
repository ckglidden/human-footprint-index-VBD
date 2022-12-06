##Spatio-temporal cross-validation for: Human footprint is associated with shifts in assemblages of major vector-borne diseases
##Author: Caroline Glidden & Eloise Skinner
##Date: Dec 6, 2022

##this code uses spatiotemporal cross validation to evaluate model performance
##time step folds are separated using an expanding window
##spatial folds are separated using k-folds clustering

###load packages
library("randomForestSRC"); library("caTools"); library("sf"); library("spatialsample"); library("tidyr")

##----------------------##
##load & set up data    ##       
##----------------------##
data <- read.csv("../data/full_dataset.csv")

#Brail muncipality coordinates for identifying 15 spatial folds
#the below shape file is too large to load in github but can be downloaded here: 
#https://www.ibge.gov.br/en/geosciences/territorial-organization/territorial-meshes/18890-municipal-mesh.html?=&t=acesso-ao-produto
mun <- read_sf("../data/br_municipios_20200807")

#identify 15 spatial folds using k-means clustering
spatial_split <- spatial_clustering_cv(mun, v = 10)

splits_df <- c()
for(i in 1:15){
  
  new_df <- assessment(spatial_split$splits[[i]])
  new_df$fold <- i
  new_df <- new_df[,c("CD_MUN", "fold")]
  
  splits_df <- rbind(splits_df, new_df)
  
}
splits_df <- st_drop_geometry(splits_df)
splits_df$CD_MUN = substr(splits_df$CD_MUN,1,nchar(splits_df$CD_MUN)-1)

#final data
data <- merge(data, splits_df, by = "CD_MUN")


#include temporal expanding window
data$year_fold_test <- ifelse(data$year == 2017, 1,
                              ifelse(data$year == 2018, 2,
                                     ifelse(data$year == 2019, 3)))

data$year_fold_train <- ifelse(data$year < 2017, 1,
                               ifelse(data$year < 2018, 2,
                                      ifelse(data$year < 2019, 3)))

#log transform variables
data$pop_log <- log(data$population)
data$cropland_log <- log(data$cropland + 0.01)
data$urban_log <- log(data$urban + 0.01)

#set up occurrence - convert incidence to 0/1 factor (1 = incidence > 0)
data$chik <- as.factor(ifelse(data$chik_incidence > 0, 1, 0))
data$zika <- as.factor(ifelse(data$zika_incidence > 0, 1, 0))
data$dengue <- as.factor(ifelse(data$denv_incidence > 0, 1, 0))
data$malaria <- as.factor(ifelse(data$malaria_incidence > 0, 1, 0))
data$cutaneous_leish <- as.factor(ifelse(data$c_leish_incidence > 0, 1, 0))
data$visceral_leish <- as.factor(ifelse(data$v_leish_incidence > 0, 1, 0))

##-------------------------------------------------------------##
##for loop for running flavi virus models through spatial_cv   ##       
##-------------------------------------------------------------##

oob_out_flavi <- data.frame(matrix(vector(), 0, 7,
                              dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type", "spatial_fold", "temporal_fold"))), 
                       stringsAsFactors=F, row.names=NULL)
 
oob_in_flavi <- data.frame(matrix(vector(), 0, 7,
                             dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type", "spatial_fold", "temporal_fold"))), 
                      stringsAsFactors=F, row.names=NULL)
 
 
flavi <- c("chik", "zika", "dengue") #based on y0
 
for(k in 1:3){
    
  #specify flavi-virus
   flavi_formula = as.formula(noquote(paste0(flavi[k],"~human_footprint+cropland_log+forest+pasture+pop_log+annual_tmp+annual_pre+annual_wet")))
 
   oob_out_patho <- data.frame(matrix(vector(),0, 7,
                                     dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type", "spatial_fold", "temporal_fold"))),
                              stringsAsFactors=F, row.names=NULL)

  oob_in_patho <- data.frame(matrix(vector(),0, 7,
                                    dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type", "spatial_fold", "temporal_fold"))),
                             stringsAsFactors=F, row.names=NULL)
  for(j in 1:3){

    #subset temporal folds here
    temporal_train <- subset(data, year_fold_train == j)
    temporal_test <- subset(data, year_fold_test == j)

    oob_out_cv <- data.frame(matrix(vector(),0, 7,
                                       dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type", "spatial_fold", "temporal_fold"))),
                                stringsAsFactors=F, row.names=NULL)

    oob_in_cv <- data.frame(matrix(vector(),0, 7,
                                      dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type", "spatial_fold", "temporal_fold"))),
                               stringsAsFactors=F, row.names=NULL)

  for(i in 1:15){

    set.seed(i*99)
    
  #subset spatial folds here
    datrain = subset(temporal_train, fold != i)
    datest  = subset(temporal_test, fold == i)

     if(table(datrain[,flavi[k]])[2] > 0 & table(datest[,flavi[k]])[2] > 0){ #run model if there are greater than 0 observations within each fold

    #tune model
    o <- randomForestSRC::tune(formula = flavi_formula, data = datrain, rfq = TRUE, perf.type = "gmean")

    #out-of-sample test
    rf0_imb_rfq=imbalanced(flavi_formula,
                                ntree=500,
                                data=datrain,
                                mtry = as.numeric(o$optimal[2]),
                                nodesize = as.numeric(o$optimal[1]),
                                method = "rfq",
                                do.trace=FALSE,
                                importance="random",
                                statistics = T)

    oob_out_0 <- predict(rf0_imb_rfq, newdata=datest)
    auc_out <- pROC::roc(response = oob_out_0$yvar, predictor= oob_out_0$predicted[,2], levels=c(0,1), auc = TRUE)
    best_threshold_out <- pROC::coords(auc_out, "best", ret = c("threshold", "sensitivity", "specificity"))
    sensitivity_out <- best_threshold_out$sensitivity
    specificity_out <- best_threshold_out$specificity
    out_error <- data.frame(Type = flavi[k], oob = oob_out_0$err.rate[500,1],
                            sens = sensitivity_out,
                            spec = specificity_out,
                            auc = auc_out$auc,
                            spatial_fold = i,
                            temporal_fold = j)

    oob_out_cv <- rbind(oob_out_cv, out_error)

    #in sample test
    in_sample_test <- imbalanced(flavi_formula,
                           ntree=500,
                           data=data,
                           mtry = as.numeric(o$optimal[2]),
                           nodesize = as.numeric(o$optimal[1]),
                           method = "rfq",
                           do.trace=FALSE,
                           importance="random",
                           statistics = T)

    oob_in_0 <- predict(in_sample_test, newdata=data)
    auc_in <- pROC::roc(response = oob_in_0$yvar, predictor= oob_in_0$predicted[,2], levels=c(0,1), auc = TRUE)
    best_threshold_in <- pROC::coords(auc_in, "best", ret = c("threshold", "sensitivity", "specificity"))
    sensitivity_in <- best_threshold_in$sensitivity
    specificity_in <- best_threshold_in$specificity

    in_error <- data.frame(Type = flavi[k], oob = oob_in_0$err.rate[500,1],
                           sens = sensitivity_in,
                           spec = specificity_in,
                           auc = auc_in$auc,
                           spatial_fold = i,
                           temporal_fold = j)
    oob_in_cv <- rbind(oob_in_cv, in_error)

    #print to track progress in case code breaks
    message <- paste0("iteration ",flavi[k],"-",j,"-",i)
    print(message)

    } else {

      out_error <- data.frame(Type = flavi[k], oob = NA,
                              sens = NA,
                              spec = NA,
                              auc = NA,
                              spatial_fold = i,
                              temporal_fold = j)

      oob_out_cv <- rbind(oob_out_cv, out_error)

      in_error <- data.frame(Type = flavi[k], oob = NA,
                             sens = NA,
                             spec = NA,
                             auc = NA,
                             spatial_fold = i,
                             temporal_fold = j)
      oob_in_cv <- rbind(oob_in_cv, in_error)

      message <- paste0("iteration ",flavi[k],"-",j,"-",i)
      print(message)

    }

  }

    ###add in temporal data binding
    oob_out_patho <- rbind(oob_out_patho, oob_out_cv)
    oob_in_patho <- rbind(oob_in_patho, oob_in_cv)

    flavi_message <- paste0(flavi[k], "done")
    print(flavi_message)

  }

  oob_out_flavi <- rbind(oob_out_flavi, oob_out_patho)
  oob_in_flavi <- rbind(oob_in_flavi, oob_in_patho)

}

#save results in case code below breaks
write.csv(oob_out_flavi, "spatial_cv_flavi_out_error.csv")
write.csv(oob_in_flavi, "spatial_cv_flavi_in_error.csv")

##---------------------------------------##
#for loop for malaria model              ##
##---------------------------------------##


ml_oob_out <- data.frame(matrix(vector(), 0, 7,
                                dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type", "spatial_fold", "temporal_fold"))), 
                         stringsAsFactors=F, row.names=NULL)

ml_oob_in <- data.frame(matrix(vector(), 0, 7,
                               dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type", "spatial_fold", "temporal_fold"))), 
                        stringsAsFactors=F, row.names=NULL)

ml_formula = malaria~human_footprint+cropland_log+forest+pasture+pop_log+annual_tmp+annual_pre+annual_wet

for(j in 1:3){
  
  temporal_train <- subset(data, year_fold_train == j)
  temporal_test <- subset(data, year_fold_test == j)
  
  oob_out_cv <- data.frame(matrix(vector(),0, 7,
                                  dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type", "spatial_fold", "temporal_fold"))), 
                           stringsAsFactors=F, row.names=NULL)
  
  oob_in_cv <- data.frame(matrix(vector(),0, 7,
                                 dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type", "spatial_fold", "temporal_fold"))), 
                          stringsAsFactors=F, row.names=NULL)
  
  for(i in 1:15){
    
    set.seed(i*99)
    
    datrain = subset(temporal_train, fold != i)
    datest  = subset(temporal_test, fold == i)
  
    if(table(datrain$malaria)[2] > 0 & table(datest$malaria)[2] > 0){
  
  #tune model
  o <- tune(ml_formula, data = datrain, rfq = TRUE, perf.type = "gmean")
  
  #out-of-sample test
  rf0_imb_rfq=imbalanced(ml_formula, 
                         ntree=500, 
                         data=datrain,
                         mtry = as.numeric(o$optimal[2]),
                         nodesize = as.numeric(o$optimal[1]),
                         method = "rfq",
                         do.trace=FALSE, 
                         importance="random", 
                         statistics = T)
  
  oob_out_0 <- predict(rf0_imb_rfq, newdata=datest)
  auc_out <- pROC::roc(response = oob_out_0$yvar, predictor= oob_out_0$predicted[,2], levels=c(0,1), auc = TRUE)
  best_threshold_out <- pROC::coords(auc_out, "best", ret = c("threshold", "sensitivity", "specificity"))
  sensitivity_out <- best_threshold_out$sensitivity
  specificity_out <- best_threshold_out$specificity
  out_error <- data.frame(Type = 'malaria', oob = oob_out_0$err.rate[500,1], 
                          sens = sensitivity_out,
                          spec = specificity_out,
                          auc = auc_out$auc,
                          spatial_fold = i,
                          temporal_fold = j)
  oob_out_cv <- rbind(oob_out_cv, out_error)
  
  #in sample test
  in_sample_test <- imbalanced(ml_formula, 
                               ntree=500, 
                               data=data,
                               mtry = as.numeric(o$optimal[2]),
                               nodesize = as.numeric(o$optimal[1]),
                               method = "rfq",
                               do.trace=FALSE, 
                               importance="random", 
                               statistics = T)
  
  oob_in_0 <- predict(in_sample_test, newdata=data)
  auc_in <- pROC::roc(response = oob_in_0$yvar, predictor= oob_in_0$predicted[,2], levels=c(0,1), auc = TRUE)
  best_threshold_in <- pROC::coords(auc_in, "best", ret = c("threshold", "sensitivity", "specificity"))
  sensitivity_in <- best_threshold_in$sensitivity
  specificity_in <- best_threshold_in$specificity
  
  in_error <- data.frame(Type = 'malaria', oob = oob_in_0$err.rate[500,1],
                         sens = sensitivity_in,
                         spec = specificity_in,
                         auc = auc_in$auc,
                         spatial_fold = i,
                         temporal_fold = j)
  oob_in_cv <- rbind(oob_in_cv, in_error)
  
  message <- paste0("iteration malaria-",j,"-",i)
  print(message)
  
    } else {
      
      out_error <- data.frame(Type = "malaria", oob = NA, 
                              sens = NA,
                              spec = NA,
                              auc = NA,
                              spatial_fold = i,
                              temporal_fold = j)
      
      oob_out_cv <- rbind(oob_out_cv, out_error)
      
      in_error <- data.frame(Type = "malaria", oob = NA, 
                             sens = NA,
                             spec = NA,
                             auc = NA,
                             spatial_fold = i,
                             temporal_fold = j)
      oob_in_cv <- rbind(oob_in_cv, in_error)
      
      message <- paste0("iteration malaria-",j,"-",i)
      print(message)
    }
      
      
      }
    
  ##add temporal data
  ml_oob_out <- rbind(ml_oob_out, oob_out_cv)
  ml_oob_in <- rbind(ml_oob_in, oob_in_cv)
    
  }
  
print("malaria done")

##----------------------------------------##
##for loop for cutaneous leish model      ##
##----------------------------------------##

cl_oob_out <- data.frame(matrix(vector(), 0, 7,
                                dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type", "spatial_fold", "temporal_fold"))), 
                         stringsAsFactors=F, row.names=NULL)

cl_oob_in <- data.frame(matrix(vector(), 0, 7,
                               dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type", "spatial_fold", "temporal_fold"))), 
                        stringsAsFactors=F, row.names=NULL)

cl_formula = cutaneous_leish~human_footprint+cropland_log+forest+pasture+pop_log+annual_tmp+annual_pre+annual_wet


for(j in 1:3){
  
  temporal_train <- subset(data, year_fold_train == j)
  temporal_test <- subset(data, year_fold_test == j)
  
  oob_out_cv <- data.frame(matrix(vector(),0, 7,
                                  dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type", "spatial_fold", "temporal_fold"))), 
                           stringsAsFactors=F, row.names=NULL)
  
  oob_in_cv <- data.frame(matrix(vector(),0, 7,
                                 dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type", "spatial_fold", "temporal_fold"))), 
                          stringsAsFactors=F, row.names=NULL)
  
  for(i in 1:15){
    
    set.seed(i*99)
    
    #sample = sample.split(data[,flavi[k]], SplitRatio = .8)
    datrain = subset(temporal_train, fold != i)
    datest  = subset(temporal_test, fold == i)
  
    if(table(datrain$cutaneous_leish)[2] > 0 & table(datest$cutaneous_leish)[2] > 0){
      
  #tune
  o <- tune(cl_formula, data = datrain, rfq = TRUE, perf.type = "gmean")
  
  #out-of-sample test
  rf0_imb_rfq=imbalanced(cl_formula, 
                         ntree=500, 
                         data=datrain,
                         mtry = as.numeric(o$optimal[2]),
                         nodesize = as.numeric(o$optimal[1]),
                         method = "rfq",
                         do.trace=FALSE, 
                         importance="random", 
                         statistics = T)
  
  oob_out_0 <- predict(rf0_imb_rfq, newdata=datest)
  auc_out <- pROC::roc(response = oob_out_0$yvar, predictor= oob_out_0$predicted[,2], levels=c(0,1), auc = TRUE)
  best_threshold_out <- pROC::coords(auc_out, "best", ret = c("threshold", "sensitivity", "specificity"))
  sensitivity_out <- best_threshold_out$sensitivity
  specificity_out <- best_threshold_out$specificity
  out_error <- data.frame(Type = 'cutaneous_leish', oob = oob_out_0$err.rate[500,1], 
                          sens = sensitivity_out,
                          spec = specificity_out,
                          auc = auc_out$auc,
                          spatial_fold = i,
                          temporal_fold = j)
  oob_out_cv <- rbind(oob_out_cv, out_error)
  
  in_sample_test <- imbalanced(cl_formula, 
                               ntree=500, 
                               data=data,
                               mtry = as.numeric(o$optimal[2]),
                               nodesize = as.numeric(o$optimal[1]),
                               method = "rfq",
                               do.trace=FALSE, 
                               importance="random", 
                               statistics = T)
  
  oob_in_0 <- predict(in_sample_test, newdata=data)
  auc_in <- pROC::roc(response = oob_in_0$yvar, predictor= oob_in_0$predicted[,2], levels=c(0,1), auc = TRUE)
  best_threshold_in <- pROC::coords(auc_in, "best", ret = c("threshold", "sensitivity", "specificity"))
  sensitivity_in <- best_threshold_in$sensitivity
  specificity_in <- best_threshold_in$specificity
  
  #in sample test
  in_error <- data.frame(Type = 'cutaneous_leish', oob = oob_in_0$err.rate[500,1], 
                         sens = sensitivity_in,
                         spec = specificity_in,
                         auc = auc_in$auc,
                         spatial_fold = i,
                         temporal_fold = j)
  oob_in_cv <- rbind(oob_in_cv, in_error)
  
  message <- paste0("iteration cl-",j,"-",i)
  print(message)
  
    } else {
      
      out_error <- data.frame(Type = "cutaneous_leish", oob = NA, 
                              sens = NA,
                              spec = NA,
                              auc = NA,
                              spatial_fold = i,
                              temporal_fold = j)
      
      oob_out_cv <- rbind(oob_out_cv, out_error)
      
      in_error <- data.frame(Type = "cutaneous_leish", oob = NA, 
                             sens = NA,
                             spec = NA,
                             auc = NA,
                             spatial_fold = i,
                             temporal_fold = j)
      oob_in_cv <- rbind(oob_in_cv, in_error)
      
      message <- paste0("iteration cl-",j,"-",i)
      print(message)
      
    }
  
}

  #add temporal folds
  cl_oob_out <- rbind(cl_oob_out, oob_out_cv)
  cl_oob_in <- rbind(cl_oob_in, oob_in_cv)
  
  
}

print("cl done")


##------------------------------------##
##for loop for visceral leish model   ##
##------------------------------------##

vl_oob_out <- data.frame(matrix(vector(), 0, 5,
                                dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type"))), 
                         stringsAsFactors=F, row.names=NULL)

vl_oob_in <- data.frame(matrix(vector(), 0, 5,
                               dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type"))), 
                        stringsAsFactors=F, row.names=NULL)

vl_formula = visceral_leish~human_footprint+cropland_log+forest+pasture+pop_log+annual_tmp+annual_pre+annual_wet

for(j in 1:3){
  
  temporal_train <- subset(data, year_fold_train == j)
  temporal_test <- subset(data, year_fold_test == j)
  
  oob_out_cv <- data.frame(matrix(vector(),0, 7,
                                  dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type", "spatial_fold", "temporal_fold"))), 
                           stringsAsFactors=F, row.names=NULL)
  
  oob_in_cv <- data.frame(matrix(vector(),0, 7,
                                 dimnames=list(c(), c("auc", "sens", "spec", "oob", "Type", "spatial_fold", "temporal_fold"))), 
                          stringsAsFactors=F, row.names=NULL)
  
  for(i in 1:15){
    
    set.seed(i*99)
    
    datrain = subset(temporal_train, fold != i)
    datest  = subset(temporal_test, fold == i)
  
    if(table(datrain$visceral_leish)[2] > 0 & table(datest$visceral_leish)[2] > 0){
  
  #tune
  o <- tune(vl_formula, data = datrain, rfq = TRUE, perf.type = "gmean")
  
  #out-of-sample test
  rf0_imb_rfq=imbalanced(vl_formula, 
                         ntree=500, 
                         data=datrain,
                         mtry = as.numeric(o$optimal[2]),
                         nodesize = as.numeric(o$optimal[1]),
                         method = "rfq",
                         do.trace=FALSE, 
                         importance="random", 
                         statistics = T)

  
  oob_out_0 <- predict(rf0_imb_rfq, newdata=datest)
  auc_out <- pROC::roc(response = oob_out_0$yvar, predictor= oob_out_0$predicted[,2], levels=c(0,1), auc = TRUE)
  best_threshold_out <- pROC::coords(auc_out, "best", ret = c("threshold", "sensitivity", "specificity"))
  sensitivity_out <- best_threshold_out$sensitivity
  specificity_out <- best_threshold_out$specificity
  out_error <- data.frame(Type = 'visceral_leish', oob = oob_out_0$err.rate[500,1], 
                          sens = sensitivity_out,
                          spec = specificity_out,
                          auc = auc_out$auc,
                          spatial_fold = i,
                          temporal_fold = j)
  oob_out_cv <- rbind(oob_out_cv, out_error)
  
  #in sample test
  in_sample_test <- imbalanced(vl_formula, 
                               ntree=500, 
                               data=data,
                               mtry = as.numeric(o$optimal[2]),
                               nodesize = as.numeric(o$optimal[1]),
                               method = "rfq",
                               do.trace=FALSE, 
                               importance="random", 
                               statistics = T)
  
  oob_in_0 <- predict(in_sample_test, newdata=data)
  auc_in <- pROC::roc(response = oob_in_0$yvar, predictor= oob_in_0$predicted[,2], levels=c(0,1), auc = TRUE)
  best_threshold_in <- pROC::coords(auc_in, "best", ret = c("threshold", "sensitivity", "specificity"))
  sensitivity_in <- best_threshold_in$sensitivity
  specificity_in <- best_threshold_in$specificity
  
  in_error <- data.frame(Type = 'visceral_leish', oob = oob_in_0$err.rate[500,1],
                         sens = sensitivity_in,
                         spec = specificity_in,
                         auc = auc_in$auc,
                         spatial_fold = i,
                         temporal_fold = j)
  oob_in_cv <- rbind(oob_in_cv, in_error)
  
  message <- paste0("iteration vl-",j,"-",i)
  print(message)
  
    } else {
      
      out_error <- data.frame(Type = "visceral_leish", oob = NA, 
                              sens = NA,
                              spec = NA,
                              auc = NA,
                              spatial_fold = i,
                              temporal_fold = j)
      
      oob_out_cv <- rbind(oob_out_cv, out_error)
      
      in_error <- data.frame(Type = "visceral_leish", oob = NA, 
                             sens = NA,
                             spec = NA,
                             auc = NA,
                             spatial_fold = i,
                             temporal_fold = j)
      oob_in_cv <- rbind(oob_in_cv, in_error)
      
      message <- paste0("iteration vl-",j,"-",i)
      print(message)
    }
  
}

  #add temporal folds
  vl_oob_out <- rbind(vl_oob_out, oob_out_cv)
  vl_oob_in <- rbind(vl_oob_in, oob_in_cv)
  
}

print("vl done")



##bind and save all performance metrics
perf_out <- rbind(ml_oob_out, cl_oob_out, vl_oob_out); write.csv(perf_out, "spatial_temporal_cv_out_performance.csv")
perf_in <- rbind(ml_oob_in, cl_oob_in, vl_oob_in); write.csv(perf_in, "spatial_temporal_cv_in_performance.csv")


##get mean table
# perf_out <- read.csv("/Users/carolineglidden/Desktop/forElle/bootstrapping_no_vector/revised code/spatial_temporal_cv_out_performance.csv")
# perf_flavi <- read.csv("/Users/carolineglidden/Desktop/forElle/bootstrapping_no_vector/revised code/spatial_temporal_cv_flavi_out_error.csv")
# perf_out <- rbind(perf_out, perf_flavi)
# 
# oob_performance <- aggregate(oob ~ Type, data = perf_out, 
#                              FUN = function(x) c(mean = Rmisc::CI(x, 0.95)[2], 
#                                                  lowerCI = Rmisc::CI(x, 0.95)[3], 
#                                                  upperCI = Rmisc::CI(x, 0.95)[1]))
# 
# sens_performance <- aggregate(sens ~ Type, data = perf_out, 
#                              FUN = function(x) c(mean = Rmisc::CI(x, 0.95)[2], 
#                                                  lowerCI = Rmisc::CI(x, 0.95)[3], 
#                                                  upperCI = Rmisc::CI(x, 0.95)[1]))
# 
# 
# spec_performance <- aggregate(spec ~ Type, data = perf_out, 
#                               FUN = function(x) c(mean = Rmisc::CI(x, 0.95)[2], 
#                                                   lowerCI = Rmisc::CI(x, 0.95)[3], 
#                                                   upperCI = Rmisc::CI(x, 0.95)[1]))
# 
# auc_performance <- aggregate(auc ~ Type, data = perf_out, 
#                               FUN = function(x) c(mean = Rmisc::CI(x, 0.95)[2], 
#                                                   lowerCI = Rmisc::CI(x, 0.95)[3], 
#                                                   upperCI = Rmisc::CI(x, 0.95)[1]))



























######get confidence intervals for inflection points
# hfi <- subset(all_data, variables == "human_footprint")
# 
# inflection_df <- data.frame(matrix(vector(), 0, 3,
#                                    dimnames=list(c(), c("Type", "iter", "inflection"))), stringsAsFactors=F,
#                             row.names=NULL)
# 
# for(i in 1:6){
#   
#   patho <- subset(hfi, Type == unique(hfi$Type)[i])
#   
#   patho_inflection_df <- data.frame(matrix(vector(), 0, 3,
#                                            dimnames=list(c(), c("Type", "iter", "inflection"))), stringsAsFactors=F,
#                                     row.names=NULL)
#   for(j in 1:50){
#     
#     iteration <- subset(patho, iter == j)
#     
#     infun <- approxfun(iteration$prob_adj, iteration$variable_value)
#     
#     df0 <- data.frame(Type = unique(hfi$Type)[i],
#                       iter = j,
#                       inflection = infun(0.5))
#     
#     patho_inflection_df <- rbind(patho_inflection_df, df0)
#     
#   }
#   
#   inflection_df  <- rbind(inflection_df, patho_inflection_df)
#   
# }
# 
# mean_inflection <- aggregate(inflection ~ Type, data = inflection_df, 
#                              FUN = function(x) c(mean = Rmisc::CI(x, 0.95)[2], 
#                                                  lowerCI = Rmisc::CI(x, 0.95)[3], 
#                                                  upperCI = Rmisc::CI(x, 0.95)[1]))
# 
# mean_inflection <- do.call(data.frame, mean_inflection)
# names(mean_inflection) <- c("Type", "inflection_mean", "lower_ci", "upper_ci")
# write.csv(mean_inflection, "mean_inflections.csv")
# write.csv(inflection_df, "bootstrap_infecltions.csv")
# 
# ########supplemental figures
# all_data <- read.csv("vbd_hfi_pdp_all_pathos_noVect.csv")
# all_data$Type[all_data$Type == "chik"] <- "Chikungunya"
# all_data$Type[all_data$Type == "zika"] <- "Zika"
# all_data$Type[all_data$Type == "dengue"] <- "Dengue"
# all_data$Type[all_data$Type == "malaria"] <- "Malaria"
# all_data$Type[all_data$Type == "cutaneous_leish"] <- "Cutaneous Leishmanisis"
# all_data$Type[all_data$Type == "visceral_leish"] <- "Visceral Leishmanisis"
# 
# colors0 <- cmocean::cmocean("phase")(25)
# colors <- c("darkorange","cyan3", "red3", "darkorchid", "#4887D5FF", "gold")
# 
# #might need to read in data here
# library(ggplot2)
# 
# human_footprint_supp <- ggplot(subset(all_data, variables == 'human_footprint'), aes(x=variable_value, y=prob_1, color=Type)) + 
#   stat_smooth(aes(group=iter), color='lightgrey', method='loess', size=0.5, se=FALSE) + 
#   stat_smooth(aes(), method='loess', size=1, se=FALSE) +
#   scale_color_manual(values = colors) + 
#   #scale_color_brewer(palette="Dark2") +
#   ylab('probability') + xlab("Human Footprint") +
#   facet_wrap(~Type, scales='free', ncol=3) +
#   theme_bw(base_size = 12) + 
#   theme(legend.position = "none")
# 
# temp_supp <- ggplot(subset(all_data, variables == 'annual_tmp'), aes(x=variable_value, y=prob_1, color=Type)) + 
#   stat_smooth(aes(group=iter), color='lightgrey', method='loess', size=0.5, se=FALSE) + 
#   stat_smooth(aes(), method='loess', size=1, se=FALSE) +
#   scale_color_manual(values = colors) + 
#   #scale_color_brewer(palette="Dark2") +
#   ylab('probability') + xlab("Annual Temperature") +
#   facet_wrap(~Type, scales='free', ncol=3) +
#   theme_bw(base_size = 12) + 
#   theme(legend.position = "none")
# 
# 
# population_supp <- ggplot(subset(all_data, variables == 'pop_log'), aes(x=variable_value, y=prob_1, color=Type)) + 
#   stat_smooth(aes(group=iter), color='lightgrey', method='loess', size=0.5, se=FALSE) + 
#   stat_smooth(aes(), method='loess', size=1, se=FALSE) +
#   scale_color_manual(values = colors) + 
#   #scale_color_brewer(palette="Dark2") +
#   ylab('probability') + xlab("log(population size)") +
#   facet_wrap(~Type, scales='free', ncol=3) +
#   theme_bw(base_size = 12) + 
#   theme(legend.position = "none")
# 
# forest_supp <- ggplot(subset(all_data, variables == 'forest'), aes(x=variable_value, y=prob_1, color=Type)) + 
#   stat_smooth(aes(group=iter), color='lightgrey', method='loess', size=0.5, se=FALSE) + 
#   stat_smooth(aes(), method='loess', size=1, se=FALSE) +
#   scale_color_manual(values = colors) + 
#   #scale_color_brewer(palette="Dark2") +
#   ylab('probability') + xlab("% forest cover") +
#   facet_wrap(~Type, scales='free', ncol=3) +
#   theme_bw(base_size = 12) + 
#   theme(legend.position = "none")
# 
# ggsave("vbd_hfi_plots/hfi_pdp_supp_figure.png", human_footprint_supp, dpi = 600, units = 'in', width=6.5, height = 5)
# ggsave("vbd_hfi_plots/tmp_pdp_supp_figure.png", temp_supp, dpi = 600, units = 'in', width=6.5, height = 5)
# ggsave("vbd_hfi_plots/pop_pdp_supp_figure.png", population_supp, dpi = 600, units = 'in', width=6.5, height = 5)
# ggsave("vbd_hfi_plots/forest_pdp_supp_figure.png", forest_supp, dpi = 600, units = 'in', width=6.5, height = 5)
# 
# 
# ###overlay plots
# #inflections <- read.csv("inflections.csv")
# 
# hfi_0 <- ggplot(subset(all_data, variables == 'human_footprint'), aes(x=variable_value, y=prob_adj, color=Type)) +
#   geom_smooth(aes(group=interaction(iter, Type)), stat="smooth", fullrange = TRUE,  method="gam", size=0.25, se=FALSE, alpha=0.25) + 
#   geom_smooth(aes(),method="gam", size=0.75, se=FALSE, fullrange = TRUE) +
#   scale_color_manual(values = colors) + 
#   #scale_color_brewer(palette="Dark2") +
#   xlab("Human footprint index")+
#   ylab("Rescaled probabilty of occurrence")+
#   coord_cartesian(ylim=c(0,1)) +
#   ylim(0, 1.2) +
#   #geom_hline(yintercept = 0.5) + 
#   theme_classic()+
#   ggtitle('A.') +
#   theme(plot.title = element_text(size=10),
#         axis.title=element_text(size=10),
#         axis.text = element_text(colour = "black", size = 10),
#         legend.position = "bottom",
#         legend.title=element_blank(),
#         strip.text.x = element_text(size = 10),
#         legend.text=element_text(size=10))
# 
# hfi <- hfi_0 + scale_linetype_manual(values = c("solid", "solid", "solid", "solid", "solid", "solid")) + 
#   #scale_color_brewer(palette="Dark2") +
#   scale_x_continuous(breaks = seq(0,40,5)) +
#   geom_segment(aes(x=0, y = 0.5, xend=13.3, yend=0.5), linetype = "dashed", size = 0.75, colour = "grey") +
#   geom_segment(data = mean_inflection, aes(x=inflection_mean, y = 0.5, xend= inflection_mean, yend = 0), linetype = "dashed", size = 0.75)
# 
# tmp <- ggplot(subset(all_data, variables == 'annual_tmp'), aes(x=variable_value, y=prob_adj, color=Type)) +
#   geom_line(aes(group=interaction(iter, Type)), stat="smooth", method="gam", size=0.25, se=FALSE, alpha=0.25) + 
#   geom_smooth(aes(),method="gam", size=0.75, se=FALSE, fullrange = FALSE) +
#   scale_color_manual(values = colors) + 
#   #scale_color_brewer(palette="Dark2") +
#   xlab("Annual temperature")+
#   ylab("Rescaled probabilty of occurrence")+
#   coord_cartesian(ylim=c(0,1)) +
#   ylim(0, 1.2) +
#   #geom_hline(yintercept = 0.5) + 
#   theme_classic()+
#   ggtitle('B.') +
#   theme(plot.title = element_text(size=10),
#         axis.title=element_text(size=10),
#         axis.text = element_text(colour = "black", size = 10),
#         legend.position = "bottom",
#         legend.title=element_blank(),
#         strip.text.x = element_text(size = 10),
#         legend.text=element_text(size=10))
# 
# occ <- gridExtra::arrangeGrob(hfi, tmp,nrow=2)
# ggsave(file="vbd_hfi_plots/occurence_pdps.pdf", occ, width = 4.5, height = 8, units = "in", dpi=600)
# #combined_pdp <- cowplot::plot_grid(hfi, tmp, align = 'v', ncol=1)
# 
# ####supplement pdp plot
# pop <- ggplot(subset(all_data, variables == 'pop_log'), aes(x=variable_value, y=prob_adj, color=Type)) +
#   geom_line(aes(group=interaction(iter, Type)), stat="smooth", method="gam", size=0.25, se=FALSE, alpha=0.25) + 
#   geom_smooth(aes(),method="gam", size=0.75, se=FALSE, fullrange = FALSE) +
#   scale_color_manual(values = colors) + 
#   #scale_color_brewer(palette="Dark2") +
#   xlab("Log(population)")+
#   ylab("Rescaled probabilty of occurrence")+
#   coord_cartesian(ylim=c(0,1)) +
#   #ylim(0, 1) +
#   #geom_hline(yintercept = 0.5) + 
#   theme_classic()+
#   theme(plot.title = element_text(size=10),
#         axis.title=element_text(size=10),
#         axis.text = element_text(colour = "black", size = 10),
#         legend.position = "bottom",
#         legend.title=element_blank(),
#         strip.text.x = element_text(size = 12),
#         legend.text=element_text(size=12))
# ggsave(file="vbd_hfi_plots/occurence_pop_supplement_pdps.png", pop, width = 6.5, height = 4, units = "in", dpi=600)
# 
# 
# ######mean around classification error
# ###errors
# oob_out <- read.csv("flavi_out_error.csv")
# mean_flavi_oob <- aggregate(oob ~ Type, data = oob_out, 
#                              FUN = function(x) c(mean = (1- Rmisc::CI(x, 0.95)[2]), 
#                                                  lowerCI = (1 - Rmisc::CI(x, 0.95)[3]), 
#                                                  upperCI = (1 - Rmisc::CI(x, 0.95)[1])))
# 
# oob_in <- read.csv("flavi_in_error.csv")
# mean_flavi_oob <- aggregate(oob ~ Type, data = oob_in, 
#                             FUN = function(x) c(mean = (1- Rmisc::CI(x, 0.95)[2]), 
#                                                 lowerCI = (1 - Rmisc::CI(x, 0.95)[3]), 
#                                                 upperCI = (1 - Rmisc::CI(x, 0.95)[1])))
# 
# cl_out <- read.csv("cl_outError_noVect.csv"); cl_in <- read.csv("cl_inError_noVect.csv")
# mean_cl_oob <- aggregate(oob ~ Type, data = cl_in, 
#                             FUN = function(x) c(mean = (1- Rmisc::CI(x, 0.95)[2]), 
#                                                 lowerCI = (1 - Rmisc::CI(x, 0.95)[3]), 
#                                                 upperCI = (1 - Rmisc::CI(x, 0.95)[1])))
# 
# vl_out <- read.csv("vl_outError_noVect.csv"); vl_in <- read.csv("cl_inError_noVect.csv")
# mean_vl_oob <- aggregate(oob ~ Type, data = vl_in, 
#                          FUN = function(x) c(mean = (1- Rmisc::CI(x, 0.95)[2]), 
#                                              lowerCI = (1 - Rmisc::CI(x, 0.95)[3]), 
#                                              upperCI = (1 - Rmisc::CI(x, 0.95)[1])))
# 
# ml_oob_in <- read.csv("ml_inError_noVect.csv")
# mean_ml_oob <- aggregate(oob ~ Type, data = ml_oob_in, 
#                          FUN = function(x) c(mean = (1- Rmisc::CI(x, 0.95)[2]), 
#                                              lowerCI = (1 - Rmisc::CI(x, 0.95)[3]), 
#                                              upperCI = (1 - Rmisc::CI(x, 0.95)[1])))
# 
# 
# 
# 

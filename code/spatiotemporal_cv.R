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
spatial_split <- spatial_clustering_cv(mun, v = 15)

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
test_year <- c(2017, 2018, 2019)


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

##----------------------------------------------------------------------##
##for loop for running flavi virus models through spatial-temporal cv   ##       
##----------------------------------------------------------------------##

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
    temporal_train <- subset(data, year < test_year[j])
    temporal_test <- subset(data, year == test_year[j])

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
  
  temporal_train <- subset(data, year < test_year[j])
  temporal_test <- subset(data, year == test_year[j])
  
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
  
  temporal_train <- subset(data, year < test_year[j])
  temporal_test <- subset(data, year == test_year[j])
  
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
  
  temporal_train <- subset(data, year < test_year[j])
  temporal_test <- subset(data, year == test_year[j])
  
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
perf_out <- rbind(flavi_oob_out, ml_oob_out, cl_oob_out, vl_oob_out); write.csv(perf_out, "spatial_temporal_cv_out_performance.csv")
perf_in <- rbind(flavi_oob_in, ml_oob_in, cl_oob_in, vl_oob_in); write.csv(perf_in, "spatial_temporal_cv_in_performance.csv")


##get mean table
# perf_out <- read.csv("spatial_temporal_cv_out_performance.csv")
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




























##Partial dependence plots for: Human footprint is associated with shifts in assemblages of major vector-borne diseases
##Author: Caroline Glidden & Eloise Skinner
##Date: Dec 6, 2022

##this code creates partial dependence plots for each covariate in the models & calculates HFI inflection points (where probability of occurrence shifts from increasing to decreasing, or vice versa)
##bootstrapping is used to estimate uncertainty and confidence intervals for the pdps and inflection points

###load packages
library("randomForestSRC"); library("caTools")

###load data
data <- read.csv("../data/full_dataset.csv")

#log transform variables
data$pop_log <- log(data$population)
data$cropland_log <- log(data$cropland + 0.01)
data$urban_log <- log(data$urban + 0.01)

#set up occurrence
data$chik <- as.factor(ifelse(data$chik_incidence > 0, 1, 0))
data$zika <- as.factor(ifelse(data$zika_incidence > 0, 1, 0))
data$dengue <- as.factor(ifelse(data$denv_incidence > 0, 1, 0))
data$malaria <- as.factor(ifelse(data$malaria_incidence > 0, 1, 0))
data$cutaneous_leish <- as.factor(ifelse(data$c_leish_incidence > 0, 1, 0))
data$visceral_leish <- as.factor(ifelse(data$v_leish_incidence > 0, 1, 0))


######for loop for running flavi virus models
flavi_pdp <- data.frame(matrix(vector(), 0, 7,
                                      dimnames=list(c(), c("prob_1", "prob_se", "prob_adj", "variable", "variable_value",
                                                           "Type", "iter"))), stringsAsFactors=F,
                               row.names=NULL)

flavi <- c("chik", "zika", "dengue") #based on y0

for(k in 1:3){
  
  #subset data here
  flavi_formula = as.formula(noquote(paste0(flavi[k],"~human_footprint+cropland_log+forest+pasture+pop_log+annual_tmp+annual_pre+annual_wet")))

  flavi_df <- data.frame(matrix(vector(), 0, 7,
                                  dimnames=list(c(), c("prob_1", "prob_se", "prob_adj", "variable", "variable_value",
                                                       "Type", "iter"))), stringsAsFactors=F,
                           row.names=NULL)
  
  for(i in 1:50){
    
    set.seed(i*99)
    
    sample = sample.split(data[,flavi[k]], SplitRatio = .8)
    datrain = subset(data, sample == TRUE)
    datest  = subset(data, sample == FALSE)
    
    #tune model
    o <- tune(flavi_formula, data = datrain, rfq = TRUE, perf.type = "gmean")
    
    #train model on 80% of the data
    rf0_imb_rfq=imbalanced(flavi_formula, 
                           ntree=500, 
                           data=datrain,
                           mtry = as.numeric(o$optimal[2]),
                           nodesize = as.numeric(o$optimal[1]),
                           method = "rfq",
                           do.trace=FALSE, 
                           importance="random", 
                           statistics = T)

    
    #####loop through all variables
    var <- c("human_footprint", "annual_tmp", "pop_log", "forest")
    
    flavi_df_1 <- data.frame(matrix(vector(), 0, 7,
                                   dimnames=list(c(), c("prob_1", "prob_se", "prob_adj", "variable", "variable_value",
                                                       "Type", "iter"))), stringsAsFactors=F,
                            row.names=NULL)
    
    #create bind data-frames with pdp output for each variable
    for(j in 1:length(var)){
      hfp <- plot.variable(rf0_imb_rfq, xvar.names = var[j], partial = TRUE,target = "1",  npts = 25, grid.resoluion = 100)
      flavi_df_0 <- data.frame(prob_1 = hfp$pData[[c(1, 2)]],
                                prob_se = hfp$pData[[c(1, 3)]], 
                                prob_adj = (hfp$pData[[c(1, 2)]] - min(hfp$pData[[c(1, 2)]]))/(max(hfp$pData[[c(1, 2)]])-min(hfp$pData[[c(1, 2)]])),
                                variables = var[j],
                                variable_value = hfp$pData[[c(1, 5)]],
                                Type =  flavi[k],
                                iter = i)
      flavi_df_1 <- rbind(flavi_df_1, flavi_df_0)
    }
    
    flavi_df <- rbind(flavi_df, flavi_df_1)
    
    
  }
  
  flavi_pdp <- rbind(flavi_pdp, flavi_df)
  
}

write.csv(flavi_pdp, "flavi_pdps.csv")


##---------------------------------------##
#for loop for malaria model              ##
##---------------------------------------##


ml_formula = malaria~human_footprint+cropland_log+forest+pasture+pop_log+annual_tmp+annual_pre+annual_wet
  
ml_df <- data.frame(matrix(vector(), 0, 7,
                                dimnames=list(c(), c("prob_1", "prob_se", "prob_adj", "variable", "variable_value",
                                                     "Type", "iter"))), stringsAsFactors=F,row.names=NULL)
for(i in 1:50){
  
  set.seed(i*99)
  
  sample = sample.split(data$malaria, SplitRatio = .8)
  datrain = subset(data, sample == TRUE)
  datest  = subset(data, sample == FALSE)
  
  #tune model
  o <- tune(ml_formula, data = datrain, rfq = TRUE, perf.type = "gmean")
  
  #train model on 80% of the data
  rf0_imb_rfq=imbalanced(ml_formula, 
                         ntree=500, 
                         data=datrain,
                         mtry = as.numeric(o$optimal[2]),
                         nodesize = as.numeric(o$optimal[1]),
                         method = "rfq",
                         do.trace=FALSE, 
                         importance="random", 
                         statistics = T)
  
  #####loop through all variables
  var <- c("human_footprint", "annual_tmp", "pop_log", "forest")
  
  ml_df_1 <- data.frame(matrix(vector(), 0, 7,
                                  dimnames=list(c(), c("prob_1", "prob_se", "prob_adj", "variable", "variable_value",
                                                       "Type", "iter"))), stringsAsFactors=F,
                           row.names=NULL)
  
  for(j in 1:length(var)){
    hfp <- plot.variable(rf0_imb_rfq, xvar.names = var[j], partial = TRUE,target = "1",  npts = 25, grid.resoluion = 100)
    ml_df_0 <- data.frame(prob_1 = hfp$pData[[c(1, 2)]],
                             prob_se = hfp$pData[[c(1, 3)]], 
                             prob_adj = (hfp$pData[[c(1, 2)]] - min(hfp$pData[[c(1, 2)]]))/(max(hfp$pData[[c(1, 2)]])-min(hfp$pData[[c(1, 2)]])),
                             variables = var[j],
                             variable_value = hfp$pData[[c(1, 5)]],
                             Type =  "malaria",
                             iter = i)
    ml_df_1 <- rbind(ml_df_1, ml_df_0)
  }
  
  ml_df <- rbind(ml_df, ml_df_1)
  
  
}

##----------------------------------------##
##for loop for cutaneous leish model      ##
##----------------------------------------##

cl_oob_out <- data.frame(matrix(vector(), 0, 2,
                                dimnames=list(c(), c("oob","Type"))), 
                         stringsAsFactors=F, row.names=NULL)

cl_oob_in <- data.frame(matrix(vector(), 0, 2,
                               dimnames=list(c(), c("oob","Type"))), 
                        stringsAsFactors=F, row.names=NULL)

cl_formula = cutaneous_leish~human_footprint+cropland_log+forest+pasture+pop_log+annual_tmp+annual_pre+annual_wet

cl_df <- data.frame(matrix(vector(), 0, 7,
                           dimnames=list(c(), c("prob_1", "prob_se", "prob_adj", "variable", "variable_value",
                                                "Type", "iter"))), stringsAsFactors=F,row.names=NULL)
for(i in 1:50){
  
  set.seed(i*99)
  
  sample = sample.split(data$cutaneous_leish, SplitRatio = .8)
  datrain = subset(data, sample == TRUE)
  datest  = subset(data, sample == FALSE)
  
  #tune model
  o <- tune(cl_formula, data = datrain, rfq = TRUE, perf.type = "gmean")
  
  #train model on 80% of the data
  rf0_imb_rfq=imbalanced(cl_formula, 
                         ntree=500, 
                         data=datrain,
                         mtry = as.numeric(o$optimal[2]),
                         nodesize = as.numeric(o$optimal[1]),
                         method = "rfq",
                         do.trace=FALSE, 
                         importance="random", 
                         statistics = T)
  
  #####loop through all variables
  var <- c("human_footprint", "annual_tmp", "pop_log", "forest")
  
  cl_df_1 <- data.frame(matrix(vector(), 0, 7,
                               dimnames=list(c(), c("prob_1", "prob_se", "prob_adj", "variable", "variable_value",
                                                    "Type", "iter"))), stringsAsFactors=F,
                        row.names=NULL)
  
  for(j in 1:length(var)){
    hfp <- plot.variable(rf0_imb_rfq, xvar.names = var[j], partial = TRUE,target = "1",  npts = 25, grid.resoluion = 100)
    cl_df_0 <- data.frame(prob_1 = hfp$pData[[c(1, 2)]],
                          prob_se = hfp$pData[[c(1, 3)]], 
                          prob_adj = (hfp$pData[[c(1, 2)]] - min(hfp$pData[[c(1, 2)]]))/(max(hfp$pData[[c(1, 2)]])-min(hfp$pData[[c(1, 2)]])),
                          variables = var[j],
                          variable_value = hfp$pData[[c(1, 5)]],
                          Type =  "cutaneous_leish",
                          iter = i)
    cl_df_1 <- rbind(cl_df_1, cl_df_0)
  }
  
  cl_df <- rbind(cl_df, cl_df_1)
  
}


##------------------------------------##
##for loop for visceral leish model   ##
##------------------------------------##

vl_formula = visceral_leish~human_footprint+cropland_log+forest+pasture+pop_log+annual_tmp+annual_pre+annual_wet

vl_df <- data.frame(matrix(vector(), 0, 7,
                           dimnames=list(c(), c("prob_1", "prob_se", "prob_adj", "variable", "variable_value",
                                                "Type", "iter"))), stringsAsFactors=F,row.names=NULL)
for(i in 1:50){
  
  set.seed(i*99)
  
  sample = sample.split(data$visceral_leish, SplitRatio = .8)
  datrain = subset(data, sample == TRUE)
  datest  = subset(data, sample == FALSE)
  
  #tune model
  o <- tune(vl_formula, data = datrain, rfq = TRUE, perf.type = "gmean")
  
  #train model on 80% of the data
  rf0_imb_rfq=imbalanced(vl_formula, 
                         ntree=500, 
                         data=datrain,
                         mtry = as.numeric(o$optimal[2]),
                         nodesize = as.numeric(o$optimal[1]),
                         method = "rfq",
                         do.trace=FALSE, 
                         importance="random", 
                         statistics = T)
  
  #####loop through all variables
  var <- c("human_footprint", "annual_tmp", "pop_log", "forest")
  
  vl_df_1 <- data.frame(matrix(vector(), 0, 7,
                               dimnames=list(c(), c("prob_1", "prob_se", "prob_adj", "variable", "variable_value",
                                                    "Type", "iter"))), stringsAsFactors=F,
                        row.names=NULL)
  
  for(j in 1:length(var)){
    hfp <- plot.variable(rf0_imb_rfq, xvar.names = var[j], partial = TRUE,target = "1",  npts = 25, grid.resoluion = 100)
    vl_df_0 <- data.frame(prob_1 = hfp$pData[[c(1, 2)]],
                          prob_se = hfp$pData[[c(1, 3)]], 
                          prob_adj = (hfp$pData[[c(1, 2)]] - min(hfp$pData[[c(1, 2)]]))/(max(hfp$pData[[c(1, 2)]])-min(hfp$pData[[c(1, 2)]])),
                          variables = var[j],
                          variable_value = hfp$pData[[c(1, 5)]],
                          Type =  "visceral_leish",
                          iter = i)
    vl_df_1 <- rbind(vl_df_1, vl_df_0)
  }
  
  vl_df <- rbind(vl_df, vl_df_1)
  
}

#bind all pdps for plotting

all_data <- rbind(flavi_pdp, ml_df, cl_df, vl_df)

##-------------------------------------------------##
##get confidence intervals for inflection points   ##
##-------------------------------------------------##

hfi <- subset(all_data, variables == "human_footprint")

inflection_df <- data.frame(matrix(vector(), 0, 3,
                                   dimnames=list(c(), c("Type", "iter", "inflection"))), stringsAsFactors=F,
                            row.names=NULL)

for(i in 1:6){
  
  patho <- subset(hfi, Type == unique(hfi$Type)[i])
  
  patho_inflection_df <- data.frame(matrix(vector(), 0, 3,
                                           dimnames=list(c(), c("Type", "iter", "inflection"))), stringsAsFactors=F,
                                    row.names=NULL)
  for(j in 1:50){
    
    iteration <- subset(patho, iter == j)
    
    infun <- approxfun(iteration$prob_adj, iteration$variable_value)
    
    df0 <- data.frame(Type = unique(hfi$Type)[i],
                      iter = j,
                      inflection = infun(0.5))
    
    patho_inflection_df <- rbind(patho_inflection_df, df0)
    
  }
  
  inflection_df  <- rbind(inflection_df, patho_inflection_df)
  
}

mean_inflection <- aggregate(inflection ~ Type, data = inflection_df, 
                             FUN = function(x) c(mean = Rmisc::CI(x, 0.95)[2], 
                                                 lowerCI = Rmisc::CI(x, 0.95)[3], 
                                                 upperCI = Rmisc::CI(x, 0.95)[1]))

mean_inflection <- do.call(data.frame, mean_inflection)
names(mean_inflection) <- c("Type", "inflection_mean", "lower_ci", "upper_ci")
write.csv(mean_inflection, "mean_inflections.csv")
write.csv(inflection_df, "bootstrap_infecltions.csv")

##----------------------------------------##
##scaled pdp plots for main text          ##
##----------------------------------------##

hfi_0 <- ggplot(subset(all_data, variables == 'human_footprint'), aes(x=variable_value, y=prob_adj, color=Type)) +
  geom_smooth(aes(group=interaction(iter, Type)), stat="smooth", fullrange = TRUE,  method="gam", size=0.25, se=FALSE, alpha=0.25) + 
  geom_smooth(aes(),method="gam", size=0.75, se=FALSE, fullrange = TRUE) +
  scale_color_manual(values = colors) + 
  #scale_color_brewer(palette="Dark2") +
  xlab("Human footprint index")+
  ylab("Rescaled probabilty of occurrence")+
  coord_cartesian(ylim=c(0,1)) +
  ylim(0, 1.2) +
  #geom_hline(yintercept = 0.5) + 
  theme_classic()+
  ggtitle('A.') +
  theme(plot.title = element_text(size=10),
        axis.title=element_text(size=10),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "bottom",
        legend.title=element_blank(),
        strip.text.x = element_text(size = 10),
        legend.text=element_text(size=10))

hfi <- hfi_0 + scale_linetype_manual(values = c("solid", "solid", "solid", "solid", "solid", "solid")) + 
  #scale_color_brewer(palette="Dark2") +
  scale_x_continuous(breaks = seq(0,40,5)) +
  geom_segment(aes(x=0, y = 0.5, xend=13.3, yend=0.5), linetype = "dashed", size = 0.75, colour = "grey") +
  geom_segment(data = mean_inflection, aes(x=inflection_mean, y = 0.5, xend= inflection_mean, yend = 0), linetype = "dashed", size = 0.75)

tmp <- ggplot(subset(all_data, variables == 'annual_tmp'), aes(x=variable_value, y=prob_adj, color=Type)) +
  geom_line(aes(group=interaction(iter, Type)), stat="smooth", method="gam", size=0.25, se=FALSE, alpha=0.25) + 
  geom_smooth(aes(),method="gam", size=0.75, se=FALSE, fullrange = FALSE) +
  scale_color_manual(values = colors) + 
  #scale_color_brewer(palette="Dark2") +
  xlab("Annual temperature")+
  ylab("Rescaled probabilty of occurrence")+
  coord_cartesian(ylim=c(0,1)) +
  ylim(0, 1.2) +
  #geom_hline(yintercept = 0.5) + 
  theme_classic()+
  ggtitle('B.') +
  theme(plot.title = element_text(size=10),
        axis.title=element_text(size=10),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "bottom",
        legend.title=element_blank(),
        strip.text.x = element_text(size = 10),
        legend.text=element_text(size=10))

occ <- gridExtra::arrangeGrob(hfi, tmp,nrow=2)
ggsave(file="vbd_hfi_plots/occurence_pdps.pdf", occ, width = 4.5, height = 8, units = "in", dpi=600)

##-----------------------##
##supplemental figures   ##
##-----------------------##

#supplemental unscaled pdp plots
all_data$Type[all_data$Type == "chik"] <- "Chikungunya"
all_data$Type[all_data$Type == "zika"] <- "Zika"
all_data$Type[all_data$Type == "dengue"] <- "Dengue"
all_data$Type[all_data$Type == "malaria"] <- "Malaria"
all_data$Type[all_data$Type == "cutaneous_leish"] <- "Cutaneous Leishmanisis"
all_data$Type[all_data$Type == "visceral_leish"] <- "Visceral Leishmanisis"

colors0 <- cmocean::cmocean("phase")(25)
colors <- c("darkorange","cyan3", "red3", "darkorchid", "#4887D5FF", "gold")

human_footprint_supp <- ggplot(subset(all_data, variables == 'human_footprint'), aes(x=variable_value, y=prob_1, color=Type)) + 
  stat_smooth(aes(group=iter), color='lightgrey', method='loess', size=0.5, se=FALSE) + 
  stat_smooth(aes(), method='loess', size=1, se=FALSE) +
  scale_color_manual(values = colors) + 
  #scale_color_brewer(palette="Dark2") +
  ylab('probability') + xlab("Human Footprint") +
  facet_wrap(~Type, scales='free', ncol=3) +
  theme_bw(base_size = 12) + 
  theme(legend.position = "none")

temp_supp <- ggplot(subset(all_data, variables == 'annual_tmp'), aes(x=variable_value, y=prob_1, color=Type)) + 
  stat_smooth(aes(group=iter), color='lightgrey', method='loess', size=0.5, se=FALSE) + 
  stat_smooth(aes(), method='loess', size=1, se=FALSE) +
  scale_color_manual(values = colors) + 
  #scale_color_brewer(palette="Dark2") +
  ylab('probability') + xlab("Annual Temperature") +
  facet_wrap(~Type, scales='free', ncol=3) +
  theme_bw(base_size = 12) + 
  theme(legend.position = "none")


population_supp <- ggplot(subset(all_data, variables == 'pop_log'), aes(x=variable_value, y=prob_1, color=Type)) + 
  stat_smooth(aes(group=iter), color='lightgrey', method='loess', size=0.5, se=FALSE) + 
  stat_smooth(aes(), method='loess', size=1, se=FALSE) +
  scale_color_manual(values = colors) + 
  #scale_color_brewer(palette="Dark2") +
  ylab('probability') + xlab("log(population size)") +
  facet_wrap(~Type, scales='free', ncol=3) +
  theme_bw(base_size = 12) + 
  theme(legend.position = "none")

forest_supp <- ggplot(subset(all_data, variables == 'forest'), aes(x=variable_value, y=prob_1, color=Type)) + 
  stat_smooth(aes(group=iter), color='lightgrey', method='loess', size=0.5, se=FALSE) + 
  stat_smooth(aes(), method='loess', size=1, se=FALSE) +
  scale_color_manual(values = colors) + 
  #scale_color_brewer(palette="Dark2") +
  ylab('probability') + xlab("% forest cover") +
  facet_wrap(~Type, scales='free', ncol=3) +
  theme_bw(base_size = 12) + 
  theme(legend.position = "none")

ggsave("vbd_hfi_plots/hfi_pdp_supp_figure.png", human_footprint_supp, dpi = 600, units = 'in', width=6.5, height = 5)
ggsave("vbd_hfi_plots/tmp_pdp_supp_figure.png", temp_supp, dpi = 600, units = 'in', width=6.5, height = 5)
ggsave("vbd_hfi_plots/pop_pdp_supp_figure.png", population_supp, dpi = 600, units = 'in', width=6.5, height = 5)
ggsave("vbd_hfi_plots/forest_pdp_supp_figure.png", forest_supp, dpi = 600, units = 'in', width=6.5, height = 5)



####supplement scaled pdp plot
pop <- ggplot(subset(all_data, variables == 'pop_log'), aes(x=variable_value, y=prob_adj, color=Type)) +
  geom_line(aes(group=interaction(iter, Type)), stat="smooth", method="gam", size=0.25, se=FALSE, alpha=0.25) + 
  geom_smooth(aes(),method="gam", size=0.75, se=FALSE, fullrange = FALSE) +
  scale_color_manual(values = colors) + 
  #scale_color_brewer(palette="Dark2") +
  xlab("Log(population)")+
  ylab("Rescaled probabilty of occurrence")+
  coord_cartesian(ylim=c(0,1)) +
  #ylim(0, 1) +
  #geom_hline(yintercept = 0.5) + 
  theme_classic()+
  theme(plot.title = element_text(size=10),
        axis.title=element_text(size=10),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "bottom",
        legend.title=element_blank(),
        strip.text.x = element_text(size = 12),
        legend.text=element_text(size=12))
ggsave(file="vbd_hfi_plots/occurence_pop_supplement_pdps.png", pop, width = 6.5, height = 4, units = "in", dpi=600)





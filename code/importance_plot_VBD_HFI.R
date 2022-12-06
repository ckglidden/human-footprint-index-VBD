##Estimating variable importance for: Human footprint is associated with shifts in assemblages of major vector-borne diseases
##Author: Caroline Glidden & Eloise Skinner
##Date: Dec 6, 2022

##this code uses permutation based methods to estimate variable importance for each model

#load packages
library(ggplot2)
library(randomForestSRC)

######data frame to describe traits
meta <- c()
meta$feature <- c("pop_log", "human_footprint","annual_tmp","annual_pre","annual_wet","forest","pasture","cropland_log"); meta <- as.data.frame(meta)
meta$feature2 <- c("log(pop)", "human footprint","annual tmp","annual pre","annual wet","forest","pasture","log(cropland)")
meta$order <- seq(1, nrow(meta), by = 1)
meta$trait_type <- c(rep("anthropogenic",2),rep("climate",3),rep("land class",3))

colors_phase <- cmocean::cmocean('phase')(20)
colors <- c(colors_phase[3], colors_phase[12], "forestgreen", "gold")

##########################
#####dengue
dengue_model <- readRDS("../saved_models/dengue_model_noVect.rds")
oo <- subsample(dengue_model, verbose = FALSE)
vimpCI <- extract.subsample(oo)$var.jk.sel.Z
vimpCI$feature <- rownames(vimpCI) #use only one vector and rename feature

vimpCI <- merge(vimpCI, meta, by = 'feature')

dengue <- ggplot(vimpCI, aes(x = reorder(feature2, -order), y = mean, color=trait_type)) +
  geom_point(size = 3) + xlab('feature') + ylab('importance') +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = "dodge", width = 0.4, size = 1) +
  scale_color_manual(values=colors) +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(legend.position="none", 
        plot.margin = unit(c(0,0,0,0), "cm"), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),) + 
  ggtitle('a. Dengue')

##########################
#####chik
chik_model <- readRDS("../saved_models/chik_model_noVect.rds")
oo <- subsample(chik_model, verbose = FALSE)
vimpCI <- extract.subsample(oo)$var.jk.sel.Z
vimpCI$feature <- rownames(vimpCI)

vimpCI <- merge(vimpCI, meta, by = 'feature')

chik <- ggplot(vimpCI, aes(x = reorder(feature, -order), y = mean, color=trait_type)) +
  geom_point(size = 3) + xlab('feature') + ylab('importance') +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = "dodge", width = 0.4, size = 1) +
  scale_color_manual(values=colors) +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(legend.position="none", 
        axis.title.y = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm"), 
        axis.title.x = element_blank(),
        axis.text.y = element_blank()) + 
  ggtitle('b.Chikungunya')

##########################
#####zika
zik_model <- readRDS("../saved_models/zika_model_noVect.rds")
oo <- subsample(zik_model, verbose = FALSE)
vimpCI <- extract.subsample(oo)$var.jk.sel.Z
vimpCI$feature <- rownames(vimpCI)

vimpCI <- merge(vimpCI, meta, by = 'feature')

zik <- ggplot(vimpCI, aes(x = reorder(feature, -order), y = mean, color=trait_type)) +
  geom_point(size = 3) + xlab('feature') + ylab('importance') +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = "dodge", width = 0.4, size = 1) +
  scale_color_manual(values=colors) +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(legend.position="none",
        axis.title.y = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm"), 
        axis.title.x = element_blank(),
        axis.text.y = element_blank()) + 
  ggtitle('c. Zika')

##########################
#####malaria
ml_model <- readRDS("../saved_models/ml_noVect_model.rds")
oo <- subsample(ml_model, verbose = FALSE)
vimpCI <- extract.subsample(oo)$var.jk.sel.Z
vimpCI$feature <- rownames(vimpCI)

vimpCI <- merge(vimpCI, meta, by = 'feature')

malaria <- ggplot(vimpCI, aes(x = reorder(feature2, -order), y = mean, color=trait_type)) +
  geom_point(size = 3) + xlab('feature') + ylab('importance') +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = "dodge", width = 0.4, size = 1) +
  scale_color_manual(values=colors) +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(legend.position="none",
        plot.margin = unit(c(0,0,0,0), "cm"), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) + 
  ggtitle('d. Malaria')

##########################
#####cutaneous leish
cl_model <- readRDS("../saved_models/cl_noVect_model.rds")
oo <- subsample(cl_model, verbose = FALSE)
vimpCI <- extract.subsample(oo)$var.jk.sel.Z
vimpCI$feature <- rownames(vimpCI)

vimpCI <- merge(vimpCI, meta, by = 'feature')

cl <- ggplot(vimpCI, aes(x = reorder(feature, -order), y = mean, color=trait_type)) +
  geom_point(size = 3) + xlab('feature') + ylab('importance') +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = "dodge", width = 0.4, size = 1) +
  scale_color_manual(values=colors) +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(legend.position="none", 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"), 
        axis.title.x = element_blank()) + 
  ggtitle('e. Cutaneous leishmaniasis')

##########################
#####visceral leish
vl_model <- readRDS("../saved_models/vl_noVect_model.rds")
oo <- subsample(vl_model, verbose = FALSE)
vimpCI <- extract.subsample(oo)$var.jk.sel.Z
vimpCI$feature <- rownames(vimpCI)

vimpCI <- merge(vimpCI, meta, by = 'feature')

vl <- ggplot(vimpCI, aes(x = reorder(feature, -order), y = mean, color=trait_type)) +
  geom_point(size = 3) + xlab('feature') + ylab('importance') +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = "dodge", width = 0.4, size = 1) +
  scale_color_manual(values=colors) +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(legend.position="none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"), 
        axis.title.x = element_blank()) + 
  ggtitle('f. Visceral leishmaniasis')

#gridExtra::grid.arrange(dengue, chik, zik, malaria, cl, vl, ncol=3)
combined <- cowplot::plot_grid(dengue, chik, zik, malaria, cl, vl, align = 'hv', ncol=3)
cowplot::save_plot('vbd_hfi_plots/hfi_variable_importance.pdf', combined, dpi=600, units = "in", base_width = 9, base_height=6)
# Manuscript: Multi-taxon inventory reveals highly consistent biodiversity responses to ecospace variation
# Author: Ane Kirstine Brunbjerg
# Date: 30-01-2020

#setup
#### load packages ####
library(MASS)
library(DHARMa)
library(dplyr)
library(mgcv)
#### load data - all explanatory variables are standardized  ####
eco_data <- read.csv2(here::here("eco_processing","eco_data.csv"), header=T)


# Models:
# For each model, we made a preliminary screening and selection of relevant variables, 
# only keeping variables with a hypothesized relationship to the species 
# group in question. We further constrained the response direction and shape to
# ecologically plausible responses (Burnham and Anderson 2002, Zuur, et al. 2010) 
# implying an exclusion of negative effects of expansion, continuity and 
# heterogeneity variables on species richness – based on the reasoning that
# more resources, more diverse resources, more environmental variation and 
# increasing temporal and spatial continuity are all hypothesized to have 
# one-sided positive effects on richness. We select variables and constrain
# responses to reduce the risk of including spurious correlations in the 
# models and thereby covering important causal relationships (Supplementary material Table A3).
# Log transformation was preferred if model improvement was indicated by Akaike’s Information
# Criterion (AIC) (Johnson and Omland 2004). The number of explanatory 
# variables were further reduced in order to avoid collinearity (VIF values
# < 3, (Zuur, et al. 2010)). A preliminary set of full models was built 
# using all remaining variables: a general linear mixed poisson model 
# (GLMM) with region as random variable and a GLM with poisson errors 
# using the log link function. We selected the best model type using 
# the ΔAIC < 2 criterium (Burnham and Anderson 2002). Negative binomial
# errors were used if overdispersion was detected (Hilbe 2011) in poisson
# models. We included a quadratic term of the abiotic position variables 
# if the full model significantly improved according to the ΔAIC < 2 
# criterium. Expansion and continuity variables having a negative effect
# in the full model after variable transformation and adding of quadratic
# terms were deleted sequentially starting with the variable with the 
# lowest z-value. The residuals of full models were checked for model 
# misfit, overdispersion and spatial autocorrelation using simulated 
# residuals and R package DHARMa (Hartig 2016). 

# Plants:----
# Model for variable selection using delta AIC < 2 criteria:
MASS::glm.nb(plant_richness ~  eco_pool_log + natural_landscape + 
               shrub_layer_variability + ph_variability_log + soil_fertility_variability_log + soil_moisture_variability_log + ph + smi +I(smi^2)+sfi +I(sfi^2)+
               airtemp +  light + litter_mass + soil_org_C + 
               temporal_continuity + geo_pool_log,
             data = eco_data) -> m_full_nb
summary(m_full_nb)

stepAIC(m_full_nb) # using stepAIC to reduce full model - and the delta AIC<2 criteria for
# manually reducing further (the stepAIC function does not account for squared variables 
# in the model, i.e. we have to make sure the linear effect is kept in the models
# if the squared effect is significant - therefore manual reduction is used here
# 
# Manually from here:
MASS::glm.nb(plant_richness ~ eco_pool_log + shrub_layer_variability + soil_moisture_variability_log+geo_pool_log+
             ph +  smi + I(smi^2) +sfi + I(sfi^2) + 
               litter_mass + soil_org_C + temporal_continuity ,
             data = eco_data) -> m_full_nb
summary(m_full_nb) #geo_pool_log, soil_moisture_variability_log deleted sequentially

# Final model:
MASS::glm.nb(plant_richness ~ eco_pool_log +  shrub_layer_variability+
               ph + smi+ I(smi^2) + sfi + I(sfi^2)  +litter_mass +soil_org_C+ temporal_continuity,
             data = eco_data) -> final_plantmodel
summary(final_plantmodel) 

# Model validation: QQ-plot, residuals vs. fitted, scale-location and residuals vs. leverage plot:
par(mfrow=c(2,2))
plot(final_plantmodel)

# Plotting residuals vs. variables in the model:
final_plantmodel %>% 
  simulateResiduals(n = 10000) -> sim_res_nb

plotResiduals(eco_data$eco_pool_log, sim_res_nb$scaledResiduals,
              xlab = "Eco_pool_log", ylab = "DHARMa simulated residuals")

plotResiduals(eco_data$shrub_layer_variability, sim_res_nb$scaledResiduals,
              xlab = "Shrub layer variability", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$ph, sim_res_nb$scaledResiduals,
              xlab = "pH", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$smi, sim_res_nb$scaledResiduals,
              xlab = "Soil moisture index", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$sfi, sim_res_nb$scaledResiduals,
              xlab = "Soil fertility index", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$temporal_continuity, sim_res_nb$scaledResiduals,
              xlab = "Temporal continuity", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$soil_org_C, sim_res_nb$scaledResiduals,
              xlab = "Soil organic C", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$litter_mass, sim_res_nb$scaledResiduals,
              xlab = "Litter mass", ylab = "DHARMa simulated residuals")


# Model performance:
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
theta=summary(final_plantmodel)$theta 
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(plant_richness ~ eco_pool_log + shrub_layer_variability + 
              ph + I(smi^2) + I(sfi^2) + sfi + smi  + temporal_continuity+soil_org_C+litter_mass,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1, type="response")
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
R2cv=1-(sum((eco_data$plant_richness-PTOT)^2)/(length(PTOT))/(var(eco_data$plant_richness)))
cor(PTOT, eco_data$plant_richness)^2  

# save predicted plants for later
pred_plants <- PTOT

# Mosses:----
# Model for variable selection using delta AIC < 2 criteria:
MASS::glm.nb(moss_richness ~ natural_landscape+ph_variability_log+soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+sfi+I(sfi^2)+ smi+light+boulderPA + litter_mass_log+soil_org_C_log+ soil_org_matter_log+dead_wood_debris_log+shrub_layer+ spatial_continuity+temporal_continuity,
             data = eco_data) -> m_full_nb
summary(m_full_nb)

stepAIC(m_full_nb) # using stepAIC to reduce full model - and the delta AIC<2 criteria for
# manually reducing further (the stepAIC function does not account for squared variables 
# in the model, i.e. we have to make sure the linear effect is kept in the models
# if the squared effect is significant - therefore manual reduction is used here
# 
# Manually from here:
MASS::glm.nb(moss_richness ~ natural_landscape + ph +sfi + I(sfi^2) + 
               smi + boulderPA + litter_mass_log + shrub_layer,
             data = eco_data) -> m_full_nb
summary(m_full_nb) #boulderPA, ph deleted sequentially

# Final model:
MASS::glm.nb(moss_richness ~ natural_landscape    +sfi +  I(sfi^2) + 
               smi  +litter_mass_log  + shrub_layer,
             data = eco_data) -> final_moss_model
summary(final_moss_model) 

# Model validation: QQ-plot, residuals vs. fitted, scale-location and residuals vs. leverage plot:
par(mfrow=c(2,2))
plot(final_moss_model)

# Plotting residuals vs. variables in the model:
final_moss_model %>% 
  simulateResiduals(n = 10000) -> sim_res_nb

plotResiduals(eco_data$natural_landscape, sim_res_nb$scaledResiduals,
              xlab = "Natural landscapes", ylab = "DHARMa simulated residuals")

plotResiduals(eco_data$smi, sim_res_nb$scaledResiduals,
              xlab = "Soil moisture index", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$sfi, sim_res_nb$scaledResiduals,
              xlab = "Soil fertility index", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$litter_mass_log, sim_res_nb$scaledResiduals,
              xlab = "Litter mass", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$shrub_layer, sim_res_nb$scaledResiduals,
              xlab = "Shrub layer", ylab = "DHARMa simulated residuals")


# Model performance:
theta=summary(final_moss_model)$theta
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(moss_richness ~ natural_landscape  + 
              I(sfi^2) + sfi + smi  + litter_mass_log + 
              shrub_layer,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1, type="response")
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
R2cv=1-(sum((eco_data$moss_richness-PTOT)^2)/(length(PTOT))/(var(eco_data$moss_richness)))
cor(PTOT, eco_data$moss_richness)^2 

# save predicted moss richness for later
pred_moss <- PTOT

# Lichens----
# Model for variable selection using delta AIC < 2 criteria:
MASS::glm.nb(lichen_richness ~   sfi + smi+I(smi^2) + boulderPA + litter_mass_log  + 
               shrub_layer  + temporal_continuity,
             data = eco_data) -> m_full_nb

stepAIC(m_full_nb) # using stepAIC to reduce full model - and the delta AIC<2 criteria for
# manually reducing further (the stepAIC function does not account for squared variables 
# in the model, i.e. we have to make sure the linear effect is kept in the models
# if the squared effect is significant - therefore manual reduction is used here
# 
# Manually from here:
MASS::glm.nb(lichen_richness ~   sfi + smi+I(smi^2) + boulderPA + litter_mass_log  + 
               shrub_layer  + temporal_continuity,
             data = eco_data) -> m_full_nb
summary(m_full_nb)# litter_mass_log, temporal_continuity deleted

# Final model:
MASS::glm.nb(lichen_richness ~   sfi + smi+I(smi^2) + boulderPA  + 
               shrub_layer  ,
             data = eco_data) -> final_lichenmodel
summary(final_lichenmodel)

# Model validation: QQ-plot, residuals vs. fitted, scale-location and residuals vs. leverage plot:
par(mfrow=c(2,2))
plot(final_lichenmodel)

# Plotting residuals vs. variables in the model:
final_lichenmodel %>% 
  simulateResiduals(n = 10000) -> sim_res_nb

plotResiduals(eco_data$smi, sim_res_nb$scaledResiduals,
              xlab = "Soil moisture index", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$sfi, sim_res_nb$scaledResiduals,
              xlab = "Soil fertility index", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$boulderPA, sim_res_nb$scaledResiduals,
              xlab = "Boulders", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$shrub_layer, sim_res_nb$scaledResiduals,
              xlab = "Shrub layer", ylab = "DHARMa simulated residuals")

# Resiuals indicate non-linear patterns for sfi^2. Cheking for significant gam on residuals:
resid(final_lichenmodel) -> eco_data$sim_res_nbtest
T1 <- gam(sim_res_nbtest ~ s(sfi), data = eco_data)
summary(T1)
# GAM not significant - we do not add sfi^2 to the model

# Model performance:
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
theta=summary(final_lichenmodel)$theta 
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(lichen_richness ~ sfi + smi+I(smi^2) + boulderPA  + 
              shrub_layer ,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1, type="response")
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
R2cv=1-(sum((eco_data$lichen_richness-PTOT)^2)/(length(PTOT))/(var(eco_data$lichen_richness)))
cor(PTOT, eco_data$lichen_richness)^2  

# save predicted lichens for later
pred_lichen <- PTOT

# Producers----
# Model for variable selection using delta AIC < 2 criteria:
MASS::glm.nb(producer_richness ~eco_pool_log + natural_landscape + ph_variability_log +  
               soil_moisture_variability_log + ph + airtemp +smi + I(smi^2) +sfi +  I(sfi^2) +  
                light + boulderPA + litter_mass_log + shrub_layer + dead_wood_debris_log +  
               soil_org_C_log +  
               temporal_continuity + geo_pool_log,
             data = eco_data) -> m_full_nb
summary(m_full_nb)

stepAIC(m_full_nb) # using stepAIC to reduce full model - and the delta AIC<2 criteria for
# manually reducing further (the stepAIC function does not account for squared variables 
# in the model, i.e. we have to make sure the linear effect is kept in the models
# if the squared effect is significant - therefore manual reduction is used here
# 
# Manually from here:
MASS::glm.nb(producer_richness ~ eco_pool_log + natural_landscape + ph_variability_log + 
               soil_moisture_variability_log + ph +  smi + I(smi^2) + sfi + I(sfi^2) +
               light + boulderPA + litter_mass_log + shrub_layer + dead_wood_debris_log + 
               soil_org_C_log + temporal_continuity + geo_pool_log,
             data = eco_data) -> m_full_nb #litter_mass_log, dead_wood_debris_log,  ph_variability_log, light,
#geo_pool_log,soil_org_C_log, natural_landscape,  soil_moisture_variability_log deleted sequentially

# Final model:
MASS::glm.nb(producer_richness ~ eco_pool_log  + ph + sfi + 
               smi + I(smi^2) + I(sfi^2) +  boulderPA  + shrub_layer + 
               temporal_continuity ,
             data = eco_data) -> final_producermodel

# Model validation: QQ-plot, residuals vs. fitted, scale-location and residuals vs. leverage plot:
par(mfrow=c(2,2))
plot(final_producermodel)

# Plotting residuals vs. variables in the model:
final_producermodel %>% 
  simulateResiduals(n = 10000) -> sim_res_nb

plotResiduals(eco_data$eco_pool_log, sim_res_nb$scaledResiduals,
              xlab = "Eco_pool_log", ylab = "DHARMa simulated residuals")

plotResiduals(eco_data$ph, sim_res_nb$scaledResiduals,
              xlab = "pH", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$smi, sim_res_nb$scaledResiduals,
              xlab = "Soil moisture index", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$sfi, sim_res_nb$scaledResiduals,
              xlab = "Soil fertility index", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$temporal_continuity, sim_res_nb$scaledResiduals,
              xlab = "Temporal continuity", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$boulderPA, sim_res_nb$scaledResiduals,
              xlab = "Boulders", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$shrub_layer, sim_res_nb$scaledResiduals,
              xlab = "Shrub layer", ylab = "DHARMa simulated residuals")


# Model performance:
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
theta=summary(final_producermodel)$theta 
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(producer_richness ~ eco_pool_log  + ph + sfi + 
              smi + I(smi^2) + I(sfi^2) +  boulderPA  + shrub_layer + 
              temporal_continuity,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1, type="response")
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
R2cv=1-(sum((eco_data$producer_richness-PTOT)^2)/(length(PTOT))/(var(eco_data$producer_richness)))
cor(PTOT, eco_data$producer_richness)^2  

# save predicted producers for later:
pred_producers <- PTOT

# Check correlations between predicted producers and the sum of predicted plants, mosses and lichens:
pred_summed <- pred_plants+pred_moss+pred_lichen
cor(pred_summed, eco_data$producer_richness)^2  

# Symbionts----
# Model for variable selection using delta AIC < 2 criteria:
MASS::glm.nb(symbiont_richness ~natural_landscape+ph_variability_log+soil_fertility_variability_log+soil_moisture_variability_log+
               ph+airtemp+ sfi+I(sfi^2)+smi+light +plant_richness_std_log+ soil_org_C_log+
               fungi_symbiont_plants_log+dungPA+soil_org_matter+ large_tree_density+shrub_layer+temporal_continuity   ,
             data = eco_data) -> m_full_nb
summary(m_full_nb)

stepAIC(m_full_nb) # using stepAIC to reduce full model - and the delta AIC<2 criteria for
# manually reducing further (the stepAIC function does not account for squared variables 
# in the model, i.e. we have to make sure the linear effect is kept in the models
# if the squared effect is significant - therefore manual reduction is used here
# 

# Manually from here:
MASS::glm.nb(symbiont_richness ~  natural_landscape + ph_variability_log + soil_moisture_variability_log + 
               airtemp + sfi + I(sfi^2) + light + plant_richness_std_log + fungi_symbiont_plants_log + 
               dungPA + shrub_layer + temporal_continuity,
             data = eco_data) -> m_full_nb
summary(m_full_nb) #tempcon,  soil_moisture_variability_log, airtemp, plant_richness_std_log, natural_landscape deleted sequentially

# Final model:
MASS::glm.nb(symbiont_richness ~   ph_variability_log + 
               sfi + I(sfi^2) + light + fungi_symbiont_plants_log + 
               dungPA + shrub_layer ,
             data = eco_data) -> final_symbiontmodel
summary(final_symbiontmodel) #


# Model validation: QQ-plot, residuals vs. fitted, scale-location and residuals vs. leverage plot:
par(mfrow=c(2,2))
plot(final_symbiontmodel)

# Plotting residuals vs. variables in the model:
final_symbiontmodel %>% 
  simulateResiduals(n = 10000) -> sim_res_nb

plotResiduals(eco_data$ph_variability_log, sim_res_nb$scaledResiduals,
              xlab = "Soil pH variability", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$sfi, sim_res_nb$scaledResiduals,
              xlab = "Soil fertility index", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$light, sim_res_nb$scaledResiduals,
              xlab = "Light intensity", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$fungi_symbiont_plants_log, sim_res_nb$scaledResiduals,
              xlab = "Fungi symbiont plants", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$dungPA, sim_res_nb$scaledResiduals,
              xlab = "Dung", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$shrub_layer, sim_res_nb$scaledResiduals,
              xlab = "Shrub layer", ylab = "DHARMa simulated residuals")

# Model performance:
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
theta=summary(final_symbiontmodel)$theta 
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(symbiont_richness ~   ph_variability_log + 
              sfi + I(sfi^2) + light + fungi_symbiont_plants_log + 
              dungPA + shrub_layer ,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1, type="response")
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
R2cv=1-(sum((eco_data$symbiont_richness-PTOT)^2)/(length(PTOT))/(var(eco_data$symbiont_richness)))
cor(PTOT, eco_data$symbiont_richness)^2  

# save predicted symbionts for later
pred_symbionts <- PTOT

# Decomposers----
# Model for variable selection using delta AIC < 2 criteria:
MASS::glm.nb(decomposer_richness ~natural_landscape+ph_variability_log+ph+airtemp+ sfi+smi+light +plant_richness_std_log+ soil_org_C_log+dungPA+litter_mass_log+ dead_wood_vol+shrub_layer,
             data = eco_data) -> m_full_nb
summary(m_full_nb)

stepAIC(m_full_nb) # using stepAIC to reduce full model - and the delta AIC<2 criteria for
# manually reducing further (the stepAIC function does not account for squared variables 
# in the model, i.e. we have to make sure the linear effect is kept in the models
# if the squared effect is significant - therefore manual reduction is used here
# 

# Manually from here:
MASS::glm.nb(decomposer_richness ~ smi + light + plant_richness_std_log + 
               soil_org_C_log + dungPA + litter_mass_log + dead_wood_vol + 
               shrub_layer, data = eco_data) -> m_full_nb
summary(m_full_nb) # delta 2 criteria:soil_org_C_log,smi deleted sequentially

# Final model:
MASS::glm.nb(decomposer_richness ~ light + plant_richness_std_log 
             + dungPA + litter_mass_log + dead_wood_vol + 
               shrub_layer, data = eco_data) -> final_decomposermodel
summary(final_decomposermodel)


# Model validation: QQ-plot, residuals vs. fitted, scale-location and residuals vs. leverage plot:
par(mfrow=c(2,2))
plot(final_decomposermodel)

# Plotting residuals vs. variables in the model:
final_decomposermodel %>% 
  simulateResiduals(n = 10000) -> sim_res_nb

plotResiduals(eco_data$light, sim_res_nb$scaledResiduals,
              xlab = "Light intensity", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$plant_richness_std_log, sim_res_nb$scaledResiduals,
              xlab = "Plant richness", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$litter_mass_log, sim_res_nb$scaledResiduals,
              xlab = "Litter mass", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$dungPA, sim_res_nb$scaledResiduals,
              xlab = "Dung", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$dead_wood_vol, sim_res_nb$scaledResiduals,
              xlab = "Dead wood volume", ylab = "DHARMa simulated residuals")


# Model performance:
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
theta=summary(final_decomposermodel)$theta 
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(decomposer_richness ~ light + plant_richness_std_log + 
              dungPA + litter_mass_log + dead_wood_vol + 
              shrub_layer ,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1, type="response")
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
R2cv=1-(sum((eco_data$decomposer_richness-PTOT)^2)/(length(PTOT))/(var(eco_data$decomposer_richness)))
cor(PTOT, eco_data$decomposer_richness)^2  

# save predicted decomposers for later
pred_decomposer <- PTOT

# Detritivores----
# Model for variable selection using delta AIC < 2 criteria:
glm(detritivore_richness ~  natural_landscape+ph_variability_log+soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+ sfi+smi+light 
    +plant_richness_std_log+ soil_org_C+litter_mass_log+ dead_wood_vol_log+large_tree_density_log+shrub_layer     ,
    family = "poisson",
    data = eco_data) -> m_full_poi
summary(m_full_poi)

stepAIC(m_full_poi) # using stepAIC to reduce full model - and the delta AIC<2 criteria for
# manually reducing further (the stepAIC function does not account for squared variables 
# in the model, i.e. we have to make sure the linear effect is kept in the models
# if the squared effect is significant - therefore manual reduction is used here
# 

# Manually from here:
glm(detritivore_richness ~  soil_moisture_variability_log + ph + sfi + smi + 
      light + plant_richness_std_log + soil_org_C + litter_mass_log + 
      large_tree_density_log + shrub_layer  ,
    family = "poisson",
    data = eco_data) -> m_full_poi
summary(m_full_poi)#delta 2 criteria: litterlog, sfi deleted sequentially

# Final model:#
glm(detritivore_richness ~  soil_moisture_variability_log + ph +  smi + 
      light + plant_richness_std_log + soil_org_C  + 
      large_tree_density_log + shrub_layer  ,
    family = "poisson",
    data = eco_data) -> final_detritivoremodel  


# Model validation: QQ-plot, residuals vs. fitted, scale-location and residuals vs. leverage plot:
par(mfrow=c(2,2))
plot(final_detritivoremodel)

# Plotting residuals vs. variables in the model:
final_detritivoremodel %>% 
  simulateResiduals(n = 10000) -> sim_res_nb

plotResiduals(eco_data$soil_moisture_variability_log, sim_res_nb$scaledResiduals,
              xlab = "Soil moisture variability", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$ph, sim_res_nb$scaledResiduals,
              xlab = "Soil ph", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$smi, sim_res_nb$scaledResiduals,
              xlab = "Soil moisture index", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$light, sim_res_nb$scaledResiduals,
              xlab = "Light intensity", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$plant_richness_std_log, sim_res_nb$scaledResiduals,
              xlab = "Plant richness", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$soil_org_C, sim_res_nb$scaledResiduals,
              xlab = "Soil organic C", ylab = "DHARMa simulated residuals")

plotResiduals(eco_data$large_tree_density_log, sim_res_nb$scaledResiduals,
              xlab = "Large tree density", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$shrub_layer, sim_res_nb$scaledResiduals,
              xlab = "Shrub layer", ylab = "DHARMa simulated residuals")

# Model performance:
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
theta=summary(final_detritivoremodel)$theta
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(detritivore_richness ~  soil_moisture_variability_log + ph +  smi + 
              light + plant_richness_std_log + soil_org_C  + large_tree_density_log + 
              shrub_layer  ,
            family = "poisson", data = DataCV)
  P1=predict(M1, DataC1, type="response")
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
R2cv=1-(sum((eco_data$detritivore_richness-PTOT)^2)/(length(PTOT))/(var(eco_data$detritivore_richness)))
cor(PTOT, eco_data$detritivore_richness)^2 

# Save predicted detritivores for later:
pred_detriti<- PTOT

# Flying insects----
# Model for variable selection using delta AIC < 2 criteria:
MASS::glm.nb(flying_insect_richness ~  natural_landscape+soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+ sfi+smi+light + 
               soil_org_C+plant_richness_std_log+flower_abundance_log+ dead_wood_debris+shrub_layer,
             data = eco_data) -> m_full_nb
summary(m_full_nb)

stepAIC(m_full_nb) # using stepAIC to reduce full model - and the delta AIC<2 criteria for
# manually reducing further (the stepAIC function does not account for squared variables 
# in the model, i.e. we have to make sure the linear effect is kept in the models
# if the squared effect is significant - therefore manual reduction is used here

# Manually from here:  
MASS::glm.nb(flying_insect_richness ~ soil_fertility_variability_log + airtemp + sfi + 
               smi + light + plant_richness_std_log ,
             data = eco_data) -> m_full_nb
summary(m_full_nb) # sfi deleted


# Final model:
MASS::glm.nb(flying_insect_richness ~ soil_fertility_variability_log + airtemp  + 
               smi + light + plant_richness_std_log ,
             data = eco_data) -> final_flyersmodel
summary(final_flyersmodel) #

# Model validation: QQ-plot, residuals vs. fitted, scale-location and residuals vs. leverage plot:
par(mfrow=c(2,2))
plot(final_flyersmodel) #looks fine!

# Testing for spatial autocorrelation:
final_flyersmodel %>% 
  simulateResiduals(n = 10000) -> sim_res_nb
plot(sim_res_nb)
testSpatialAutocorrelation(sim_res_nb, x=eco_data$utm_x,y=eco_data$utm_y) #  significant , p<0.05

# We use backwards selection based on five-fold cross-validation to reduce to final model instead:
MASS::glm.nb(flying_insect_richness ~  natural_landscape+soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+ sfi+smi+light + 
               soil_org_C+plant_richness_std_log+flower_abundance_log+ dead_wood_debris+shrub_layer,
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 953.55

theta=summary(m_full_nb)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(flying_insect_richness ~ natural_landscape+soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+sfi+ smi+light+soil_org_C + 
              plant_richness_std_log+flower_abundance_log+dead_wood_debris+ shrub_layer,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$flying_insect_richness)^2 

# Run the above code for all models - leaving out one variable at a time and calculate cross-validated R2 values
# full CV R2: 0.01676973 
# - shrub_layer: 0.04428795
# - dead_wood_debris:  0.09653129   biggest increase in CV R2 - take out, 
# - flower_abundance_log:0.03242412
# - plant_richness_std_log : 0.004365274
# - soil_org_C: 0.02007331
# - light: 0.006869675
# - smi: 0.04701313
# - sfi: 0.01608653
# - airtemp: 0.008810747
# - ph: 0.0141552
# - soil_moisture_variability_log : 0.01643176
# - soil_fertility_variability_log:0.01769044
# - natural_landscape: 0.03069723

#backwards selection 2:
MASS::glm.nb(flying_insect_richness ~  natural_landscape+soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+ sfi+smi+light + 
               soil_org_C+plant_richness_std_log+flower_abundance_log+shrub_layer,
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 951.55

theta=summary(m_full_nb)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(flying_insect_richness ~  natural_landscape+soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+ sfi+smi+light + 
              soil_org_C+plant_richness_std_log+flower_abundance_log+shrub_layer,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$flying_insect_richness)^2 
# Run the above code for all models - leaving out one variable at a time
# full CV R2:0.09653149  # AIC=951.55
# - shrub_layer: 0.116453
# - flower_abundance_log: 0.1048288
# - soil_org_C: 0.1094681
# - plant_richness_std_log : 0.07512146
# - light: 0.09107267
# - smi: 0.1268702
# - sfi: 0.09754931
# - airtemp: 0.09681744
# - ph:0.09358857
# - soil_moisture_variability_log : 0.09662606
# - soil_fertility_variability_log: 0.1073907
# - natural_landscape:0.1387325  biggest increase in CV R2 - take out, 


#backwards selection 3:
MASS::glm.nb(flying_insect_richness ~  soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+ sfi+smi+light + 
               soil_org_C+plant_richness_std_log+flower_abundance_log+shrub_layer,
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 950.2

theta=summary(m_full_nb)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(flying_insect_richness ~  soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+ sfi+smi+light + 
              soil_org_C+plant_richness_std_log+flower_abundance_log+shrub_layer,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$flying_insect_richness)^2 
# Run the above code for all models - leaving out one variable at a time
# full CV R2: 0.1387472 
# - shrub_layer: 0.155209
# - flower_abundance_log: 0.1488059
# - soil_org_C: 0.1607174
# - plants_log : 0.09954409
# - light: 0.1302147
# - smi: 0.1737952 biggest increase in CV R2 - take out, 
# - sfi: 0.1262469
# - airtemp: 0.1385553
# - ph: 0.1427723
# - soil_moisture_variability_log : 0.1398875
# - soil_fertility_variability_log: 0.138177

#backwards selection 4:
MASS::glm.nb(flying_insect_richness ~  soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+ sfi+light + 
               soil_org_C+plant_richness_std_log+flower_abundance_log+shrub_layer,
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 951.19

theta=summary(m_full_nb)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(flying_insect_richness ~  soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+ sfi+light + 
              soil_org_C+plant_richness_std_log+flower_abundance_log+shrub_layer,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$flying_insect_richness)^2 
# Run the above code for all models - leaving out one variable at a time
# full CV R2: 0.1738404  
# - shrub_layer:0.1820299
# - flower_abundance_log: 0.1867288 biggest increase in CV R2 - take out, 
# - soil_org_C: 0.1322052
# - plant_richness_std_log : 0.1110702
# - light: 0.1539386
# - sfi: 0.1634398
# - airtemp: 0.1636072
# - ph: 0.1829944 
# - soilm_sd_may : 0.1674087
# - soil_fertility_variability_log: 0.1760462

#backwards selection 5:
MASS::glm.nb(flying_insect_richness ~  soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+ sfi+light + 
               soil_org_C+plant_richness_std_log+shrub_layer,
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 950.26

theta=summary(m_full_nb)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(flying_insect_richness ~   soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+ sfi+light + 
              soil_org_C+plant_richness_std_log+shrub_layer,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$flying_insect_richness)^2
# Run the above code for all models - leaving out one variable at a time
# full CV R2: 0.1867487  
# - shrub_layer:0.1962114
# - soil_org_C: 0.151106
# - plant_richness_std_log :0.1051168 
# - light: 0.1605509
# - sfi: 0.1815061
# - airtemp: 0.1840217
# - ph:  0.1975754  biggest increase in CV R2 - take out, 
# - soil_moisture_variability_log : 0.1810386
# - soil_fertility_variability_log: 0.1829705

#backwards selection 6:
MASS::glm.nb(flying_insect_richness ~  soil_fertility_variability_log+soil_moisture_variability_log+airtemp+ sfi+light + 
               soil_org_C+plant_richness_std_log+shrub_layer,
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 949.78

theta=summary(m_full_nb)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(flying_insect_richness ~   soil_fertility_variability_log+soil_moisture_variability_log+airtemp+sfi+ light+soil_org_C + 
              plant_richness_std_log+shrub_layer,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$flying_insect_richness)^2
# Run the above code for all models - leaving out one variable at a time
# full CV R2: 0.1976169  
# - shrub_layer:0.2085814  biggest increase in CV R2 - take out, 
# - soil_org_C: 0.1713677
# - plant_richness_std_log : 0.07339184
# - light: 0.1688393
# - sfi: 0.1983921
# - airtemp: 0.197681
# - soil_moisture_variability_log : 0.1905566
# - soil_fertility_variability_log: 0.1880746

#backwards selection 7:
MASS::glm.nb(flying_insect_richness ~  soil_fertility_variability_log+soil_moisture_variability_log+airtemp+ sfi+light + 
               soil_org_C+plant_richness_std_log,
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 947.9

theta=summary(m_full_nb)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(flying_insect_richness ~   soil_fertility_variability_log+soil_moisture_variability_log+airtemp+ sfi+light + 
              soil_org_C+plant_richness_std_log,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$flying_insect_richness)^2
# Run the above code for all models - leaving out one variable at a time
# full CV R2: 0.2085868  
# - soil_org_C: 0.1753869
# - plant_richness_std_log : 0.08927768
# - light: 0.1813645
# - sfi: 0.2087571 biggest increase in CV R2 - take out,
# - airtemp: 0.2082473
# - soilm_sd_may : 0.1986692
# - soil_fertility_variability_log: 0.1972336

#backwards selection 8:
MASS::glm.nb(flying_insect_richness ~  soil_fertility_variability_log+soil_moisture_variability_log+airtemp+ light + 
               soil_org_C+plant_richness_std_log,
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 947.52

theta=summary(m_full_nb)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(flying_insect_richness ~   soil_fertility_variability_log+soil_moisture_variability_log+airtemp+ light + 
              soil_org_C+plant_richness_std_log,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$flying_insect_richness)^2
# Run the above code for all models - leaving out one variable at a time
# full CV R2: 0.2087952  
# - soil_org_C: 0.1643846
# - plant_richness_std_log : 0.09038673
# - light: 0.1551485
# - airtemp: 0.2112657 biggest increase in CV R2 - take out,
# - soil_moisture_variability_log : 0.199245
# - soil_fertility_variability_log: 0.2068985


#backwards selection 9:
MASS::glm.nb(flying_insect_richness ~  soil_fertility_variability_log+soil_moisture_variability_log+ light + 
               soil_org_C+plant_richness_std_log,
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 951.25 -  deltaAIC>2  from previous model - we stick to the above as final model


MASS::glm.nb(flying_insect_richness ~  soil_fertility_variability_log+soil_moisture_variability_log+ airtemp+light + 
               soil_org_C+plant_richness_std_log,
             data = eco_data) -> final_flyersmodel  


# Model validation: QQ-plot, residuals vs. fitted, scale-location and residuals vs. leverage plot:
par(mfrow=c(2,2))
plot(final_flyersmodel)

# Plotting residuals vs. variables in the model:
final_flyersmodel %>% 
  simulateResiduals(n = 10000) -> sim_res_nb

plotResiduals(eco_data$soil_fertility_variability_log, sim_res_nb$scaledResiduals,
              xlab = "Soil fertility variability", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$soil_moisture_variability_log, sim_res_nb$scaledResiduals,
              xlab = "Soil moisture variability", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$airtemp, sim_res_nb$scaledResiduals,
              xlab = "Air temperature", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$light, sim_res_nb$scaledResiduals,
              xlab = "Light intensity", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$plant_richness_std_log, sim_res_nb$scaledResiduals,
              xlab = "Plant richness", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$soil_org_C, sim_res_nb$scaledResiduals,
              xlab = "Soil organic C", ylab = "DHARMa simulated residuals")

# Model performance:
theta=summary(final_flyersmodel)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(flying_insect_richness ~   soil_fertility_variability_log+soil_moisture_variability_log+airtemp+ light + 
              soil_org_C+plant_richness_std_log,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$flying_insect_richness)^2

# Save predicted flying insects for later:
pred_flyers<- PTOT

# Predatory arthropods----
# Model for variable selection using delta AIC < 2 criteria:
MASS::glm.nb(predator_richness ~natural_landscape+soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+ sfi+smi+light + I(ph^2)
             +soil_org_C_log+insect_host_plants+dungPA+ shrub_layer,
             data = eco_data) -> m_full_nb
summary(m_full_nb) 

stepAIC(m_full_nb)
# using stepAIC to reduce full model - and the delta AIC<2 criteria for
# manually reducing further (the stepAIC function does not account for squared variables 
# in the model, i.e. we have to make sure the linear effect is kept in the models
# if the squared effect is significant - therefore manual reduction is used here
# 
# Manually from here: 
MASS::glm.nb(predator_richness ~ natural_landscape + soil_fertility_variability_log + 
               soil_moisture_variability_log + ph +I(ph^2) +  airtemp + sfi + smi + light + 
               insect_host_plants + shrub_layer,
             data = eco_data) -> m_full_nb# delta 2 criteria:sfi,soil_fertility_variability_log deleted sequentially

#final model:
MASS::glm.nb(predator_richness ~ natural_landscape  + 
               soil_moisture_variability_log + ph + I(ph^2) +airtemp  + smi + light + insect_host_plants + 
                shrub_layer,
             data = eco_data) -> final_predatormodel
summary(final_predatormodel) #


# Model validation: QQ-plot, residuals vs. fitted, scale-location and residuals vs. leverage plot:
par(mfrow=c(2,2))
plot(final_predatormodel)

# Plotting residuals vs. variables in the model:
final_predatormodel %>% 
  simulateResiduals(n = 10000) -> sim_res_nb

plotResiduals(eco_data$natural_landscape, sim_res_nb$scaledResiduals,
              xlab = "Natural landscape", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$soil_moisture_variability_log, sim_res_nb$scaledResiduals,
              xlab = "Soil moisture variability", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$ph, sim_res_nb$scaledResiduals,
              xlab = "Soil ph", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$smi, sim_res_nb$scaledResiduals,
              xlab = "Soil moisture index", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$airtemp, sim_res_nb$scaledResiduals,
              xlab = "Air temperature", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$light, sim_res_nb$scaledResiduals,
              xlab = "Light intensity", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$insect_host_plants, sim_res_nb$scaledResiduals,
              xlab = "Insect host plant avaliability", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$shrub_layer, sim_res_nb$scaledResiduals,
              xlab = "Shrub layer", ylab = "DHARMa simulated residuals")

# Model performance:
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
theta=summary(final_predatormodel)$theta
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(predator_richness ~ natural_landscape  + 
              soil_moisture_variability_log + ph + airtemp  + smi + light + insect_host_plants + 
              I(ph^2) + shrub_layer ,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1, type="response")
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
R2cv=1-(sum((eco_data$predator_richness-PTOT)^2)/(length(PTOT))/(var(eco_data$predator_richness)))
cor(PTOT, eco_data$predator_richness)^2 

# Save predicted predators for later:
pred_pred<- PTOT

# Herbivores----
# Model for variable selection using delta AIC < 2 criteria:
MASS::glm.nb(herbivore_richness ~natural_landscape+ph_variability_log+soil_fertility_variability_log+ph+airtemp+ sfi+smi+light + plant_richness_std_log+
               soil_org_C+flower_abundance_log+ dead_wood_vol_log+shrub_layer ,
             data = eco_data) -> m_full_nb
summary(m_full_nb) 

stepAIC(m_full_nb)
# using stepAIC to reduce full model - and the delta AIC<2 criteria for
# manually reducing further (the stepAIC function does not account for squared variables 
# in the model, i.e. we have to make sure the linear effect is kept in the models
# if the squared effect is significant - therefore manual reduction is used here

# Manually from here: 
MASS::glm.nb(herbivore_richness ~  natural_landscape +  airtemp +ph+
               plant_richness_std_log  + shrub_layer +ph + flower_abundance_log,
             data = eco_data) -> m_full_nb # delta AIC criteria: ph,flower_abundance_log deleted sequentially
summary(m_full_nb) #

#final model:
MASS::glm.nb(herbivore_richness ~  natural_landscape +  airtemp +
               plant_richness_std_log  + shrub_layer ,
             data = eco_data)  -> final_herbivoremodel
summary(final_herbivoremodel) #


# Model validation: QQ-plot, residuals vs. fitted, scale-location and residuals vs. leverage plot:
par(mfrow=c(2,2))
plot(final_herbivoremodel)

# Plotting residuals vs. variables in the model:
final_herbivoremodel %>% 
  simulateResiduals(n = 10000) -> sim_res_nb

plotResiduals(eco_data$natural_landscape, sim_res_nb$scaledResiduals,
              xlab = "Natural landscape", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$airtemp, sim_res_nb$scaledResiduals,
              xlab = "Air temperature", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$plant_richness_std_log, sim_res_nb$scaledResiduals,
              xlab = "Plant richness", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$shrub_layer, sim_res_nb$scaledResiduals,
              xlab = "Shrub layer", ylab = "DHARMa simulated residuals")

# Model performance:
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
theta=summary(final_herbivoremodel)$theta
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(herbivore_richness ~ natural_landscape +  airtemp +
              plant_richness_std_log  + shrub_layer ,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1, type="response")
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
R2cv=1-(sum((eco_data$herbivore_richness-PTOT)^2)/(length(PTOT))/(var(eco_data$herbivore_richness)))
cor(PTOT, eco_data$herbivore_richness)^2 

# Saving predicted herbivores for later:
pred_herb<- PTOT

# Consumers----
# Model for variable selection using delta AIC < 2 criteria:
MASS::glm.nb(consumer_richness ~  natural_landscape+ph_variability_log+soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+
               sfi+smi+light+litter_mass_log+flower_abundance_log+shrub_layer+dead_wood_vol_log+large_tree_density+
               soil_org_C+ dungPA+insect_host_plants_log,
             data = eco_data) -> m_full_nb
summary(m_full_nb) 

stepAIC(m_full_nb)
# using stepAIC to reduce full model - and the delta AIC<2 criteria for
# manually reducing further (the stepAIC function does not account for squared variables 
# in the model, i.e. we have to make sure the linear effect is kept in the models
# if the squared effect is significant - therefore manual reduction is used here
# 
# Manually from here: 
MASS::glm.nb(consumer_richness ~ ph_variability_log + soil_fertility_variability_log + 
               airtemp + smi + light + shrub_layer + dead_wood_vol_log + soil_org_C + 
               insect_host_plants_log,
             data = eco_data) -> m_full_nb #delta AIC criteria: soil_fertility_variability_log, soil_org_C,
# smi, dead_wood_vol_log,  ph_variability_log deleted sequentially

#final model:
MASS::glm.nb(consumer_richness ~  airtemp +  light + shrub_layer + 
               insect_host_plants_log,
             data = eco_data) -> final_consumermodel
summary(final_consumermodel) #


# Model validation: QQ-plot, residuals vs. fitted, scale-location and residuals vs. leverage plot:
par(mfrow=c(2,2))
plot(final_consumermodel)

# Plotting residuals vs. variables in the model:
final_consumermodel %>% 
  simulateResiduals(n = 10000) -> sim_res_nb

plotResiduals(eco_data$airtemp, sim_res_nb$scaledResiduals,
              xlab = "Air temperature", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$light, sim_res_nb$scaledResiduals,
              xlab = "Light intensity", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$shrub_layer, sim_res_nb$scaledResiduals,
              xlab = "Shrub layer", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$insect_host_plants_log, sim_res_nb$scaledResiduals,
              xlab = "Insect host plant availability", ylab = "DHARMa simulated residuals")

# Model performance:
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
theta=summary(final_consumermodel)$theta
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(consumer_richness ~  airtemp +  light + shrub_layer + 
              insect_host_plants_log, family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1, type="response")
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
R2cv=1-(sum((eco_data$consumer_richness-PTOT)^2)/(length(PTOT))/(var(eco_data$consumer_richness)))
cor(PTOT, eco_data$consumer_richness)^2 

# save predicted consumers for later:
pred_consumers <- PTOT

# Check correlations between predicted consumers and the sum of predicted symbiotic fungi, decomposing fungi, 
# detritivores, flying insects, predators and herbivores:
pred_cons_summed <- pred_symbionts+pred_decomposer+pred_detriti+pred_flyers+pred_pred+pred_herb

cor(pred_cons_summed, eco_data$consumer_richness)^2  

# All groups ----
# Model for variable selection using delta AIC < 2 criteria:
MASS::glm.nb(all_group_richness ~natural_landscape+ph_variability_log+soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+ sfi+smi+
               light  +soil_org_C_log+litter_mass_log+flower_abundance_log+boulderPA+dungPA+I(sfi^2)+ dead_wood_debris_log+shrub_layer+temporal_continuity,
             data = eco_data) -> m_full_nb
summary(m_full_nb) 

stepAIC(m_full_nb)
# using stepAIC to reduce full model - and the delta AIC<2 criteria for
# manually reducing further (the stepAIC function does not account for squared variables 
# in the model, i.e. we have to make sure the linear effect is kept in the models
# if the squared effect is significant - therefore manual reduction is used here
# 
# Manually from here: 
MASS::glm.nb(all_group_richness ~  soil_fertility_variability_log + ph + airtemp + sfi +I(sfi^2) + smi + light + 
               soil_org_C_log + flower_abundance_log + boulderPA + dungPA + 
                shrub_layer,
             data = eco_data) -> m_full_nb #DeltaAIC criteria: soil_fertility_variability_log,boulderPA, airtemp,light, 
# smi deleted sequentially


#final model:
MASS::glm.nb(all_group_richness ~ ph + sfi+ I(sfi^2) +soil_org_C_log + flower_abundance_log +  dungPA + 
                shrub_layer,
             data = eco_data) -> final_allgroupsmodel
summary(final_allgroupsmodel) #



# Model validation: QQ-plot, residuals vs. fitted, scale-location and residuals vs. leverage plot:
par(mfrow=c(2,2))
plot(final_allgroupsmodel)

# Plotting residuals vs. variables in the model:
final_allgroupsmodel %>% 
  simulateResiduals(n = 10000) -> sim_res_nb


plotResiduals(eco_data$ph, sim_res_nb$scaledResiduals,
              xlab = "Soil ph", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$sfi, sim_res_nb$scaledResiduals,
              xlab = "Soil fertility index", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$flower_abundance_log, sim_res_nb$scaledResiduals,
              xlab = "Flower abundance", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$dungPA, sim_res_nb$scaledResiduals,
              xlab = "Dung", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$soil_org_C, sim_res_nb$scaledResiduals,
              xlab = "Soil organic C", ylab = "DHARMa simulated residuals")

plotResiduals(eco_data$shrub_layer, sim_res_nb$scaledResiduals,
              xlab = "Shrub layer", ylab = "DHARMa simulated residuals")

# Model performance:
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
theta=summary(final_allgroupsmodel)$theta
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(all_group_richness ~ ph + sfi+ soil_org_C_log + flower_abundance_log +  dungPA + 
              I(sfi^2) + shrub_layer,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1, type="response")
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
R2cv=1-(sum((eco_data$all_group_richness-PTOT)^2)/(length(PTOT))/(var(eco_data$all_group_richness)))
cor(PTOT, eco_data$all_group_richness)^2 

#Saving predicted all group richness:
pred_allgroups<- PTOT

# Check correlations between predicted 'all groups' and the sum of predicted plants, mosses, lichens, symbiotic fungi, decomposing fungi, 
# detritivores, flying insects, predators and herbivores:
pred_summed <- pred_plants+pred_moss+pred_lichen+pred_symbionts+pred_decomposer+pred_detriti+pred_flyers+pred_pred+pred_herb
cor(pred_summed, eco_data$all_group_richness)^2 

# Check correlations between the sum of predicted producers and consumers and richness of all groups:
pred_cons_prod <- pred_consumers+pred_producers
cor(pred_cons_prod, eco_data$all_group_richness)^2  



# eDNA fungi----
# Model for variable selection using delta AIC < 2 criteria:
MASS::glm.nb(fungi_OTU_rich ~  natural_landscape+ph_variability_log+soil_moisture_variability_log+ph+airtemp+ sfi+smi+light+
               +shrub_layer+dead_wood_vol_log+large_tree_density+ plant_richness_std_log+dungPA+ fungi_symbiont_plants,
             data = eco_data) -> m_full_nb
summary(m_full_nb) 
stepAIC(m_full_nb)


# using stepAIC to reduce full model - and the delta AIC<2 criteria for
# manually reducing further (the stepAIC function does not account for squared variables 
# in the model, i.e. we have to make sure the linear effect is kept in the models
# if the squared effect is significant - therefore manual reduction is used here
# 
# Manually from here:  
MASS::glm.nb(fungi_OTU_rich ~ natural_landscape + ph_variability_log + 
               ph + sfi + smi + dead_wood_vol_log + plant_richness_std_log + fungi_symbiont_plants,
             data = eco_data) -> m_full_nb # DeltaAIC criteria: sfi deleted

#final model:
MASS::glm.nb(fungi_OTU_rich ~ natural_landscape + ph_variability_log + 
               ph +  smi + dead_wood_vol_log + plant_richness_std_log + fungi_symbiont_plants,
             data = eco_data) -> final_ednafungimodel
summary(final_ednafungimodel) #

# Model validation: QQ-plot, residuals vs. fitted, scale-location and residuals vs. leverage plot:
par(mfrow=c(2,2))
plot(final_ednafungimodel)

# Plotting residuals vs. variables in the model:
final_ednafungimodel %>% 
  simulateResiduals(n = 10000) -> sim_res_nb

plotResiduals(eco_data$ph_variability_log, sim_res_nb$scaledResiduals,
              xlab = "Soil ph variability", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$ph, sim_res_nb$scaledResiduals,
              xlab = "Soil ph", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$smi, sim_res_nb$scaledResiduals,
              xlab = "Soil moisture index", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$natural_landscape, sim_res_nb$scaledResiduals,
              xlab = "Natural landscape", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$plant_richness_std_log, sim_res_nb$scaledResiduals,
              xlab = "Plant richness", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$dead_wood_vol_log, sim_res_nb$scaledResiduals,
              xlab = "Dead wood volume", ylab = "DHARMa simulated residuals")

plotResiduals(eco_data$fungi_symbiont_plants, sim_res_nb$scaledResiduals,
              xlab = "Fungi symbiont plants", ylab = "DHARMa simulated residuals")

# Model performance:
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
theta=summary(final_ednafungimodel)$theta
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(fungi_OTU_rich ~ natural_landscape + ph_variability_log + 
              ph +  smi + dead_wood_vol_log + plant_richness_std_log + fungi_symbiont_plants,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1, type="response")
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
R2cv=1-(sum((eco_data$fungi_OTU_rich-PTOT)^2)/(length(PTOT))/(var(eco_data$fungi_OTU_rich)))
cor(PTOT, eco_data$fungi_OTU_rich)^2 

# DNA flying insects----
# Model for variable selection using delta AIC < 2 criteria:
MASS::glm.nb(insect_OTU_rich ~natural_landscape+ph_variability_log+soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+ sfi+smi+light+plant_richness_std_log+flower_abundance_log +I(sfi^2),
             data = eco_data) -> m_full_nb
summary(m_full_nb)

stepAIC(m_full_nb)
# using stepAIC to reduce full model - and the delta AIC<2 criteria for
# manually reducing further (the stepAIC function does not account for squared variables 
# in the model, i.e. we have to make sure the linear effect is kept in the models
# if the squared effect is significant - therefore manual reduction is used here
# 
# Manually from here: 
MASS::glm.nb(insect_OTU_rich ~ soil_fertility_variability_log + soil_moisture_variability_log + 
               sfi + light + plant_richness_std_log + I(sfi^2) ,
             data = eco_data) -> m_full_nb # Delta AIC criteria: I(sfi^2)+soil_fertility_variability_log deleted sequentially

#final model:
MASS::glm.nb(insect_OTU_rich ~  soil_moisture_variability_log + 
               sfi + light + plant_richness_std_log  ,
             data = eco_data) -> m_full_nb
summary(m_full_nb) #

# Model validation: QQ-plot, residuals vs. fitted, scale-location and residuals vs. leverage plot:

par(mfrow=c(2,2))
plot(m_full_nb) #looks fine!

m_full_nb %>% 
  simulateResiduals(n = 10000) -> sim_res_nb
plot(sim_res_nb) # looks ok - deviation NS

testResiduals(sim_res_nb) # no overdispersion, 

# Testing for spatial autocorrelation:
testSpatialAutocorrelation(sim_res_nb, x=eco_data$utm_x,y=eco_data$utm_y) # significant , p<0.05

# We use backwards selection based on five-fold cross-validation to reduce to final model instead:
MASS::glm.nb(insect_OTU_rich ~ natural_landscape+ph_variability_log+soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+ 
               sfi+smi+light+plant_richness_std_log+flower_abundance_log +I(sfi^2),
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 1189.8

theta=summary(m_full_nb)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(insect_OTU_rich ~ natural_landscape+ph_variability_log+soil_fertility_variability_log+soil_moisture_variability_log+ph+airtemp+ 
              sfi+I(sfi^2)+smi+light+plant_richness_std_log+flower_abundance_log ,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$insect_OTU_rich)^2 

# Run the above code for all models - taking out one variable at a time
# full CV R2: 0.2148022 
# - flower_abundance_log: 0.2237306
# - plants_log : 0.2073788
# - light: 0.06007104
# - smi: 0.2169403
# - sfi: 
# - sfi^2:0.1981472
# - airtemp: 0.2215688
# - ph: 0.2291165 biggest increase in CV R2 - take out
# - soilm_sd_may : 0.1980968
# - soil_fertility_variability_log: 0.1871068
# - ph_variability_log: 0.2158924 
# - natural_landscape: 0.2186753

# Backwards 2:
MASS::glm.nb(insect_OTU_rich ~ natural_landscape+ph_variability_log+soil_fertility_variability_log+soil_moisture_variability_log+airtemp+ 
               sfi+smi+light+plant_richness_std_log+flower_abundance_log +I(sfi^2),
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 1187.8

theta=summary(m_full_nb)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(insect_OTU_rich ~ natural_landscape+ph_variability_log+soil_fertility_variability_log+soil_moisture_variability_log+airtemp+ 
              sfi+I(sfi^2)+smi+light+plant_richness_std_log+flower_abundance_log ,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$insect_OTU_rich)^2 
# Run the above code for all models - taking out one variable at a time
# full CV R2: 0.2291165 # 
# - flower_abundance_log: 0.2389478  biggest increase in CV R2 - take out
# - plant_richness_std_log : 0.2132893
# - light: 0.0634009
# - smi: 0.2314609
# - sfi: 
# - sfi^2:0.2157656
# - airtemp: 0.2360301
# - soil_moisture_variability_log : 0.2114804
# - soil_fertility_variability_log: 0.2050575
# - ph_variability_log: 0.2300251 
# - natural_landscape: 0.2318445

# Backwards 3:
MASS::glm.nb(insect_OTU_rich ~ natural_landscape+ph_variability_log+soil_fertility_variability_log+soil_moisture_variability_log+airtemp+ 
               sfi+smi+light+plant_richness_std_log +I(sfi^2),
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 1186.3

theta=summary(m_full_nb)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(insect_OTU_rich ~ natural_landscape+ph_variability_log+soil_fertility_variability_log+soil_moisture_variability_log+airtemp+ 
              sfi+I(sfi^2)+smi+light+plant_richness_std_log ,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$insect_OTU_rich)^2 
# Run the above code for all models - taking out one variable at a time
# full CV R2:  0.2389467
# - plant_richness_std_log : 0.2106768
# - light: 0.06469624
# - smi: 0.240302
# - sfi: 
# - sfi^2:0.2249307
# - airtemp: 0.2447366 biggest increase in CV R2 - take out
# - soil_moisture_variability_log : 0.2154837
# - soil_fertility_variability_log: 0.2162413
# - ph_variability_log: 0.2397099 
# - natural_landscape: 0.2417309

# Backwards 4:
MASS::glm.nb(insect_OTU_rich ~ natural_landscape+ph_variability_log+soil_fertility_variability_log+soil_moisture_variability_log+ 
               sfi+smi+light+plant_richness_std_log +I(sfi^2),
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 1184.5

theta=summary(m_full_nb)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(insect_OTU_rich ~ natural_landscape+ph_variability_log+soil_fertility_variability_log+soil_moisture_variability_log+ 
              sfi+I(sfi^2)+smi+light+plant_richness_std_log ,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$insect_OTU_rich)^2 
# Run the above code for all models - taking out one variable at a time
# full CV R2: 0.2447363 
# - plant_richness_std_log : 0.2155872
# - light: 0.06174184
# - smi: 0.2462553
# - sfi: 
# - sfi^2:0.2329698
# - soil_moisture_variability_log : 0.2239139
# - soil_fertility_variability_log: 0.2211403
# - ph_variability_log: 0.2452247 
# - natural_landscape: 0.2477612 biggest increase in CV R2 - take out

# Backwards 5:
MASS::glm.nb(insect_OTU_rich ~ ph_variability_log+soil_fertility_variability_log+soil_moisture_variability_log+ 
               sfi+smi+light+plant_richness_std_log +I(sfi^2),
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 1182.6

theta=summary(m_full_nb)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(insect_OTU_rich ~ ph_variability_log+soil_fertility_variability_log+soil_moisture_variability_log+ 
              sfi+I(sfi^2)+smi+light+plant_richness_std_log ,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$insect_OTU_rich)^2 
# Run the above code for all models - taking out one variable at a time
# full CV R2:  0.247761 
# - plant_richness_std_log : 0.2211346
# - light: 0.06451802
# - smi: 0.2506747 biggest increase in CV R2 - take out
# - sfi: 
# - sfi^2:0.2352292
# - soil_moisture_variability_log : 0.2294701
# - soil_fertility_variability_log: 0.2293346
# - ph_variability_log: 0.2482972 

# Backwards 6:
MASS::glm.nb(insect_OTU_rich ~ soil_fertility_variability_log+soil_moisture_variability_log+ ph_variability_log+
               sfi+light+plant_richness_std_log +I(sfi^2),
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 1180.7

theta=summary(m_full_nb)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(insect_OTU_rich ~ soil_fertility_variability_log+soil_moisture_variability_log+ ph_variability_log+
              sfi+light+plant_richness_std_log +I(sfi^2)  ,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$insect_OTU_rich)^2 
# Run the above code for all models - taking out one variable at a time
# full CV R2:   0.2506744 
# - plant_richness_std_log : 0.2303178
# - light: 0.06602455
# - sfi: 
# - sfi^2:0.2376822
# - soil_moisture_variability_log : 0.2342475
# - soil_fertility_variability_log: 0.2333979
# - ph_variability_log: 0.2536961 biggest increase in CV R2 - take out


# Backwards 7:
MASS::glm.nb(insect_OTU_rich ~ soil_fertility_variability_log+soil_moisture_variability_log+
               sfi+light+plant_richness_std_log +I(sfi^2),
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 1179.1

theta=summary(m_full_nb)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(insect_OTU_rich ~ soil_fertility_variability_log+soil_moisture_variability_log+ 
              sfi+I(sfi^2)+light+plant_richness_std_log ,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$insect_OTU_rich)^2 
# Run the above code for all models - taking out one variable at a time
# full CV R2:0.2536942   
# - plant_richness_std_log : 0.2293785 
# - light: 0.06983434
# sfi: 
# sfi^2:0.2416447 biggest increase in CV R2 - take out
# soil_moisture_variability_log : 0.2384677
# soil_fertility_variability_log: 0.2395862

# Backwards 8:
MASS::glm.nb(insect_OTU_rich ~ soil_fertility_variability_log+soil_moisture_variability_log+ 
               sfi+light+plant_richness_std_log ,
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 1180.7

theta=summary(m_full_nb)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(insect_OTU_rich ~ soil_fertility_variability_log+soil_moisture_variability_log+ 
              sfi+light+plant_richness_std_log ,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$insect_OTU_rich)^2 
# Run the above code for all models - taking out one variable at a time
# full CV R2:0.2416181  
# - plant_richness_std_log: 0.1879101
# - light: 0.06582965
# - sfi: 0.2591986 biggest increase in CV R2 - take out
# - soil_moisture_variability_log : 0.2264395
# - soil_fertility_variability_log: 0.2339878

# Backwards 9:
MASS::glm.nb(insect_OTU_rich ~ soil_fertility_variability_log+soil_moisture_variability_log+ 
               light+plant_richness_std_log ,
             data = eco_data) -> m_full_nb
summary(m_full_nb)# AIC= 1183.7 deltaAIC>2  from previous model - we stick to the above as final model


MASS::glm.nb(insect_OTU_rich ~ soil_fertility_variability_log+soil_moisture_variability_log+ 
               sfi+light+plant_richness_std_log ,
             data = eco_data) -> final_ednaflyermodel
summary(final_ednaflyermodel)


# Model validation: QQ-plot, residuals vs. fitted, scale-location and residuals vs. leverage plot:
par(mfrow=c(2,2))
plot(final_ednaflyermodel)



# Plotting residuals vs. variables in the model:
final_ednaflyermodel %>% 
  simulateResiduals(n = 10000) -> sim_res_nb

plotResiduals(eco_data$soil_moisture_variability_log, sim_res_nb$scaledResiduals,
              xlab = "Soil moisture variability", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$soil_fertility_variability_log, sim_res_nb$scaledResiduals,
              xlab = "Soil fertility variability", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$sfi, sim_res_nb$scaledResiduals,
              xlab = "Soil fertility index", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$light, sim_res_nb$scaledResiduals,
              xlab = "Light intensity", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$plant_richness_std_log, sim_res_nb$scaledResiduals,
              xlab = "Plant richness", ylab = "DHARMa simulated residuals")


# Model performance:
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
theta=summary(final_ednaflyermodel)$theta
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(insect_OTU_rich ~ soil_fertility_variability_log+soil_moisture_variability_log+ 
              sfi+light+plant_richness_std_log ,family=negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1)
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
cor(PTOT, eco_data$insect_OTU_rich)^2


# eDNA eukaryotes----
# Model for variable selection using delta AIC < 2 criteria:
glm.nb(eukaryote_OTU_rich ~natural_landscape+ph_variability_log+ph+airtemp+ sfi+smi+light+dead_wood_vol+
         soil_org_C_log+dungPA+boulderPA+soil_org_matter+insect_host_plants +spatial_continuity_log+temporal_continuity+
         I(sfi^2),
       data = eco_data) -> m_full_nb
summary(m_full_nb)

stepAIC(m_full_nb) # using stepAIC to reduce full model - and the delta AIC<2 criteria for
# manually reducing further (the stepAIC function does not account for squared variables 
# in the model, i.e. we have to make sure the linear effect is kept in the models
# if the squared effect is significant - therefore manual reduction is used here
# 

# Manually from here:
glm.nb(eukaryote_OTU_rich ~ natural_landscape + 
         ph + sfi + smi + light + soil_org_C_log + I(sfi^2),
       data = eco_data) -> m_full_nb
summary(m_full_nb)# DeltaAIC criteria: soil_org_C_log, natural_landscape, ph ,light deleted sequentially

# Final model:
glm.nb(eukaryote_OTU_rich ~ sfi + smi +  I(sfi^2),
       data = eco_data) -> final_eukaryotemodel 

summary(final_eukaryotemodel)
# Model validation: QQ-plot, residuals vs. fitted, scale-location and residuals vs. leverage plot:
par(mfrow=c(2,2))
plot(final_eukaryotemodel)

# Plotting residuals vs. variables in the model:
final_eukaryotemodel %>% 
  simulateResiduals(n = 10000) -> sim_res_nb


plotResiduals(eco_data$smi, sim_res_nb$scaledResiduals,
              xlab = "Soil moisture index", ylab = "DHARMa simulated residuals")
plotResiduals(eco_data$sfi, sim_res_nb$scaledResiduals,
              xlab = "Soil fertility index", ylab = "DHARMa simulated residuals")

# Model performance:
region <- c("Njut","Wjut","Ejut",  "Zeal" ,"FLM"  )
theta=summary(final_eukaryotemodel)$theta
PTOT=NULL
for (i in region)
{
  ##Data that will be predicted
  DataC1=eco_data[eco_data$region==i,]
  ###To train the model
  DataCV=eco_data[!eco_data$region==i,]
  M1 <- glm(eukaryote_OTU_rich ~ sfi + smi +  I(sfi^2),
            family = negative.binomial(theta=theta, link=log), data = DataCV)
  P1=predict(M1, DataC1, type="response")
  names(P1)=NULL
  P1
  PTOT= c(PTOT, P1)
}
R2cv=1-(sum((eco_data$eukaryote_OTU_rich-PTOT)^2)/(length(PTOT))/(var(eco_data$eukaryote_OTU_rich)))
cor(PTOT, eco_data$eukaryote_OTU_rich)^2 


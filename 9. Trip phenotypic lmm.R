#########################################################
# This code: 
# --- For chinstraps, Trip behaviour = Max distance ---
# Does model selection using lme4
# 
# Predicts effects from best model
# Plots figures of model output
# 
# chinstraps: Inc vs. Chick rearing
# 
# Leandri
# August 2022
#############################################################

# load packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(MuMIn)
library(ggplot2)
library(nlme)
library(effects)
library(sjPlot)
library(lattice)
library(ggeffects)
library(broom.mixed)
library(MetBrewer)

#---------------------------------------------------------------------------------------------
# read in trip metrics 
# List clean data files from the directory that you want
cleandata <- list.files(path = "C:/Users/cwoosthuizen/Desktop/krilltokt penguins/outputs/summary",
                        pattern = '.*trip_metrics_chinstraps.*\\.rds$',
                        full.names = TRUE)

print(cleandata)

# make empty list
cleanlist = list()

for (i in cleandata) {             # to run on all individuals in the folder
  # Import tdr data 
  cleanlist[[i]] = readRDS(i)
}

trip = as_tibble(bind_rows(cleanlist))  # 742 trips


head(trip)
names(trip)

# Change Bro and Cre to chick-rearing phase:
#--------------------------------------------------------------
# rename 
trip$phase = if_else(trip$group == "Bro", "Chick", trip$group)
trip$phase = if_else(trip$phase == "Cre", "Chick", trip$phase)

unique(trip$phase)
unique(trip$group)
#--------------------------------------------------------------------


# dist.max
trip$trip.duration <- as.numeric(trip$trip.duration.hr)

# quick correlation test
# library(corrplot)
# trip.cor = trip %>%
#   dplyr::select(trip.duration, dist.max, col.dist, dist.max, trip.dist)
# 
# n = cor(trip.cor)
# corrplot(n, method = 'number')
# trip dist and dist max are highly correlated
# trip duration and trip distance are highly correlated
# and trip duration and dist max are highy correlated, so can use just one response variable to describe the trip
# behaviours between chinstraps 

# Do LMM models for max distance:

# response: dist.max:
plot(trip$trip.duration, trip$dist.max)
ggplot(trip, aes(x=trip.duration, col = island))+geom_density() # normally distributed

trip = trip %>%
filter(!dist.max >= 400) 

ggplot(trip, aes(x=dist.max, col = island))+geom_density() # normally distributed

ggplot(trip, aes(x=dist.max, col = island))+geom_density() # normally distributed
ggplot(trip, aes(x=dist.max, col = solar.ele.start.day))+geom_density()
ggplot(trip, aes(y=dist.max, x = solar.ele.trip.start, col = island))+geom_point()
unique(trip$solar.ele.start.day)

# models: 
# 1. null model:
# Testing the random effects of id (set REML as T) with the most complex model: 
ran.mod1 <- lmer(dist.max ~ phase*island + (1|id), trip, REML = T)
ran.mod2 <- gls(dist.max ~ phase*island, trip, method = 'REML')
#mod1 <- glmer(dist.max ~ 1 + (1|id), trip, family = Gamma(link = "inverse"))
AIC(ran.mod1, ran.mod2)

summary(ran.mod1)

# 2. Model selection: Fit mixed-effects models with id as random effect
# Baseline explanatory variables: 
#  Phase ( 2 Levels = 'Inc' or 'Chick')
#  Group: (2 levels = 'Inc' or 'Bro' or 'Cre')
# * I want to see if it will be better to split breeding season into only two phases or into the three different breeding
# phases as recorded in the field
#  Island (2 levels = 'Kopaitic' or 'Nelson')

# null
nullmod <- lmer(dist.max ~ 1 + (1|id), trip, REML = F)

# island effects
mod1 <- lmer(dist.max ~ island + (1|id), trip, REML = F)

# Breeding stage effects
mod2 <- lmer(dist.max ~ phase + (1|id), trip, REML = F)
mod3 <- lmer(dist.max ~ group + (1|id), trip, REML = F)

AIC(nullmod, mod1, mod2, mod3)

# retain PHASE as explanatory variable

# Breeding stage and island effects
mod4 <- lmer(dist.max ~ phase + island + (1|id), trip, REML = F)
mod5 <- lmer(dist.max ~ phase*island + (1|id), trip, REML = F)
mod6 <- lmer(dist.max ~ group + island + (1|id), trip, REML = F)
mod7 <- lmer(dist.max ~ group*island + (1|id), trip, REML = F)
AIC(nullmod, mod1, mod2, mod3, mod4, mod5, mod6, mod7)

summary(mod5)

models <- c('nullmod', "mod1","mod2", "mod3","mod4", "mod5", 'mod6', 'mod7')
aics <- AIC(nullmod, mod1,mod2,mod3,mod4,mod5, mod6, mod7)
delta.aics <- aics$AIC - min(aics$AIC) # To work out the change in delta AIC #see definitions in book (from Anderson (2008))
exp.delta <- exp(-0.5*delta.aics)
wi <- exp.delta/sum(exp.delta)# these are the Akaike weights for each model #The probability that model is the actual (fitted) k-l best model in the set (Anderson 2008) (See Burnham et al (2011) paper) 
(modtable <- data.frame(models, numpar=aics$df, aics$AIC,  delta.aics,wi))
# Clearly mod 5 is the best model to use 

#-----
# Add in deviance? 

# can't get deviance from glmmTMB model object: 
# https://github.com/glmmTMB/glmmTMB/issues/621
# HOWEVER: there are some fairly deep issues relating to the definition of 'deviance' for GLMMS: 
# see here (https://github.com/lme4/lme4/blob/master/man/merMod-class.Rd#L201-L224). In terms of that 
# discussion, the function given above returns the "absolute unconditional" deviance - which won't, 
# for example, match the deviance() value of the corresponding result from glm()

# https://stackoverflow.com/questions/46879494/using-broomglance-in-conjunction-with-glmmtmb

#deviance.lme <- function(object, ...) {
#   -2*c(logLik(object))
#}
#
#dev.m1 <- deviance.lme(m1)

# extract deviance on its own for every model:
#-- output from model
dev.null <- glance(nullmod)
dev.m1 <- glance(mod1)
dev.m2 <- glance(mod2)
dev.m3 <- glance(mod3)
dev.m4 <- glance(mod4)
dev.m5 <- glance(mod5)
dev.m6 <- glance(mod6)
dev.m7 <- glance(mod7)


#-- real deviance:
dev.null <- dev.null$deviance
dev.m1 <- dev.m1$deviance
dev.m2 <- dev.m2$deviance
dev.m3 <- dev.m3$deviance
dev.m4 <- dev.m4$deviance
dev.m5 <- dev.m5$deviance
dev.m6 <- dev.m6$deviance
dev.m7 <- dev.m7$deviance


dev = c(dev.null, dev.m1, dev.m2, dev.m3, dev.m4, dev.m5, dev.m6, dev.m7)

(modtable <- data.frame(models,  numpar=aics$df, aics$AIC, delta.aics,wi, dev))

write.csv(modtable, './outputs/phenotypic lmms/trip/lme4_max dist_chinstraps_model selection table.csv')

# use mod 5
summary(mod5)
summary(nullmod)

# Now refit the most parsimonious model with REML = T and plot it:
mod.r <- lmer(dist.max ~ phase*island + (1|id), trip, REML = T)
summary(mod.r)

#nullmod.r <- lmer(dist.max ~ 1 + (1|id), trip, REML = T)
#summary(nullmod.r)

# quick visualisation of results: 
plot(ggpredict(mod.r, terms = c("phase", "island")))

#----------------------------------------------------------------------
# 9. Model checking: 
library(performance)
library(see)
# check model assumptions
check_model(mod.r)

plot(mod.r)
acf(resid(mod.r)) # looks good
#qqplot(mod1)

# #Check for linearity
# linearity<-plot(resid(mod.r),#extract the residuals
#                 trip$dist.max) #specify original y variable
# #Check for normality
# #normality <- qqmath(mod4) # not working
# qqnorm(resid(mod.r)); qqline(resid(mod.r))
# 
# plot(mod.r, trip.duration ~ resid(., scaled=F))
# 
#Try with model_plot (argument for type can be varied)
#plot_model(mod.r, type='diag') # you can ajust type (see package info: ?plot_model)


# Examine random effects:
dotplot(ranef(mod.r, condVar = T))

performance::icc(mod.r) # repeatability estimate
performance::r2(mod.r)

#----------------------------------------------------------------------------------------


# 6. R-squared values
r.squaredGLMM(mod.r)
# r.squaredGLMM(nullmod.r)
# The marginal R-squared considers only the variance of the fixed effects, while the conditional R-squared 
# takes both the fixed and random effects into account.

# 8.  Summary of model output in an HTML table

#-----------------------------------------------------------
# Now do it the easy way:
# https://cran.r-project.org/web/packages/sjPlot/vignettes/tab_mixed.html
# Summary of Mixed Models as HTML Table
# Daniel L?decke  2020-05-23
# library(sjPlot)
#----------------------------------------------------------

#output for one model
tab_model(mod.r, CSS = css_theme("cells"), 
          file = paste0('./outputs/phenotypic lmms/trip/lme4_max dist_chinstraps_best model summary.doc'))
#tab_model(m1, best.mod, CSS = css_theme("cells")

tab_model(mod.r, CSS = css_theme("cells"))
# more precise p-value can be computed based on conditional F-tests with Kenward-Roger approximation
# tab_model(basemodREML, p.val = "kr", show.df = TRUE)  #  takes a very long time to run
# variance = residual variance
# tau00 ID = Variance from ID as random effect
# ICC (Inter-correlation coefficient) = repeatability of individuals
# N ID = Number of individuals
# Observations: Number of trip duration observations
# Marginal R-squared: Variance of fixed effects
# Conditional R-squared: Variance of fixed and random effects


# Customizing HTML tables - Daniel L?decke
# see https://cran.r-project.org/web/packages/sjPlot/vignettes/table_css.html


##---------------------- Plot ------------------------------------------------------------##
#library(ggeffects)
summary(mod.r)

# the trips for chinstraps at Nelson are much longer compared to the trips during creche. Trip durations for 
# chinstraps at Kopaitic are more similar between breeding phases

p = ggpredict(mod.r, terms = c('phase', 'island'))
p
plot(p)

pred <- as.data.frame(p)
head(pred)


ggplot(pred, aes(x, predicted, col = group)) +
  geom_point(position = position_dodge(.5)) +
  #facet_wrap(~facet, ncol = 3) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    position = position_dodge(.5)
  ) #+
#scale_x_discrete(breaks = 1:3, labels = get_x_labels(pred))

# plot dist.max output from lmm model: 
pred = pred %>%
  dplyr::rename(phase = x, 
                dist.max = predicted, 
                island = group)

head(pred)

ggplot() + 
  geom_point(data = pred, aes(phase, dist.max, col = island))#+ 
  #facet_wrap(~phase, ncol = 3) #+
#geom_line(pred, aes(ymin = conf.low, ymax = conf.high, col = phase), alpha = .1) 
#geom_line(aes(pred, aes(y = conf.low, col = phase)))

#ggplot(pred, aes(x = x, y = predicted, colour = phase)) +
#  stat_smooth(method = "lm", se = FALSE) +
#  facet_wrap(~facet, ncol = 2)

# specify facet plot order - predicted
pred$phase_f = factor(pred$phase, levels=c('Inc','Chick'))
levels(pred$phase_f) <- list(Incubation = 'Inc', Chick = 'Chick') 
pred$phase_f

# specify facet plot order - observed
trip$phase_f = factor(trip$phase, levels=c('Inc','Chick'))
levels(trip$phase_f) <- list(Incubation = 'Inc', Chick = 'Chick') 
trip$phase_f

# specify facet plot order - predicted
pred$island_f = factor(pred$island, levels=c('Kopaitic','Nelson'))
levels(pred$island_f) <- list(Nelson = 'Nelson', Kopaitic = 'Kopaitic') 
pred$island_f

# specify facet plot order - observed
trip$island_f = factor(trip$island, levels=c('Kopaitic','Nelson'))
levels(trip$island_f) <- list(Nelson = 'Nelson', Kopaitic = 'Kopaitic') 
trip$island_f


library(MetBrewer)
met.brewer("Derain",n=2,type="discrete")

s = ggplot(data = pred, aes(phase_f, dist.max, col = island)) +
  geom_point(position = position_dodge(.8), size =5) +   #facet_wrap(~phase, ncol = 3) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    position = position_dodge(.8), size = 1.2, width = 0.6
  )+ 
  scale_colour_manual(values=met.brewer("Derain", 2), name = 'Island')
#geom_point(trip, mapping = aes(island, dist.max, col = species), position = position_dodge(.5), 
#shape = 1) #+ facet_wrap(~phase, ncol =3)
s

# try adding observed points to the predicted estimates
s <- s + geom_point(trip, mapping = aes(x = phase_f, y = dist.max, col = island), 
                    position = position_jitterdodge(.3),  shape = 4, alpha = 0.3)+
 # facet_wrap(~phase_f) + 
  coord_cartesian(ylim=c(NA, 200))                        # cut the y scale at 200: removes about 3/4 points
s  
#scale_x_discrete(breaks = 1:3, labels = get_x_labels(pred))
s <- s + theme_bw() + 
  theme(panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              strip.background = element_blank(),
                              panel.border = element_rect(colour = "black"),
                              strip.text.x = element_text(size = 14), 
                              strip.text.y = element_text(size = 14), 
                             #panel.grid.major.x = element_text(size = 14), 
                              # panel.grid.major.y = element_text(size = 14),
                               legend.text = element_text(size = 14), 
                              legend.title = element_text(size = 16),
                              legend.position = c(0.90, 0.85),
                              axis.text.x = element_text(size = 12, colour = "black"), 
                              axis.text.y = element_text(size = 12, colour = "black"), 
                                axis.title = element_text(size = 16))+
                      ggtitle("Chinstrap penguins")+
                    xlab("Breeding stage") +
                      ylab("Maximum distance (km)")
s

# save image 
rstudioapi::savePlotAsImage(paste0("./outputs/phenotypic lmms/trip/plots/", trip$species[1], "_max distance phenotypic lmm.png"), width=900,height=600)

# add in violin plots
s = ggplot(data =trip, mapping = aes(phase_f, dist.max, col =island_f))+
  geom_point(data = trip, position = position_jitterdodge(dodge.width = .9, jitter.width = 0.3),  shape = 4, alpha = 0.6)+
  geom_violin(trim = T, alpha = 0.5, size = 0.6) +  #facet_wrap(group_f~island, ncol = 6) +
 # coord_cartesian(ylim= c(0.0, 0.5))+ 
  geom_point(data= pred, aes(phase_f, dist.max, group= island_f), 
                  position = position_dodge(.9), size = 5, alpha = 3, col = 'grey20')+
  geom_errorbar(data= pred, aes(phase_f, dist.max, ymin = conf.low, ymax = conf.high, group = island_f), 
                position = position_dodge(.9), size = 1.2, width = 0.6, col = 'grey20')+
  scale_color_manual(values=met.brewer("Derain", 2, direction = -1), name = 'Island')+
   coord_cartesian(ylim=c(NA, 150)) 
  #scale_colour_manual(values = c('day'= 'orange', 'twilight' = 'purple'))
  #geom_point(size = 1.5, alpha = 3, position = position_dodge(.9))+
  #geom_boxplot(width=0.1, position = position_dodge(.9), size = 0.1)
  #stat_summary(fun.data=mean_sdl, mult=1, 
   #              geom="pointrange", color="red")
  s

  s <- s + theme_bw() + 
  theme(panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              strip.background = element_blank(),
                              panel.border = element_rect(colour = "black"),
                              strip.text.x = element_text(size = 14), 
                              strip.text.y = element_text(size = 14), 
                             #panel.grid.major.x = element_text(size = 14), 
                              # panel.grid.major.y = element_text(size = 14),
                               legend.text = element_text(size = 14), 
                              legend.title = element_text(size = 16),
                              legend.position = c(0.90, 0.85),
                              axis.text.x = element_text(size = 12, colour = "black"), 
                              axis.text.y = element_text(size = 12, colour = "black"), 
                                axis.title = element_text(size = 16))+
                      ggtitle("Chinstrap penguins")+
                    xlab("Breeding stage") +
                      ylab("Maximum distance (km)")
s

# save image 
rstudioapi::savePlotAsImage(paste0("./outputs/phenotypic lmms/trip/plots/", trip$species[1], "_max distance phenotypic lmm_VIOLIN.png"), width=900,height=600)


# Adding in marginal histogram?
# https://r-graph-gallery.com/277-marginal-histogram-for-ggplot2.html
library(ggExtra)

# normal plot
s1 = ggplot(data = pred, aes(phase_f, dist.max, col = island)) +
  geom_point(position = position_dodge(.8), size =5) +   #facet_wrap(~phase, ncol = 3) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    position = position_dodge(.8), size = 1.2
  ) + 
  scale_colour_manual(values=met.brewer("Derain", 2))#+
#geom_point(trip, mapping = aes(island, trip.duration, col = species), position = position_dodge(.5), 
#shape = 1) #+ facet_wrap(~phase, ncol =3)
s1

# try adding observed points to the predicted estimates
s1 <- s1 + geom_point(trip, mapping = aes(x = phase_f, y = dist.max, col = island), 
                      position = position_jitterdodge(.35),  shape = 4, alpha = 0.2)+
  #facet_wrap(~phase, ncol = 3) + 
  coord_cartesian(ylim=c(NA, 150))                        # cut the y scale at 150: removes 3 points
s1  
#scale_x_discrete(breaks = 1:3, labels = get_x_labels(pred))
s1 <- s1 + theme_bw() + 
  theme(legend.position = 'bottom', 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        #strip.text.x = element_text(size = 20), 
        #strip.text.y = element_text(size = 20), 
        #panel.grid.major.x = element_text(size = 14), 
        #panel.grid.major.y = element_text(size = 14),
        legend.text = element_text(size = 9), 
        legend.title = element_text(size = 9), 
        axis.text.x = element_text(size = 10, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"), 
        axis.title = element_text(size = 12))+
  ggtitle('Chinstraps')+
  xlab("Breeding phase") +
  ylab("Maximum distance (km)")
s1

s2 <- ggMarginal(s1, type = 'density', groupColour = TRUE, groupFill = TRUE, margins = 'y')
s2

# save image 
#rstudioapi::savePlotAsImage(paste0("./plots/", trip$species[1], "_max distance phenotypic lmm_marginal.png"), width=900,height=600)

 
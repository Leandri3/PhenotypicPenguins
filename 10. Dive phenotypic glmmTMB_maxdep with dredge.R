#########################################################
# This code: 
# --- For gentoos, Dive behaviour = Max depth ---
# 
# Do model selection using - glmmTMB
# Includes autocorrelation structures
# family = Gamma(link = "log"), 
# 
# Predicts effects from best model
# Plots figures of model output
# 
# Leandri
# October 2022
#############################################################

# load packages
library(tidyverse)
#library(lme4)
library(lmerTest)
library(car)
library(MuMIn)
library(ggplot2)
#library(nlme)
library(effects)
library(sjPlot)
library(lattice)
library(ggeffects)
library(broom.mixed)
library(arm)
library(mgcv)
library(glmmTMB)
library(DHARMa)
library(MetBrewer)


#-----------------------------------------------------------------
# read in dive metrics
gent <- readRDS('./outputs/environmental covariates/gentoos_foraging_5m_dive metrics_and_extracted_environmental_covariates.rds')

# include all dive data from both populations
dive <- as_tibble(bind_rows(gent))

head(dive)
names(dive)

# No need to change Bro and Cre to chick-rearing phase: already done in previous step
#--------------------------------------------------------------
# rename 
#dive$phase = if_else(dive$group == "Bro", "Chick", dive$group)
#dive$phase = if_else(dive$phase == "Cre", "Chick", dive$phase)

unique(dive$phase)
unique(dive$group)
#--------------------------------------------------------------------

dive = dive %>%
  rename(feed= bathy)# %>%
  #rename(bathy = depth) %>%
dive = dive %>%
  rename(bathy = depth)  

names(dive)
str(dive)

# make sure date.time is set to GMT (or UTC)
attr(dive$begdesc, "tzone") # check time zone

# how many dive observations per island?
dive.isl = dive %>%
  group_by(island)%>%
  summarise(n = n())
dive.isl

# response variable = max depth
hist(dive$maxdep)
ggplot(dive, aes(x=maxdep))+geom_density(adjust = 5)
ggplot(dive, aes(x=maxdep, col = island))+geom_density()
#ggplot(dive, aes(x=island, y = maxdep, col = island))+geom_hex()
ggplot(dive, aes(x=maxdep, col = phase))+geom_density()
ggplot(dive, aes(x=solar.ele.begdesc, y = maxdep, col = island))+geom_point()


plot(dive$maxdep)

# more variables to add in: 
dive <- dive %>%
  rename(diurnal = solar.ele.begdesc.day)
names(dive)

# diurnal patterns: 
ggplot(dive, aes(diurnal, maxdep, fill = island)) + geom_boxplot() + facet_wrap(island~.)
#ggplot(dive, aes(diurnal, maxdep)) + geom_hex(bins = 10) + facet_wrap(island~.)
ggplot(dive, aes(diurnal, maxdep, fill = island)) + geom_boxplot() + facet_wrap(island~group)


# solar elevation
ggplot(dive, aes(solar.ele.begdesc, maxdep, col = island)) +geom_point(shape = '.')+ facet_grid(island~phase)
ggplot(dive, aes(solar.ele.begdesc, maxdep, col = island)) +geom_point(shape = '.')+ facet_grid(island~group)


# add in bathymetry?/ proportion of use of water column
# 
# dive$bathy.std <- scale(dive$bathy)
# ggplot(dive, aes(bathy, maxdep, col = island))+ geom_point()
# ggplot(dive, aes(bathy.std, maxdep, col = island))+ geom_point()
# ggplot(dive, aes(bathy.prop.100, maxdep, col = island))+ geom_point()

## quick correlation test  
#library(corrplot)
#dive.cor = dive %>%
#  dplyr::select(maxdep, divetim, botttim, dive.res, solar.ele.begdesc)
#str(dive.cor)
#n = cor(dive.cor)
#corrplot(n, method = 'number')
# Max dep and dive time are highly correlated: can use either max depth or dive time
# Proportion of water column use might be interesting to look at, but quite highly correlated with max depth
# Use either dive efficiency (Kokubun et al 2010, Zimmer et al 2010) 
# OR dive residuals  (Bestley et al 2015, Lowther et al 2015)


# Calculate 'time' for correlation structure:
g = dive %>% 
  group_by(ID) %>%   
  dplyr::mutate(t = difftime(begdesc, begdesc[1], units='hours'))  

g$t2 = plyr::round_any(as.numeric(g$t), 0.5, f = ceiling) 
length(unique(g$t2))

g$times = glmmTMB::numFactor(as.numeric(g$t2))

# model selection: -using: 
# 1. glmmTMB()
# Gamma distribution with log link
# ou (Ornstein-Uhlenbeck) autocorrelation structure
# with method = 'ML'

## dredge model selection:
# ---------------------------------------------------------------------------------------
# Fit all possible models with MUMIn::dredge()
# ---------------------------------------------------------------------------------------
# library(MuMIn)
# options(na.action = "na.fail")
#  
#  fit1 <- glmmTMB(maxdep ~solar.ele.begdesc*group*island+(1|ID)+
#               ou(times + 0|ID),
#               data= g, 
#              family = Gamma(link = "log"), 
#               REML = F)
#  
#  summary(fit1)
#  
#  m1 = dredge(fit1)# , m.lim = c(1, 4))
#  m1
#  
#  ms1 <- subset(m1, delta < 5)
#  ms1
#  
# mod1 = glmmTMB(maxdep ~ solar.ele.begdesc*group*island+(1|ID)+
#              ou(times + 0|ID),
#              data= g, 
#              family = Gamma(link = "log"), 
#              REML = F)
# summary(mod1)
# # or: 
# #mod1 = glmmTMB(maxdep ~ solar.ele.begdesc*group*island+(1|ID)+
# #             ou(times + 0|ID),
# #             data= g, 
# #             family = binomial(link = 'logit'),
# #             REML = T)
# #-------------------------------------------------------------------------------
# plot_model(mod1)
# 
# performance::icc(mod1) 
# performance::r2(mod1) 
# 
# plot(ggpredict(mod1, terms = c('solar.ele.begdesc','island','group')))
# 
# summary(mod1)   
# tab_model(mod1, CSS = css_theme("cells"), show.re.var = F)
# 
# 
# #-------------------------------------------------------------------------------------------
#  fit2 <- glmmTMB(maxdep ~solar.ele.begdesc*phase*island+(1|ID)+
#               ou(times + 0|ID),
#               data= g, 
#               family = Gamma(link = "log"), 
#               REML = F)
#  
#  summary(fit2)
#  
#  
#  m2 = dredge(fit2)# , m.lim = c(1, 4))
#  m2
#  
#  ms2 <- subset(m2, delta < 5)
#  ms2
#  
# mod2 = glmmTMB(maxdep ~ solar.ele.begdesc*phase*island+(1|ID)+
#             ou(times + 0|ID),
#             data= g, 
#             family = Gamma(link = "log"), 
#             REML = F)
# 
# plot_model(mod2)
# 
# performance::icc(mod2) 
# performance::r2(mod2) 
# 
# plot(ggpredict(mod2, terms = c('solar.ele.begdesc','island')))
# 
# summary(mod2)   
# tab_model(mod2, CSS = css_theme("cells"), show.re.var = F)
# 
#  #--------------------------------------------------------------------------------------------------------------
#  
#   fit3 <- glmmTMB(maxdep ~diurnal*phase*island+(1|ID)+
#               ou(times + 0|ID),
#               data= g, 
#               family = Gamma(link = "log"), 
#               REML = F)
#  
#  summary(fit3)
#  
#  
#  m3 = dredge(fit3)# , m.lim = c(1, 4))
#  m3
#  
#  ms3 <- subset(m3, delta < 5)
#  ms3
#  
# mod3 = glmmTMB(maxdep ~ diurnal+(1|ID)+
#             ou(times + 0|ID),
#             data= g, 
#             family = Gamma(link = "log"), 
#             REML = F)
# 
# plot_model(mod3)
# 
# performance::icc(mod3) 
# performance::r2(mod3) 
# 
# plot(ggpredict(mod3, terms = c('diurnal')))
# 
# summary(mod3)   
# tab_model(mod3, CSS = css_theme("cells"), show.re.var = F)
#  #--------------------------------------------------------------------------------------------------------------
#  
#   fit4 <- glmmTMB(maxdep ~diurnal*group*island+(1|ID)+
#               ou(times + 0|ID),
#               data= g, 
#               family = Gamma(link = "log"), 
#               REML = F)
#  
#  summary(fit4)
#  
#  
#  m4 = dredge(fit4)# , m.lim = c(1, 4))
#  m4
#  
#  ms4 <- subset(m4, delta < 5)
#  ms4
#  
# mod4 = glmmTMB(maxdep ~ diurnal + group*island+(1|ID)+
#             ou(times + 0|ID),
#             data= g, 
#             family = Gamma(link = "log"), 
#             REML = F)
# 
# plot_model(mod4)
# 
# performance::icc(mod4) 
# performance::r2(mod4) 
# 
# plot(ggpredict(mod4, terms = c('diurnal', 'group', 'island')))
# 
# summary(mod4)   
# tab_model(mod4, CSS = css_theme("cells"), show.re.var = F)
# 
# AIC(mod1, mod2, mod3, mod4)
# 
# so mod1 is the best. 
# rerun GLMER With REML = T

best.mod = glmmTMB(maxdep ~ solar.ele.begdesc*group*island+(1|ID)+
            ou(times + 0|ID),
            data= g, 
            family = Gamma(link = "log"), 
            REML = T)

plot_model(best.mod)

performance::icc(best.mod) 
performance::r2(best.mod) 

plot(ggpredict(best.mod, terms = c('solar.ele.begdesc','island','group')))

summary(best.mod)   
tab_model(best.mod, CSS = css_theme("cells"))
acf(resid(best.mod))

#------------------------------------------
# check model performance: 
library(DHARMa)

res = simulateResiduals(best.mod)
plot(res, asFactor = T)


# # add in dispersion formula? - to improve homogeniety in variance
#  best.mod2 =  glmmTMB(maxdep ~  diurnal*group*island+ (1|ID)+
#                ou(times + 0|ID),
#                dispformula = ~island+diurnal+group,
#                data= g, 
#                family = Gamma(link = "log"), 
#                REML = T)
#  summary(best.mod2)
#  plot(ggpredict(best.mod2, terms = c( 'group','diurnal', 'island')))
# r.squaredGLMM(best.mod2) # can't get r-squared when dispersion model is included
# tab_model( best.mod,best.mod2)

# res = simulateResiduals(best.mod2)
# plot(res, asFactor = T)

# icc(best.mod2) 
# 
# AIC(best.mod, best.mod2)
# 
#g$bathy.std <- scale(g$bathy)
#
#ggplot(g, aes(bathy, maxdep, col = island))+ geom_point()
#ggplot(g, aes(bathy.std, maxdep, col = island))+ geom_point()
# 
# best.mod3 =  glmmTMB(maxdep ~  diurnal+island+ phase+bathy.std+ (1|ID)+
#               ou(times + 0|ID),
#               dispformula = ~diurnal+phase+island,
#               data= g, 
#               family = Gamma(link = "log"), 
#               REML = T)
#summary(best.mod3)
# #
# plot(ggpredict(best.mod3, terms = c('bathy.std')))
## res2 = simulateResiduals(best.mod2)
# plot(res2)  # Looks very similar but the edges of boxplots are a bit smaller. 

#--------------------------------------------------------------------------
# https://easystats.github.io/see/articles/performance.html
library(performance)
library(see)
c <-check_collinearity(best.mod)
#check_collinearity(m9) 
plot(c)

icc(best.mod)
testDispersion(best.mod2)

#----------------------------------------------------------------------
# check model summary 

# R-squared values: 
par(mar=c(2, 2, 2, 2));qqnorm(residuals(best.mod),main=NULL)
r.squaredGLMM(best.mod)

# Variance partitioning: 
#library(broom.mixed)
v<- tidy(best.mod)

print(v, n =15)
# now pull out the variances
var.mod = tidy(best.mod, effects="ran_pars")
(var.modid = (var.mod[1,5])^2)  #Gives standard deviance => ^2 : to get variance
#(var.modResidual = (var.mod[2,4])^2)

#Total phenotypic variance:
#(var.mod.total = var.modid + var.modResidual)

# proportion of variance accounted for by random effects:
#var.modid /var.mod.total  # proportion of var accounted for by ID random effect
#var.modResidual/var.mod.total   # proportion of var accounted for by Residual random effect

# repeatability
print(VarCorr(best.mod),comp = c("Variance","Std.Dev."))
VarCorr(best.mod)$"ID"

# 8.  Summary of model output in an HTML table

#-----------------------------------------------------------
# Now do it the easy way:
# https://cran.r-project.org/web/packages/sjPlot/vignettes/tab_mixed.html
# Summary of Mixed Models as HTML Table
# Daniel L?decke  2020-05-23
# library(sjPlot)
#----------------------------------------------------------

#output for one model
tab_model(m1, best.mod, CSS = css_theme("cells"), 
          file = paste0('./outputs/phenotypic lmms/dive/glmmTMB_maxdep_gentoos_best model summary.doc'))
tab_model(m1, best.mod, CSS = css_theme("cells"))

# don't show all the random effects
tab_model(best.mod, CSS = css_theme("cells"), show.re.var = F)
# more precise p-value can be computed based on conditional F-tests with Kenward-Roger approximation
# tab_model(basemodREML, p.val = "kr", show.df = TRUE)  #  takes a very long time to run
# variance = residual variance
# tau00 ID = Variance from ID as random effect intercept
# tau11 ID.sst.std = variance from SST as random slope
# p01 ID: Correlation structure between ID and sst.std
# ICC (Inter-correlation coefficient) = repeatability of individuals; captures the within-class similarity of the covariate adjusted data values
# N ID = Number of individuals
# Observations: Number of max dive observations
# Marginal R-squared: Variance of fixed effects
# Conditional R-squared: Variance of fixed and random effects


# Customizing HTML tables - Daniel L?decke
# see https://cran.r-project.org/web/packages/sjPlot/vignettes/table_css.html


#-----------------------------
#dotplot(ranef(best.mod))
plot_model(best.mod)

#---------------------- plot ---------------------------------------------------# 
#Try with model_plot (argument for type can be varied)
#plot_model(best.mod, type='diag') # you can adjust type (see package info: ?plot_model)

summary(best.mod)
p = ggpredict(best.mod, terms = c('solar.ele.begdesc','island',  'group'))
p
plot(p)

#observed data
ggplot() +
  geom_point(dive, mapping = aes(solar.ele.begdesc, maxdep, col = island), 
             position = position_dodge(.5), shape =1, alpha = 0.5) + 
  facet_wrap(~group)


pred <- as.data.frame(p)
head(pred)

pred = pred %>%
  dplyr::rename(solar.ele.begdesc = x, 
                maxdep = predicted, 
                island= group, 
                group = facet)

head(pred)

ggplot() + 
  geom_line(data = pred, aes(x=solar.ele.begdesc, y=maxdep, col = island))+ 
  facet_wrap(~group) +
geom_ribbon(pred, mapping = aes(x=solar.ele.begdesc, y=maxdep, ymin = conf.low, ymax = conf.high, col = island), alpha = 1) 
#geom_line(aes(pred, aes(y = conf.low, col = group)))

#ggplot(pred, aes(x = x, y = predicted, colour = group)) +
#  stat_smooth(method = "lm", se = FALSE) +
#  facet_wrap(~facet, ncol = 2)

# # specify facet plot order - predicted
# pred$phase_f = factor(pred$phase, levels=c('Inc','Chick'))
# levels(pred$phase_f) <- list(Incubation = 'Inc', Chick = 'Chick') 
# pred$phase_f
# 
# # specify facet plot order - observed
# dive$phase_f = factor(dive$phase, levels=c('Inc','Chick'))
# levels(dive$phase_f) <- list(Incubation = 'Inc', Chick = 'Chick') 
# dive$phase_f

# specify facet plot order - predicted
pred$group_f = factor(pred$group, levels=c('Inc','Bro', 'Cre'))
levels(pred$group_f) <- list(Incubation = 'Inc', Brood = 'Bro', Creche = 'Cre') 
pred$group_f

# specify facet plot order - observed
dive$group_f = factor(dive$group, levels=c('Inc','Bro', 'Cre'))
levels(dive$group_f) <- list(Incubation = 'Inc', Brood = 'Bro', Creche = 'Cre') 
dive$group_f

library(MetBrewer)

s = ggplot(data =pred, aes(solar.ele.begdesc, maxdep, col = island)) +
  geom_line(size = 1.2, alpha = 1, position = position_dodge(.5)) +   facet_wrap(group_f~.) +
  #stat_smooth(method = "lm", se = FALSE)
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    position = position_dodge(.5), 
    alpha = 0.2, size = 0.3
  ) +
  scale_y_reverse()+ 
  scale_color_manual(values=met.brewer("Derain", 2, direction = 1), name = 'Island')
  #geom_point(trip, mapping = aes(island, trip.duration, col = species), position = position_dodge(.5), 
             #shape = 1) #+ facet_wrap(~group, ncol =3)
s

# try adding observed points to the predicted estimates
s <- s + geom_point(dive, mapping = aes(x = solar.ele.begdesc, y = maxdep, col = island), 
             position = position_jitterdodge(.7), shape = '.', alpha = 0.7) +
  facet_wrap(group_f~.) +
  coord_cartesian(ylim = c(150, NA)) 
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
                              axis.text.x = element_text(size = 12, colour = "black"), 
                              axis.text.y = element_text(size = 12, colour = "black"), 
                                axis.title = element_text(size = 16))+
                      ggtitle("Gentoo penhuins")+
                    xlab("Solar elevation (degrees above the horizon)") +
                      ylab("Maximum depth (m)")
s

# save image 
 rstudioapi::savePlotAsImage(paste0("./outputs/phenotypic lmms/dive/plots/", dive$species[1], "_maxdep_population phenotypic lmm2.png"), width=900,height=600)

# violin plots? 
library(MetBrewer)
#met.brewer(name="Renoir", n=2, type="discrete")


s = ggplot(data =dive, mapping = aes(solar.ele.begdesc, maxdep, col =island))+
  geom_point(data = dive, position = position_jitterdodge(0.8), shape = 20, alpha = 0.1)+
  geom_violin(trim = T, alpha = 0.2, size = 0.8) +  facet_wrap(island~group_f) +
 # coord_cartesian(ylim= c(0.0, 0.5))+ 
  geom_line(data= pred, aes(solar.ele.begdesc, maxdep), 
                  size = 1, alpha = 1, position = position_dodge(.5), col = 'grey20')+#
  geom_ribbon(data= pred, aes(solar.ele.begdesc, maxdep, ymin = conf.low, ymax = conf.high), 
                position = position_dodge(.5),  alpha = 0.2, size = 0.4, col = 'grey20')+#
  scale_color_manual(values=met.brewer("Derain", 2, direction = 1), name = 'Island')+
  scale_y_reverse()+ coord_cartesian(xlim = c(-7, 52))#+ # coord_cartesian(ylim = c(150, NA)) 
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
                              axis.text.x = element_text(size = 12, colour = "black"), 
                              axis.text.y = element_text(size = 12, colour = "black"), 
                                axis.title = element_text(size = 16))+
                      ggtitle("Gentoo penguins")+
                    xlab("Solar elevation (degrees above the horizon)") +
                      ylab("Maximum depth (m)")
s

# save image 
rstudioapi::savePlotAsImage(paste0("./outputs/phenotypic lmms/dive/plots/", dive$species[1], "_maxdep_population phenotypic lmm_violin2.png"), width=900,height=600)

#-----------------------------------------------------------------------------------------------------------------------------------------
# Summarise data to create histogram counts
histo = dive %>% 
        group_by(maxdep, island, group_f) %>%
        mutate(breaks = cut(solar.ele.begdesc, 
                        breaks=seq(-6, 52, 1),  # range of solar.ele.begdesc
                        labels=seq(-6, 51, 1),  # range of solar.ele.begdesc, but for some reason must be a little smaller range than above
                        include.lowest=TRUE),
         breaks = as.numeric(as.character(breaks))) %>%
  group_by(maxdep, breaks, island, group_f) %>% 
  summarise(n = n()) %>%
  # If you want to change the absolute heights of the bars, just multiply n/sum(n) by a scaling factor 
  mutate(pct = ifelse(maxdep==0, n/sum(n)*0.3, 1 - n/sum(n)*0.3))   # the * 0.3 is to scale how tall the histograms are. 

histo

ggplot() +
  geom_segment(data = histo, 
               aes(x = breaks, xend = breaks, y = maxdep, yend = pct, colour=factor(island)),
               size= 1, show.legend=FALSE) +
   facet_wrap(island~group_f) + scale_y_reverse()+
   theme_bw(base_size=12)


#------------------
s = ggplot(data =pred, aes(solar.ele.begdesc, maxdep, col = island)) +
  geom_line(size = 1, alpha = 1, position = position_dodge(.5)) +   facet_wrap(island~group_f) +
  #stat_smooth(method = "lm", se = FALSE)
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    position = position_dodge(.5), 
    alpha = 0.1#, size = 0.3
  ) + #scale_y_reverse() +
  scale_color_manual(values=met.brewer("Derain", 2, direction = -1), name = 'Island')
s

# try adding observed points to the predicted estimates
#s <- s + geom_point(dive, mapping = aes(x = solar.ele.begdesc, y = type, col = island), 
#                    position = position_jitter(height = .15), shape = '.', alpha = 0.6) +
#  facet_wrap(group_f~island, ncol = 6) +
#  geom_hline(yintercept = 1.0)+
#  geom_hline(yintercept = 0.0)
##coord_cartesian(ylim = c(120, NA)) 
#s  

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
        axis.text.x = element_text(size = 12, colour = "black"), 
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 16))+
  ggtitle("Gentoo penguins")+
  xlab("Solar elevation (degrees above the horizon)") +
  ylab("Maximum depth (m)")
s

s = s +
  geom_segment(data = histo,  
               aes(x = breaks, xend = breaks, y = maxdep, yend = pct),  # remove the colour function here as it is already given in s.
               size= 1, show.legend=FALSE,  # adjust the size parameter to change the "bar" widths in the histogram
               alpha = 0.05
               ) + scale_y_reverse()+
  facet_wrap(island~group_f)

s




# try adding observed points to the predicted estimates
#s <- s + geom_point(dive, mapping = aes(x = solar.ele.begdesc, y = type, col = island), 
#                    position = position_jitter(height = .15), shape = '.', alpha = 0.6) +
#  facet_wrap(group_f~island, ncol = 6) +
#  geom_hline(yintercept = 1.0)+
#  geom_hline(yintercept = 0.0)
##coord_cartesian(ylim = c(120, NA)) 
#s  

# other way around
#---------------------------------------------------------------------------------------------------------------------

s = ggplot()+
  geom_segment(data = histo,  
               aes(x = breaks, xend = breaks, y = maxdep, yend = pct, col = island),  # remove the colour function here as it is already given in s.
               size= 1, show.legend=FALSE,  # adjust the size parameter to change the "bar" widths in the histogram
               alpha = 0.05
               ) + scale_y_reverse()+
  facet_wrap(island~group_f)

s


s = s +
  geom_line(data =pred, aes(solar.ele.begdesc, maxdep), 
            size = 1, alpha = 1, position = position_dodge(.5), col = 'grey20') +   facet_wrap(island~group_f) +
  #stat_smooth(method = "lm", se = FALSE)
  geom_ribbon(data =pred, aes(solar.ele.begdesc, maxdep, ymin = conf.low, ymax = conf.high),
    position = position_dodge(.5), 
    alpha = 0.1, col = 'grey20'#, size = 0.3
  ) + coord_cartesian(xlim = c(-7, 52))+ #scale_y_reverse() +
  scale_color_manual(values=met.brewer("Derain", 2, direction = 1), name = 'Island')
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
        axis.text.x = element_text(size = 12, colour = "black"), 
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 16))+
  ggtitle("Gentoo penguins")+
  xlab("Solar elevation (degrees above the horizon)") +
  ylab("Maximum depth (m)")
s

# save image 
rstudioapi::savePlotAsImage(paste0("./outputs/phenotypic lmms/dive/plots/", dive$species[1], "_maxdep_population phenotypic lmm_histo.png"), width=900,height=600)

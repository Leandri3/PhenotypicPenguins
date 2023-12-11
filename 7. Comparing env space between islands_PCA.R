# Compare environmental spaces of all CS together
library(tidyverse)
library(MetBrewer)
library(ggplot2)
library(ggfortify)

dat1 <- readRDS('./outputs/environmental space covariates/95_UD_CS_Nel_Chick_with_updated_resolution and env variables.rds')

dat1 <- dat1 %>%
  mutate(island = 'Nelson')%>%
  mutate(phase = 'chick')%>%
  mutate(species = "chinstraps")


dat2 <- readRDS('./outputs/environmental space covariates/95_UD_CS_Kop_Chick_with_updated_resolution and env variables.rds')

dat2 <- dat2 %>%
  mutate(island = 'Kopaitic')%>%
  mutate(phase = 'chick')%>%
  mutate(species = "chinstraps")

dat3 <- readRDS('./outputs/environmental space covariates/95_UD_CS_Nel_Inc_with_updated_resolution and env variables.rds')
dat3 <- dat3 %>%
  mutate(island = 'Nelson')%>%
  mutate(phase = 'inc')%>%
  mutate(species = "chinstraps")

dat4 <- readRDS('./outputs/environmental space covariates/95_UD_CS_Kop_Inc_with_updated_resolution and env variables.rds')

dat4 <- dat4 %>%
  mutate(island = 'Kopaitic')%>%
  mutate(phase = 'inc')%>%
  mutate(species = "chinstraps")

dat <- as_tibble(bind_rows(dat3,dat4)) 
head(dat)
dat$depth2 <- (dat$depth$IBCSO_v2_bed*-1)
dat$depth2

# Do a principal component analysis 
# An introduction to statistical learning 
# Book by James et al 2017
# Chapter 10: Unsupervised learning
# Section 10.2 Principal Component analysis

names(dat)
dat = na.omit(dat)

split <- dat %>%
 dplyr::select(depth2, slope, sst,shelf_posneg, u_current, v_current,  sal, island, phase) %>%   #fsle,
  dplyr::filter(phase == 'inc')

split = na.omit(split)
head(split)

divetemps2.isl <- split %>%
  group_by(island) %>%
  dplyr::summarise(n = n())
divetemps2.isl    

# Kop = 2233
# Nel = 36611

split <- split %>%
  dplyr::select(sst, sal, depth2, slope, shelf_posneg, u_current, v_current)  #fsle,

# We first briefly examine the data. We notice that the variables have vastly different means.
apply(split , 2, mean) 
apply(split , 2, var)

# We now perform principal components analysis using the prcomp() function, which is one of several functions in R that perform PCA.
pr.out = prcomp(split, scale =TRUE)

names(pr.out)
summary(pr.out)
# The center and scale components correspond to the means and standard deviations of the variables that were used for 
# scaling prior to implementing PCA.
pr.out$center
pr.out$scale
# The rotation matrix provides the principal component loadings; each column of pr.out$rotation contains the 
# corresponding principal component loading vector.
pr.out$rotation
# We see that there are six distinct principal components. This is to be expected because there are in general 
# min(n ??? 1, p) informative principal components in a data set with n observations and p variables.

# Using the prcomp() function, we do not need to explicitly multiply the data by the principal component loading 
# vectors in order to obtain the principal component score vectors. Rather the 50 ?4 matrix x has as its columns 
# the principal component score vectors. That is, the kth column is the kth principal component score vector.
dim(pr.out$x)
# We can plot the first two principal components as follows:
#biplot(pr.out, scale = T)

#library(devtools)
#install_github("ggbiplot", "vqv")
#pr.comp <- as.data.frame(pr.out)
#ggbiplot(pr.out)  #geom_point(aes(shape = 21), size = 1))

p <-autoplot(pr.out, data = dat, colour = "island", loadings = T, loadings.colour = 'black', 
             loadings.label = T, loadings.label.size = 5, loadings.label.colour = 'black', 
             frame = T ) #frame.type = 't'
print(p)

p+ labs(title = "PCA plot for Chinstraps_Incubation")

#save image: 
rstudioapi::savePlotAsImage(paste0("./outputs/plots/chinstrap/phenotypic environmental space/", dat$species[1], "_", dat$phase[1], "_PCA plot env space2.png"), width=900,height=600)

#---------------------------------------------------------------------------------------
# now for chick rearing:
dat <- as_tibble(bind_rows(dat1,dat2)) 
head(dat)
dat$depth2 <- (dat$depth$IBCSO_v2_bed*-1)
dat$depth2

names(dat)
dat = na.omit(dat)

split <- dat %>%
  dplyr::select(depth2, slope, sst,shelf_posneg, u_current, v_current,  sal, island, phase) %>% #fsle,
  dplyr::filter(phase == 'chick')

split = na.omit(split)
head(split)

divetemps2.isl <- split %>%
  group_by(island) %>%
  dplyr::summarise(n = n())
divetemps2.isl    

# Kop = 697
# Nel = 2525

split <- split %>%
  dplyr::select(sst,sal, depth2, slope, shelf_posneg,  u_current, v_current)   #fsle,

# We first briefly examine the data. We notice that the variables have vastly different means.
apply(split , 2, mean) 
apply(split , 2, var)

# We now perform principal components analysis using the prcomp() function, which is one of several functions in R that perform PCA.
pr.out = prcomp(split, scale =TRUE)

names(pr.out)
summary(pr.out)
# The center and scale components correspond to the means and standard deviations of the variables that were used for 
# scaling prior to implementing PCA.
pr.out$center
pr.out$scale
# The rotation matrix provides the principal component loadings; each column of pr.out$rotation contains the 
# corresponding principal component loading vector.
pr.out$rotation
# We see that there are six distinct principal components. This is to be expected because there are in general 
# min(n ??? 1, p) informative principal components in a data set with n observations and p variables.

# Using the prcomp() function, we do not need to explicitly multiply the data by the principal component loading 
# vectors in order to obtain the principal component score vectors. Rather the 50 ?4 matrix x has as its columns 
# the principal component score vectors. That is, the kth column is the kth principal component score vector.
dim(pr.out$x)
# We can plot the first two principal components as follows:
#biplot(pr.out, scale = T)

#library(devtools)
#install_github("ggbiplot", "vqv")
#pr.comp <- as.data.frame(pr.out)
#ggbiplot(pr.out)  #geom_point(aes(shape = 21), size = 1))

p <-autoplot(pr.out, data = dat, colour = "island", loadings = T, loadings.colour = 'black', 
             loadings.label = T, loadings.label.size = 5, loadings.label.colour = 'black', 
             frame = T) #, frame.type = 't'
print(p)

p+ labs(title = "PCA plot for Chinstraps_Chick-rearing")

#save image: 
rstudioapi::savePlotAsImage(paste0("./outputs/plots/chinstrap/phenotypic environmental space/", dat$species[1], "_", dat$phase[1], "_PCA plot env space2.png"), width=900,height=600)

rm(list=ls())
#----------------------------------------------------------------------------------------------------------
# now do the same for GENTOOS

dat1 <- readRDS('./outputs/environmental space covariates/95_UD_GT_Nel_Chick_with_updated_resolution and env variables.rds')

dat1 <- dat1 %>%
  mutate(island = 'Nelson')%>%
  mutate(phase = 'chick')%>%
  mutate(species = "gentoos")


dat2 <- readRDS('./outputs/environmental space covariates/95_UD_GT_Kop_Chick_with_updated_resolution and env variables.rds')

dat2 <- dat2 %>%
  mutate(island = 'Kopaitic')%>%
  mutate(phase = 'chick')%>%
  mutate(species = "gentoos")

dat3 <- readRDS('./outputs/environmental space covariates/95_UD_GT_Nel_Inc_with_updated_resolution and env variables.rds')
dat3 <- dat3 %>%
  mutate(island = 'Nelson')%>%
  mutate(phase = 'inc')%>%
  mutate(species = "gentoos")

dat4 <- readRDS('./outputs/environmental space covariates/95_UD_GT_Kop_Inc_with_updated_resolution and env variables.rds')

dat4 <- dat4 %>%
  mutate(island = 'Kopaitic')%>%
  mutate(phase = 'inc')%>%
  mutate(species = "gentoos")

dat <- as_tibble(bind_rows(dat3,dat4)) 
head(dat)
dat$depth2 <- (dat$depth$IBCSO_v2_bed*-1)
dat$depth2

# Do a principal component analysis 
# An introduction to statistical learning 
# Book by James et al 2017
# Chapter 10: Unsupervised learning
# Section 10.2 Principal Component analysis

names(dat)
dat = na.omit(dat)

split <- dat %>%
  dplyr::select(depth2, slope, sst,shelf_posneg, u_current, v_current,sal, island, phase) %>%   #fsle, 
  dplyr::filter(phase == 'inc')

split = na.omit(split)
head(split)

divetemps2.isl <- split %>%
  group_by(island) %>%
  dplyr::summarise(n = n())
divetemps2.isl    

# Kop =  1084
# Nel =  1140

split <- split %>%
  dplyr::select(sst, sal, depth2, slope, shelf_posneg, u_current, v_current) #fsle,

# We first briefly examine the data. We notice that the variables have vastly different means.
apply(split , 2, mean) 
apply(split , 2, var)

# We now perform principal components analysis using the prcomp() function, which is one of several functions in R that perform PCA.
pr.out = prcomp(split, scale =TRUE)

names(pr.out)
summary(pr.out)
# The center and scale components correspond to the means and standard deviations of the variables that were used for 
# scaling prior to implementing PCA.
pr.out$center
pr.out$scale
# The rotation matrix provides the principal component loadings; each column of pr.out$rotation contains the 
# corresponding principal component loading vector.
pr.out$rotation
# We see that there are six distinct principal components. This is to be expected because there are in general 
# min(n ??? 1, p) informative principal components in a data set with n observations and p variables.

# Using the prcomp() function, we do not need to explicitly multiply the data by the principal component loading 
# vectors in order to obtain the principal component score vectors. Rather the 50 ?4 matrix x has as its columns 
# the principal component score vectors. That is, the kth column is the kth principal component score vector.
dim(pr.out$x)
# We can plot the first two principal components as follows:
#biplot(pr.out, scale = T)

#library(devtools)
#install_github("ggbiplot", "vqv")
#pr.comp <- as.data.frame(pr.out)
#ggbiplot(pr.out)  #geom_point(aes(shape = 21), size = 1))

p <-autoplot(pr.out, data = dat, colour = "island", loadings = T, loadings.colour = 'black', 
             loadings.label = T, loadings.label.size = 5, loadings.label.colour = 'black', 
             frame = T) # , frame.type = 't'
print(p)

p+ labs(title = "PCA plot for Gentoos_Incubation")

#save image: 
rstudioapi::savePlotAsImage(paste0("./outputs/plots/gentoo/phenotypic environmental space/", dat$species[1], "_", dat$phase[1], "_PCA plot env space2.png"), width=900,height=600)

#split.island <- c(rep("Kopaitic", 202770), rep("Nelson", 2611423) )
#ggplot2::ggbiplot(pr.out, ellipse=TRUE, groups = split.island, alpha = 0.2)
#
#ggbiplot(pr.out, ellipse = TRUE, choices=c(3,4), groups = divetemps2.island)
#ggbiplot(pr.out, ellipse = TRUE, choices=c(5,6), groups = divetemps2.species)

#------------------------------------------------------------------------------------
# now for chick rearing:
dat <- as_tibble(bind_rows(dat1,dat2)) 
head(dat)
dat$depth2 <- (dat$depth$IBCSO_v2_bed*-1)
dat$depth2

names(dat)
dat = na.omit(dat)

split <- dat %>%
  dplyr::select(depth2, slope, sst,shelf_posneg, u_current, v_current,  sal, island, phase) %>% #fsle,
  dplyr::filter(phase == 'chick')

split = na.omit(split)
head(split)

divetemps2.isl <- split %>%
  group_by(island) %>%
  dplyr::summarise(n = n())
divetemps2.isl    

# Kop chick =  1349
# Nel  chick = 1197

split <- split %>%
  dplyr::select(sst, sal, depth2, slope, shelf_posneg,  u_current, v_current)  #fsle,

# We first briefly examine the data. We notice that the variables have vastly different means.
apply(split , 2, mean) 
apply(split , 2, var)

# We now perform principal components analysis using the prcomp() function, which is one of several functions in R that perform PCA.
pr.out = prcomp(split, scale =TRUE)

names(pr.out)
summary(pr.out)
# The center and scale components correspond to the means and standard deviations of the variables that were used for 
# scaling prior to implementing PCA.
pr.out$center
pr.out$scale
# The rotation matrix provides the principal component loadings; each column of pr.out$rotation contains the 
# corresponding principal component loading vector.
pr.out$rotation
# We see that there are six distinct principal components. This is to be expected because there are in general 
# min(n ??? 1, p) informative principal components in a data set with n observations and p variables.

# Using the prcomp() function, we do not need to explicitly multiply the data by the principal component loading 
# vectors in order to obtain the principal component score vectors. Rather the 50 ?4 matrix x has as its columns 
# the principal component score vectors. That is, the kth column is the kth principal component score vector.
dim(pr.out$x)
# We can plot the first two principal components as follows:
#biplot(pr.out, scale = T)

#library(devtools)
#install_github("ggbiplot", "vqv")
#pr.comp <- as.data.frame(pr.out)
#ggbiplot(pr.out)  #geom_point(aes(shape = 21), size = 1))

p <-autoplot(pr.out, data = dat, colour = "island", loadings = T, loadings.colour = 'black', 
             loadings.label = T, loadings.label.size = 5, loadings.label.colour = 'black', 
             frame = T)#,frame.type = 't'
print(p)

p+ labs(title = "PCA plot for Gentoos_Chick-rearing")

#save image: 
rstudioapi::savePlotAsImage(paste0("./outputs/plots/gentoo/phenotypic environmental space/", dat$species[1], "_", dat$phase[1], "_PCA plot env space2.png"), width=900,height=600)

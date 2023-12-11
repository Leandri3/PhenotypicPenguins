# PhenotypicPenguins
This repository contains examples of the important R scripts that I used in my data analysis for my MSc dissertation, Chapter 2: Phenotypic plasticity of foraging behaviours in chinstrap and gentoo penguins in the Antarctic Peninsula. I did not include scripts where I plotted the figures. 

1) GPS tracks to trips
   -Divide GPS track into separate trips to sea based on diving data
   -Calculate maximum and cumulative trip distances and trip durations based on GPS
2) DiveMove Analysis
   - Zero-offset correction
   - Calculate dive statistics for all individuals in a deployment round
  
3) CrawlWrap for dive locations
  - speed filter GPS data (remove unreliable locations) with trip::sda
  - crawlWrap dive locations of multiple animals with a function (momentuHMM)

4) Utilisation distributions _ 'AdehabitatHR'
   - Creates species-specific utilisations distribtuions (UDs) for each island
   - In order to characterize an environmental space
   - Calculates size of area used by species at both islands
   - Plots the UDs against bathy background

5) EM classification of 5m dives
   - Estimates dive residuals (using linear mixed-effects models)
   - Applies an Expectation-Maximization algorithm to Classify diving types into 2 clusters (foraging and non-foraging) per species, per island
 Using the following predictors: Bottom time
                                 Dive residuals 
                                 Maximum depth
   - Separates into three different Clusters: 1 = Foraging

6) Extract solar elevation for dive locations
   - Calculate solar elevation for every foraging dive location
  
7) Extract environmental covariates for 95% UDs
   - Imports remote-sensed environmental variables from various sources
   - Then import 95% species-specific UDs (VUD rasters are first converted to polygons from which information will be extracted)
   - Extracts remote-sensed environmental information from 95% species-specific UD polygons and for each dive location
  
8) Comparing environmental space between islands using PCA
 - Uses Principal Component Analysis to compare environmental spaces of the two islands

9) Trip behaviours linear-mixed effects models
   For Trip behaviours e.g. maximum distance, trip durations
   - Model selection using lme4 and AIC
   - Model checking
   - Repeatability estimation
   - Predicts effects from best model
   - Plots figures of model output

10) Dive behaviours with generalised linear-mixed effects models
    For dive behaviours e.g. maximum dive depth and foraging dive type
   -Do model selection using - glmmTMB
   - Includes autocorrelation structures
   - family = Gamma(link = "log"),
   - repeatability estimation
   - Predicts effects from best model
   - Plots figures of model output
  
   

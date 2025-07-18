#This project looks at the effects of trophic interactions of foundation species on dependent communities
#this script is to load packages needed for this project

#This code says if the package needs to be installed, install it, if not, load from library
if(!require(gsheet))install.packages("gsheet");library(gsheet)
if(!require(tidyverse))install.packages("tidyverse");library(tidyverse)
if(!require(vegan))install.packages("vegan");library(vegan)
if(!require(ggrepel))install.packages("ggrepel");library(ggrepel)
if(!require(labdsv))install.packages("labdsv");library(labdsv)
if(!require(patchwork))install.packages("patchwork");library(patchwork)
if(!require(ggplot2))install.packages("ggplot2");library(ggplot2)
if(!require(lubridate))install.packages("lubridate");library(lubridate)
if(!require(viridis))install.packages("viridis");library(viridis)
if(!require(GpGp))install.packages("GpGp");library(GpGp)
if(!require(openxlsx))install.packages("openxlsx");library(openxlsx)
if(!require(viridis))install.packages("viridis");library(viridis)
if(!require(glmmTMB))install.packages("glmmTMB");library(glmmTMB)
if(!require(DHARMa))install.packages("DHARMa");library(DHARMa)
if(!require(car))install.packages("car");library(car)

#This project looks at the effects of trophic interactions of foundation species on dependent communities
#this is is to load packages needed for this project
#
install.packages("viridis")
#These are for uploading the data from google sheets
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



#load packages script
lp<-function(pck){
  if(!require(pck,character.only = TRUE))install.packages(pck);library(pck,character.only = TRUE)
}

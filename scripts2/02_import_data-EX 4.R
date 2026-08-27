#This project looks at the effects of trophic interactions of foundation species on dependent communities
#this script is to import original data from google
#We assigned each google sheet link to url and then imported the data into R.

#####Community Data-------
###Biomass
url<-"https://docs.google.com/spreadsheets/d/10NUI80qSWpe6QwEbiq7zBRt5cDi67r5Z/edit?gid=1086801212#gid=1086801212"

biomass<-read.csv(text=gsheet2text(url,format='csv'),stringsAsFactors = FALSE)  

###Amphipod and Isopod subsampling
url<-"https://docs.google.com/spreadsheets/d/10NUI80qSWpe6QwEbiq7zBRt5cDi67r5Z/edit?gid=1930529029#gid=1930529029"

allaisub<-read.csv(text=gsheet2text(url,format='csv'),stringsAsFactors = FALSE)  

###Body Size
url<-"https://docs.google.com/spreadsheets/d/10NUI80qSWpe6QwEbiq7zBRt5cDi67r5Z/edit?gid=328129538#gid=328129538"

bodysize<-read.csv(text=gsheet2text(url,format='csv'),stringsAsFactors = FALSE)

#Taxa ID
url<-"https://docs.google.com/spreadsheets/d/1UDc0NGUatjNf8FFfXCRDwhPHyDxvcg5J/edit?gid=246528797#gid=246528797"

taxaid<-read.csv(text=gsheet2text(url,format='csv'),stringsAsFactors = FALSE)

#####Oyster reef data-----
#shell height
url<-"https://docs.google.com/spreadsheets/d/1GNmZZDJD0dkYHPJuaQxv-KaGDl5QqoVJ/edit?gid=899021306#gid=899021306"

shellheight<-read.csv(text=gsheet2text(url,format='csv'),stringsAsFactors = FALSE)

#reef characteristics from quadrat data
url<-"https://docs.google.com/spreadsheets/d/1GNmZZDJD0dkYHPJuaQxv-KaGDl5QqoVJ/edit?gid=2110957063#gid=2110957063"

reefquads<-read.csv(text=gsheet2text(url,format='csv'),stringsAsFactors = FALSE)

#Sandra paired shell height and volume 
url<-"https://docs.google.com/spreadsheets/d/1rc124Wh4t2nmqy_cJ1I0OuZ3xl11JUNx/edit?gid=584157151#gid=584157151"

paired.sh.v<-read.csv(text=gsheet2text(url,format='csv'),stringsAsFactors = FALSE)

##### real taxa id's
url<-"https://docs.google.com/spreadsheets/d/1MsEpvFw1IqiymmitmpGCQEM3CACVttO6/edit?gid=1179566381#gid=1179566381"

real.id<-read.csv(text=gsheet2text(url,format='csv'),stringsAsFactors = FALSE)

#functional feeding groups
url<- "https://docs.google.com/spreadsheets/d/1vDwu5Ru_1HSK4Gw5T_27yKZuEAKogwp-lCTi5d6y3n8/edit?usp=sharing"

ffg<-read.csv(text=gsheet2text(url,format='csv'),stringsAsFactors = FALSE)


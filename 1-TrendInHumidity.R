suppressMessages(library(fields))
suppressMessages(library(ncdf4))
suppressMessages(library(RColorBrewer))
suppressMessages(library(abind))
suppressMessages(library(base))
suppressMessages(library(Rfit))
#This reads in specific humidity (shum) values calculated outside of this project 
#where the each grid cell for each date is either the humidity measurement or NA,
#with NA values on dates that were not determined to be 'heat waves'. Thus the majority
#of the shum.hw file below is NA's, since heat waves occur relatively rarely. I do not
#include the scripts required to calculate heat wave days and extract the humidity values
#because it is quite long and takes several days to run on a single Sherlock node.
#File below is too large (150mb) to include in the zipfile; Therefore you will not 
#be able to compute anything below. The output is saved to .rds files and read in in the
#other scripts, however, so those should run on a personal computer.
shum.hw=readRDS("/scratch/users/tballard/shum/percentile.threshold/shum.hw")

## Define month indices ##
  jan=1:124; feb=125:236; mar=237:360; apr=361:480;may=481:604; jun=605:724;
  jul=725:848; aug=849:972; sep=973:1092; oct=1093:1216; nov=1217:1336; dec=1337:1460;

## Extract the shum values for the month of interest. If I were better at matrix operations
## I might be able to avoid this, but running regression on 192x94x1460x36 is tricky when
## you want 3 dimensions instead
month.extract=function(month){
    NA.matrix=matrix(rep(NA,192*94),nrow=192) #used to initialize
    shum.all.values=array(NA.matrix,c(192,94,1)) #192x94x1 array of NA's to initialize
    for (i in 1:36){	
      shum.values=shum.hw[,,month,i] 
      shum.all.values=abind(shum.all.values,shum.values)
    }	
    shum.all.values=shum.all.values[,,-1] #remove the NA.matrix used to initialize
    return(shum.all.values)
  }
  
  daily.jan.values=month.extract(jan) #192lat x 94lon x (31day*4val/day*36yr=4464time)
  daily.jul.values=month.extract(jul)

##### Read in ocean/land mask file 1's are land #####
  fileName="/scratch/PI/omramom/reanalysis/ncep-doe-r2/4xdaily/ocean.mask/land.sfc.gauss.nc"
  land = ncvar_get(nc_open(fileName), "land") #192 x 94
  land[land==0]=NA #set ocean values to NA instead of 0
  mask=function(data,land){
    new.data=land*data
    return(new.data)
  }
  daily.jan.values.mask=apply(daily.jan.values, c(3), mask, land=land) #dim=18048x4464
  daily.jan.values.mask=array(daily.jan.values.mask, dim=c(192,94,4464))
  daily.jul.values.mask=apply(daily.jul.values, c(3), mask, land=land)
  daily.jul.values.mask=array(daily.jul.values.mask,dim=c(192,94,4464))
    
##### Compute them OLS trends! #####
  #Function below computes OLS regression of shum vs year and outputs the slope value and p-value
  #Then apply this function to every pixel using 'apply' command
  #Note the lm function automatically skips over NA's
  fit.lm=function(dataset, month){
    years=rep(c(1979:2014),each=length(month)) #1979, 1979, ... 2014, 2014 
    
    a=tryCatch(summary(lm(dataset~years))$coefficient[2,c(1,4)], error=function(e) c(NA,NA)) #slope for 'year' and p-value
    return(a)
  }
 
 
  
  #Array of lon,lat,month,results; results is 2D of the slope and its SE from running the regression
  #'aperm' rearranges order of arrays. the aperm below switches the apply output from dim=2,192,94 to dim=192,94,2
  lm.trends=array(rep(NA,192*94*2*2),c(192,94,2,2)) #initialize
  lm.trends[,,1,]=aperm(apply(daily.jan.values.mask, c(1,2), fit.lm, month=jan), c(2,3,1))
  lm.trends[,,2,]=aperm(apply(daily.jul.values.mask, c(1,2), fit.lm, month=jul), c(2,3,1))

  
  saveRDS(lm.trends,"/scratch/users/tballard/shum/class.project/shum.hw.trend.rds") #194x94x2monthsx2variables
  
  
  
##### Compute rank based regression trends #####
   rfit.lm=function(dataset, month){
    years=rep(c(1979:2014),each=length(month)) #1979, 1979, ... 2014, 2014 
    
    a=tryCatch(summary(rfit(dataset~years))$coefficient[2,c(1,4)], error=function(e) c(NA,NA)) #slope for 'year' and p-value
    return(a)
  }
 
 
  
  #Array of lon,lat,month,results; results is 2D of the slope and its SE from running the regression
  #'aperm' rearranges order of arrays. the aperm below switches the apply output from dim=2,192,94 to dim=192,94,2
  rank.trends=array(rep(NA,192*94*2*2),c(192,94,2,2)) #initialize
  rank.trends[,,1,]=aperm(apply(daily.jan.values.mask, c(1,2), rfit.lm, month=jan), c(2,3,1))
  rank.trends[,,2,]=aperm(apply(daily.jul.values.mask, c(1,2), rfit.lm, month=jul), c(2,3,1))

  saveRDS(rank.trends,"/scratch/users/tballard/shum/class.project/shum.hw.rank.trend.rds") #194x94x2monthsx2variables
  
  
  
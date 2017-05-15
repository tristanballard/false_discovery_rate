suppressMessages(library(fields))
suppressMessages(library(ncdf4))
suppressMessages(library(RColorBrewer))
suppressMessages(library(abind))
suppressMessages(library(base))

#This reads in the trends calculated in TrendInHumidity.R and then applies the 
#multiple comparison correction methods. It outputs an .rds file with adjusted p-values
#for each of the tests on the 192x94 spatial grid for both January and July.
lm.trend=readRDS("/scratch/users/tballard/shum/class.project/shum.hw.trend.rds")
#dim = 194 lat x 92 lon x 2 months x 2 variables (trend and p-value)

source("/scratch/users/tballard/shum/class.project/jcli-3199-fdr.R")
#loads the fdr function coded by Ventura for the function in Ventura et al. (2004)

  trend.jan=lm.trend[,,1,1] #192x94 lat and lon
  trend.jul=lm.trend[,,2,1] #192x94
  pval.jan=lm.trend[,,1,2] #192x94
  pval.jul=lm.trend[,,2,2] #192x94
  rm(lm.trend)


##### Perform Multiple Comparisons Corrections #####

#Benjamini & Hochberg (1995)
  pval.jan2=p.adjust(pval.jan, method="BH") #vector of length 192*94=18048
  pval.jan2=array(pval.jan2, dim=c(192,94)) #put back into 2dim array
  pval.jul2=p.adjust(pval.jul, method="BH") 
  pval.jul2=array(pval.jul2, dim=c(192,94)) 

#Benjamini & Yekutieli (2001)
  pval.jan3=p.adjust(pval.jan, method="BY") #vector of length 192*94=18048
  pval.jan3=array(pval.jan3, dim=c(192,94)) #put back into 2dim array
  pval.jul3=p.adjust(pval.jul, method="BY") 
  pval.jul3=array(pval.jul3, dim=c(192,94)) 

#Ventura et al. (2004) with q=.05
  ##Unlike p.adjust, the fdr() function does not take inputs with NA's, so first remove them
  pval.jan4.i=as.vector(pval.jan)
  pval.jan4.index=!is.na(pval.jan4.i) #True if there is a p-val, False if it's NA; may need later
  pval.jan4.ii=pval.jan4.i[!is.na(pval.jan4.i)] #length=5914 instead of original 18048
  ##Now apply the fdr() function. method='original' follows the BH (1995) method, adjustment.method="mean" 
  ##performs the adjustment Ventura et al. suggest to that method if using spatial data.
  ##Output is a vector of the indices of the significant tests
  pval.jan4.fdr=fdr(pval.jan4.ii, qlevel=0.05, method="original", adjustment.method="mean") #length=2436
  length(pval.jan4.fdr) #length=2436
  #Now set the indices not present in that output to NA
  pval.jan4.ii.ind=rep(0,length(pval.jan4.ii))
  pval.jan4.ii.ind[pval.jan4.fdr]=1 #Vector of 0's and 1's
  pval.jan4.ii.adj=pval.jan4.ii.ind*pval.jan4.ii #vector length 5914 of original p-vals and 0's
  
  pval.jan4.i[!is.na(pval.jan4.i)]=pval.jan4.ii.adj #Take original vector, replace the non-NA list of p-values with this new list
  pval.jan4.i[pval.jan4.i==0]=.99 #Set the 0's added in to pval's of .99 instead; in Plot.MCC.R all pval's >=.05 are set to NA
  pval.jan4=pval.jan4.i
  pval.jan4=array(pval.jan4, dim=c(192,94))
  
  ##Repeat for July
  pval.jul4.i=as.vector(pval.jul)
  pval.jul4.index=!is.na(pval.jul4.i) #True if there is a p-val, False if it's NA; may need later
  pval.jul4.ii=pval.jul4.i[!is.na(pval.jul4.i)]
  pval.jul4.fdr=fdr(pval.jul4.ii, qlevel=0.05, method="original", adjustment.method="mean") #length=3045
  length(pval.jul4.fdr) #length=3045
  #Now set the indices not present in that output to NA
  pval.jul4.ii.ind=rep(0,length(pval.jul4.ii))
  pval.jul4.ii.ind[pval.jul4.fdr]=1 #Vector of 0's and 1's
  pval.jul4.ii.adj=pval.jul4.ii.ind*pval.jul4.ii #vector length 4886 of original p-vals and 0's
  
  pval.jul4.i[!is.na(pval.jul4.i)]=pval.jul4.ii.adj #Take original vector, replace the non-NA list of p-values with this new list
  pval.jul4.i[pval.jul4.i==0]=.99 #Set the 0's added in to pval's of .99 instead; in Plot.MCC.R all pval's >=.05 are set to NA
  pval.jul4=pval.jul4.i
  pval.jul4=array(pval.jul4, dim=c(192,94))
  

#Ventura et al. (2004) with q=.01
  ##Unlike p.adjust, the fdr() function does not take inputs with NA's, so first remove them
  pval.jan5.i=as.vector(pval.jan)
  pval.jan5.index=!is.na(pval.jan5.i) #True if there is a p-val, False if it's NA; may need later
  pval.jan5.ii=pval.jan5.i[!is.na(pval.jan5.i)] #length=5914 instead of original 18048
  ##Now apply the fdr() function. method='original' follows the BH (1995) method, adjustment.method="mean" 
  ##performs the adjustment Ventura et al. suggest to that method if using spatial data.
  ##Output is a vector of the indices of the significant tests
  pval.jan5.fdr=fdr(pval.jan5.ii, qlevel=0.01, method="original", adjustment.method="mean") #length=2053
  length(pval.jan5.fdr) #length=2053
  #Now set the indices not present in that output to NA
  pval.jan5.ii.ind=rep(0,length(pval.jan5.ii))
  pval.jan5.ii.ind[pval.jan5.fdr]=1 #Vector of 0's and 1's
  pval.jan5.ii.adj=pval.jan5.ii.ind*pval.jan5.ii #vector length 5914 of original p-vals and 0's
  
  pval.jan5.i[!is.na(pval.jan5.i)]=pval.jan5.ii.adj #Take original vector, replace the non-NA list of p-values with this new list
  pval.jan5.i[pval.jan5.i==0]=.99 #Set the 0's added in to pval's of .99 instead; in Plot.MCC.R all pval's >=.05 are set to NA
  pval.jan5=pval.jan5.i
  pval.jan5=array(pval.jan5, dim=c(192,94))
  
  ##Repeat for July
  pval.jul5.i=as.vector(pval.jul)
  pval.jul5.index=!is.na(pval.jul5.i) #True if there is a p-val, False if it's NA; may need later
  pval.jul5.ii=pval.jul5.i[!is.na(pval.jul5.i)]
  pval.jul5.fdr=fdr(pval.jul5.ii, qlevel=0.01, method="original", adjustment.method="mean") #length=2767
  length(pval.jul5.fdr) #length=2767
  #Now set the indices not present in that output to NA
  pval.jul5.ii.ind=rep(0,length(pval.jul5.ii))
  pval.jul5.ii.ind[pval.jul5.fdr]=1 #Vector of 0's and 1's
  pval.jul5.ii.adj=pval.jul5.ii.ind*pval.jul5.ii #vector length 4886 of original p-vals and 0's
  
  pval.jul5.i[!is.na(pval.jul5.i)]=pval.jul5.ii.adj #Take original vector, replace the non-NA list of p-values with this new list
  pval.jul5.i[pval.jul5.i==0]=.99 #Set the 0's added in to pval's of .99 instead; in Plot.MCC.R all pval's >=.05 are set to NA
  pval.jul5=pval.jul5.i
  pval.jul5=array(pval.jul5, dim=c(192,94))
  

#Bonferroni correction for kicks 
  pval.jan6=p.adjust(pval.jan, method="bonferroni") #vector of length 192*94=18048
  pval.jan6=array(pval.jan6, dim=c(192,94)) #put back into 2dim array
  pval.jul6=p.adjust(pval.jul, method="bonferroni") 
  pval.jul6=array(pval.jul6, dim=c(192,94)) 
  

  
##### Save adjusted p-values to combined array #####
  pvals=aperm(abind(pval.jan, pval.jul, pval.jan2, pval.jul2, pval.jan3, pval.jul3, pval.jan4, pval.jul4, pval.jan5, pval.jul5, pval.jan6, pval.jul6, along=.5), c(2,3,1))
  #dim(pvals)= 192 x 94 x 12; aperm above just reorders dimensions the way I prefer
  saveRDS(pvals,"/scratch/users/tballard/shum/class.project/pvals.rds")
 
  
  
  
  
##### Some summary stats on the correction methods results #####
  sum(pval.jan<.05, na.rm=T) #2438 significant pval's originally
  sum(pval.jan2<.05, na.rm=T)/sum(pval.jan<.05, na.rm=T) #75% remain significant
  sum(pval.jan3<.05, na.rm=T)/sum(pval.jan<.05, na.rm=T) #37%
  sum(pval.jan4<.05, na.rm=T)/sum(pval.jan<.05, na.rm=T) #99.9%
  sum(pval.jan5<.05, na.rm=T)/sum(pval.jan<.05, na.rm=T) #59%
  sum(pval.jan6<.05, na.rm=T)/sum(pval.jan<.05, na.rm=T) #14%
  
  
  sum(pval.jul<.05, na.rm=T) #2724 significant pval's originally
  sum(pval.jul2<.05, na.rm=T)/sum(pval.jul<.05, na.rm=T) #90% remain
  sum(pval.jul3<.05, na.rm=T)/sum(pval.jul<.05, na.rm=T) #64%
  sum(pval.jul4<.05, na.rm=T)/sum(pval.jul<.05, na.rm=T) #100%
  sum(pval.jul5<.05, na.rm=T)/sum(pval.jul<.05, na.rm=T) #83%
  sum(pval.jul6<.05, na.rm=T)/sum(pval.jul<.05, na.rm=T) #36%

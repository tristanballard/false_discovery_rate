suppressMessages(library(fields))
suppressMessages(library(ncdf4))
suppressMessages(library(RColorBrewer))
suppressMessages(library(abind))
suppressMessages(library(base))
### Only difference is read in the p-values and trend values from the rank regression insetad of OLS
pvals=readRDS("/scratch/users/tballard/shum/class.project/pvals.rank.rds") #192x94x8
lm.trend=readRDS("/scratch/users/tballard/shum/class.project/shum.hw.rank.trend.rds")
  trend.jan=lm.trend[,,1,1] #192x94
  trend.jul=lm.trend[,,2,1] #192x94
  rm(lm.trend)
  
##### Set non-significant trend values to NA #####
  
### Original LM data ###
  trend.jan.mod1=trend.jan
  trend.jan.mod1[pvals[,,1]>=.05]=NA
  trend.jul.mod1=trend.jul
  trend.jul.mod1[pvals[,,2]>=.05]=NA
  

### Benjamini & Hochberg 1995 ###
  trend.jan.mod2=trend.jan
  trend.jan.mod2[pvals[,,3]>=.05]=NA
  trend.jul.mod2=trend.jul
  trend.jul.mod2[pvals[,,4]>=.05]=NA
  
  
### Benjamini & Yekutieli (2001) ###
  trend.jan.mod3=trend.jan
  trend.jan.mod3[pvals[,,5]>=.05]=NA
  trend.jul.mod3=trend.jul
  trend.jul.mod3[pvals[,,6]>=.05]=NA
  
### Ventura et al. (2004) q=.05 ###
  trend.jan.mod4=trend.jan
  trend.jan.mod4[pvals[,,7]>=.05]=NA
  trend.jul.mod4=trend.jul
  trend.jul.mod4[pvals[,,8]>=.05]=NA
  
### Ventura et al. (2004) q=.01 ###
  trend.jan.mod5=trend.jan
  trend.jan.mod5[pvals[,,9]>=.05]=NA
  trend.jul.mod5=trend.jul
  trend.jul.mod5[pvals[,,10]>=.05]=NA
  
### Bonferroni ###
  trend.jan.mod6=trend.jan
  trend.jan.mod6[pvals[,,11]>=.05]=NA
  trend.jul.mod6=trend.jul
  trend.jul.mod6[pvals[,,12]>=.05]=NA
  


  
#############################################################
###########           Plotting Time           ###############
#############################################################

### Read in (cleaned up a bit for the specific Rpackage settings) lat and lon coordinates 
lat=readRDS("/scratch/users/tballard/plots/global.attributes/shum/lat") 
lon=readRDS("/scratch/users/tballard/plots/global.attributes/shum/lon")


dir="/scratch/users/tballard/shum/class.project/" #where to save plots
    plot.name=paste(dir,"plot.MCC.rank.png",sep="")
    plot.width=6000 #units are pixels
    plot.height=1500
    plot.res=200*.9
    colorbar=rev(brewer.pal(11,"RdYlBu"))
    colorbar[6]="purple"
    legend.lim=c(-.35,.35)
  
##### Make the plots #####
png(plot.name, units="px", width=plot.width, height=plot.height, res=plot.res)
par(mfrow=c(2,6)) #2 x 6 plot

quilt.plot(lon, lat, as.vector(trend.jan.mod1),nx = 192, ny = 94,zlim=legend.lim, las=1,
             #main="Median Number of Heat Events per Year",col=colorbar)
             main="Original Trend in Humidity (g/kg/yr)",col=colorbar,bty='n', axes=F, add.legend=F, cex.main=2)
  par(font=2); legend(x=-210,y=15, "January", bty='n', cex=2)
  par(fg=NA); 
  image.plot(legend.only=T, zlim=legend.lim, col=colorbar, horizontal=T, bty='n'
             ,legend.shrink = .7)
  world(add = TRUE, col = "black",lwd=.7)
  
  
quilt.plot(lon, lat, as.vector(trend.jan.mod2),nx = 192, ny = 94,zlim=legend.lim, las=1,
             #main="Median Number of Heat Events per Year",col=colorbar)
             main="Benjamini & Hochberg Correction",col=colorbar,bty='n', axes=F, add.legend=F, cex.main=2)
  par(fg=NA); 
  image.plot(legend.only=T, zlim=legend.lim, col=colorbar, horizontal=T, bty='n'
             ,legend.shrink = .7)
  world(add = TRUE, col = "black",lwd=.7)
  
quilt.plot(lon, lat, as.vector(trend.jan.mod3),nx = 192, ny = 94,zlim=legend.lim, las=1,
             #main="Median Number of Heat Events per Year",col=colorbar)
             main="Benjamini & Yekutieli Correction",col=colorbar,bty='n', axes=F, add.legend=F, cex.main=2)
  par(fg=NA); 
  image.plot(legend.only=T, zlim=legend.lim, col=colorbar, horizontal=T, bty='n'
             ,legend.shrink = .7)
  world(add = TRUE, col = "black",lwd=.7)
  
quilt.plot(lon, lat, as.vector(trend.jan.mod4),nx = 192, ny = 94,zlim=legend.lim, las=1,
             #main="Median Number of Heat Events per Year",col=colorbar)
             main="Ventura Correction q=.05",col=colorbar,bty='n', axes=F, add.legend=F, cex.main=2)
  par(fg="black") #reset par parameter otherwise the next line doesn't work
  par(font=2); legend(x=-210, y=15, "January", bty='n', cex=2)
  par(fg=NA); 
  image.plot(legend.only=T, zlim=legend.lim, col=colorbar, horizontal=T, bty='n'
             ,legend.shrink = .7)
  world(add = TRUE, col = "black",lwd=.7)
  
quilt.plot(lon, lat, as.vector(trend.jan.mod5),nx = 192, ny = 94,zlim=legend.lim, las=1,
             #main="Median Number of Heat Events per Year",col=colorbar)
             main="Ventura Correction q=.01",col=colorbar,bty='n', axes=F, add.legend=F, cex.main=2)
  par(fg=NA); 
  image.plot(legend.only=T, zlim=legend.lim, col=colorbar, horizontal=T, bty='n'
             ,legend.shrink = .7)
  world(add = TRUE, col = "black",lwd=.7)
  
quilt.plot(lon, lat, as.vector(trend.jan.mod6),nx = 192, ny = 94,zlim=legend.lim, las=1,
             #main="Median Number of Heat Events per Year",col=colorbar)
             main="Bonferroni Correction",col=colorbar,bty='n', axes=F, add.legend=F, cex.main=2)
  par(fg=NA); 
  image.plot(legend.only=T, zlim=legend.lim, col=colorbar, horizontal=T, bty='n'
             ,legend.shrink = .7)
  world(add = TRUE, col = "black",lwd=.7)
  
quilt.plot(lon, lat, as.vector(trend.jul.mod1),nx = 192, ny = 94,zlim=legend.lim, las=1,
             #main="Median Number of Heat Events per Year",col=colorbar)
             main="",col=colorbar,bty='n', axes=F, add.legend=F, cex.main=2)
  par(fg="black") #reset par parameter otherwise the next line doesn't work
  par(font=2); legend(x=-210, y=15, "July", bty='n', cex=2)
  par(fg=NA); 
  image.plot(legend.only=T, zlim=legend.lim, col=colorbar, horizontal=T, bty='n'
             ,legend.shrink = .7)
  world(add = TRUE, col = "black",lwd=.7)
  
  
quilt.plot(lon, lat, as.vector(trend.jul.mod2),nx = 192, ny = 94,zlim=legend.lim, las=1,
             #main="Median Number of Heat Events per Year",col=colorbar)
             main="",col=colorbar,bty='n', axes=F, add.legend=F, cex.main=2)
  par(fg=NA); 
  image.plot(legend.only=T, zlim=legend.lim, col=colorbar, horizontal=T, bty='n'
             ,legend.shrink = .7)
  world(add = TRUE, col = "black",lwd=.7)


quilt.plot(lon, lat, as.vector(trend.jul.mod3),nx = 192, ny = 94,zlim=legend.lim, las=1,
             #main="Median Number of Heat Events per Year",col=colorbar)
             main="",col=colorbar,bty='n', axes=F, add.legend=F, cex.main=2)
  par(fg=NA); 
  image.plot(legend.only=T, zlim=legend.lim, col=colorbar, horizontal=T, bty='n'
             ,legend.shrink = .7)
  world(add = TRUE, col = "black",lwd=.7)
  
quilt.plot(lon, lat, as.vector(trend.jul.mod4),nx = 192, ny = 94,zlim=legend.lim, las=1,
             #main="Median Number of Heat Events per Year",col=colorbar)
             main="",col=colorbar,bty='n', axes=F, add.legend=F, cex.main=2)
  par(fg="black") #reset par parameter otherwise the next line doesn't work
  par(font=2); legend(x=-210, y=15, "July", bty='n', cex=2)
  par(fg=NA); 
  image.plot(legend.only=T, zlim=legend.lim, col=colorbar, horizontal=T, bty='n'
             ,legend.shrink = .7)
  world(add = TRUE, col = "black",lwd=.7)
  
quilt.plot(lon, lat, as.vector(trend.jul.mod5),nx = 192, ny = 94,zlim=legend.lim, las=1,
             #main="Median Number of Heat Events per Year",col=colorbar)
             main="",col=colorbar,bty='n', axes=F, add.legend=F, cex.main=2)
  par(fg=NA); 
  image.plot(legend.only=T, zlim=legend.lim, col=colorbar, horizontal=T, bty='n'
             ,legend.shrink = .7)
  world(add = TRUE, col = "black",lwd=.7)
  
quilt.plot(lon, lat, as.vector(trend.jul.mod6),nx = 192, ny = 94,zlim=legend.lim, las=1,
             #main="Median Number of Heat Events per Year",col=colorbar)
             main="",col=colorbar,bty='n', axes=F, add.legend=F, cex.main=2)
  par(fg=NA); 
  image.plot(legend.only=T, zlim=legend.lim, col=colorbar, horizontal=T, bty='n'
             ,legend.shrink = .7)
  world(add = TRUE, col = "black",lwd=.7)

dev.off()
  

  
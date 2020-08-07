#################################################
## Paper Figure 1 ###############################
#################################################

library(plyr)
library(dplyr)
library(zoo)
library(TeachingDemos)
library(nest)
library(ncdf4)

#################################################

load("Data/Timeseries.RData")

#################################################

# Age limits
xlimz<-c(-1200,100)

# Number of timeseries in the plot
num_ts <- 6

# Set layout:

# Set standard colors for the timeseries
col_hadcm3_temp <- "darkblue"
col_pages2k <- "maroon"#"lawngreen"

col_hadcm3_d18O <- "#59a900"
col_d18O_240 <- "#1f78b4"
col_d18O_242 <- "#a6cee3"

col_co2 <- "firebrick3"
col_volc <- "gray50"
col_solar <- "orange"


# Settings for size of axes, units, and names
axisnumscex <- 1.1 #0.75
axslinno <- 1.1 #0.75
unitscex <- 1.1 #0.75
unitslinno <- 5 #2.5
namlin <- 1 #0.5
namcex <- 1.1 #1


## 4) Step-by-step plotting
for(plot in 1:1){
  cairo_pdf(width=8,height=10,file="Paper_Plots/Fig1_Timeseries.pdf")
  # Raster graphic:
  #png(width=12,height=8,units = "cm",file=paste(outputdir,"/Fig1_talk.png",sep=""),res=150)
  
  # Set margins, out margin, axes styles, standard text size, and line width
  par(mar=c(0,0,0,0),oma=c(0,16,6,1),xaxs="i",yaxs="i",cex=1,lwd=2)
  
  # The first row denotes the header (e.g. names of geologic periods), the remaining rows correspond to the individual timeseries plots
  layout(matrix(rep(1:(num_ts+1),each=2),num_ts+1,2,byrow=TRUE),heights=c(0.2,1,1,1))
  
  #HEADER
  # Start plot
  plot(xlimz,c(0,1),axes=FALSE,type="n",xlab="",ylab="")
  # Lines to mark LIA and MCA
  arrows(x0=-100,x1=-400,y0=0.01,y1=0.01,code=3,angle=90,length=0.025,col="black")
  arrows(x0=-700,x1=-1000,y0=0.01,y1=0.01,code=3,angle=90,length=0.025,col="black")
  # Text for MIS3, Holocene, and LGM
  text(mean(c(-100,-400)),0.5,"LIA",col="black",cex=namcex)
  text(mean(c(-700,-1000)),0.5,"MCA",col="black",cex=namcex)
  # Add top axis
  mtext(side=3,line=3,"time [y CE]",cex=unitscex)
  tmp.y <- TeachingDemos::cnvrt.coords(x=NA,y=par('usr')[4])$dev$y
  axis(3,at=seq(-1050,50,by=200),labels=FALSE)
  mtext(side=3,at=seq(-1050,50,by=200),seq(900,1900,by=200),line=0.5,cex=axisnumscex)
  
  # 1: Temperature
  # y-axis extent:
  range_temp <- c(-1.5,1.0)
  # Start plot
  plot(xlimz,range_temp,axes=FALSE,type="n",xlab="",ylab="")
  # Add full and filtered timeseries
  #polygon(c(rev(Timeseries$pages2k$time[850:2000]), Timeseries$pages2k$time[850:2000]), 
  #        c(rev(Timeseries$pages2k$q_2.5th[850:2000]), Timeseries$pages2k$q_97.5th[850:2000])+as.numeric(Timeseries$HadCM3_TAM[index(Timeseries$HadCM3_TAM) == 0]), 
  #        col = adjustcolor(col_pages2k,0.3), border = NA)
  
  lines(window(Timeseries$HadCM3_GMST_a, start = -1140, end = -100),type="l",col=adjustcolor(col_hadcm3_temp,0.3))
  lines(gaussdetr(window(Timeseries$HadCM3_GMST_a, start = -1140, end = -100),tsc.in=100)$Xsmooth,type="l",col=col_hadcm3_temp)
  
  lines(Timeseries$pages2k$time[850:2000], 
        Timeseries$pages2k$value[850:2000], 
        type = "l", col = adjustcolor(col_pages2k,0.3))
  lines(gaussdetr(zoo(x = Timeseries$pages2k$value[850:2000],
                      order.by = Timeseries$pages2k$time[850:2000]), tsc.in = 100)$Xsmooth,type = "l", col = col_pages2k)
  
  lines(Timeseries$HadCRUT4, type = "l", col = adjustcolor("black", 0.3))
  lines(gaussdetr(Timeseries$HadCRUT4, tsc.in = 100)$Xsmooth, type = "l", col = "black")
  
  
  # Add proxy name (on left side)
  mtext(side=2,"GMST anomaly",        cex = unitscex, line = unitslinno, las = 1, col = "black",         at = 1)
  mtext(side=2,"[° C]",       cex = unitscex, line = unitslinno, las = 1, col = "black",         at = 0.5)
  #mtext(side=2,"",     cex = unitscex, line = unitslinno, las = 1, col = "black",         at = 0)
  mtext(side=2,"GMST HadCM3", cex=namcex,     line=namlin+4,     las = 1, col = col_hadcm3_temp, at = -0)
  mtext(side=2,"PAGES2k",     cex=namcex,     line=namlin+4,     las = 1, col = col_pages2k,     at = -0.5)
  mtext(side=2,"HadCRUT4",    cex=namcex,     line=namlin+4,     las = 1, col = "black",         at = -1)
  # Add units (on right side)
  
  

  
  # Add axis (on right side)
  axis(2,at=seq(-1,1,by=0.5),labels=FALSE,col="black")
  mtext(side=2,at=seq(-1,1,by=1),seq(-1,1,by=1),line=axslinno,las=1,cex=axisnumscex,col="black")
  
  text(-1100, 0.6, "(a)", cex = namcex+0.5, col = "black")
  
  # 2: d18O at Bunker (dripwater equivalent)
  
  range_d18O <- c(-9, -5)
  plot(xlimz,range_d18O,axes=FALSE,type="n",xlab="",ylab="")
  lines(Timeseries$HadCM3_Bunker_a,col=adjustcolor(col_hadcm3_d18O,0.3))
  lines(gaussdetr(Timeseries$HadCM3_Bunker_a,tsc.in=100)$Xsmooth,col=col_hadcm3_d18O)
  lines(Timeseries$SISAL_Bunker_240, col = adjustcolor(col_d18O_240, 0.3))
  lines(gaussdetr(Timeseries$SISAL_Bunker_240,tsc.in=100)$Xsmooth,col=col_d18O_240)
  lines(Timeseries$SISAL_Bunker_242, col = adjustcolor(col_d18O_242, 0.3))
  lines(gaussdetr(na.omit(Timeseries$SISAL_Bunker_242),tsc.in=100)$Xsmooth, col=col_d18O_242)
  
  #xtkz <-axTicks(2)[seq(1,by=2,length(axTicks(2))-1)]
  axis(2,col=col_hadcm3_d18O,at=seq(-9,-6),labels=FALSE)
  mtext(side=2,col=col_hadcm3_d18O,at=c(-9,-7),c(-9,-7),line=axslinno,las=1,cex=axisnumscex)
  #mtext(side=2,expression(paste("\u0394",bgroup("[","\u2030" ,"]"))),cex=unitscex,col="black",line=unitslinno,las=1)
  mtext(side=2,expression(paste("[",delta^{plain(18)}, plain(O),"]")),cex=unitscex,col=col_hadcm3_d18O,line=unitslinno,las=1, at = -7)
  mtext(side=2,"Bunker cave", cex = unitscex, col = "black",         las=1, line = unitslinno,  at = -6)
  mtext(side=2,"HadCM3",      cex = namcex,       col = col_hadcm3_d18O, las=1, line = unitslinno+4,at = -7)
  mtext(side=2,"eID240",      cex = namcex,       col = col_d18O_240,    las=1, line = unitslinno+4,at = -8)
  mtext(side=2,"eID242",      cex = namcex,       col = col_d18O_242,    las=1, line = unitslinno+4,at = -9)
  
  text(-1100, -6.3 , "(b)", cex = namcex+0.5, col = col_hadcm3_d18O)
  
  
  # FORCINGS
  
  # 3: CO2
  range_co2 = range(Timeseries$co2[index(Timeseries$co2)>-1100])
  range_co2[1] = 250
  #rgngrp<-range(ngrip[which(index(ngrip) <= 0 & index(ngrip) >= -67000)])
  plot(xlimz,range_co2,axes=FALSE,type="n",xlab="",ylab="")
  lines(Timeseries$co2_100Gsmooth,col=col_co2)
  
  axis(2,col=col_co2,at=c(250,300,350),labels=FALSE)
  mtext(side=2,col=col_co2,at=c(250,300,350),c(250,300,350),line=axslinno,las=1,cex=axisnumscex)
  mtext(side=2,expression(paste("[",plain(ppm),"]")) ,cex=unitscex,col=col_co2,line=unitslinno,las=1, at = 310)
  mtext(side=2,expression(plain(CO[plain(2)])),       cex=namcex,  col=col_co2,line=unitslinno,    las=1, at = 350)
  
  text(-1100, 350 , "(c)", cex = namcex+0.5, col = col_co2)
  
  # 4: volcanic forcing
  range_volc = range(Timeseries$volcanic)
  #rgngrp<-range(ngrip[which(index(ngrip) <= 0 & index(ngrip) >= -67000)])
  plot(xlimz,range_volc,axes=FALSE,type="n",xlab="",ylab="")
  lines(Timeseries$volcanic,col=col_volc)
  axis(2,col=col_volc,at=c(0,0.25,0.5),labels=FALSE)
  #text(x = 1980-1950-20, y = 0.25, "Mt Pinatubo", col = "black", cex = namcex)
  text(x = -100, y = 0.5, "Tambora", col = "black", cex = namcex*1.2)
  text(x = -620, y = 0.55, "Samalas", col = "black", cex = namcex*1.2)
  #text(x = -485, y = 0.4, "1465 mystery", col = "black", cex = namcex)
  mtext(side=2,col=col_volc,at=c(0,0.25,0.5),c(0.0,0.25,0.5),line=axslinno,las=1,cex=axisnumscex)
  mtext(side=2,expression(paste("volcanic forcing")),cex=namcex,  col=col_volc,line=unitslinno,las=1, at = 0.5)
  mtext(side=2,expression(paste("[",plain(AOD),"]")),cex=unitscex,col=col_volc,line=unitslinno,las=1, at = 0.375)
  
  
  text(-1100, 0.5, "(d)", cex = namcex+0.5, col = col_volc)
  
  # 4: solar forcing
  range_solar = c(1364, 1368)
  #rgngrp<-range(ngrip[which(index(ngrip) <= 0 & index(ngrip) >= -67000)])
  plot(xlimz,range_solar,axes=FALSE,type="n",xlab="",ylab="")
  lines(Timeseries$solar,col=col_solar)
  
  
  axis(2,col=col_solar,at=c(1365, 1366, 1367),labels=FALSE)
  mtext(side=2,col=col_solar,at=c(1365,1366,1367),c(1365,1366,1367),line=axslinno,las=1,cex=axisnumscex)
  mtext(side=2,expression(paste("TSI")),cex=namcex,col=col_solar,line=unitslinno,las=1, at = 1367)
  mtext(side=2,expression(paste("[",plain("W m"^"-2"),"]")),cex=unitscex,col=col_solar,line=unitslinno,las=1, at = 1366)
  
  
  text(-1100, 1367 , "(e)", cex = namcex+0.5, col = col_solar)
  
  # #FORCINGS
  
  par(xpd=NA)
  tmp.y2 <- TeachingDemos::cnvrt.coords(x=NA, y=tmp.y, input='dev')$usr$y
  segments(c(min(xlimz),max(xlimz),min(xlimz),min(xlimz)),c(par('usr')[3],par('usr')[3],par('usr')[3],tmp.y2),c(par('usr')[1],par('usr')[2],max(xlimz),max(xlimz)),c(tmp.y2,tmp.y2,par('usr')[3],tmp.y2))
  segments(x0=c(-1000, -500, 0),x1=c(-1000, -500, 0),y1=par('usr')[3],y0=1*tmp.y2,lty=2,col="grey")
  #segments(x0=c(-1050, -950, -850, -750, -650, -550, -450, -350, -250, -150, -50),x1=c(-1050, -950, -850, -750, -650, -550, -450, -350, -250, -150, -50),y1=par('usr')[3],y0=1*tmp.y2,lty=2,col="black")
  axis(1,at=seq(-1100,50,by=100),labels=FALSE);
  mtext(side=1,at=seq(-1000,0,by=500),seq(1000,0,by=-500),line=1,cex=axisnumscex)
  mtext(side=1,line=3,"time [y BP]",cex=unitscex)
  
  # mtext(side=3,line=3.25,"time [yrs BP]",cex=unitscex)
  # axis(3,at=seq(-1100,50,by=100),labels=FALSE);
  # mtext(side=3,at=seq(-1000,0,by=500),seq(1000,0,by=-500),line=0.5,cex=axisnumscex)
  
  
  dev.off()
}

rm(namcex, namlin, num_ts, range_co2, range_d18O, range_solar, range_temp, range_volc, tmp.y, tmp.y2, unitscex, unitslinno, xlimz)
rm(col_co2, col_d18O_240, col_d18O_242, col_hadcm3_d18O, col_hadcm3_temp, col_pages2k, col_solar, col_volc, axisnumscex, axslinno, plot)


##############################################
##############################################
########## R functions, own library ##########
##############################################
##############################################


# update 28 - 06 - 2016

# update 215 - 02 - 2017
#  Adding North function


##################
##### REMARK #####
##################

# save2, png2, not working
# improve heatmap_correlation

####################
##### OVERVIEW #####
####################

### SAVE
# save2(data, file, save)
#        -> safe save in case a file is already present
# png2(data, file, width, height, save) 
#        -> idem for picture


## GRAPHS
# Map(lon_min, lon_max, lat_min, lat_max, no_cores)
#        -> Layer of the world borders for maps
# heatmap_correlation(pval, sign)
#        -> Heatmap with colors showing the significancy
# North(x, y, shape, size)
#        -> Plot North arrow on ggplot and ggmap plots


## CYCLONE
# plot_track(data, nplot, low, mid, high, title, legend, mp, force_long, bg_map, savePlot, name, return_plot)
#        -> Create a map showing tracks
# point_density(data, lat = seq(-90, 90, 2.5), lon = seq(-180, 180, 2.5), dmax = 400000, factor = 1, bg_map = NA, xstep = 40, ystep = 20, title = "", legend_name = "",  max = NA, savePlot = FALSE, saveMatrix = FALSE, return_plot = FALSE, name = NA, no_cores = 1, verbose = TRUE) 
# track_density(data, lat, lon, dmax, method, n.spline, factor, bg_map, savePlot, saveMatrix, return_plot, name, no_cores, verbose)
#        -> Create a map showing a density of tracks
# plot_cyclonePhaseDiagram(data, savePlot, name)
#        -> Create a plot of the diagramme phase for a specific storm, under development


###################
### Preliminary ###
###################

# To load pictures from the repository on demand
.base_fun.dir <- dirname(sys.frame(1)$ofile)

#################
###### Save #####
#################

## To not replace a file with a same name

save2 <- function(data, name, save = TRUE) {

  if (save) {
    if (file.exists(paste(name, ".dat", sep = ""))) {
      i <- 1
      while (file.exists(paste(name, i, ".dat", sep = ""))) {
        i <- i + 1 
      }
      save(data, file = paste(name, i, ".dat", sep = "")) 
    }

    save(data, file = paste(name, ".dat", sep = "")) 
  }
}

png2 <- function(plot, name, width = 1200, height = 600, save = TRUE) {

  if (save) {
    if (file.exists(paste(name, ".png", sep = ""))) {
      i <- 1
      while (file.exists(paste(name, i, ".png", sep = ""))) {
        i <- i + 1 
      }
      png(file = paste(name, i, ".png", sep = ""), width = 1200, height = 600) 
    }

    png(file = paste(name, ".png", sep = ""), width = 1200, height = 600) 
    print(plot)
    graphics.off()
  }
}



####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


###################
###### Graphs #####
###################

# function Map
##############

# give a SpatialPolygons object of the world coasts & boundaries
# Nicely shaped depending on the area selected

Map <- function(lon_min, lon_max, lat_min, lat_max, no_cores = 1) {

  # longitude may vary between -360 and 360, so we can choose longitude break for the worlmap

  ## Verifying the data provided
  if (lat_min >= lat_max) {
    warning("latitude min >= latitude max")
    toto <- lat_min
    lat_min <- lat_max
    lat_max <- toto
  } else if (lon_min >= lon_max) {
    warning("Error : longitude min >= longitude max")
    toto <- lon_min
    lon_min <- lon_max
    lon_max <- toto
  } else if (lon_min < -360 | lon_max > 360) {
    warning("longitude min < -180 ou longitude max > 540")
    lon_min <- max(-360, lon_min)
    lon_max <- min(360, lon_max)
  } else if (lat_min < -90 | lat_max > 90) {
    warning("Error : latitude min < -90 ou latitude max > 90")
    lat_min <- max(-90, lat_min)
    lat_max <- min(90, lat_max)
  } else if (lon_max - lon_min > 360) {
    warning("The longitude overlaps")
  }

  if (no_cores < 0 ) {
    print("Error : no_cores is negatif")
  }

  require(maps)
  require(maptools)
  require(raster)
  require(rgeos)

  # Downloadind of the map and improving of some features
  w_map <- map("world", fill = TRUE, plot = FALSE)

  w_map <- map2SpatialPolygons(w_map, IDs = w_map$names, 
             proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  w_map <- gSimplify(w_map, tol = 0.00001)
  w_map <- suppressWarnings(gBuffer(w_map, byid = TRUE, width = 0))

  # Extreme Siberie
  #if (lon_min < -160) {
  clip.extent <- as(extent(180, 200, -90, 90), "SpatialPolygons")
  proj4string(clip.extent) <- CRS(proj4string(w_map))
  gI <- gIntersects(w_map , clip.extent , byid = TRUE)
  out <- lapply(which(gI) , function(x){ gIntersection(w_map[x, ], clip.extent) })
  keep <- sapply(out, class)
  out <- out[keep == "SpatialPolygons"]
  sib_map <- lapply(1:length(out), 
                           function(i) { 
                             Pol <- slot(out[[i]], "polygons")[[1]]
                             Pol <- Polygons(list(Polygon(cbind(slot(Pol@Polygons[[1]],"coords")[,1]-360, slot(Pol@Polygons[[1]],"coords")[,2]))), ID = as.character(i))
                             return(Pol)
                           })
  w_map <- SpatialPolygons(c(w_map@polygons,sib_map), proj4string=CRS(proj4string(w_map)))
  #}
  
  fun <- function(lon_min, lon_max, lat_min, lat_max, w_map) {

    # Creation of the window
    clip.extent <- as(extent(lon_min, lon_max, lat_min, lat_max), "SpatialPolygons")
    proj4string(clip.extent) <- CRS(proj4string(w_map))

    # Making the new map
    gI <- gIntersects(w_map , clip.extent , byid = TRUE)

    if (no_cores > 1) { # to speed up the program
      require(parallel)
      cl <- makeCluster(no_cores, type="FORK")
      out <- parLapply(cl, which(gI) , function(x){ gIntersection(w_map[x, ], clip.extent) })
      stopCluster(cl)
    } else {
      out <- lapply(which(gI) , function(x){ gIntersection(w_map[x, ], clip.extent) })
    }

    keep <- sapply(out, class)
    out <- out[keep == "SpatialPolygons"]

    return(out)

  }


  if (lon_max > 180) {
    out1 <- fun(lon_min, 180, lat_min, lat_max, w_map)
    out2 <- fun(180, lon_max, lat_min, lat_max, shift(w_map,360))
    out <- c(out1,out2)
  } else {
    out <- fun(lon_min, lon_max, lat_min, lat_max, w_map)
  }

  w_map <- SpatialPolygons(lapply(1:length(out), 
                           function(i) { 
                             Pol <- slot(out[[i]], "polygons")[[1]]
                             slot(Pol, "ID") <- as.character(i)
                             return(Pol)
                           }))
  return(w_map)

}


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


# function heatmap_correlation
##############################

# Plot a heatmap, the color show how significant is the correlation.

heatmap_correlation <- function(data, threshold = c(0.1, 0.05, 0.01), title = "", savePlot = FALSE, name = NA, return_plot = FALSE) {

  # data is a data.frame with the following colomn :
  #       data$X : values to be plotted on the X-axis (must be factor)
  #       data$Y : values to be plotted on the Y-axis (must be factor)
  #   data$group : possible values that group each data$X values
  #    data$pval : p-value of each correlation
  #    data$sign : sign of the correlation

  # threshold : pvalue thresholds that indicate when to change the color
  # title : tilte of the plot

  require(ggplot2)

##### computing the threshold
  data$threshold <- 0
  for (i in 1:length(threshold)) {
    data$threshold[which(data$pval < threshold[i])] <- i
  }

  data$threshold = factor(data$threshold * data$sign, levels = -3:3)

##### some parameters :
  col <- colorRampPalette(c("blue4","white","red4"))(n = length(threshold)*2 + 1)

  theme <- theme(
             plot.title   = element_text(size = 24, face = "bold", margin = margin(0, 0, 20, 0)),
             plot.margin  = unit(c(1, 0, 1, 1), units = "cm"),
             panel.background = element_blank(),
             legend.title = element_text(size = 16),
             legend.text = element_text(size = 12),
             legend.key = element_blank(),
             axis.title   = element_blank(),
             axis.text    = element_text(size = 16),
             axis.ticks = element_blank(),
             strip.background = element_blank(),
             strip.text = element_text(size = 22) )

##### The legend
  p <- ggplot(data, aes(y = Y, x = X)) +
         geom_point(aes(y = Y, x = X, colour = threshold), shape = 15, size = 8) + 
         scale_colour_manual(drop=FALSE, name = "LINK :\n\npositive           negative", 
                             values = col, breaks = c(3,2,1,-3,-2,-1), 
                             labels = c("pvalue < 0.01  ", "pvalue < 0.05", "pvalue < 0.1",
                                        "pvalue < 0.01", "pvalue < 0.05", "pvalue < 0.1"),
                             guide = guide_legend(ncol = 2, order = 1)) +
         geom_point(aes(y = Y, x = X, shape = threshold), size = 8, fill = "white") + 
         scale_shape_manual(name = "", values = rep(22,7), labels = "Not significant", breaks = 0, 
                            guide = guide_legend(order = 2))

##### The plot
  p <- p +
         geom_tile(aes(fill = threshold), colour = "black") +
         scale_fill_manual(drop=FALSE, values = col, guide = FALSE) +
         facet_grid(. ~ group) + # To split the plot along a variable
         labs(title = title) + theme


  if (savePlot) {
    png2(p,name, width = 1200, height = 600, savePlot)
  }

  if (return_plot) {
    return(p)
  } else if (!savePlot) {
    print(p)
  }

}


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


# functions geom_North, scale_North and associates
###################################################

#Intro
library(png)
library(grid)
library(ggplot2)

#Some function for grid
NorthGrob <- function(x, y, country, size=1, alpha=1){
  grob(x=x, y=y, country=country, size=size, cl = "North")
}

drawDetails.North <- function(x, recording=FALSE){
  grid.raster(x$x, x$y, 
              width = x$size*unit(4,"mm"), height = x$size*unit(6.35,"mm"),
              image = North_img, interpolate=FALSE)
}

scale_North <- function(..., guide = "legend") {
  sc <- discrete_scale("shape", "identity", scales::identity_pal(), ..., guide = guide)

  sc$super <- ScaleDiscreteIdentity
  class(sc) <- class(ScaleDiscreteIdentity)
  sc
}

GeomNorth <- ggproto("GeomNorth", Geom,
               required_aes = c("x", "y"),
               draw_key = function (data, params, size) 
               {
                 NorthGrob(0.5,0.5, shape=data$shape,  size=data$size)
               },

               draw_group = function(data, panel_scales, coord) {
                 coords <- coord$transform(data, panel_scales)     
                 NorthGrob(coords$x, coords$y, coords$shape, coords$size)
               }
)

geom_North <- function(mapping = NULL, data = NULL, stat = "identity",
                position = "identity", na.rm = FALSE, show.legend = NA, 
                inherit.aes = TRUE, ...) {

  fig_file <- paste0(.base_fun.dir, "/North_symbols/", 
                formatC(as.numeric(as.character(data$shape)), width=2, flag="0"), ".png")

    if (file.exists(fig_file)) {
      North_img <<- readPNG(fig_file) 
    } else {
      stop(paste("There is no North symbol of shape", shape,
             "it should be picked between 1 and 19"))
    }

  layer(
    geom = GeomNorth, mapping = mapping,  data = data, stat = stat, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


# function scalebar
###################

scalebar = function(x, y, w, lat, len = 100){
  # x,y = lower left coordinate of bar
  # w = width of bar
  # len = first tick in km, then 2*len, 4*len and 8*len
  # lat where to approximate the distance

  len_deg <- len / (111.30*cos(lat*2*pi/360))

  bar = data.frame( 
    xmin = c(x, x + len_deg/2, x + len_deg, x + 2*len_deg, x + 4*len_deg),
    xmax = c(x + len_deg/2, x + len_deg, x + 2*len_deg, x + 4*len_deg, x + 8*len_deg),
    ymin = y,
    ymax = y+w,
    fill.col = c("black", "white", "black", "white", "black"),
    z = c(1, 0, 1, 0, 1)
  )

  labs = data.frame(
    xlab = c(bar$xmin[-2], x + 8*len_deg, x + 9.2*len_deg),
    ylab = c(rep(y-w*1, 5), y+w/2),
    text = c(as.character(c(0, len, 2*len, 4*len, 8*len)), "kilometres")
    )

  #geom_rect(data=bar, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
   #         show.legend = F,  color = "black", fill = bar$fill.col) +
  #geom_text(data=labs, aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F) 

  list(bar, labs)
}

#sb = scalebar(62,19.5,.4,28)

# Plot map

#p10 +
#  geom_rect(data=sb[[1]], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
#            show.legend = F,  color = "black", fill = sb[[1]]$fill.col) +
#  geom_text(data=sb[[2]], aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F) 



####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


# function GeomRaster for coord_map()
#####################################

###




####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################



###################
##### Cyclone #####
###################

# function plot_track
#####################

# Plot the given tracks

plot_track <- function(data, type, scale = NA, nplot = NA, plotnum = FALSE, limits = c(-180, 180, -90, 90), bg_map = NA, title = NULL, legend = NULL, savePlot = FALSE, name = NA, return_plot = FALSE, pt = TRUE) {

  # Dimension of data : variables x storms x timesteps
  # 	   variable 1 : latitude
  # 	   variable 2 : longitude
  # 	   variable 3 : color variable (ex : wind)
  # Missing value (not ploted data) = NA

  # type : "wind" -> the color variable provided is continuous (ex : wind)
  #        "fac"  -> the color variable provided is discrete   (ex : warm core) 
  #        "cat"  -> the color variable provided is the wind in knots,
  #                  and the color depend on the S-S category

  # nplot : number of tracks to plot. NA : all ploted

  # limits : limits of the graph in the order : minimum longitude
  #                                             maximum longitude
  #                                             minimum latitude
  #                                             maximum latitude

  # plotnum : plot the numero of each tracks

  # bg_map : map in the background, of type "S4", might ne given by : Map(-180, 180, -90, 90)
  # title, legend : title of the plot and of the legend

  # savePlot : flag to save the plot
  # name : name of the save

  # return_plot : return the ggplot to further modification

  require(ggplot2)


  # verifing the data provided
  #####

  if (type == "cat") {
    if (scale == "1min") {
      scales <- c(34,64,83,96,113,137)
    } else  if (scale == "10min") {
      scales <- c(30,56,73,84,99,120)
    } else  if (scale == "1hour") {
      scales <- c(27,51,66,77,90,110)
    } else  if (scale == "3hours") {
      scales <- c(25,46,60,70,82,100)
    } else {
      stop("error scale not provided when type = 'cat'")
    }
  }

  if (length(dim(data)) == 3) {
    nstorm <- dim(data)[2]
    ntime  <- dim(data)[3]
  } else {
    stop("Wrong dimension for data")
  }

  if (min(data[1, , ], na.rm = TRUE) < -90 | max(data[1, , ], na.rm = TRUE) > 90) {
    stop(paste("latitude minimale :", min(data[1, , ], na.rm = TRUE), 
               "latitude minimale :", max(data[1, , ], na.rm = TRUE)))
  }

  if (max(data[2, , ], na.rm = TRUE) - min(data[2, , ], na.rm = TRUE) > 360) {
    stop(paste("longitude minimale :", min(data[2, , ], na.rm = TRUE), 
               "longitude maximale :", max(data[2, , ], na.rm = TRUE)))
  } else if (min(data[2, , ], na.rm = TRUE) < -360 | max(data[2, , ], na.rm = TRUE) > 360 ) {
    warning(paste("longitude minimale :", min(data[2, , ], na.rm = TRUE),
                  "longitude maximale :", max(data[2, , ], na.rm = TRUE)))
  }

  ## if there is empty rows
  s_empty <- c() 
  for (s in 1:nstorm) {
    if (length(which(is.na(data[1, s, ]))) != length(which(is.na(data[2, s, ])))) {
      stop(paste("for storm :", s,
                  "\n number of latitude provided =",  ntime - length(which(is.na(data[1, s, ]))),
                  "\n number of longitude provided =", ntime - length(which(is.na(data[2, s, ])))))
    } else {
      if (length(which(is.na(data[1, s, ]))) != length(which(is.na(data[3, s, ])))) {
        warning(paste("Not the same number of variable color and points for storm :", s,
                    "\n number of points provided =",       ntime - length(which(is.na(data[1, s, ]))),
                    "\n number of third values provided =", ntime - length(which(is.na(data[3, s, ])))))
      }

      if (length(which(is.na(data[1, s, ]))) %in% c(ntime,ntime -1 )) {
        warning(paste("the storm", s,"do not contain any values, or only one"))
        s_empty=c(s_empty,s)
      }
    }

  }

  # if there is 0 instead of NA
  test <- which(data[1, , ] == 0 & data[1, , ] == 0 & data[1, , ] == 0, arr.ind = T)
  data[cbind(1, test)] <- NA
  data[cbind(2, test)] <- NA
  data[cbind(3, test)] <- NA


  if ( (length(limits) != 4) |
       (limits[1] > limits[2]) |
       (limits[3] > limits[4]) |
       (limits[1] < -360) |
       (limits[2] > 360) |
       (limits[3] < -90) |
       (limits[4] > 90)) {
    stop("Wrong values for limits, the formats is : c(minimum longitude, maximum longitude, minimum latitude, maximum latitude)")
  }


  # name not given
  if (is.na(name) & (savePlot)) {
    stop("saves activated but no name provided")
  }

  # bg_map not given
  if (typeof(bg_map) != "S4") {
    bg_map <- Map(limits[1], limits[2], limits[3], limits[4], 4)
  }

  # in the case of empty row
  if (length(s_empty) == 0) {  
    samp <- 1:nstorm
  } else {
    samp <- (1:nstorm)[-s_empty]
    nstorm <- nstorm - length(s_empty)
  }

  # Selecting the tracks
  if (is.na(nplot)) {
    nplot <- nstorm
    test <- samp
  } else {
    if (nstorm > nplot) {
      test <- sample(samp,nplot)
      title <- paste(title," (", nplot,"/", nstorm,")",sep="")
    } else {
      nplot <- nstorm
      test <- samp
    }
  }

  

  # Defining the plot
  #####
  p <- ggplot() +
       geom_polygon(data = bg_map, aes(x = long, y = lat, group = group), 
                    colour = "black", fill = "#009E73", size = .4) +
       labs(title = title, y = "latitude (in degree N)", x = "longitude (in degree W)") + 

       theme(plot.title       = element_text(size = 24, face = "bold", margin = margin(0, 0, 20, 0)),
             panel.background = element_rect(fill = "white", colour= "black", size = 1),
             plot.margin      = unit(c(1, 0, 1, 1), units = "cm"),
             panel.grid.major = element_line(colour = "grey70"),
             panel.grid.minor = element_line(colour = "grey90"),
             axis.title       = element_text(size = 18),             
             axis.title.x     = element_text(vjust = 0),
             axis.title.y     = element_text(vjust = 1),
             axis.text        = element_text(size=12),
             legend.title     = element_text(size = 16),
             legend.text      = element_text(size = 12),
             legend.key       = element_rect(fill = "white")) + 

       scale_y_continuous(breaks = pretty(c(limits[3], limits[4]), n = 10),
                          limits = c(limits[3], limits[4]),
                          expand = c(0,0)) + 
       scale_x_continuous(breaks = pretty(c(limits[1], limits[2]), n = 10),
                          limits = c(limits[1], limits[2]),
                          expand = c(0,0))

  ## The scales
  if (type == "wind") {
    p <- p + scale_colour_gradientn(name = legend, breaks = seq(0,80,20), limits = c(0,80), colours=c("darkred","red","violet","blue","darkblue"), values = c(1,0.7,0.5,0.3,0), oob = scales::squish)

  } else if (type == "fac") {
    if (length(unique(data[3,,]) == 2)) {
      p <- p + scale_colour_manual(values = c("0" = "darkblue", "1" = "firebrick3"), name = legend)
    } else {
      p <- p + scale_colour_brewer(type = "qual", palette = "Set1", name = legend)
    }

  } else if (type == "cat") {
    p <- p + scale_colour_manual(
               values = c("#5ebaff", "#00faf4", "#ffffcc", "#ffe775", "#ffc140", "#ff8f20", "#ff6060"),
               label = c("Tropical Depression", "Tropical Cyclone", "Hurricane 1", "Hurricane 2", "Major Hurricane 3", "Major Hurricane 4", "Major Hurricane 5"),
               name = legend, guide = guide_legend(reverse = TRUE))   
  } else {
    stop("Wrong name for type, choose among 'wind', 'fac' and 'cat'")
  }


  # Formatting the data
  ####

  data2=c()
  for (i in 1:nplot) {
    data2 <- rbind(data2,rep(NA,3),t(data[, test[i], ]))
  }

  data2<- data.frame(data2)
  colnames(data2) <- c("lat", "lon", "wind")

  if (type == "fac") {
    data2$wind <- as.factor(data2$wind)
  }

  if (type == "cat") {

    data3 <- rbind(data.frame(data2, group = -1), # Tropical depression
                   data.frame(data2, group =  0), # Tropical cyclone
                   data.frame(data2, group =  1), # Hurricane category
                   data.frame(data2, group =  2),
                   data.frame(data2, group =  3),
                   data.frame(data2, group =  4),
                   data.frame(data2, group =  5))

    data3$lat[data3$group == 0 & data3$wind < scales[1] ] = NA
    data3$lon[data3$group == 0 & data3$wind < scales[1] ] = NA

    data3$lat[data3$group == 1 & data3$wind < scales[2] ] = NA
    data3$lon[data3$group == 1 & data3$wind < scales[2] ] = NA

    data3$lat[data3$group == 2 & data3$wind < scales[3] ] = NA
    data3$lon[data3$group == 2 & data3$wind < scales[3] ] = NA

    data3$lat[data3$group == 3 & data3$wind < scales[4] ] = NA
    data3$lon[data3$group == 3 & data3$wind < scales[4] ] = NA

    data3$lat[data3$group == 4 & data3$wind < scales[5] ] = NA
    data3$lon[data3$group == 4 & data3$wind < scales[5] ] = NA

    data3$lat[data3$group == 5 & data3$wind < scales[6] ] = NA
    data3$lon[data3$group == 5 & data3$wind < scales[6] ] = NA

    data2 <- data.frame(lat = data3$lat, lon = data3$lon, wind = data3$group)
  }
 
  ## Formatting the longitude, splitting when necessary
  data2$lon[which(data2$lon > limits[2])] <- data2$lon[which(data2$lon > limits[2])] - 360
  data2$lon[which(data2$lon < limits[1])] <- data2$lon[which(data2$lon < limits[1])] + 360

  split = c(1, which(abs(data2$lon[2:length(data2$lon)] - data2$lon[2:length(data2$lon)-1]) > 300) + 1, length(data2$lon) + 1)

  data3 <- c()
  for ( i in 2:length(split)) {
    if (split[i-1] != 1 & split[i] != length(data2$lon) + 1) {
      data3 <- rbind(data3,  
                     c(NA, NA, data2[(split[i-1] - 1),3]),
                     c(data2[(split[i-1] - 1),1], 
                       data2[(split[i-1] - 1),2] + sign(data2[(split[i-1]) + 1,2])*360,
                       data2[(split[i-1] - 1),3]), 
                     data2[split[i-1]:(split[i]-1),],
                     c(data2[(split[i]-1) + 1,1],
                       data2[(split[i]-1) + 1,2] - sign(data2[(split[i]-1) + 1,2])*360, 
                       data2[(split[i]-1) + 1,3]))
    } else if (split[i-1] == 1 & split[i] == length(data2$lon) + 1) {
      data3 = data2
    } else if (split[i-1] == 1) {
      data3 <- rbind(data3,  
                     c(NA, NA, data2[(split[i-1]) + 1,3]), 
                     data2[split[i-1]:(split[i]-1),],
                     c(data2[(split[i]-1) + 1,1],
                       data2[(split[i]-1) + 1,2] - sign(data2[(split[i]-1) + 1,2])*360, 
                       data2[(split[i]-1) + 1,3]))
    } else if (split[i] == length(data2$lon) + 1) {
      data3 <- rbind(data3,  
                     c(NA, NA, data2[(split[i-1] - 1),3]),
                     c(data2[(split[i-1] - 1),1], 
                       data2[(split[i-1] - 1),2] + sign(data2[(split[i-1]) + 1,2])*360,
                       data2[(split[i-1] - 1),3]), 
                     data2[split[i-1]:(split[i]-1),])
    }
  }

  data2 <- data3
  data2$lon <- scales::squish(data2$lon, range = c(limits[1],limits[2]))
  data2$lat <- scales::squish(data2$lat, range = c(limits[3],limits[4]))

  if (type == "cat" | type == "fac") {
    data2$wind <- as.factor(data2$wind)
  }

  # Filling the plot
  ####

  p <- p + geom_path(data = data2, aes(lon, lat, col = wind), size = 1.5)

  p <- p + geom_polygon(data = bg_map, aes(x = long, y = lat, group = group), 
                    colour = "black", fill = "transparent", size = .4)

  if (pt) {
    data.point <- data.frame(t(data[, test, 1]))
    names(data.point) <- c("lat","lon", "wind")
    p <- p + geom_point(data = data.point, aes(x = lon, y = lat), colour = "black", size = 1.5)
  }

  if (plotnum) {
    p <- p + geom_text(data = data.frame(x = data[2, test, 1], y = data[1, test, 1], label = test), aes(x, y, label = label), hjust = -0.5, vjust = -0.5)
  }

  
  if (savePlot) {
    setEPS()
    postscript(name, width = 16, height = 8)
    print(p)
    dev.off()
  }

  if (return_plot) {
    return(p)
  } else if (!savePlot) {
    print(p)
  }

}


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


# function point_density
#######################

# Plot a map of point density (ex : cyclogenesis)


point_density <- function(data, lat = seq(-90, 90, 2.5), lon = seq(-180, 180, 2.5), dmax = 400000, factor = 1, bg_map = NA, xstep = 40, ystep = 20, title = "", legend_name = "",  max = NA, savePlot = FALSE, saveMatrix = FALSE, return_plot = FALSE, name = NA, no_cores = 1, verbose = TRUE) {


  # Dimension of data : variables x storms
  # 	   variable 1 : latitude
  # 	   variable 2 : longitude

  # lat, lon : grid point where to compute the density

  # dmax : maximum distance between a strom point and a grid point to count it. (in meter)

  # factor : number by which the density is multipied (to get a more friendly number)

  # bg_map : map in the background
  # savePlot : flag to save the plot
  # saveMatrix : flag to save the matrix of the density
  # name : name of the saves

  # no_cores : number of cores used to parallelize (for methods "timestep" and "tracks1")

  # verbose : to print information


  # verifing the data provided
  if (length(dim(data)) == 2) {
    nstorm <- dim(data)[2]
  } else {
    stop("Wrong dimension for data")
  }

  if (min(data[1, ], na.rm = TRUE) < -90 | max(data[1, ], na.rm = TRUE) > 90) {
    stop(paste("latitude minimale :", min(data[1, ], na.rm = TRUE), 
               "latitude minimale :", max(data[1, ], na.rm = TRUE)))
  }

  if (max(data[2, ], na.rm = TRUE) - min(data[2, ], na.rm = TRUE) > 360) {
    stop(paste("longitude minimale :", min(data[2, ], na.rm = TRUE), 
               "longitude maximale :", max(data[2, ], na.rm = TRUE)))
  } else if (min(data[2, ], na.rm = TRUE) < -360 | max(data[2, ], na.rm = TRUE) > 360 ) {
    stop(paste("longitude minimale :", min(data[2, ], na.rm = TRUE),
               "longitude maximale :", max(data[2, ], na.rm = TRUE)))
  }


  if (is.na(name) & (savePlot | saveMatrix)) {
    stop("saves activated but no name provided")
  }


  #libs
  require(ggplot2)
  require(geosphere)
  require(parallel)

  # parameters
  nstorm <- dim(data)[2]
  ntime <- dim(data)[3]



  # list of storm point :
  lat_storm <- data[1, ]
  lat_storm <- lat_storm[which( ! is.na(lat_storm), arr.ind = TRUE)]

  lon_storm <- data[2, ]
  lon_storm <- lon_storm[which( ! is.na(lon_storm), arr.ind = TRUE)]


  # function for the main loop
  gridPointLoop_timestep <- function(i,j) {
    # i is the longitude
    # j is the latitude
      
    d <- distCosine(c(lon[i], lat[j]), cbind(lon_storm, lat_storm))

    return( length(which(d < dmax)))
  }

  # parallelise loop
#  cl <- makeCluster(no_cores, type="FORK")
#  n <- clusterMap(cl, gridPointLoop_timestep, rep(1:length(lon), each = length(lat)), rep(1:length(lat), length(lon)), SIMPLIFY = TRUE)
#  stopCluster(cl)
  n <- mapply(gridPointLoop_timestep, rep(1:length(lon), each = length(lat)), rep(1:length(lat), length(lon)), SIMPLIFY = TRUE)

  n2 <- n/nstorm*factor

  if (is.na(max)) {
    max = max(n2)
  }

  data.plot <- data.frame(lat = rep(lat, length(lon)), lon = rep(lon, each = length(lat)), density = n2)

  breaks <- pretty(c(0, max))


  p <- ggplot(data = data.plot, aes(x = lon, y = lat)) +
      geom_polygon(data = bg_map, aes(x = long, y = lat, group = group), 
                   colour = "black", fill = "#009E73", size = .4) +
      geom_raster(aes(alpha = density), fill = "red", interpolate = TRUE) +
      geom_contour(aes(z = density), colour = "darkred", breaks = c(0.01, breaks)) + 
      scale_alpha_continuous(range = c(0.01, 1), name = legend_name, breaks = breaks, limits = c(0,max)) +
      theme(axis.text    = element_text(size = 12),
        legend.text  = element_text(size = 12),
        axis.title   = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.key.size = unit(1, "cm"),
        legend.key   = element_rect(fill = "white", colour = "darkred", size = 1),
        plot.title   = element_text(size = 24, face = "bold", margin = margin(0,0,20,0)),
        plot.margin  = unit(c(1,0,1,1), units = "cm"),
        axis.title.x = element_text(margin = margin(20,0,0,0)),
        axis.title.y = element_text(margin = margin(0,20,0,0)),
        panel.background = element_rect(fill = "transparent", colour= "black"),
        panel.grid.major = element_line(colour = "grey50"),
        panel.grid.minor = element_line(colour = "grey70")) +
      scale_y_continuous(expand = c(0,0), breaks = seq(-80, 80, ystep), limits = c(min(lat), max(lat))) + 
      scale_x_continuous(expand = c(0,0), breaks = seq(-360, 360, xstep), limits = c(min(lon), max(lon))) +
      labs(title = title, y = "latitude (in degree N)", x = "longitude (in degree W)")
        


  if(saveMatrix) {
    save2(n2, file = paste(name, ".dat", sep = ""), saveMatrix) 
  }

  if (savePlot){
    png2(p, name, width = 1200, height = 600)
  }

  if (return_plot) {
    return(p)
  } else if (!savePlot)  {
    print(p)
  }

}





####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################



# function track_density
#######################

# Plot a map of the storm density




track_density <- function(data, lat = seq(-90, 90, 2.5), lon = seq(-180, 180, 2.5), dmax = 400000, method = "timestep", n.spline = 500, factor = 1, bg_map = NA, xstep = 40, ystep = 20, title = "", legend_name = "",  maxi = NA, savePlot = FALSE, saveMatrix = FALSE, return_plot = FALSE, name = NA, no_cores = 1, verbose = TRUE) {

  # Dimension of data : variables x storms x timesteps
  # 	   variable 1 : latitude
  # 	   variable 2 : longitude
  # Missing value (not ploted data) = NA
  # Rq : if there is NA in the middle of the trajectory, we do not split it

  # lat, lon : grid point where to compute the density

  # dmax : maximum distance between a strom point and a grid point to count it. (in meter)

  # type of method
  # The density at a grid point is :
  #   for "timestep", the number of timesteps with a low at less than dmax of the grid point
  #   for "tracks1", the number of tracks with one interpolated point inside the circle of
  #             radius dmax around the grid point. We do not count NA in the middle of
  #             the track. 
  #        n.spline is the number of points interpolated for each trajectories
  #   for "tracks2", the number of tracks which intersect a square of side 2*dmax 
  #             at the equator. For higher latitude, the square is reshaped in rectangle 
  #             of same area. The values are then interpolated on the lat/lon grid provided


  # factor : number by which the density is multipied (to get a more friendly number)

  # bg_map : map in the background
  # savePlot : flag to save the plot
  # saveMatrix : flag to save the matrix of the density
  # name : name of the saves

  # no_cores : number of cores used to parallelize (for methods "timestep" and "tracks1")

  # verbose : to print information



  # verifing the data provided
  if (length(dim(data)) == 3) {
    nstorm <- dim(data)[2]
    ntime <- dim(data)[3]
  } else {
    stop("Wrong dimension for data")
  }

  if (min(data[1, , ], na.rm = TRUE) < -90 | max(data[1, , ], na.rm = TRUE) > 90) {
    stop(paste("latitude minimale :", min(data[1, , ], na.rm = TRUE), 
               "latitude minimale :", max(data[1, , ], na.rm = TRUE)))
  }

  if (max(data[2, , ], na.rm = TRUE) - min(data[2, , ], na.rm = TRUE) > 360) {
    stop(paste("longitude minimale :", min(data[2, , ], na.rm = TRUE), 
               "longitude maximale :", max(data[2, , ], na.rm = TRUE)))
  } else if (min(data[2, , ], na.rm = TRUE) < -360 | max(data[2, , ], na.rm = TRUE) > 360 ) {
    stop(paste("longitude minimale :", min(data[2, , ], na.rm = TRUE),
               "longitude maximale :", max(data[2, , ], na.rm = TRUE)))
  }


  s_empty <- c() # if there is empty row
  for (s in 1:nstorm) {
    if (length(which(is.na(data[1, s, ]))) != length(which(is.na(data[2, s, ])))) {
      stop(paste("for storm :", s,
                  "\n number of latitude provided =",  ntime - length(which(is.na(data[1, s, ]))),
                  "\n number of longitude provided =", ntime - length(which(is.na(data[2, s, ])))))
    } else {
      if (length(which(is.na(data[1, s, ]))) %in% c(ntime,ntime -1 )) {
        warning(paste("the storm", s,"do not contain any values, or only one"))
        s_empty=c(s_empty,s)
      }
    }

  }

  
  if (length(s_empty) > 0) {
    data <- data[,-s_empty, ]
  }

  if (is.na(name) & (savePlot | saveMatrix)) {
    stop("saves activated but no name provided")
  }



  #libs
  require(ggplot2)


  # parameters
  nstorm <- dim(data)[2]
  ntime <- dim(data)[3]


  # case method timestep
  ######################
  if (method == "timestep") {

    require(geosphere)
    require(parallel)

    # list of storm point :
    lat_storm <- data[1, , ]
    lat_storm <- lat_storm[which( ! is.na(lat_storm), arr.ind = TRUE)]

    lon_storm <- data[2, , ]
    lon_storm <- lon_storm[which( ! is.na(lon_storm), arr.ind = TRUE)]


    # function for the main loop
    gridPointLoop_timestep <- function(i,j) {
      # i is the longitude
      # j is the latitude
      
      d <- distCosine(c(lon[i], lat[j]), cbind(lon_storm, lat_storm))

      return( length(which(d < dmax)))
    }

    # parallelise loop
    cl <- makeCluster(no_cores, type="FORK")
    n <- clusterMap(cl, gridPointLoop_timestep, rep(1:length(lon), each = length(lat)), rep(1:length(lat), length(lon)), SIMPLIFY = TRUE)
    stopCluster(cl)



  # Case method tracks1
  #####################
  } else if (method == "tracks1") {

    # interpolation or the tracks : Each tracks get the same number of points
    n.spline <- 500
    spline.lat <- array(0, dim = c(n.spline, nstorm))
    spline.lon <- array(0, dim = c(n.spline, nstorm))
    for (s in 1:nstorm) {
      nb_ts <- seq(min(which( ! is.na(data[1, s, ]))), max(which( ! is.na(data[1, s, ]))))
      spline.lat[, s] <- spline(nb_ts, data[1, s, nb_ts], n = n.spline)$y
      spline.lon[, s] <- spline(nb_ts, data[2, s, nb_ts], n = n.spline)$y
    }

    # list of storm point :
    lon_storm <- c(spline.lon)
    lat_storm <- c(spline.lat)


    # function for the main loop
    gridPointLoop_tracks1 <- function(i, j) {
      # i is the longitude
      # j is the latitude
    
      d <- distCosine(cbind(lon[i], lat[j]), cbind(lon_storm, lat_storm))
    
      #then test if one of the point is in the circle for each strom
      ns <- 0
      for (s in 1:nstorm) {
        dist <- d[1:n.spline + (s - 1) * n.spline]
        if (any(dist < dmax)) {
          ns <- ns + 1
      }}
      return(ns)
    }


    # parallelise loop
    cl<-makeCluster(no_cores, type="FORK")
    system.time(n <- clusterMap(cl, gridPointLoop_tracks1, rep(1:length(lon), each = length(lat)), rep(1:length(lat), length(lon)), SIMPLIFY = TRUE))
    stopCluster(cl)




  } else if (method == "tracks2") {

    warning("not working well yet")

    require(maps)
    require(maptools)
    require(raster)
    require(rgeos)
    require(geosphere)
    require(akima)



    # Creating the lines tracks
    tracks <- list()
    for (s in 1:nstorm) {
      ts_max <- max(which(! is.na(data[1,s,])))
      tracks <- c(tracks, list(Lines(Line(t(data[2:1,s,1:ts_max])), ID = as.character(s))))
    }

    all_tracks <- SpatialLines(tracks, proj4string = CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
      

    ### Creating the Polygons

    # Defining the grid
    dx = abs(mean(lon[2:length(lon)-1]-lon[2:length(lon)]))
    #dx = 2 / sqrt(2) * dmax  / 40076000 * 360 # divide by the length of the equator to get the º of long
    dx = 90 / ceiling(90 / dx)

    # grid point for the polygons
    lon2 = seq(-180,180,dx)
    lat2 = seq(-90,90,dx)

    # Creating the polygons
    make_poly <- function(j,id) {
      P <- list()
      for (i in 2:length(lon2)) {
        id <- id + 1
        lonx <- lon2[i-1]
        lonn <- lon2[i]
        latx <- lat2[j-1]
        latn <- lat2[j]
        P1 <- Polygon(cbind(c(lonn, lonx, lonx, lonn, lonn),c(latn,latn,latx,latx,latn)))
        P <- c(P,list(Polygons(list(P1), ID = as.character(id))))
      }
        return(P)
    }

    cl <- makeCluster(no_cores, type="FORK")
    P <- clusterMap(cl,make_poly,2:length(lat2),0:length(lat2)*length(lon2), SIMPLIFY = TRUE)
    stopCluster(cl)

    P <- SpatialPolygons(P[1:length(P)], proj4string = CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
    area <- areaPolygon(P)


    # grid point for the center of the polygons
    lon3 <- (lon2[2:length(lon2)] + lon2[2:length(lon2) - 1])/2
    lat3 <- (lat2[2:length(lat2)] + lat2[2:length(lat2) - 1])/2
    nlon3 <- length(lon3)
    nlat3 <- length(lat3)

    # Count the intersection :
    gI <- gIntersects(all_tracks, P , byid = TRUE)
    dens <- rowSums(gI) / area * pi * dmax * dmax
    dens2 <- matrix(dens, nlon3, nlat3)


    # Smoothing over the squares included in the circle of radius dmax
    # hypothesis :
    #      - meridien are parallel
    #      - the circle of radius dmax cover at most 3 squares on the same longitude 
    #        (true by the definition of dx, lat2/lon2, except at the poles)

    # To compute the are of a circular segment
    #    h is the closest distance between the center and the segment
    area <- function(h) {
      if (h > dmax) {
        return(0)
      }
      A <- dmax**2 * acos(h / dmax)
      B <- h * sqrt(dmax**2 - h**2)
      if (A > B) {
        return(A-B)
      } else {
        warning("small angle for circular segment  : acos return a two small value")
        return(0)
      }
    }

    smoothDens <- matrix(0, nlon3, nlat3) # the smoothed values

    # length of a square, in m
    a <- distCosine(c(0,0),c(0,dx))

    # Loop over each latitude to update the caracteristique of the width of the square
    for ( j in (nlat3/2 + 1):nlat3 ) { # We begin at the equator
      b <- a * cos(lat3[j] * pi / 180) # the width of the square in m; a is the length in m

      # We suppose the circle cover 2n-1 square in a row 
      n <- round(dmax / b) + 1

      # We suppose the circle cover 2m-1 square in a column
      m <- round(dmax / a) + 1

      # We removed the most poleward points
      if (j + m > nlat3) {
        break
      }

      # the unknown fraction of square covered by the circle :
      #      the first element is the square at the center of the circle,
      #      then we go left until reaching the nth square
      #      then we do the same for the row below
      # and so one untli reaching th mth lines
      #      the fraction for the other squares may be computed from the first ones
      
      # We will form a system of equation :
      M <- matrix(0, m*n, m*n)
      r <- rep(0, m*n)      

      i <- 1 # next equation line

      # equations 1 : to find where the circle cut the longitude between two row
      #    n1 is the number of square completely covered by the circle in the kth row
      #    n - n1 - 1 is the number of square completely outside of the circle
      #      in both lower and higher row
      #    n1 is biggest number so that the diagonal of the square formed by the 2n1 -1 smaller squares is inferior to the diameter of the circle

      if (m > 1) {
        for (k in 1:(m-1)) {

          n1 <- 1
          while( sqrt( ((2*(n1-1) + 1)*b)**2 + ((2*(k-1) + 1)*a)**2 )  < dmax*2 ) {
            n1 <- n1 + 1
            if (n1 > n) {
              stop(paste("we didn't find n1 for the line ",k))
            }
          }
          n1 <- n1 - 1
  
          if (n1 == 1) {
            M[i, n1 + n*(k-1)] <- 1
            r[i] <- 1
            i <- i+1
          } else {
            diag(M[i:(i + n1 -1), 1:n1 + n*(k-1)]) <- 1
            r[i:(i + n1 -1)] <- 1
            i <- i + n1
          }

          if (n1 == n - 1) {
#           nothing
          } else if (n1 == n - 2) {
            M[i,n*(k+1)] <- 1
            i <- i + 1
          } else {
            diag(M[i + 0:(n-n1-2),n*k + (n1+2):n]) <- 1
            i <- i + n-n1-1
          }
        }
      }

      # equations 2 : to compute the area of circle covering each column of square

      # we begin with the central column
      r[i] <- (pi * dmax**2 - 2*area(b / 2)) / ((2*m - 1)*a*b)
      M[i,1] <- 1/(2*m - 1)
      M[i,(1:(m - 1))*n + 1] <- 2/(2*m - 1)
      i <- i + 1

      # then we loop over the other columns
      if (n > 2) {
        for ( k in 2:(n-1) ) {
          r[i] <- (area(b*(k - 3/2)) - area(b*(k - 1/2))) / ((2*m - 1)*a*b)
          M[i,k] <- 1/(2*m - 1)
          M[i,(1:(m - 1))*n + k] <- 2/(2*m - 1)
          i <- i + 1
        }
      }

      # we finish with the lefter column
      if (n > 1) {
        r[i] <- area(b*(n - 3/2)) / ((2*m - 1)*a*b)
        M[i,n] <- 1/(2*m - 1)
        M[i,(2:m)*n] <- 2/(2*m - 1)
        i <- i + 1
      }

      # equation 3 : same as equation 2 but for each row (except the central one, not needed)


      # we do nothing with the central row

      # then we loop over the other rows
      if (m > 2) {
        for ( k in 2:(m-1) ) {
          r[i] <- (area(a*(k - 3/2)) - area(a*(k - 1/2))) / ((2*n - 1)*a*b)
          M[i,(k-1)*n + 1] <- 1/(2*n - 1)
          M[i,(2:n) + (k - 1)*n] <- 2/(2*n - 1)
          i <- i + 1
        }
      }

      # we finish with the lowest row
      if (m > 1) {
        r[i] <- area(a*(m - 3/2)) / ((2*n - 1)*a*b)
        M[i,(m - 1)*n + 1] <- 1/(2*n - 1)
        M[i,(2:n) + (m - 1)*n] <- 2/(2*n - 1)
        i <- i + 1
      }

      if (i != m*m+1) {
        warning("there is not enough equations to solve it")
      }

      # Resolution :
      f <- solve(M,r) # fraction of each squares covered by the circle

      f <- f * a * b / (pi * dmax**2) # fraction of the circle covering each squares

      ftot <- 4*sum(f) - 2*sum(f[1:n]) -2*sum(f[(0:(m-1))*n + 1])

      if ( signif(ftot,5) != 1) {
        warning(paste("we do not exactly covert a circle",ftot))
      }

      # storing the results
      # !!!! cas ou i-n / i+n out of bound
      for ( i in 1:nlon3 ) {
        I1 = i:(i-n+1)
        I2 = i:(i+n-1)     
        I1[which(I1 <= 0)] <- I1[which(I1 <= 0)] + nlon3
        I2[which(I2 > nlon3)] <- I2[which(I2 > nlon3)] - nlon3

        j1 <- j
        J1 <- 0:(1-m) + j1
        J2 <- 0:(m-1) + j1

        smoothDens[i,j1] <- sum((dens2[I1,J1] + dens2[I1,J2] + dens2[I2,J1] + dens2[I2,J2])*f) - sum((dens2[I1,j1] + dens2[I2,j1])*f[1:n]) - sum((dens2[i,J1] + dens2[i,J2])*f[(0:(m-1))*n + 1])


        j2 <- nlat3 - j + 1
        J1 <- 0:(1-m) + j2
        J2 <- 0:(m-1) + j2

        smoothDens[i,j1] <- sum((dens2[I1,J1] + dens2[I1,J2] + dens2[I2,J1] + dens2[I2,J2])*f) - sum((dens2[I1,j2] + dens2[I2,j2])*f[1:n]) - sum((dens2[i,J1] + dens2[i,J2])*f[(0:(m-1))*n + 1])
      }
    }


    # interpolation
    n <- bicubic(lon3, lat3, smoothDens, rep(lon, each = length(lat)), rep(lat, length(lon)))$z
    n [which(n < 0, arr.ind = TRUE)] <- 0

  } else {
    print("ERROR : wrong value for parameter method")
    return(NA)
  }

  n2 <- n/nstorm*factor

  if (is.na(maxi)) {
    maxi = max(n2)
  }

  data.plot <- data.frame(lat = rep(lat, length(lon)), lon = rep(lon, each = length(lat)), density = n2)

  breaks <- pretty(c(0, maxi))


  p <- ggplot(data = data.plot, aes(x = lon, y = lat)) +
      geom_polygon(data = bg_map, aes(x = long, y = lat, group = group), 
                   colour = "black", fill = "#009E73", size = .4) +
      geom_raster(aes(alpha = density), fill = "red", interpolate = TRUE) +
      geom_contour(aes(z = density), colour = "darkred", breaks = c(0.01, breaks)) + 
      scale_alpha_continuous(range = c(0.01, 1), name = legend_name, breaks = breaks) +
      theme(axis.text    = element_text(size = 12),
        legend.text  = element_text(size = 12),
        axis.title   = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.key.size = unit(1, "cm"),
        legend.key   = element_rect(fill = "white", colour = "darkred", size = 1),
        plot.title   = element_text(size = 24, face = "bold", margin = margin(0,0,20,0)),
        plot.margin  = unit(c(1,0,1,1), units = "cm"),
        axis.title.x = element_text(margin = margin(20,0,0,0)),
        axis.title.y = element_text(margin = margin(0,20,0,0)),
        panel.background = element_rect(fill = "transparent", colour= "black"),
        panel.grid.major = element_line(colour = "grey50"),
        panel.grid.minor = element_line(colour = "grey70")) +
      scale_y_continuous(expand = c(-0.001,-0.001), breaks = seq(-80, 80, ystep)) + 
      scale_x_continuous(expand = c(-0.001,-0.001), breaks = seq(-360, 360, xstep)) +
      labs(title = title, y = "latitude (in degree N)", x = "longitude (in degree W)")
        



  if(saveMatrix) {
    save2(n2, file = paste(name, ".dat", sep = ""), saveMatrix) 
  }

  if (savePlot){
    png2(p, name, width = 1200, height = 600)
  }

  if (return_plot) {
    return(p)
  } else if (!savePlot) {
    print(p)
  }



}


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


plot_cyclonesPhaseDiagram <- function(data, savePlot = FALSE, name = NA) {

  # data is an array of dimension : variables x storms x timestep
  #   - variable 1 : parameter B, 
  #   - variable 2 : thermal wind bottom troposphere
  #   - variable 3 : thermal wind top troposphere
  #   - variable 4 : (Minimal pressure) / Wind
  #   - variable 5 : Mean radius of 34 knots wind

  # verification of the data
  if (length(dim(data)) == 2) {
    nstorm <- 1
    ntime <- dim(data)[2]
    warning("lack one dimension : may crash")
  } else if (length(dim(data)) == 3) {
    nstorm <- dim(data)[2]
    ntime <- dim(data)[3]
  } else {
    stop("Wrong dimension for the data")
  }

  # trunc some outstanding values
  param <- c("B","Vl","Vt","Max Wind","mean radius")
  extr <- rbind(c(-25,125),c(-600,300),c(-600,300),c(0,80),c(0,500))
  for (i in 1:dim(data)[1]) {
    test <- which(data[i,,] < extr[i,1], arr.ind=TRUE)
    if (length(test) > 0) {
      warning(paste("some value exceed the maximal limit for",param[i]))
      data[cbind(i,test)] <- extr[i,1]
    }
    
    test <- which(data[i,,] > extr[i,2], arr.ind=TRUE)
    if (length(test) > 0) {
      warning(paste("some value exceed the minimal limit for",param[i]))
      data[cbind(i,test)] <- extr[i,2]
    }
  }


  s_empty <- c() # if there is empty row
  for (s in 1:nstorm) {
    if (length(which(is.na(data[1, s, ]))) != length(which(is.na(data[2, s, ]))) | length(which(is.na(data[1, s, ]))) != length(which(is.na(data[3, s, ])))) {
      warnings(paste("for storm :", s,
                     "\n number of varaible 1 provided =",  ntime - length(which(is.na(data[1, s, ]))),
                     "\n number of varaible 2 provided =",  ntime - length(which(is.na(data[2, s, ]))),
                     "\n number of variable 3 provided =",  ntime - length(which(is.na(data[3, s, ])))))
    } else {
      if (length(which(is.na(data[1, s, ]))) %in% c(ntime,ntime -1 )) {
        warnings(paste("the storm", s,"do not contain any values, or only one"))
        s_empty=c(s_empty,s)
      }
    }

  }

  if (length(s_empty) > 0) {
    data <- data[,-s_empty,]
    nstorm <- dim(data)[2]
  }

  # name not given
  if (is.na(name) & (savePlot)) {
    stop("saves activated but no name provided")
  }

  theme <- theme(
             plot.title   = element_text(size = 24, face = "bold", margin = margin(0, 0, 20, 0)),
             plot.margin  = unit(c(1, 0, 1, 1), units = "cm"),
             legend.title = element_text(size = 16),
             axis.title   = element_text(size = 18),
#             axis.title.x = element_text(vjust = 0),
#             axis.title.y = element_text(vjust = 1),
             axis.text    = element_text(size = 14) )

#  scalecol <- scale_colour_gradientn(name = "Minimum Sea\nLevel Pressure\n(in hPa)", breaks = seq(970,1020,10), limits = c(950,1020), colours=c("darkred","red","violet","blue","darkblue"), values = c(0,0.3,0.5,0.7,1))

  scalecol <- scale_colour_gradientn(name = "Maximum surface\nwind (in km/h)", breaks = seq(0,80,20), limits = c(0,80), colours=c("darkred","red","violet","blue","darkblue"), values = c(1,0.7,0.5,0.3,0))

  scalesiz <- scale_size(name = "Mean radius\nof the 34 knots\nsurface wind\n(in km)", breaks = seq(0,500,100), limits = c(0,500))

  p1 <- ggplot() + theme + scalecol + scalesiz +
        scale_y_continuous(expand = c(0,0), limits = c(-25,125), breaks = seq(-20,120,20)) + 
        scale_x_continuous(expand = c(0,0), limits = c(-600,300), breaks = seq(-600,300,100)) +
        geom_vline(xintercept = 0, size = 3, color = "grey70") + 
        geom_hline(yintercept = 10, size = 3, color = "grey70") +
        geom_text(aes(x = -400,y = 75, label = "Asymetric cold-core"), size = 8, color = "grey50") +
        geom_text(aes(x = 200,y = 75, label = "Asymetric warm-core"), size = 8, color = "grey50") +
        geom_text(aes(x = -400,y = -10, label = "Symetric cold-core"), size = 8, color = "grey50") +
        geom_text(aes(x = 200,y = -10, label = "Symetric warm-core"), size = 8, color = "grey50") +
        labs(x = bquote(paste(-V[T]^L, " (Thermal wind 900-600 hPa, in .1 m/s)"))) +
        labs(y = bquote(paste(B[], "(900-600hPa storm relative thickness symetry)"))) +
        theme(axis.title.y = element_text(margin = margin(0,11.5,0,7)))


  p2 <- ggplot() + theme + scalecol + scalesiz +
        scale_y_continuous(expand = c(0,0), limits = c(-600,300), breaks = seq(-600,300,100)) + 
        scale_x_continuous(expand = c(0,0), limits = c(-600,300), breaks = seq(-600,300,100)) + 
        geom_vline(xintercept = 0, size = 3, color = "grey70") +
        geom_hline(yintercept = 0, size = 3, color = "grey70") +
        geom_text(aes(x = 200,y = 200, label = "Deep warm-core"), size = 8, color = "grey50") +
        geom_text(aes(x = 200,y = -100, label = "Moderate warm-core"), size = 8, color = "grey50") +
        geom_text(aes(x = 200,y = -400, label = "Shallow warm-core"), size = 8, color = "grey50") +
        geom_text(aes(x = -300,y = -300, label = "Deep cold-core"), size = 8, color = "grey50") +
        labs(x = bquote(paste({-V[T]^L}," (Thermal wind 900-600 hPa, in .1 m/s)"))) +
        labs(y = bquote(paste(-V[T]^U, " (Thermal wind 600-300 hPa, in .1 m/s)")))



  for ( s in 1:nstorm ) {

    nb_ts <- which(! is.na(data[1, s, ]))

    data.plot <- data.frame(B = data[1, s, nb_ts], Vl = data[2, s, nb_ts], Vu = data[3, s, nb_ts], P = data[4, s, nb_ts], Wind = data[5, s, nb_ts])


    title1 <- paste("Phase Diagram 1 for the storm",s)
    title2 <- paste("Phase Diagram 2 for the storm",s)

    p1s <- p1 +
           labs(title = title1) +
           geom_path(data = data.plot,aes(x = Vl, y = B)) + 
           geom_point(data = data.plot,aes(x = Vl, y = B, col = P, size = Wind)) +
           geom_text(data = head(data.plot,1), aes(x = Vl, y = B, label = "A"), size = 8, colour = "grey30", hjust=-0.3, vjust=-0.3) + 
           geom_text(data = tail(data.plot,1), aes(x = Vl, y = B, label = "Z"), size = 8, colour = "grey30", hjust=-0.3, vjust=-0.3)
           

    p2s <- p2 + 
           labs(title = title2) +
           geom_path(data = data.plot,aes(x = Vl, y = Vu)) + 
           geom_point(data = data.plot,aes(x = Vl, y = Vu, col = P, size = Wind)) +
           geom_text(data = head(data.plot,1), aes(x = Vl, y = Vu, label = "A"), size = 8, colour = "grey30", hjust=-0.3, vjust=-0.3) + 
           geom_text(data = tail(data.plot,1), aes(x = Vl, y = Vu, label = "Z"), size = 8, colour = "grey30", hjust=-.3, vjust=-.3)


    if (savePlot){
      png2(p1s,paste(name,"_",s,"_phase1.png", sep = ""), width = 1200, height = 600, savePlot)
      png2(p2s,paste(name,"_",s,"_phase2.png", sep = ""), width = 1200, height = 600, savePlot)
    } else {
      X11()
      print(p1s)
      X11()
      print(p2s)
    }

  }
}


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


plot_cyclonesDiagram <- function(data, name = NA) {

  # data is an array of dimension : variables x timestep
  #   - variable 1 : latitude
  #   - variable 2 : longitude
  #   - variable 3 : Wind
  #   - variable 4 : WarmCore
  #   - variable 5 : parameter B, 
  #   - variable 6 : thermal wind bottom troposphere
  #   - variable 7 : thermal wind top troposphere
  #   - variable 8 : Mean radius of 34 knots wind
  #   - variable 9 : Minimal pressure
  #   - variable 10 : IKE
  #   - variable 11 : date
  #   - variable 12 : latitude of landfall(s)
  #   - variable 13 : longitude of landfall(s)

#########################################################
#########################################################

 # Function from devtools::install_github("baptiste/egg") 
 # Used to have the same width for the phase diagrams

  gtable_frame <- function (g, width = unit(1, "null"), height = unit(1, "null"), 
      debug = FALSE) {
    panels <- g[["layout"]][grepl("panel", g[["layout"]][["name"]]), 
        ]
    ll <- unique(panels$l)
    tt <- unique(panels$t)
    fixed_ar <- g$respect
    if (fixed_ar) {
        ar <- as.numeric(g$heights[tt[1]])/as.numeric(g$widths[ll[1]])
        height <- width * ar
        g$respect <- FALSE
    }
    core <- g[seq(min(tt), max(tt)), seq(min(ll), max(ll))]
    top <- g[seq(1, min(tt) - 1), seq(min(ll), max(ll))]
    bottom <- g[seq(max(tt) + 1, nrow(g)), seq(min(ll), max(ll))]
    left <- g[seq(min(tt), max(tt)), seq(1, min(ll) - 1)]
    right <- g[seq(min(tt), max(tt)), seq(max(ll) + 1, ncol(g))]
    fg <- nullGrob()
    if (length(left)) {
        lg <- gtable::gtable_add_cols(left, unit(1, "null"), 
            0)
        lg <- gtable::gtable_add_grob(lg, fg, 1, l = 1)
    }
    else {
        lg <- fg
    }
    if (length(right)) {
        rg <- gtable_add_cols(right, unit(1, "null"))
        rg <- gtable_add_grob(rg, fg, 1, l = ncol(rg))
    }
    else {
        rg <- fg
    }
    if (length(top)) {
        tg <- gtable_add_rows(top, unit(1, "null"), 0)
        tg <- gtable_add_grob(tg, fg, t = 1, l = 1)
    }
    else {
        tg <- fg
    }
    if (length(bottom)) {
        bg <- gtable_add_rows(bottom, unit(1, "null"), -1)
        bg <- gtable_add_grob(bg, fg, t = nrow(bg), l = 1)
    }
    else {
        bg <- fg
    }
    grobs = list(fg, tg, fg, lg, core, rg, fg, bg, fg)
    widths <- unit.c(sum(left$widths), width, sum(right$widths))
    heights <- unit.c(sum(top$heights), height, sum(bottom$heights))
    all <- gtable_matrix("all", grobs = matrix(grobs, ncol = 3, 
        nrow = 3, byrow = TRUE), widths = widths, heights = heights)
    if (debug) {
        hints <- rectGrob(gp = gpar(fill = NA, lty = 2, lwd = 0.2))
        tl <- expand.grid(t = 1:3, l = 1:3)
        all <- gtable::gtable_add_grob(all, replicate(9, hints, 
            simplify = FALSE), t = tl$t, l = tl$l, z = Inf, name = "debug")
    }
    all[["layout"]][5, "name"] <- "panel"
    if (fixed_ar) 
        all$respect <- TRUE
    all
  }


#########################################################
#########################################################



  # verification of the data
  if (length(dim(data)) == 2) {
    ntime <- dim(data)[2]
  } else {
    stop("Wrong dimension for the data")
  }

  # Some formatiing
  date <- as.POSIXct(data[11,], origin = as.Date("1970-01-01"))
  date2 <- (data[11,] - data[11,1]) / 86400

  date_lf <- (data[14,] - data[11,1]) / 86400

  if( (max(data[2,], na.rm = TRUE) - min(data[2,], na.rm = TRUE)) > 300) {
    nb_neg <- which(data[2,] < 0)
    data[2, nb_neg] <- data[2, nb_neg] + 360

    nb_neg_lf <- which(data[13, ] < 0)
    data[13, nb_neg_lf] <- data[13, nb_neg_lf] +360
  }


  rangeB <- c(-3,10)
  rangeV <- c(-40,40)
  rangeR <- c(0,750)


  data[5,] <- scales::squish(data[5,], range = rangeB)
  data[6,] <- scales::squish(data[6,], range = rangeV)
  data[7,] <- scales::squish(data[7,], range = rangeV)
  data[8,] <- scales::squish(data[8,], range = rangeR)
  
# trunc some outstanding values  !!! +Wind (twice) + Minimal pressure + IKE
#  param <- c("B","Vl","Vt","Max Wind","mean radius")
#  extr <- rbind(c(-25,125),c(-600,300),c(-600,300),c(0,80),c(0,500))
#  for (i in 1:dim(data)[1]) {
#    test <- which(data[i,,] < extr[i,1], arr.ind=TRUE)
#    if (length(test) > 0) {
#      warning(paste("some value exceed the maximal limit for",param[i]))
#      data[cbind(i,test)] <- extr[i,1]
#    }
#    
#    test <- which(data[i,,] > extr[i,2], arr.ind=TRUE)
#    if (length(test) > 0) {
#      warning(paste("some value exceed the minimal limit for",param[i]))
#      data[cbind(i,test)] <- extr[i,2]
#    }
#  }


  # name not given
  #if (is.na(name) & (savePlot)) {
  #  stop("saves activated but no name provided")
  #}


 ############
 ## GGPLOT ##
 ############

  # Main Theme

  theme <- theme(
             plot.title   = element_text(size = 12, face = "bold", margin = margin(0, 0, 20, 0)),
             plot.margin  = unit(c(1, 0, 1, 1), units = "cm"),
             legend.title = element_text(size = 10),
             legend.key   = element_rect(fill = "white"),
             legend.text  = element_text(size = 9),
             axis.title   = element_text(size = 10),
             axis.text    = element_text(size = 9),
             panel.background = element_rect(fill = "white", colour= "black", size = 1),
             panel.grid.major = element_line(colour = "grey70"),
             panel.grid.minor = element_line(colour = "grey90") )


 ## Phase Diagramme
 ##################

  scalecol <- scale_fill_gradientn(name = "Maximum surface\nwind (in knots)", breaks = seq(0,80,20), limits = c(0,80), colours=c("darkred","red","violet","blue","darkblue"), values = c(1,0.7,0.5,0.3,0), oob = scales::squish)
  scalecol_F <- scale_fill_gradientn(name = "Maximum surface\nwind (in knots)", breaks = seq(0,80,20), limits = c(0,80), colours=c("darkred","red","violet","blue","darkblue"), values = c(1,0.7,0.5,0.3,0), guide = FALSE, oob = scales::squish)


  scalesiz <- scale_radius(name = "Mean radius\nof the 34 knots  \nsurface wind\n(in km)", breaks = c(100,200,300,500,750), range = c(.1,10), limits = rangeR)
  scalesiz_F <- scale_radius(name = "Mean radius\nof the 34 knots\nsurface wind\n(in km)", breaks = c(100,200,300,500,750), range = c(.1,10), guide = FALSE, limits = rangeR)

  p1 <- ggplot() + theme + scalecol + scalesiz_F +
        scale_y_continuous(expand = c(0,0), limits = rangeB, breaks = pretty(rangeB)) + 
        scale_x_continuous(expand = c(0,0), limits = rangeV, breaks = pretty(rangeV, n = 10)) +
        geom_vline(xintercept = 0, size = 2.5, color = "grey70") + 
        geom_hline(yintercept = 1, size = 2.5, color = "grey70") +
        geom_text(aes(x = -20,y = 7, label = "Asymetric cold-core"), size = 3, color = "grey30") +
        geom_text(aes(x = 20,y = 7, label = "Asymetric warm-core"), size = 3, color = "grey30") +
        geom_text(aes(x = -20,y = -1, label = "Symetric cold-core"), size = 3, color = "grey30") +
        geom_text(aes(x = 20,y = -1, label = "Symetric warm-core"), size = 3, color = "grey30") +
        labs(x = bquote(paste(-V[T]^L, " (Thermal wind 900-600 hPa, in m/s)"))) +
        labs(y = "Beta Parameter 900-600hPa\n(storm relative thickness symetry)") +
        theme(axis.title.y = element_text(margin = margin(0,11.5,0,7)))


  p2 <- ggplot() + theme + scalecol_F + scalesiz +
        scale_y_continuous(expand = c(0,0), limits = rangeV, breaks = pretty(rangeV, n = 10)) + 
        scale_x_continuous(expand = c(0,0), limits = rangeV, breaks = pretty(rangeV, n = 10)) + 
        geom_vline(xintercept = 0, size = 2.5, color = "grey70") +
        geom_hline(yintercept = 0, size = 2.5, color = "grey70") +
        geom_text(aes(x = 15,y = 25, label = "Deep warm-core"), size = 3, color = "grey30") +
        geom_text(aes(x = 15,y = -5, label = "Moderate warm-core"), size = 3, color = "grey30") +
        geom_text(aes(x = 15,y = -35, label = "Shallow warm-core"), size = 3, color = "grey30") +
        geom_text(aes(x = -25,y = -25, label = "Deep cold-core"), size = 3, color = "grey30") +
        labs(x = bquote(paste({-V[T]^L}," 900-600 hPa (Thermal wind, in m/s)"))) +
        labs(y = bquote(paste(-V[T]^U, " 600-300 hPa (Thermal wind, in m/s)")))

  nb_ts <- which(! is.na(data[1, ]))

  data.plot <- data.frame(B = data[5, nb_ts], Vb = data[6, nb_ts], Vt = data[7, nb_ts], Wind = data[3, nb_ts], Rad = data[8, nb_ts])

    p1 <- p1 +
           labs(title = NULL) +
           geom_path(data = data.plot,aes(x = Vb, y = B)) + 
           geom_point(data = data.plot,aes(x = Vb, y = B, fill = Wind, size = Rad), shape = 21, col = "black") +
           geom_text(data = head(data.plot,2), aes(x = Vb, y = B, label = "A"), size = 6, colour = "black", hjust="outward", vjust="outward") + 
           geom_text(data = tail(data.plot,1), aes(x = Vb, y = B, label = "Z"), size = 6, colour = "black", hjust="outward", vjust="outward")
           

    p2 <- p2 + 
           labs(title = NULL) +
           geom_path(data = data.plot,aes(x = Vb, y = Vt)) + 
           geom_point(data = data.plot,aes(x = Vb, y = Vt, fill = Wind, size = Rad), shape = 21, col = "black") +
           geom_text(data = head(data.plot,2), aes(x = Vb, y = Vt, label = "A"), size = 6, colour = "black", hjust="outward", vjust="outward") + 
           geom_text(data = tail(data.plot,1), aes(x = Vb, y = Vt, label = "Z"), size = 6, colour = "black", hjust="outward", vjust="outward")



  ## Track map
  ############

    lon_min <- round(min(data[2,], na.rm = T) - 10, -1)
    lon_max <- round(max(data[2,], na.rm = T) + 10, -1)
    lat_min <- min(80,round(min(data[1,], na.rm = T) - 10, -1))
    lat_max <- min(80,round(max(data[1,], na.rm = T) + 10, -1))
    map <- Map(lon_min, lon_max, lat_min, lat_max, 5)

    p3 <- plot_track(data = array(data[1:3,], dim= c(3,1,ntime)), type = "wind", plotnum = FALSE, limits = c(lon_min, lon_max, lat_min, lat_max), title = NULL, legend = "Category", return_plot = TRUE, bg_map = map) +
          
           coord_equal() + theme +
           geom_path(data = data.frame(lon = c(lon_min, lon_min, lon_max, lon_max, lon_min),
                                       lat = c(lat_min, lat_max, lat_max, lat_min, lat_min)),
                     aes(x = lon, y = lat)) +

           geom_point(data = data.frame(lon = data[13, ], lat = data[12, ]),
                      aes(x = lon, y = lat), size = 3, shape = 20) +

           geom_point(data = data.frame(lon = data[13, ], lat = data[12, ]),
                      aes(x = lon, y = lat), size = 5, shape = 4, stroke = 1) +

           geom_text(data = data.frame(lon = c(data[2,1], tail(data[2,nb_ts], 1)),
                                       lat = c(data[1,1], tail(data[1,nb_ts], 1))),
                     aes(x = lon, y = lat), label = c("A", "Z"),
                     size = 6, colour = "black", hjust=-0.3, vjust=-0.3)


  ## White map for the merging

    p4 <- ggplot() + coord_equal() +
           geom_path(data = data.frame(lon = c(lon_min, lon_min, lon_max, lon_max, lon_min),
                                       lat = c(lat_min, lat_max, lat_max, lat_min, lat_min)),
                     aes(x = lon, y = lat), col = "white") +
           theme(panel.background = element_rect(fill = "white"),
                 axis.text = element_text (colour = "white"),
                 axis.title = element_text (colour = "white"),
                 axis.ticks = element_line (colour = "white"))



  ## Merging (to get same width for both phase diagram)
  ##########

  # For map (needed for grid.arrange)
  grobs <- lapply(list(p3, p4) , ggplotGrob)
  
  widths <- lapply(c(1, 1), unit, "null")
  heights <- lapply(c(1, 1), unit, "null")
 
  fg <- mapply(gtable_frame, g = grobs, width = widths, height = heights, SIMPLIFY = FALSE)
  splits <- cut(seq_along(fg), 2, labels = seq_len(2))
  spl <- split(fg, splits)
  rows <- lapply(spl, function(r) do.call(cbind.gtable, r))
  col1 <- do.call(rbind.gtable, rows)
     

  # For phase diagram
  grobs <- lapply(list(p1, p2), ggplotGrob)
  
  widths <- lapply(c(1, 1), unit, "null")
  heights <- lapply(c(1, 1), unit, "null")
 
  fg <- mapply(gtable_frame, g = grobs, width = widths, height = heights, SIMPLIFY = FALSE)
  splits <- cut(seq_along(fg), 2, labels = seq_len(2))
  spl <- split(fg, splits)
  rows <- lapply(spl, function(r) do.call(cbind.gtable, r))
  col2 <- do.call(rbind.gtable, rows)
    


 #####################
 ## Saving the Plot ##
 #####################

  setEPS()
  postscript(name, width = 16, height = 8)
  par(mfrow=c(2,2))

  plot(1)
  plot(1)

  grid.arrange(col1, col2, ncol = 2,
             top = textGrob(name_sel[i], gp = gpar(fontsize = 30, font = 2), vjust = 1))


 ########## => The time series plot (not ggplot because of multiple axis)
 ############## ############## ############## ############## ############## 

 ######################
 # Defined axis ticks #
 ######################

  at1 <- pretty(c(0, max(data[10,nb_ts])))
  at2 <- pretty(c(min(data[9,nb_ts]), max(data[9,nb_ts])))
  at3 <- pretty(c(0, max(data[3,nb_ts])))
  
  # To get the ticks for each axis at the same level
  while (length(at1) < max(length(at2), length(at3))) {
    at1 <- c(at1, tail(at1,1) + at1[2] - at1[1])
  }

  while (length(at2) < max(length(at1), length(at3))) {
    at2 <- c(head(at2,1) + at2[1] - at2[2], at2)
  }

  while (length(at3) < max(length(at1), length(at2))) {
    at3 <- c(at3, tail(at3,1) + at3[2] - at3[1])
  }

  # X axis
  atx <- pretty(c(min(date2[nb_ts]), max(date2[nb_ts])))

  # Preparinf the data for ploting WarmCore
  wc_T <- rep(NA, length(data[4,nb_ts]))
  wc_T[which(data[4,nb_ts] == 1)] <- 3*at3[1]/2 - at3[2]/2


  wc_F <- rep(NA, length(data[4,nb_ts]))
  wc_F[which(data[4,nb_ts] == 0)] <- 3*at3[1]/2 - at3[2]/2



 ##################
 # Before ploting #
 ##################

  par(xpd=TRUE) # To allow drawing outside of the plot (for legend and Warmcore)
  par(mar=c(6, 9, 2, 9) + 0.1) # Change the Margin (for double axis and legend)

  # White plot of the first plot to get the dimension
  plot(date2[nb_ts], data[10,nb_ts], col = "white",
       ylim = c(min(at1), max(at1)), xlim = c(min(atx), max(atx)),
       main = NULL, xlab = "", ylab = "", yaxt = 'n', xaxt = 'n',  bty="n")

  # Plot the background grid
  for (i in 1:length(atx)) {
    lines(c(atx[i], atx[i]), c(at1[1],tail(at1, 1)), col = "grey")
  }

  for (i in 1:length(at2)) {
    lines(c(atx[1], tail(atx, 1)), c(at1[i],at1[i]), col = "grey")
  }

 ##############
 # First Plot #
 ##############
  # new plot
  par(new = T)

  plot(date2[nb_ts], data[10,nb_ts],
       ylim = c(min(at1), max(at1)), xlim = c(min(atx), max(atx)),
       main = NULL, xlab = "", ylab = "",
       col = "darkgreen", type = "l", lwd = 2,
       yaxt = 'n', xaxt = 'n', # Remove both axis
       bty="n") # Remove the box around the plot

  # First Y axis
  axis(side = 4, at = at1, col = "darkgreen", lwd = 2, cex.axis = .9)
  mtext(side = 4, line = 2.25, 'IKE (TJ)')

 ###############
 # Second Plot #
 ###############
  # new plot
  par(new = T)

  plot(date2[nb_ts], data[9,nb_ts],
  ylim = c(min(at2), max (at2)), xlim = c(min(atx), max(atx)),
  col = "orange2",  type = "l", lwd = 2, bty="n",
  main = NULL, xlab = "", ylab = "",
  xaxt = 'n', yaxt = 'n')

  # Second Y axis
  axis(2, line = 4.5, at = at2, col = "orange2", lwd = 2, cex.axis = .9)
  mtext(side = 2, line = 6.75, 'Minimal Pressure (hPa)')

 ###############
 # Third Plot #
 ###############
  # new plot
  par(new = T)

  plot(date2[nb_ts], data[3,nb_ts],
  ylim = c(min(at3), max (at3)), xlim = c(min(atx), max(atx)),
  col = "mediumpurple",  type = "l", lwd = 2, bty="n",
  main = NULL, xlab = "", ylab = "",
  xaxt = 'n', yaxt = 'n')

  # Third Y axis
  axis(2, line = 0, at = at3, col = "mediumpurple", lwd = 2, cex.axis = .9)
  mtext(side = 2, line = 2.25, 'Maximum Wind (knots)')

  # X axis
  axis(1, line = 2, at = atx, lwd = 2, cex.axis = .9)
  mtext(side = 1, line = 4, 'Days after first detection')

 #################
 # Plot Warmcore #
 #################

  lines(date2[nb_ts], wc_T, col = "firebrick3", lwd = 8)
  lines(date2[nb_ts], wc_F, col = "dodgerblue3", lwd = 8)

  points(date2[nb_ts], wc_T, col = "firebrick3", lwd = 3, pch = 20)
  points(date2[nb_ts], wc_F, col = "dodgerblue3", lwd = 3, pch = 20)

  legend(tail(atx,1)*1.15, max(at3), c("True", "False", "Landfall"), lty = c(1,1,0), pch= c(NA,NA,"|"), lwd = c(8,8,8), col = c("firebrick3", "dodgerblue3", "black"), title = "Warmcore :", bty = "n", yjust = 1)

 ##################
 # Landfall lines #
 ##################

  for (i in 1:length(date_lf)) {
  if ( ! is.na(date_lf[i])) {
      points(date_lf[i], 3/2*at3[1] - at3[2]/2, col = "black", lwd = 8, pch = "|", cex = 2)
  }}

dev.off()


}






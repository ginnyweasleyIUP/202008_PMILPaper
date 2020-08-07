##############################################################
## extract Data for Cave Site from surrounding grid boxes ####
##############################################################

extract_gridboxes <- function(lon_cave, lat_cave, d.lon = 3.75, d.lat = 2.5){
  result_list <- list()
  
  if(lon_cave<0){lon_cave = 360+lon_cave}
  
  #corners of the square around lon_cave, lat_cave with dimensions of d.lon*d.lat
  corner <- list(
    E1_lon = lon_cave+d.lon/2, E1_lat = lat_cave+d.lat/2,
    E2_lon = lon_cave-d.lon/2, E2_lat = lat_cave+d.lat/2,
    E3_lon = lon_cave+d.lon/2, E3_lat = lat_cave-d.lat/2,
    E4_lon = lon_cave-d.lon/2, E4_lat = lat_cave-d.lat/2
  )
  
  #assuming that longitudes are in 0 to 360 and not from -180 to 180
  if(corner$E1_lon<0){corner$E1_lon = 360 + corner$E1_lon}
  if(corner$E2_lon<0){corner$E2_lon = 360 + corner$E2_lon}
  if(corner$E3_lon<0){corner$E3_lon = 360 + corner$E3_lon}
  if(corner$E4_lon<0){corner$E4_lon = 360 + corner$E4_lon}
  
  #we desire that lon_real>long_grid and lat_real<lat_grid
  corner$E1_lon_pos <- which.min(abs(DATA_past1000$SIM_mean$lon-corner$E1_lon))
  if(DATA_past1000$SIM_mean$lon[corner$E1_lon_pos]>corner$E1_lon){corner$E1_lon_pos = corner$E1_lon_pos -1}
  corner$E1_lat_pos <- which.min(abs(DATA_past1000$SIM_mean$lat-corner$E1_lat))
  if(DATA_past1000$SIM_mean$lat[corner$E1_lat_pos]>corner$E1_lat){corner$E1_lat_pos = corner$E1_lat_pos +1}
  
  corner$E2_lon_pos <- which.min(abs(DATA_past1000$SIM_mean$lon-corner$E2_lon))
  if(DATA_past1000$SIM_mean$lon[corner$E2_lon_pos]>corner$E2_lon){corner$E2_lon_pos = corner$E2_lon_pos -1}
  corner$E2_lat_pos <- which.min(abs(DATA_past1000$SIM_mean$lat-corner$E2_lat))
  if(DATA_past1000$SIM_mean$lat[corner$E2_lat_pos]>corner$E2_lat){corner$E2_lat_pos = corner$E2_lat_pos +1}

  corner$E3_lon_pos <- which.min(abs(DATA_past1000$SIM_mean$lon-corner$E3_lon))
  if(DATA_past1000$SIM_mean$lon[corner$E3_lon_pos]>corner$E3_lon){corner$E3_lon_pos = corner$E3_lon_pos -1}
  corner$E3_lat_pos <- which.min(abs(DATA_past1000$SIM_mean$lat-corner$E3_lat))
  if(DATA_past1000$SIM_mean$lat[corner$E3_lat_pos]>corner$E3_lat){corner$E3_lat_pos = corner$E3_lat_pos +1}

  corner$E4_lon_pos <- which.min(abs(DATA_past1000$SIM_mean$lon-corner$E4_lon))
  if(DATA_past1000$SIM_mean$lon[corner$E4_lon_pos]>corner$E4_lon){corner$E4_lon_pos = corner$E4_lon_pos -1}
  corner$E4_lat_pos <- which.min(abs(DATA_past1000$SIM_mean$lat-corner$E4_lat))
  if(DATA_past1000$SIM_mean$lat[corner$E4_lat_pos]>corner$E4_lat){corner$E4_lat_pos = corner$E4_lat_pos +1}

  ratio <- list()
  
  ratio$E1 <- (corner$E1_lon - DATA_past1000$SIM_mean$lon[corner$E1_lon_pos])*(corner$E1_lat - DATA_past1000$SIM_mean$lat[corner$E1_lat_pos])/(d.lon*d.lat)
  ratio$E2 <- (DATA_past1000$SIM_mean$lon[corner$E2_lon_pos]+d.lon-corner$E2_lon)*(corner$E2_lat - DATA_past1000$SIM_mean$lat[corner$E2_lat_pos])/(d.lon*d.lat)
  ratio$E3 <- (corner$E3_lon - DATA_past1000$SIM_mean$lon[corner$E3_lon_pos])*(DATA_past1000$SIM_mean$lat[corner$E3_lat_pos-1]-corner$E3_lat)/(d.lon*d.lat)
  ratio$E4 <- (DATA_past1000$SIM_mean$lon[corner$E4_lon_pos]+d.lon-corner$E4_lon)*(DATA_past1000$SIM_mean$lat[corner$E4_lat_pos-1]-corner$E4_lat)/(d.lon*d.lat)
  #ratio$E1+ratio$E2+ratio$E3+ratio$E4

    result_list <- c(corner, ratio)
  
  return(result_list)
}
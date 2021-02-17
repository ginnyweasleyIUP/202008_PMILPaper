##############################################################
## extract Data for Cave Site from surrounding grid boxes ####
##############################################################

extract_gridboxes <- function(lon_cave, lat_cave, d.lon = 3.75, d.lat = 2.5){
  result_list <- list()
  lon_list <- seq(from = 0, to = 360-d.lon, by = d.lon)
  lat_list <- seq(from = 90, to = -90, by = -d.lat)
  
  if(lon_cave<0){lon_cave = 360+lon_cave}
  
  # 1) search for top right node -> it's the one, that we'll orient by. It is node Q22
  # lon:
  pos = which.min(abs(lon_list-lon_cave))
  if(lon_cave - lon_list[pos] > 0){
    pos = pos + 1
  }
  if(pos>length(lon_list)){
    lon_ref = 1
  }else{
    lon_ref = pos
  }
  
  
  # lat:
  pos = which.min(abs(lat_list-lat_cave))
  if(lat_cave - lat_list[pos] > 0){
    pos = pos - 1
  }
  lat_ref = pos
  
  #Q22 is now marked as lon_list[lon_ref], lat_list[lat_ref] --> a value at that point would be e.g. TEMP[lon_ref, lat_ref, ...]
  #Q11 with lon_list[lon_ref-1], lat_list[lat_ref-1]
  #Q12 with lon_list[lon_ref-1], lat_list[lat_ref]
  #Q21 with lon_list[lon_ref], lat_list[lat_ref]
  
  #x1 = lon_cave - d.lon/2
  #x2 = lon_cave + d.lon/2
  #y1 = lat_cave - d.lat/2
  #y2 = lat_cave + d.lat/2
  
  #ratios should work even with the 360->0 switch because the ratios stay the same!
  # (x,y) is middle between the 4 neighboring cells. So it is always at
  y = lat_list[lat_ref] - d.lat/2
  if(lon_ref ==1){
    x = lon_list[length(lon_list)]+d.lon-d.lon/2
  }else{
    x = lon_list[lon_ref] - d.lon/2
  }
  
  ratio <- list()
  #             x2                    x   y2                      y
  ratio$Q22 = ((lon_cave + d.lon/2) - x)*((lat_cave + d.lat/2) - y)/(d.lon*d.lat)
  #             x2                    x   y   y1
  ratio$Q21 = ((lon_cave + d.lon/2) - x)*(y - (lat_cave - d.lat/2))/(d.lon*d.lat)
  #            x    x1                      y2                   y
  ratio$Q12 = (x - (lon_cave - d.lon/2))*((lat_cave + d.lat/2) - y)/(d.lon*d.lat)
  #            x    x1                    y   y1
  ratio$Q11 = (x - (lon_cave - d.lon/2))*(y - (lat_cave - d.lat/2))/(d.lon*d.lat)
  
  #ratios
  
  
  if(lon_ref == 1){
    corners <- list(Q11_lat = lat_ref+1, Q12_lat = lat_ref+1, Q21_lat = lat_ref, Q22_lat = lat_ref,
                    Q11_lon = length(lon_list), Q12_lon = lon_ref, Q21_lon = length(lon_list), Q22_lon = lon_ref)
  }else{
    corners <- list(Q11_lat = lat_ref+1, Q12_lat = lat_ref+1, Q21_lat = lat_ref, Q22_lat = lat_ref,
                    Q11_lon = lon_ref-1, Q12_lon = lon_ref, Q21_lon = lon_ref-1, Q22_lon = lon_ref)
  }
  
  
  result_list <- c(corners, ratio)
  
  return(result_list)
  
}





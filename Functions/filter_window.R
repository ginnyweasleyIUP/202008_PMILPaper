filter_window <- function(ts_filtered){
  
  mdn = median(ts_filtered, na.rm = T)
  for(ii in 1:length(ts_filtered/2)){
    if(ts_filtered[ii]-mdn<0.5){break}
  }
  for(jj in 1:length(ts_filtered/2)){
    if(rev(ts_filtered)[jj]-mdn<0.5){break}
  }

  return(window(ts_filtered,start = index(ts_filtered)[ii], end = index(ts_filtered)[length(ts_filtered)-(jj+1)]))
  
}

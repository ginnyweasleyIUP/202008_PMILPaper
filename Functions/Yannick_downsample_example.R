#################################################
## Down-sample Cave Data

source("Functions/SubsampleTimeseriesBlock_highresNA.R")

for(ii in 1:length(entity_info$entity_id)){
  print(ii)
  
  #pass hier auf, wierum deine Daten sind. Start und End muss leider immer in positiver Richtung sein. 
  #Je nachdem wierum deine proxy Daten und deren Zeitachse sotriert ist musst du sie vielleicht mit rev(data) umdrehen. 
  #Kontrolliere am Besten ob du es richtig gemacht hast mit einem Beispiel Record und plotte dann sim, sim_ds und record als Zeitreihe, dann kannst du es abschätzen
  ds <- SubsampleTimeseriesBlock_highresNA(ts(data = proxy_data[[ii]]$TEMP,
                                                     start = 850,
                                                     end = 1850),
                                            proxy_data[[ii]]$interp_age)
  data <- matrix(c(proxy_data[[ii]]$interp_age, ds, ncol= 2))
  colnames(data) = c("interp_age", 
                     "TEMP_ds")
}

remove(SubsampleTimeseriesBlock_highresNA)
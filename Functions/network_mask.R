network_mask <-function(list, method){
  mask = matrix(logical(length = sum(mask_spec)*sum(mask_spec)), ncol = sum(mask_spec))
  
  if(method == "site"){
    for(site in list$site_id){
      entity_total = DATA_past1000$CAVES$entity_info$entity_id[mask_spec]
      entity_list = DATA_past1000$CAVES$entity_info$entity_id[DATA_past1000$CAVES$entity_info$site_id == site]
      
      index.c = length(entity_list)
      for(ii in 1:length(entity_list)){
        index.c[ii] = which(entity_list[ii] == entity_total)
      }
      
      for(ii in 1:(length(index.c)-1)){
        for(jj in ii:length(index.c)){
          mask[index.c[ii],index.c[jj]] = mask[index.c[jj],index.c[ii]] = mask[index.c[ii],index.c[ii]] = mask[index.c[jj],index.c[jj]] = T
        }
      }
    }
  }
  if(method == "gridbox"){
    for(gridbox in list$gridbox_id){
      entity_total = DATA_past1000$CAVES$entity_info$entity_id[mask_spec]
      entity_list = DATA_past1000$CAVES$gridbox_list$entity_id[DATA_past1000$CAVES$gridbox_list$gridbox_id == gridbox]
      
      index.c = length(entity_list)
      for(ii in 1:length(entity_list)){
        index.c[ii] = which(entity_list[ii] == entity_total)
      }
      
      for(ii in 1:(length(index.c)-1)){
        for(jj in ii:length(index.c)){
          mask[index.c[ii],index.c[jj]] = mask[index.c[jj],index.c[ii]] = mask[index.c[ii],index.c[ii]] = mask[index.c[jj],index.c[jj]] = T
        }
      }
    }
  }
  if(method == "cluster"){
    entity_total = DATA_past1000$CAVES$entity_info$entity_id[mask_spec]
    entity_list = list$entity_id
      
    index.c = length(entity_list)
    for(ii in 1:length(entity_list)){
      index.c[ii] = which(entity_list[ii] == entity_total)
    }
      
    for(ii in 1:(length(index.c)-1)){
      for(jj in ii:length(index.c)){
        mask[index.c[ii],index.c[jj]] = mask[index.c[jj],index.c[ii]] = mask[index.c[ii],index.c[ii]] = mask[index.c[jj],index.c[jj]] = T
      }
    }
  }
  
  return(mask)
}


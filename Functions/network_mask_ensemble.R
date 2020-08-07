network_mask_ensemble <-function(list, method){
  entity_total = c(14,21,33,48,51,76,85,93,94,95,97,111,113,117,118,123,128,137,147,165,172,178,179,185,187,
                   201,202,209,212,222,226,238,240,242,253,278,286,289,294,296,305,326,329,330,335,351,361,378,390,395,399,
                   420,422,430,435,437,442,443,447,448,461,464,496,498,506,514,522,528,538,539,541,544,546,547,548,560,564,573,577,588,589,595,598,
                   613,620,621,623,672,673)
  
  mask = matrix(logical(length = length(entity_total)*length(entity_total)), ncol = length(entity_total))
  
  
  if(method == "site"){
    for(site in list$site_id){
      entity_list = DATA_past1000$CAVES$entity_info$entity_id[DATA_past1000$CAVES$entity_info$site_id == site]
      entity_list = entity_list[entity_list %in% entity_total]
      if(is_empty(entity_list)){next}
      
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
      entity_list = DATA_past1000$CAVES$gridbox_list$entity_id[DATA_past1000$CAVES$gridbox_list$gridbox_id == gridbox]
      entity_list = entity_list[entity_list %in% entity_total]
      if(is_empty(entity_list)){next}
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
    entity_list = list$entity_id
    entity_list = entity_list[entity_list %in% entity_total]
    if(is_empty(entity_list)){next}  
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


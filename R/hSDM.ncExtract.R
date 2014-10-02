hSDM.ncExtract <- function(files,what=c("eval","coef","autocor","data")){
  ## collect various metadata from a list of hSDM netCDF output
  foreach(f=files,.combine=rbind.data.frame)%do%{
    nc=nc_open(f,write=F)
    ## metadata
    model=ncatt_get(nc, 0, attname="model", verbose=FALSE )$value
    modelname=ncatt_get(nc, 0, attname="modelname", verbose=FALSE )$value
    species=ncatt_get(nc, 0, attname="species", verbose=FALSE )$value        
    if(what=="eval"){
      evaluation=data.frame(t(ncvar_get(nc, "evaluation")))
      colnames(evaluation)=strsplit(ncatt_get(nc, "evaluation", attname="colnames", verbose=FALSE )$value,",")[[1]]
      rownames(evaluation)=strsplit(ncatt_get(nc, "evaluation", attname="rownames", verbose=FALSE )$value,",")[[1]]
      return(cbind.data.frame(species,modelname,model,evaluation))
    }
    if(what=="coef"){
      parameters=data.frame(t(ncvar_get(nc, "parameters")))
      colnames(parameters)=strsplit(ncatt_get(nc, "parameters", attname="colnames", verbose=FALSE )$value,",")[[1]]
      rownames(parameters)=strsplit(ncatt_get(nc, "parameters", attname="rownames", verbose=FALSE )$value,",")[[1]]
      return(cbind.data.frame(species,modelname,model,parameters))
    }
    if(what=="autocor"){
      ac=data.frame(t(ncvar_get(nc, "ac")))
      colnames(ac)=strsplit(ncatt_get(nc, "ac", attname="colnames", verbose=FALSE )$value,",")[[1]]
      return(cbind.data.frame(species,modelname,model,ac))
    }
    if(what=="data"){
      data=data.frame(t(ncvar_get(nc, "data")))
      colnames(data)=strsplit(ncatt_get(nc, "data", attname="colnames", verbose=FALSE )$value,",")[[1]]
      ## Remove this!
      data=data[,!grepl("fit",colnames(data))]
      return(cbind.data.frame(species,modelname,model,data))
    }
    nc_close(nc)
  }
}
hSDM.ncExtract <- function(files,what=c("eval","coef","autocor","predictions","envdata")){
  ## collect various metadata from a list of hSDM netCDF output
  foreach(f=files,.combine=rbind.data.frame)%do%{
    nc=nc_open(f,write=F)
    ## metadata
    modeltype=ncatt_get(nc, 0, attname="modeltype", verbose=FALSE )$value
    modelname=ncatt_get(nc, 0, attname="modelname", verbose=FALSE )$value
    species=ncatt_get(nc, 0, attname="species", verbose=FALSE )$value
    suitability=ncatt_get(nc,0,attname="suitability",verbose=FALSE)$value
    observability=ncatt_get(nc,0,attname="observability",verbose=FALSE)$value
    meta=c(species=species,modeltype=modeltype,modelname=modelname,suitability=suitability,observability=observability)
    
    if(what=="eval"){
      evaluation=data.frame(t(ncvar_get(nc, "evaluation")))
      colnames(evaluation)=strsplit(ncatt_get(nc, "evaluation", attname="colnames", verbose=FALSE )$value,",")[[1]]
      rownames(evaluation)=strsplit(ncatt_get(nc, "evaluation", attname="rownames", verbose=FALSE )$value,",")[[1]]
      return(cbind.data.frame(t(meta),evaluation))
    }
    if(what=="coef"){
      parameters=data.frame(t(ncvar_get(nc, "parameters")))
      colnames(parameters)=strsplit(ncatt_get(nc, "parameters", attname="colnames", verbose=FALSE )$value,",")[[1]]
      rownames(parameters)=strsplit(ncatt_get(nc, "parameters", attname="rownames", verbose=FALSE )$value,",")[[1]]
      return(cbind.data.frame(t(meta),parameters))
    }
    if(what=="autocor"){
      ac=data.frame(t(ncvar_get(nc, "ac")))
      colnames(ac)=strsplit(ncatt_get(nc, "ac", attname="colnames", verbose=FALSE )$value,",")[[1]]
      return(cbind.data.frame(t(meta),ac))
    }
    if(what=="data"){
      data=data.frame(t(ncvar_get(nc, "data")))
      colnames(data)=strsplit(ncatt_get(nc, "data", attname="colnames", verbose=FALSE )$value,",")[[1]]
      ## Remove this!
      data=data[,!grepl("fit",colnames(data))]
      return(cbind.data.frame(t(meta),data))
    }
    if(what=="predictions"){
      r=raster(f,varname="p")
      cell=raster(f,varname="cell")
      rd=cbind(coordinates(r),pred=values(r),cell=values(cell))
      return(cbind.data.frame(t(meta),rd))
    }
    if(what=="envdata"){
      ## extract data from an data "input" object 
      r=stack(f,varname="envdata")
      names(r)=strsplit(ncatt_get(nc, "envdata", attname="names", verbose=FALSE )$value,",")[[1]]
      cell=raster(f,varname="cell")
      rd=cbind(coordinates(r),pred=values(r),cell=values(cell))
      return(cbind.data.frame(species,rd))
    }   
    nc_close(nc)
  }
}
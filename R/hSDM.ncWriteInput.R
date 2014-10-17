hSDM.ncWriteInput<-function(envdata,points,ncfile,meta=NULL,overwrite=F,verbose=T){
  ## check if file already exists
  if(file.exists(ncfile)&!overwrite) stop("File exists.  Specify overwrite=T to proceed.") 
  if(file.exists(ncfile)&overwrite) file.remove(ncfile)
  
  
  ## scale environmental data
  cmeans=cellStats(envdata,"mean")
  csd=cellStats(envdata,"sd")
  senv=raster::scale(envdata)
  
  ## add cell id to facilitate linking points to raster
  cell=envdata[[1]]
  raster::values(cell)=1:ncell(cell)
  names(cell)="cell"
  
  ## rasterize points
  presences=rasterize(points,envdata,fun="sum",field="presence",background=0)
  trials=rasterize(points,envdata,fun="count",field="presence",background=0)
  
  ## Set dimentions
  d_lat=ncdim_def("lat",units="degrees_north",vals=unique(coordinates(envdata)[,2]),longname="latitude")
  d_lon=ncdim_def("lon",units="degrees_east",vals=unique(coordinates(envdata)[,1]),longname="longitude")
  d_var=ncdim_def("variable",units="",create_dimvar=F,vals=1:nlayers(envdata))  

    ## Define variables
  comp=9  #define compression level: 9 is the highest
  v_var_var=ncvar_def("envdata",units="standardized",dim=list(d_lon,d_lat,d_var),missval=-999,
                       longname="environmental covariates",compress=comp,prec="float")
  v_var_obs=ncvar_def("obs",units="presence",dim=list(d_lon,d_lat),missval=-999,
                             longname="Observations",compress=comp,prec="integer")
  v_var_trials=ncvar_def("trials",units="trial",dim=list(d_lon,d_lat),missval=-999,
                      longname="Trials",compress=comp,prec="integer")
  v_var_cell=ncvar_def("cell",units="cell",dim=list(d_lon,d_lat),missval=-999,
                         longname="Unique gridcell ID",compress=comp,prec="integer")
  
  ## set up nc file
  nc_create(ncfile,vars=list(v_var_var,v_var_obs,v_var_trials,v_var_cell),verbose=F) 
  
  nc=nc_open(ncfile,write=T)
  if(verbose) print("NetCDF file created, adding data")
  
  
  ## Add data
  lapply(1:nlayers(senv), function(i) ncvar_put(nc,"envdata",vals=t(raster::as.matrix(senv[[i]])),
                                                start=c(1,1,i),c(-1,-1,1),verb=F))
  ncatt_put(nc,varid="envdata", "scales",csd,prec="float")
  ncatt_put(nc,varid="envdata", "offsets",cmeans,prec="float")
  ncatt_put(nc,varid="envdata", "names", paste(names(envdata),collapse=","))
  
  ncvar_put(nc,"obs",vals=t(raster::as.matrix(presences)),start=c(1,1),c(-1,-1),verb=F)
  ncvar_put(nc,"trials",vals=t(raster::as.matrix(trials)),start=c(1,1),c(-1,-1),verb=F)
  ncvar_put(nc,"cell",vals=t(raster::as.matrix(cell)),start=c(1,1),c(-1,-1),verb=F)
  
    
  if(verbose) print("Data added, updating attributes")
  ################################
  ## Attributes
  ## Global Attributes
  ncatt_put(nc,varid=0, "Conventions","Cf-1.4",prec="character")
  ncatt_put(nc,varid=0, "title","Model Input Data",prec="character")
  ncatt_put(nc,varid=0, "date",as.character(Sys.time()),prec="character")
  ncatt_put(nc,varid=0, "projection",projection(envdata),prec="character")
  ncatt_put(nc,varid=0, "projection_format","PROJ.4",prec="character")
  
  ## add all metadata listed in the meta object
  if(!is.null(meta))
    for(i in 1:length(meta)) ncatt_put(nc,varid=0, names(meta)[i],unlist(meta[i]),prec="character")
  
  ## Close the file
  nc_sync(nc)
  nc_close(nc)
  
  if(verbose) print("Finished...")
  
  
}
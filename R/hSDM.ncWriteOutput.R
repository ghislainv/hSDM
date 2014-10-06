###################################################################
##
## hSDM.ncWriteOutput.R
##
####################################################################
##
## Original code by Adam M. Wilson, October 2014
## YALE
## adam.wilson@yale.edu
##
####################################################################
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991. See the package LICENSE
## file for more information.
##
## Copyright (C) 2014 Adam M. Wilson
## 
####################################################################

hSDM.ncWriteOutput<-function(results,file,overwrite=T,autocor=F,keepall=F,meta=NULL,verbose=T){

  today=format(Sys.Date(),format="%Y%m%d")
  time=format(Sys.time(),format="%H%M%S")

  if(verbose) writeLines(paste("Preparing to write",file))
  
  ## assess convergence for each parameter
  if(verbose) writeLines("Calculating convergence metrics")
  parameters_list=mcmc.list(lapply(results,FUN=function(x) x$mcmc))
  c1=gelman.diag(parameters_list, confidence = 0.95, transform=FALSE, autoburnin=FALSE, multivariate=TRUE)
  colnames(c1$psrf)=c("GelmanPSRF","GelmanPSRF.CI")
  c2=geweke.diag(parameters_list[[1]], frac1=0.1, frac2=0.5)
  c3=heidel.diag(parameters_list[[1]], eps=0.1, pvalue=0.05)
  colnames(c3)=c("Heidel.Stationarity","Heidel.Start","Heidel.P","Heidel.HalfwidthTest","Heidel.Mean","Heidel.Halfwidth")
  class(c3)="matrix"  #convert to matrix for easier cbinding
  
  ## summarize posterior distributions  
  if(verbose) writeLines("Calculating posterior summaries")
  parameters=data.frame(mean=summary(parameters_list)$statistics[,"Mean"],
                        sd=summary(parameters_list)$statistics[,"SD"],
                        median=summary(parameters_list)$quantiles[,"50%"],
                        HPDinterval(mcmc(as.matrix(parameters_list))),
                        RejectionRate=rejectionRate(parameters_list),
                        c1$psrf,GewekeZ=c2$z,c3)

  ## predictions for each cell
  ## todo: add option when full posteriors are saved to export quantiles, etc.
  if(verbose) writeLines("Summarizing pixel-level posteriors")
  
  pred=data.frame(x=results[[1]]$model$predcoords$x,
                  y=results[[1]]$model$predcoords$y,
                  cell=results[[1]]$model$predcoords$cell,
                  pred=rowMeans(do.call(cbind,lapply(results,FUN=function(x) x$prob.p.pred))))
  ## convert to raster
  predr=pred
  coordinates(predr)=c("x","y")
  predr=SpatialPixelsDataFrame(predr,tolerance=0.005,data=data.frame(predr[,c("pred","cell")]))
  fullgrid(predr)=T
  predr=stack(predr)
  projection(predr)='+proj=longlat'


  ##### Autocorrelation
  if(autocor){
    if(verbose) writeLines("Calculating autocorrelation of output")
    ac=acor_table(predr[["pred"]],verbose=F)
    ## Global Spatial Autocorrelation
    spac=c(MoransI=Moran(predr[["pred"]],w=matrix(1,11,11)),
           GearyC=Geary(predr[["pred"]],w=matrix(1,11,11)))
  }
  
  ## AUC
  if(verbose) writeLines("Calculating model evaluation metrics")
    aucdat=results[[1]]$model$data[,c("presences","trials","cell")]
    aucdat$pred=pred$pred[match(aucdat$cell,pred$cell)]
    e=evaluate(p=aucdat$pred[aucdat$presences>0],a=aucdat$pred[aucdat$presences==0])
  
  evaluation=data.frame(nPresence=e@np,nTrials=e@na,auc=e@auc,cor=e@cor)
  evaluation$Deviance=parameters$mean[grepl("Deviance",rownames(parameters))]
  evaluation$Pd=(parameters$sd[grepl("Deviance",rownames(parameters))]^2)/2 #Gelman BDA pg 182
  evaluation$DIC=evaluation$Deviance+evaluation$Pd
  evaluation$nChains=nchains
  evaluation$nBurnin=mcpar(results[[1]]$mcmc)[1]
  evaluation$nIter=mcpar(results[[1]]$mcmc)[2]
  evaluation$thin=mcpar(results[[1]]$mcmc)[3]
  evaluation$GlobalGelmanMPSRF=c1$mpsrf
  if(autocor) evaluation[,c("MoransI","GearyC")]=spac
  
  ############################################################
  ## Write results to netcdf file
  if(verbose) writeLines("Setting up netCDF file")
  comp=9  #define compression level: 9 is the highest
  
  ## Set dimentions
      d_lat=ncdim_def("lat",units="degrees_north",vals=sort(unique(coordinates(predr)[,2]),decreasing=F),longname="latitude")
      d_lon=ncdim_def("lon",units="degrees_east",vals=sort(unique(coordinates(predr)[,1])),longname="longitude")
      d_iter=ncdim_def("iter",units="iterations",longname="Posterior Iterations",vals=1:nrow(results[[1]]$mcmc),unlim=TRUE)

      v_var_mean=ncvar_def("p",units="probability",dim=list(d_lon,d_lat),missval=-999,
                           longname="p(occurrence|data)",compression=comp,prec="integer")
      v_var_cell=ncvar_def("cell",units="cell",dim=list(d_lon,d_lat),missval=-999,
                       longname="Unique gridcell ID",compression=comp,prec="integer")
  
      #autocorrelation summaries
      if(autocor){
        d_ac1=ncdim_def("autocorrelation1",units="",create_dimvar=F,vals=1:ncol(ac))  
        d_ac2=ncdim_def("autocorrelation2",units="",create_dimvar=F,vals=1:nrow(ac))  
          v_var_ac=ncvar_def("ac",units="values",dim=list(d_ac1,d_ac2),missval=-999,
                             longname="Autocorrelation",compression=comp)
      }
  
  
  #posterior parameter summaries
  d_params1=ncdim_def("parameters1",units="",create_dimvar=F,vals=1:ncol(parameters))  
  d_params2=ncdim_def("parameters2",units="",create_dimvar=F,vals=1:nrow(parameters)) 
  #evaluation summaries
  d_eval1=ncdim_def("evaluation1",units="",create_dimvar=F,vals=1:ncol(evaluation)) 
  d_eval2=ncdim_def("evaluation2",units="",create_dimvar=F,vals=1:nrow(evaluation)) 
  
  ## Define variables
  v_var_parameters=ncvar_def("parameters",units="values",dim=list(d_params1,d_params2),missval=-999,
                             longname="Posterior Parameter Values",compression=comp)
  v_var_evaluate=ncvar_def("evaluation",units="values",dim=list(d_eval1,d_eval2),missval=-999,
                           longname="Evaluation Metrics",compression=comp)
    
  ## set up nc file
  if(!overwrite&file.exists(file)) stop("File exists, set overwrite=T to overwrite")
  if(overwrite&file.exists(file)) file.remove(file)

  if(!keepall&autocor) nc_create(file,vars=list(v_var_mean,v_var_cell,v_var_parameters,v_var_evaluate,v_var_ac),verbose=F)   #save every iteration
  if(!keepall&!autocor) nc_create(file,vars=list(v_var_mean,v_var_cell,v_var_parameters,v_var_evaluate),verbose=F)   #save every iteration

  nc=nc_open(file,write=T)
  if(verbose) writeLines("NetCDF file created, adding data")

  ncatt_put(nc,"parameters", "colnames", paste(colnames(parameters),collapse=","))
  ncatt_put(nc,"parameters", "rownames", paste(rownames(parameters),collapse=","))
  
  ncatt_put(nc,"evaluation","colnames",paste(colnames(evaluation),collapse=","),prec="char")
  ncatt_put(nc,"evaluation","rownames",paste(rownames(evaluation),collapse=","),prec="char")
  
  ## Add data
  ncvar_put(nc,"parameters",vals=t(as.matrix(parameters)),start=c(1,1),c(-1,-1),verbose=F)
  ncvar_put(nc,"evaluation",vals=t(as.matrix(evaluation)),start=c(1,1),c(-1,-1),verbose=F)
  
  ## Add map data
  predr2=t(raster::as.matrix(predr[["pred"]]))[,nrow(predr):1]
  ncvar_put(nc,"p",vals=predr2*1000,start=c(1,1),c(-1,-1),verbose=F)
  ncatt_put(nc,varid="p", "projection",projection(predr),prec="character")
  ncatt_put(nc,varid="p", "projection_format","PROJ.4",prec="character")
  ncatt_put(nc,varid="p", "scale_factor",.001,prec="double")

  cell2=t(raster::as.matrix(predr[["cell"]]))[,nrow(predr):1]
  ncvar_put(nc,"cell",vals=cell2,start=c(1,1),c(-1,-1),verbose=F)
  ncatt_put(nc,varid="cell", "projection",projection(predr),prec="character")
  ncatt_put(nc,varid="cell", "projection_format","PROJ.4",prec="character")
  
  if(autocor){
    ncatt_put(nc,"ac","colnames",paste(colnames(ac),collapse=","),prec="char")
    ncvar_put(nc,"ac",vals=t(raster::as.matrix(ac)),start=c(1,1),c(-1,-1),verbose=F)
  }
  
  if(verbose) writeLines("Data added, updating attributes")
  ################################
  ## Attributes
  ## Global Attributes
  ncatt_put(nc,varid=0, "Conventions","Cf-1.4",prec="character")
  ncatt_put(nc,varid=0, "title",paste0("Predicted p(occurrence)"),prec="character")
  ncatt_put(nc,varid=0, "modeltype",as.character(results[[1]]$meta$modeltype),prec="character")
  ncatt_put(nc,varid=0, "suitability",as.character(results[[1]]$meta$suitability),prec="character")
  ncatt_put(nc,varid=0, "observability",as.character(results[[1]]$meta$observability),prec="character")
  ncatt_put(nc,varid=0, "date",as.character(today),prec="character")
  ## add any additional metadata listed in the meta object
  if(!is.null(meta))
      for(i in 1:length(meta)) ncatt_put(nc,varid=0, as.character(names(meta)[i]),as.character(unlist(meta[i])),prec="character")
  
  ## Close the file
  nc_sync(nc)
  nc_close(nc)
  
  if(verbose) print("Finished...")
}
###################################################################
##
## hSDM.nc.R
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

hSDM.nc<-function(results,species,modelname,model,data,fdata,outputdir,autocor=T,keepall=F,verbose=T){

  today=format(Sys.Date(),format="%Y%m%d")
  time=format(Sys.time(),format="%H%M%S")
  #ncfile=paste0(outputdir,"/",species,"_",modelname,"_",today,"_",time,".nc",sep="")
  ncfile=paste0(outputdir,"/",species,"_",modelname,".nc",sep="")
  
  if(verbose) writeLines(paste("Preparing to write",ncfile))
  
  ## assess convergence for each parameter
  if(verbose) writeLines("Calculating convergence metrics")
  parameters_list=mcmc.list(lapply(results,FUN=function(x) x$mcmc))
  c1=gelman.diag(parameters_list, confidence = 0.95, transform=FALSE, autoburnin=FALSE, multivariate=TRUE)
  colnames(c1$psrf)=c("GelmanPSRF","GelmanPSRF.CI")
  c2=geweke.diag(parameters_list[[1]], frac1=0.1, frac2=0.5)
  c3=heidel.diag(parameters_list[[1]], eps=0.1, pvalue=0.05)
  colnames(c3)=c("Heidel.Stationarity","Heidel.Start","Heidel.P","Heidel.HalfwidthTest","Heidel.Mean","Heidel.Halfwidth")
  class(c3)="matrix"  #convert class for easier cbinding
  
  ## summarize posterior distributions  
  if(verbose) writeLines("Calculating posterior summaries")
  parameters=data.frame(mean=summary(parameters_list)$statistics[,"Mean"],
                        sd=summary(parameters_list)$statistics[,"SD"],
                        median=summary(parameters_list)$quantiles[,"50%"],
                        HPDinterval(mcmc(as.matrix(parameters_list))),
                        RejectionRate=rejectionRate(parameters_list),
                        c1$psrf,GewekeZ=c2$z,c3)
  
  ## predictions for each cell
  if(verbose) writeLines("Summarizing pixel-level posteriors")
  
  pred=data.frame(x=data$x,y=data$y,cell=data$cell,
                  pred=rowMeans(do.call(cbind,lapply(results,FUN=function(x) x$prob.p.pred))))
  predr=rasterFromXYZ(xyz=pred[,c("x","y","pred")])
  projection(predr)='+proj=longlat'
  
  ##### Autocorrelation
  if(autocor){
    if(verbose) writeLines("Calculating autocorrelation of output")
    ac=acor_table(predr,verbose=F)
    ## Global Spatial Autocorrelation
    spac=c(MoransI=Moran(predr,w=matrix(1,11,11)),
           GearyC=Geary(predr,w=matrix(1,11,11)))
  }
  
  
  ## AUC
  if(verbose) writeLines("Calculating model evaluation metrics")
    aucdat=merge(fdata[,c("presences","trials","cell")],pred,by=c("cell"))
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
  
  ## Set dimentions
  d_lat=ncdim_def("lat",units="degrees_north",vals=sort(unique(pred$y),decreasing=F),longname="latitude")
  d_lon=ncdim_def("lon",units="degrees_east",vals=sort(unique(pred$x)),longname="longitude")
  d_iter=ncdim_def("iter",units="iterations",longname="Posterior Iterations",vals=1:nrow(results[[1]]$mcmc),unlim=TRUE)
  
  ## model data
  d_var1=ncdim_def("var1",units="",create_dimvar=F,vals=1:ncol(data)) 
  d_var2=ncdim_def("var2",units="",create_dimvar=F,vals=1:nrow(data)) 
  #posterior parameter summaries
  d_params1=ncdim_def("parameters1",units="",create_dimvar=F,vals=1:ncol(parameters))  
  d_params2=ncdim_def("parameters2",units="",create_dimvar=F,vals=1:nrow(parameters)) 
  #evaluation summaries
  d_eval1=ncdim_def("evaluation1",units="",create_dimvar=F,vals=1:ncol(evaluation)) 
  d_eval2=ncdim_def("evaluation2",units="",create_dimvar=F,vals=1:nrow(evaluation)) 
  #autocorrelation summaries
  if(autocor){
    d_ac1=ncdim_def("autocorrelation1",units="",create_dimvar=F,vals=1:ncol(ac))  
    d_ac2=ncdim_def("autocorrelation2",units="",create_dimvar=F,vals=1:nrow(ac))  
  }
  
  ## Define variables
  comp=9  #define compression level: 9 is the highest
  v_var_mean=ncvar_def("p",units="probability",dim=list(d_lon,d_lat),missval=-999,
                       longname="p(occurrence|data)",compress=comp,prec="integer")
  v_var_data=ncvar_def("data",units="fitting_data",dim=list(d_var1,d_var2),missval=-999,
                             longname="data used in model fitting",compress=comp)
  v_var_parameters=ncvar_def("parameters",units="values",dim=list(d_params1,d_params2),missval=-999,
                             longname="Posterior Parameter Values",compress=comp)
  v_var_evaluate=ncvar_def("evaluation",units="values",dim=list(d_eval1,d_eval2),missval=-999,
                           longname="Evaluation Metrics",compress=comp)
  
  if(autocor){
    v_var_ac=ncvar_def("ac",units="values",dim=list(d_ac1,d_ac2),missval=-999,
                       longname="Autocorrelation",compress=comp)
  }
  
  ## set up nc file
  if(file.exists(ncfile)) file.remove(ncfile)
  if(keepall) nc_create(ncfile,vars=list(v_var,v_var_mean,v_var_sd,v_krige,v_var_valid),verbose=F)   #save every iteration
  if(!keepall) nc_create(ncfile,vars=list(v_var_mean,v_var_parameters,v_var_evaluate,v_var_ac,v_var_data),verbose=F) 
  
  nc=nc_open(ncfile,write=T)
  if(verbose) print("NetCDF file created, adding data")
  
  ncatt_put(nc,"data", "colnames", paste(colnames(data),collapse=","))
  
  ncatt_put(nc,"parameters", "colnames", paste(colnames(parameters),collapse=","))
  ncatt_put(nc,"parameters", "rownames", paste(rownames(parameters),collapse=","))
  
  ncatt_put(nc,"evaluation","colnames",paste(colnames(evaluation),collapse=","),prec="char")
  ncatt_put(nc,"evaluation","rownames",paste(rownames(evaluation),collapse=","),prec="char")
  
  if(autocor) ncatt_put(nc,"ac","colnames",paste(colnames(ac),collapse=","),prec="char")
  
  ## Add data
  ncvar_put(nc,"data",vals=t(as.matrix(data)),start=c(1,1),c(-1,-1),verb=F)
  ncvar_put(nc,"parameters",vals=t(as.matrix(parameters)),start=c(1,1),c(-1,-1),verb=F)
  ncvar_put(nc,"evaluation",vals=t(as.matrix(evaluation)),start=c(1,1),c(-1,-1),verb=F)
  if(autocor) ncvar_put(nc,"ac",vals=t(as.matrix(ac)),start=c(1,1),c(-1,-1),verb=F)
  
  ## Add map data
  ncvar_put(nc,"p",vals=1000*t(as.matrix(predr))[,nrow(predr):1],start=c(1,1),c(-1,-1),verb=F)
  ncatt_put(nc,varid="p", "projection",projection(predr),prec="character")
  ncatt_put(nc,varid="p", "projection_format","PROJ.4",prec="character")
  ncatt_put(nc,varid="p", "scale_factor",.001,prec="double")
  
  
  if(verbose) print("Data added, updating attributes")
  ################################
  ## Attributes
  ## Global Attributes
  ncatt_put(nc,varid=0, "Conventions","Cf-1.4",prec="character")
  ncatt_put(nc,varid=0, "title",paste0("Predicted p(occurrence) for ",sp2),prec="character")
  ncatt_put(nc,varid=0, "institution","Map of Life, Yale University, New Haven, CT",prec="character")
  ncatt_put(nc,varid=0, "source","Modeled Species Distributions",prec="character")
  ncatt_put(nc,varid=0, "comment","Adam M. Wilson (adam.wilson@yale.edu)",prec="character")
  ncatt_put(nc,varid=0, "model",as.character(model),prec="character")
  ncatt_put(nc,varid=0, "modelname",as.character(modelname),prec="character")
  ncatt_put(nc,varid=0, "species",as.character(sp2),prec="character")
  ncatt_put(nc,varid=0, "date",as.character(today),prec="character")
  
  ## Close the file
  nc_sync(nc)
  nc_close(nc)
  
  if(verbose) print("Finished...")
}
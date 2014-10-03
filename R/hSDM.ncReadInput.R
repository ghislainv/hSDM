hSDM.ncReadInput<-function(ncfile){
      nc=nc_open(ncfile,write=F)
      ## metadata
      species=ncatt_get(nc, 0, attname="species", verbose=FALSE )$value
      ## envdata
      ed=stack(ncfile,varname="envdata")
      names(ed)=strsplit(ncatt_get(nc, "envdata", attname="names", verbose=FALSE )$value,",")[[1]]
      ## trials, presences, and cell ID
      trials=raster(ncfile,varname="trials")
      presences=raster(ncfile,varname="obs")
      cell=raster(ncfile,varname="cell")
      ## combine to a single data.frame
      d=cbind.data.frame(species=species,coordinates(ed),values(ed),trials=values(trials),presences=values(presences),cell=values(cell))
      ## omit rows with missing data
      d2=na.omit(d)
      if(nrow(d)!=nrow(d2)) warning(paste(nrow(d)-nrow(d2)," rows deleted due to missing data"))
      ## clean up
      nc_close(nc)
      return(d2)
    }
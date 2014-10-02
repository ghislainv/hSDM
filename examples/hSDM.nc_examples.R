
library(doMC)
registerDoMC(3)
nchains=3

library(hSDM)
library(raster)
library(sp)

#===================================
#== Multivariate normal distribution
rmvn <- function(n, mu = 0, V = matrix(1), seed=1234) {
  p <- length(mu)
  if (any(is.na(match(dim(V), p)))) {
    stop("Dimension problem!")
  }
  D <- chol(V)
  set.seed(seed)
  t(matrix(rnorm(n*p),ncol=p)%*%D+rep(mu,rep(n,p)))
}

#==================
#== Data simulation

#= Set seed for repeatability
seed <- 1234

#= Landscape
xLand <- 30
yLand <- 30 
Landscape <- raster(ncol=xLand,nrow=yLand,crs='+proj=utm +zone=1')
Landscape[] <- 0
extent(Landscape) <- c(0,xLand,0,yLand)
coords <- coordinates(Landscape)
ncells <- ncell(Landscape)

#= Neighbors
neighbors.mat <- adjacent(Landscape, cells=c(1:ncells), directions=8, pairs=TRUE, sorted=TRUE)
n.neighbors <- as.data.frame(table(as.factor(neighbors.mat[,1])))[,2]
adj <- neighbors.mat[,2]

#= Generate symmetric adjacency matrix, A
A <- matrix(0,ncells,ncells)
index.start <- 1
for (i in 1:ncells) {
  index.end <- index.start+n.neighbors[i]-1
  A[i,adj[c(index.start:index.end)]] <- 1
  index.start <- index.end+1
}

#= Spatial effects
Vrho.target <- 5
d <- 1  # Spatial dependence parameter = 1 for intrinsic CAR
Q <- diag(n.neighbors)-d*A + diag(.0001,ncells) # Add small constant to make Q non-singular
covrho <- Vrho.target*solve(Q) # Covariance of rhos
set.seed(seed)
rho <- c(rmvn(1,mu=rep(0,ncells),V=covrho,seed=seed)) # Spatial Random Effects
rho <- rho-mean(rho) # Centering rhos on zero

#= Raster and plot spatial effects
r.rho <- rasterFromXYZ(cbind(coords,rho))
plot(r.rho)

## make some spatially correlated environmental data
set.seed(seed+3)
X1 <- c(rmvn(1,mu=rep(0,ncells),V=covrho,seed=seed+1)) # Spatial Random Effects
X2 <- c(rmvn(1,mu=rep(0,ncells),V=covrho,seed=seed+2)) # Spatial Random Effects
r.X1 <- rasterFromXYZ(cbind(coords,X1))
r.X2 <- rasterFromXYZ(cbind(coords,X2))
plot(stack(r.X1,r.X2))

#= Sample the observation sites in the landscape
nsite <- 250
set.seed(seed)
x.coord <- runif(nsite,0,xLand)
set.seed(2*seed)
y.coord <- runif(nsite,0,yLand)
sites.sp <- SpatialPoints(coords=cbind(x.coord,y.coord))
cells <- extract(Landscape,sites.sp,cell=TRUE)[,1]

#= Number of visits associated to each observation point
set.seed(seed)
visits <- rpois(nsite,3)
visits[visits==0] <- 1

#= Ecological process (suitability)
set.seed(seed)
x1 <- rnorm(nsite,0,1)
set.seed(2*seed)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
beta.target <- c(-1,1,-1)
logit.theta <- X %*% beta.target + rho[cells]
theta <- inv.logit(logit.theta)
set.seed(seed)
Y <- rbinom(nsite,visits,theta)

#= Relative importance of spatial random effects
RImp <- mean(abs(rho[cells])/abs(X %*% beta.target))
RImp

#= Data-sets
data.obs <- data.frame(Y,visits,x1,x2,cell=cells)

#==================================
#== Site-occupancy model

res=foreach(ch=1:nchains) %dopar% {
  mod.hSDM.binomial.iCAR <- hSDM.binomial.iCAR(presences=data.obs$Y,
                                               trials=data.obs$visits,
                                               suitability=~x1+x2,
                                               spatial.entity=data.obs$cell,
                                               data=data.obs,
                                               n.neighbors=n.neighbors,
                                               neighbors=adj,
                                               suitability.pred=NULL,
                                               spatial.entity.pred=NULL,
                                               burnin=5000, mcmc=5000, thin=5,
                                               beta.start=0,
                                               Vrho.start=1,
                                               mubeta=0, Vbeta=1.0E6,
                                               priorVrho="1/Gamma",
                                               shape=0.5, rate=0.0005,
                                               seed=1234, verbose=1,
                                               save.rho=1, save.p=0)
}



hSDM.nc(results=res,species="test",fdata=data.obs,file="test.nc",overwrite=T,meta=NULL,autocor=T)


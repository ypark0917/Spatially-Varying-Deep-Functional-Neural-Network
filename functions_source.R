#########################################
gen.groups <- function(n, nfold){
  leave.out <- trunc(n/nfold)
  o <- sample(1:n)
  groups <- vector("list", nfold)
  for (j in (1:(nfold-1))){
    jj <- (1+(j-1)*leave.out)
    groups[[j]] <- (o[jj:(jj+leave.out-1)])
  }
  groups[[nfold]] <- o[(1+(nfold-1)*leave.out):n]
  return(groups=groups)
}
###############################

grid_finder<-function(lon, lat, full.lon, full.lat,plot=TRUE){
  
  
  library(fields)
  log.grid <- seq.default(0,360,, full.lon)
  lat.grid <- seq.default(-90,90,, full.lat)
  xy.grid <- expand.grid(log.grid,lat.grid)
  
  
  start = which.min(abs(xy.grid[,1]-lon[1] ) +   abs(xy.grid[,2]-lat[1] )    )
  mid1 = which.min(abs(xy.grid[,1]-lon[2] ) +   abs(xy.grid[,2]-lat[1] )    )
  mid2 = which.min(abs(xy.grid[,1]-lon[1] ) +   abs(xy.grid[,2]-lat[2] )    )
  final=  which.min(abs(xy.grid[,1]-lon[2] ) +   abs(xy.grid[,2]-lat[2] )    )
  
  
  st=ASIA= c(start:mid1 )
  for(m in 1:( (mid2-start)/full.lon)){
    ASIA<-  cbind( ASIA, st+ full.lon*m   )}
  
  
  if( plot==TRUE){
    A<-matrix(5,nrow=1, ncol= full.lon* full.lat)
    A[1,ASIA]<- 1
    image.fit <- as.image(A,x=xy.grid, grid=list(x=log.grid, y=lat.grid))
    image.plot(image.fit, main='Region', zlim=c(1:2))
    map("world2", ylim=c(-90,90), xlim = c(0,360), add = TRUE)
    
  }
  
  print( paste('number of lon : '  ,length(which(xy.grid[ASIA,]==xy.grid[ASIA[1],1])) ) )
  print(paste('number of lat : '  , length(which(xy.grid[ASIA,]==xy.grid[ASIA,2][1]))))
  return (as.vector(ASIA))
  
}


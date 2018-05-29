########################################################### Helper Functions ###########################################################

######### rescaleEnviron #########
# Function rescales the environmental layer between a new range allowing for consistency across simulated environments
rescaleEnviro <- function(env,  new.min = 0, new.max = 1) {
  new.min + (env - env.min) * ((new.max - new.min) / (env.max - env.min))
}
######### replaceNonNAValues #########
# function replaces values in raster (x) that are not NA with new values (tr)
replaceNonNAValues <- function(x, tr) {
  x[which(!is.na(x@data@values))] = tr
  x
}
######### replaceNonNAValues #########
# function replaces values in raster (x) that are NA with new values (tr)
replaceNAValues <- function(x, tr) {
  if(class(x) == "RasterLayer"){
  x[which(is.na(x@data@values))] = tr
  } else {
    x <- unstack(x)
    if(length(x) > 1) { break }
    x <- x[[1]]
    x[which(is.na(x@data@values))] = tr
  }
  x
}
######### rasterOverlap #########
# function asks if two rasters (species rasters specificially) overlap
rastersOverlap <- function (x1, x2) {
  if(any(Which(x1 == x2, cells=TRUE))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

######### getMinimumDistanceBetweenPoints #########
# returns minimum distance between two species ranges
getMinimumDistanceBetweenPoints <- function(coords.list, coords.df) {
  min(pointDistance(coords.list[,1:2], coords.df[which(coords.df[,4] != coords.list[1,4]),1:2], longlat=F, all.pairs=T), na.rm=T)
}
######### buildPhyInSim #########
# builds phylogeny within the simulator function
buildPhyInSim <- function (edgetable, n.layers, extinct, tips){
  edgetable2 <- edgetable 
  edgetable <- edgetable[-which(is.na(edgetable[, 1])),]
  rownames(edgetable) <- NULL
  colnames(edgetable) <- NULL
  tips.table <- matrix(cbind(matrix(edgetable[which(!edgetable[, 2] %in% edgetable[, 1])[!which(!edgetable[, 2] %in% edgetable[, 1]) %in% which(extinct==TRUE)],], ncol=10), c(1:length(n.layers))), ncol=11)
  edgetable <- merge(edgetable,tips.table, all.x=T)
  edgetable<-as.matrix(edgetable, ncol=11)
  rownames(edgetable) <- NULL
  colnames(edgetable) <- NULL
  simtrtable <- edgetable
  simtrtable2 <- simtrtable[2:nrow(simtrtable), ]
  tips <- length(as.numeric(simtrtable2[, 2])[which(!as.numeric(simtrtable2[, 2]) %in% as.numeric(simtrtable2[, 1]))])
  fixmat <- rbind(c(0:max(as.numeric(simtrtable2[,2]))), c(max(as.numeric(simtrtable2[,2])):0))
  for(i in 1:length(as.numeric(simtrtable2[,1]))) {
    initval <- as.numeric(simtrtable2[i, 1])
    simtrtable2[i, 1] <- fixmat[2,][which(fixmat[1,] == initval)]
    initval <- as.numeric(simtrtable2[i, 2])
    simtrtable2[i, 2] <- fixmat[2,][which(fixmat[1,] == initval)]
  }
  simtrtabord <- simtrtable2[order(simtrtable2[,1], decreasing = T), ]
  simtrtabord[, 1] <- as.numeric(simtrtabord[, 1]) + 1
  simtrtabord[, 2] <- as.numeric(simtrtabord[, 2]) + 1
  currtips <- as.numeric(simtrtabord[,2][which(!simtrtabord[,2] %in% simtrtabord[,1])])
  currtips <- currtips[which(currtips > tips)]
  wrongbrs <- as.numeric(unique(simtrtabord[, 1][which(as.numeric(simtrtabord[, 1]) <= tips)]))
  for(i in 1:length(currtips)){
    wrongparloc <- which(simtrtabord[, 1] == wrongbrs[i])
    wrongdauloc <- which(simtrtabord[, 2] == currtips[i])
    simtrtabord[, 1][wrongparloc] <- currtips[i]
    simtrtabord[, 2][which(simtrtabord[, 2] == wrongbrs[i])] <- currtips[i]
    simtrtabord[, 2][wrongdauloc] <- wrongbrs[i]
  }
  simtrtabord <- simtrtabord[order(simtrtabord[, 1], decreasing = T), ]
  simtredge <- cbind(as.numeric(simtrtabord[, 1]), as.numeric(simtrtabord[,2]))
  brlentime <- as.numeric(simtrtabord[, 3])
  timephylo <- reorder.phylo(rtree(tips), "postorder")
  timephylo$edge <- simtredge
  timephylo$edge.length <- as.numeric(simtrtabord[, 3])
  timephylo$tip.label <- 1:tips
  timephylo <- read.tree(text = write.tree(timephylo))
  colnames(simtrtabord)<- c("parent", "daughter", "brlen", "specMode", "area", "birthDate", "nichePos", "nicheBreadth", "birthMode", "alive", "rasterMatch")
  return(list(simtrtabord, timephylo))
}
######### generateEnv #########
# simulates an environmental layer with a degree of spatial autocorreltion
generateEnv <- function(grid.size=100, psill=0.005, range=15, rescale.max=25, rescale.min=0, plot=F, original=F) {
  xy <- expand.grid(1:grid.size, 1:grid.size)
  names(xy) <- c('x','y')
  #range is the degree of autocorrelation -> bigger is more autocorrelated
  #Psill is the partial sill of the variogram model - basically at 0 it is a perfect gradient, the higher you getthe more heterogenous it is
  g.dummy <- gstat(formula=z~x+y, locations=~x+y, dummy=T, beta=c(1,0,0.005), model=vgm(psill=0.005, range=100, model='Exp'), nmax=20)
  if(original==TRUE){ g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, model=vgm(psill=10,model="Exp", range=50), nmax=30)}
  yy <- predict(g.dummy, newdata=xy, nsim=1)
  env<-rasterFromXYZ(yy)  
  env@data@values <- rescaleEnviro(env@data@values, new.min=0, new.max=25)
  env@data@max <- max(env@data@values)
  env@data@min <- min(env@data@values)
  if(plot == TRUE) {  plot(env)}
  return(env=env)
}
######### envrasterFromDistribution #########
# Turns species range from presence/absence to the values of the environmental layer
envRasterFromDistribution <- function (species.raster, env=env) {
  cells <- Which(is.na(species.raster), cells=T)
  env[cells]<-NA
  return(env)
}


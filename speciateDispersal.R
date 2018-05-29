######### speciateDispersal #########

# Function divides a species into two daughter species based on an dispersal (founder) model of speciation

############# Arguments ###############

# env = current environmental layer
# position = niche position of species 
# breadth = niche breadth of species
# phylo.sig = phylogenetic signal of species niches at speciation
# species.raster is the current species range in raster format (presence/absence)

speciateDispersal <- function(position, breadth, species.raster, env, phylo.sig, dispersal){
  ras <- species.raster 
  # Get x and y coordinates of centroid of species range
  E<-max(rasterToPoints(ras)[,1])
  W<-min(rasterToPoints(ras)[,1])
  N<-max(rasterToPoints(ras)[,2])
  S<-min(rasterToPoints(ras)[,2])
  x <- E - ((E - W) / 2)
  y <- N - ((N -S) / 2)
  # upper and lower values of species niche
  upper <- position + breadth
  lower <- position - breadth
  # raster with environmental distance of each cell from niche position of species
  new.space <- env-position
  #make values of current species range NA
  new.space[ras] <- NA
  # get distance from range centroid to unoccupied cells
  distance.ras <- distanceFromPoints(new.space, c(x,y))
  # make current species range NA in distance raster
  distance.ras <- mask(distance.ras, new.space)
  # condition that if species occupies full domain can't undergo dispersal speciation so breaks loop
  if(all(is.na(distance.ras@data@values))){
    result <- "no_speciation"
    return(result)
  } else {
    # turn distances into probabilities (close cells have a greater probability of dispersal as cells further away)
    distance.probability <- as.numeric(na.omit(1 - ((distance.ras@data@values/distance.ras@data@max) + abs(new.space@data@values/max(abs(new.space@data@values), na.rm=T))/2))+0.001)
    distance.probability[which(distance.probability < 0)] <- 0
    # sample dispersal cell based on distance probabilities
    if(length(which(!is.na(new.space@data@values))) == 1) { dispersal.cell <- Which(new.space, cells=T) } else {
      dispersal.cell <- sample(Which(distance.ras >=0, cells=T), 1,  prob = distance.probability) 
    } 
  }
  # create range boundaries around selected dispersal cell 
  coords <- xyFromCell(new.space, dispersal.cell)
  extent.dists<-runif(4, min=1, max=dispersal)
  sp2 <- crop(env, extent(c(coords[[1]] - extent.dists[1], coords[[1]] + extent.dists[2], coords[[2]] - extent.dists[3], coords[[2]] + extent.dists[4])))
  sp2[sp2@data@values %in% ras@data@values] <- NA
  sp2@data@values[!is.na(sp2@data@values)] <- 1
  sp2<-extend(sp2, env)
  # conditions that if range of new species and parent species overlap to resample cell and boundaries
  # condition breaks if after 200 resamples it is unable to emeet these conditions (may happen if a species occupies almost the entire domain - a condition that is prevented by capping a species niche niche breadth)
  both.ras<-species.raster+sp2
  loop = 1
  while(any(!is.na(both.ras@data@values) == TRUE) & loop <200) {
    if(length(!is.na(new.space@data@values)) == 1) { dispersal.cell <- Which(new.space, cells=T) } else {
      dispersal.cell <- sample(Which(distance.ras >=0, cells=T), 1,  prob = distance.probability) 
    }
    coords <- xyFromCell(new.space, dispersal.cell)
    extent.dists<-runif(2, min=1, max=dispersal)
    sp2 <- crop(env, extent(c(coords[[1]] - extent.dists[1], coords[[1]] + extent.dists[1], coords[[2]] - extent.dists[2], coords[[2]] + extent.dists[2])))
    sp2[sp2@data@values %in% ras@data@values] <- NA
    sp2@data@values[!is.na(sp2@data@values)] <- 1
    sp2 <- extend(sp2, env)
    both.ras<-species.raster+sp2
    loop<-loop+1
  }
  sp2 <- extend(sp2, env)
  # recenter new daughter species niche and adjust range to accomodate new niche breadth
  recenteredsp2 <- nicheRecenter(position, breadth, sp2, env, phylo.sig)
  # returns new daughter species niche position, niche breadth, and range. second daughter species keeps the attributes of the parent species
  result <- list(positions = recenteredsp2[[1]], breadths = recenteredsp2[[2]], species.rasters=recenteredsp2[[3]])
  return(result)
}





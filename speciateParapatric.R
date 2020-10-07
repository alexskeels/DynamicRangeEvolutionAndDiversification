######### speciateParapatric #########

# Function divides a species into two daughter species based on an parapatric model of speciation

############# Arguments ###############

# env = current environmental layer
# position = niche position of species 
# breadth = niche breadth of species
# phylo.sig = phylogenetic signal of species niches at speciation
# species.raster is the current species range in raster format (presence/absence)

speciateParapatric <- function(position, breadth, species.raster, env, phylo.sig, dispersal){
  # set niche bound
  upper <- position + breadth
  lower <- position - breadth
  # find grid cell numbers of species range and cells within the dispersal capacity of the species
  cells <- Which(!is.na(species.raster), cells=TRUE)
  dispersal.range.cells <- cells
  for ( i in 1:dispersal) {
    dispersal.range.cells <- adjacent(species.raster, unique(dispersal.range.cells), directions=8, pairs=F)
  }
  # select cell within dispersal kernal of species
  new.range.ras <- species.raster
  new.range.ras[dispersal.range.cells][which(!dispersal.range.cells %in% cells)]<-1000
  cell.select <- sample(Which(new.range.ras == 1000, cells=T), 1)
  # make sure that the cell is within the domain
  while(is.na(env[cell.select])){cell.select <- sample(Which(new.range.ras == 1000, cells=T), 1)}
  # get cell coordinates and sample new range boundaries as distances from that cell
  coords <- xyFromCell(new.range.ras, cell.select)
  extent.dists <- runif(4, min=1, max=dispersal)
  # create new species range from these range boundaries
  sp2 <- crop(env, extent(c(coords[[1]] - extent.dists[1], coords[[1]] + extent.dists[2], coords[[2]] - extent.dists[3], coords[[2]] + extent.dists[4])))
  sp2@data@values[!is.na(sp2@data@values)] <- 1
  sp2 <- extend(sp2, env)
  # recenter species niche (phylogenetic signal in species niche value) and readjust range
  recenteredsp2 <- nicheRecenter(position, breadth, sp2, env, phylo.sig)
  # return new species range, niche position, and niche breadth
  result <- list(positions = recenteredsp2[[1]], breadths = recenteredsp2[[2]], species.rasters=recenteredsp2[[3]])
  return(result)
  }





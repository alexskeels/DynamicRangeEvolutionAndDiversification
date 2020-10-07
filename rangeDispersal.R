######### rangeDispersal #########

# Function extends a species range based on a dispersal kernal and the inherent niche breadth of the species

############# Arguments ###############

# env = current environmental layer
# position = niche position of species 
# breadth = niche breadth of species
# species.ras is the current species range in raster format (presence/absence)

rangeDispersal <- function (position, breadth, species.ras, env, dispersal.range){
  # niche bounds
  upper <- position + breadth
  lower <- position - breadth
  # identify cells that the species currently occupies
  cells<-Which(species.ras == 1, cells=T)
  # find adjacent cells to those currently occpuied and repeat for length of the dispersal kernal
  # this produces a vector of all cells that fall within the species dispersal kernal
  dispersal.range.cells <- cells
  for (i in 1:dispersal.range) {
    dispersal.range.cells <- unique(adj(env, dispersal.range.cells, directions=8, pairs=F, include=T, numCol=100, numCell=10000))
  }
  # find the environmental conditions of the new cells
  env.dat <- env[dispersal.range.cells]
  # remove cells beyond the species niche bounds
  dispersal.range.cells <- dispersal.range.cells[which(env.dat > lower & env.dat < upper)]
  # creater new range raster
  new.range <- replaceNonNAValues(species.ras, NA)
  new.range[dispersal.range.cells] <- 1
  if(sum(new.range@data@values, na.rm=T) == 0) { 
    return("extinct") 
  } else { 
    return(new.range) 
  }
}

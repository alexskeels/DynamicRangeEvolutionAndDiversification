######### Dynamic Range Evolution and Diversification model (DREaD) #########

# main simulator function of geographic diversification model

############# Arguments ###############

# totaltips = clade size to generate
# dispersal = dispersal kernal parameter
# enviro.mode = environmental change mode - character "sine" or "linear"
# amp = amplitude of environmental change sine wave
# freq = frequency of environmnental change sine wave
# slope = slope of environmental change linear model
# breadth.ev.rate = rate of niche breadth evolution
# niche.ev.rate = rate of niche position evolution
# phylo.sig = phylogenetic signal of species niches at speciation
# Me = range size threshold at which a species will not go extinct stochastically
# enviro.hetero = biniary variable. whether environmental change model varies spatially or operates in synchrony across the domain
# geo.mode = geographic mode of species
# plot = whether to plot species ranges and phylogeny at end of simulation (plot results)
# animate = whethe to save a snapshot of species range at each time step (png) to turn into a GIF using otehr software (visualisation tool)
# stepsize = size of each time step
# generateSummaryStatistics = logical. whether to generate the 32 summary statistics
# lambda = speciation rates for each speciation mode lambda1 = predominate mode, lambda2 = non-predominate modes, lambda3= mixed speciation


DREaD <- function (totaltips, dispersal, amp, freq, slope, niche.ev.rate, breadth.ev.rate, phylo.sig, Me,  enviro.hetero, geo.mode, enviro.mode,
                          plot = FALSE, stepsize = 0.1, generateSummaryStatistics=TRUE,
                          lambda1 =0.04,lambda2 =0.0005,lambda3 =0.010375) {
  ##### required libraries
  require(raster)
  require(gstat)
  require(ape)
  require(phytools)
  require(geiger)
  require(moments)
  require(phyloclim)
  require(data.table)
  require(fossil)
  require(apTreeshape)
  # require(ENMTools)
  # extinction rate constant
  ext.rate = 0.02
  # matrix descrbing the degree of environmental change across the landscape (for enviro.hetero=T)
  env.change.matrix <- matrix(rep(seq(from=0.01, to =1, by=1/100), 100), ncol=100, nrow=100, byrow=F)
  # vector of parameters
  params<-data.frame(totaltips=totaltips, dispersal=dispersal, amp=amp,
                    freq=freq, niche.ev.rate=niche.ev.rate, breath.ev.rate=breadth.ev.rate, phylo.sig=phylo.sig, Me=Me,
                    geo.mode=geo.mode,  slope=slope, enviro.mode=enviro.mode, enviro.hetero=enviro.hetero)
  # set maximum runtime
  maxtime = 200
  # Keep track of number of Mass Extinctions
  extinction <- 1
  # range size at which speciation probability peaks
  B <- 5000
  # dispersaion of speciation rate around the peak
  C <- 1750
  # set rate of each geographic mode of speciation
  if (geo.mode == "mixed") {allo.rate = lambda3; sym.rate = lambda3; para.rate = lambda3; disp.rate = lambda3}
  if (geo.mode == "allopatric") {allo.rate = lambda1; sym.rate = lambda2; para.rate = lambda2; disp.rate = lambda2}
  if (geo.mode == "parapatric") {allo.rate = lambda2; sym.rate = lambda2; para.rate = lambda1; disp.rate = lambda2}
  if (geo.mode == "sympatric") {allo.rate = lambda2; sym.rate = lambda1; para.rate = lambda2; disp.rate = lambda2}
  if (geo.mode == "dispersal") {allo.rate = lambda2; sym.rate = lambda2; para.rate = lambda2; disp.rate = lambda1}

  ########################### 1. Generate background environment ##########################

  env <- generateEnv(original = T)
  starting.env <- env

  ###########################  2. Seed initial species  ###################################

  initial.species <- seedSpecies(env, dispersal = dispersal)
  # generate edgetable
  # edgetable is a matrix that stores information on each species' phylogeny, niche position, niche breadth, speciation modes, range size
  # each row is a species

  edgetable <- makeEdgeTable(10000)

  edgetable[1,5] <- sum(initial.species[[1]]@data@values)
  edgetable[1,5] <- sum(initial.species[[1]]@data@values)
  edgetable[1,7] <- initial.species[[2]]
  edgetable[1,8] <- initial.species[[3]]

  # species rasters is a list that hangs onto each species geographic range in the form of a raster
  species.rasters <- vector('list', 10000)
  species.rasters[[1]] <- initial.species[[1]]
  time <- 1
  stepsize <- stepsize
  tips <- 1
  extinct <- vector("logical", 10000)
  extinct.number <- 0

  # while loop propels the simulation. iterations repeat until the condition (number of species generated) is met

  while((tips - extinct.number) < totaltips) {

  ############################# 3. Environmental change ###########################

    # time changes
    time <- round(time + stepsize, 3)
    # if simulation runs too long restart
    if(time >= maxtime){
      initial.species <- seedSpecies(env)
      edgetable <- matrix(ncol=10, nrow=10000)
      edgetable[1,] <- c(0, 1, 0, NA, NA, 1, NA, NA, NA, "X")
      edgetable[1,5] <- sum(initial.species[[1]]@data@values)
      edgetable[1,7] <- initial.species[[2]]
      edgetable[1,8] <- initial.species[[3]]
      species.rasters <- vector('list', 10000)
      species.rasters[[1]] <- initial.species[[1]]
      time <- 1
      stepsize <- stepsize
      tips <- 1
      extinct.number=0
      extinct <- vector("logical", 10000)
    }

    # environment changes
    env <- enviroChange(start.env=starting.env, env=env, time=time, amp=amp, freq=freq, slope=slope, model= enviro.mode, hetero=enviro.hetero, env.change.matrix=env.change.matrix)
    # species.tips is the row index of non-extinct lineages (rows in the edgetable)
    if(any(extinct==TRUE)){ species.tips <- seq_along(which(!is.na(edgetable[,10])))[-which(extinct == TRUE)] } else { species.tips <- seq_along(which(!is.na(edgetable[,10]))) }

    # iterate through extant species
    for (i in species.tips) {

      #skips rows that have speciated (i.e., is an ancestral branch not an extant lineage)
      not.na.rows <- which(!is.na(edgetable[,10]))
      if (length(not.na.rows) > 1) {
        if (edgetable[i, 2] %in% edgetable[not.na.rows, 1] ) { next }
      }

      # select species rater for current iteration and niche values
      current.species <- species.rasters[[i]]
      position <- as.numeric(edgetable[i, 7])
      breadth <- as.numeric(edgetable[i, 8])

  ############################# 4. Dispersal ###########################

      # disperse species' range
      current.species <- disperseRange(position, breadth, current.species,env, starting.env, dispersal)

      # species can become extinct if no habitat within dispersal range (driven extinct by environmental change)
      # if goes extinct record the extinction and move onto next iteration
      # otherwise keep new species range
      if (class(current.species) == "character") {
        extinct[i] <- TRUE
        edgetable[i, 10] <- "extinct"
        next
      } else {
        species.rasters[[i]] <- current.species
      }

  ###########################  5. Event selection ######################

        #range size of current species
        area <- sum(current.species@data@values, na.rm=T)
        edgetable[i, 5] <- area
        # extinction probability
        prob.ext <- -log(area / Me)/(1 / ext.rate)^2
        if (prob.ext <= 0) prob.ext <- 0
        # can only speciate allopatrically and sympatrically is range size > 1 (geometric constraint)
        if(area > 1) {
        # allopatric speciation probability
        prob.allo <- allo.rate*exp(-((area-B)^2)/(2*C^2))
        if (prob.allo <= 0) prob.allo <- 0
        # sympatric speciation probability
        prob.sym <- sym.rate*exp(-((area-B)^2)/(2*C^2))
        if (prob.sym <= 0) prob.sym <- 0
        } else {prob.sym <- 0;  prob.allo <- 0}
        # parapatric speciation probability
        prob.para <- para.rate*exp(-((area-B)^2)/(2*C^2))
        if (prob.para <= 0) prob.para <- 0
        # dispersal speciation probability
        prob.disp <- disp.rate*exp(-((area-B)^2)/(2*C^2))
        if (prob.disp <= 0) prob.disp <- 0
        # no-event probability (1 - the probability of other events)
        prob.no.event <- 1 - (prob.ext + prob.allo + prob.sym + prob.disp + prob.para)
        # select event
        event <- sample(c(1,2,3,4,5,6), 1, prob =c(prob.allo, prob.sym, prob.para, prob.disp, prob.no.event, prob.ext))

  ######################### 6. Events ################################
        #specify which is the next empty row in the edgetable
        next.free.edgetable <- max(which(!is.na(edgetable[,10])), na.rm=T)+1
        # Allopatric speciation
        if(event == 1) {
          edgetable[i, 4] <- "allopatric"
          # is species range fragmented of a single contiguous area
          clumped <- clump(current.species)
          # if only a single contiguous range exists, split it via bisection
          if (clumped@data@max == 1) {
            allopatric.speciation <- speciateAllopatric(position, breadth, current.species, env, phylo.sig)
            # add new species ranges to the species.rasters list
            species.rasters <- replace(species.rasters, next.free.edgetable, allopatric.speciation$species.rasters[[1]])
            species.rasters <- replace(species.rasters, next.free.edgetable+1, allopatric.speciation$species.rasters[[2]])
            # add information on two new daughter species to the edgetable keeping track of their parent species, niche position, niche breadth, range size, and speciation mode
            newbr1 <- c(edgetable[i, 2], next.free.edgetable, stepsize, NA,
                        sum(allopatric.speciation$species.rasters[[1]]@data@values, na.rm=T), time, allopatric.speciation$pos[[1]],
                        allopatric.speciation$breadth[[1]], "allopatric", "X")

            newbr2 <- c(edgetable[i, 2], next.free.edgetable+1, stepsize, NA,
                        sum(allopatric.speciation$species.rasters[[2]]@data@values, na.rm=T), time, allopatric.speciation$pos[[2]],
                        allopatric.speciation$breadth[[2]], "allopatric", "X")
            edgetable[c(next.free.edgetable,next.free.edgetable+1),] <- rbind(newbr1, newbr2)
          }
          # if species' range is divided in two these two range fragements will become the new daughter species
          if (clumped@data@max == 2) {
            sp1 <- clumped
            sp1[sp1 == 2] <- NA
            sp2 <- clumped
            sp2[sp2 == 1] <- NA
            recenteredsp1 <- nicheRecenter(position, breadth, sp1, env, phylo.sig)
            recenteredsp2 <- nicheRecenter(position, breadth, sp2, env, phylo.sig)
            species.rasters <- replace(species.rasters, next.free.edgetable, recenteredsp1[[3]])
            species.rasters <- replace(species.rasters, next.free.edgetable+1, recenteredsp2[[3]])
            newbr1 <- c(edgetable[i, 2], next.free.edgetable, stepsize, NA,
                        sum(recenteredsp1[[3]]@data@values, na.rm=T), time, recenteredsp1[[1]],
                        recenteredsp1[[2]], "allopatric", "X")
            newbr2 <- c(edgetable[i, 2], next.free.edgetable+1, stepsize, NA,
                        sum(recenteredsp2[[3]]@data@values, na.rm=T), time, recenteredsp2[[1]],
                        recenteredsp2[[2]], "allopatric", "X")
            edgetable[c(next.free.edgetable,next.free.edgetable+1),] <- rbind(newbr1, newbr2)
          }
          # if the species' range is fragmented into more than two fragments  the two daughter species range are determined by a spatial clustering (k-means) algorithm
          if (clumped@data@max >= 3) {
            ras.values <- clumped@data@values
            ras.coords <- coordinates(current.species)
            ras.df <- data.frame(values=ras.values, x=ras.coords[,1], y=ras.coords[,2], cluster=NA)
            NAs <- which(is.na(ras.df[,1]))
            ras.df[NAs,1:3]<-NA
            k.ras <- kmeans(na.omit(ras.df[,1:3]), 2)
            ras.df$cluster[which(!is.na(ras.df$values))] <- k.ras$cluster
            k.ras.clusters <- current.species
            k.ras.clusters@data@values <- ras.df$cluster
            sp1 <- k.ras.clusters
            sp2 <- k.ras.clusters
            sp1[sp1 == 2] <- NA
            sp2[sp2 == 1] <- NA
            both.ras<-sp1+sp2
            while(any(!is.na(both.ras@data@values) == TRUE)) {
              sp1 <- sp1 - (both.ras/2)
              sp2 <- sp2 - (both.ras/2)
              both.ras<-sp1+sp2
            }
            recenteredsp1 <- nicheRecenter(position, breadth, sp1, env, phylo.sig)
            recenteredsp2 <- nicheRecenter(position, breadth, sp2, env, phylo.sig)
            species.rasters <- replace(species.rasters, next.free.edgetable, recenteredsp1[[3]])
            species.rasters <- replace(species.rasters, next.free.edgetable+1, recenteredsp2[[3]])
            newbr1 <- c(edgetable[i, 2], next.free.edgetable, stepsize, NA,
                        sum(recenteredsp1[[3]]@data@values, na.rm=T), time, recenteredsp1[[1]],
                        recenteredsp1[[2]], "allopatric", "X")
            newbr2 <- c(edgetable[i, 2], next.free.edgetable+1, stepsize, NA,
                        sum(recenteredsp2[[3]]@data@values, na.rm=T), time, recenteredsp2[[1]],
                        recenteredsp2[[2]], "allopatric", "X")
            edgetable[c(next.free.edgetable,next.free.edgetable+1),] <- rbind(newbr1, newbr2)
          }
        }
        # Sympatric speciation
        if(event == 2) {
          edgetable[i, 4] <- "sympatric"
          sympatric.speciation <- speciateSympatric(position, breadth, species.rasters[[i]], env, phylo.sig)
          species.rasters <- replace(species.rasters, next.free.edgetable, species.rasters[[i]])
          species.rasters <- replace(species.rasters, next.free.edgetable+1, sympatric.speciation$species.rasters[[1]])
          newbr1 <- c(edgetable[i, 2], next.free.edgetable, stepsize, NA,
                      sum(species.rasters[[i]]@data@values, na.rm=T), time, position,
                      breadth, "sympatric", "X")
          newbr2 <- c(edgetable[i, 2], next.free.edgetable+1, stepsize, NA,
                      sum(sympatric.speciation$species.rasters[[1]]@data@values, na.rm=T), time, sympatric.speciation$pos[[1]],
                      sympatric.speciation$breadth[[1]], "sympatric", "X")
          edgetable[c(next.free.edgetable,next.free.edgetable+1),] <- rbind(newbr1, newbr2)
        }

        # Parapatric Speciation
        if(event == 3) {
          edgetable[i, 4] <- "parapatric"
          position <- as.numeric(edgetable[i,7])
          breadth <- as.numeric(edgetable[i,8])
          parapatric.speciation <- speciateParapatric(position, breadth, species.rasters[[i]], env, phylo.sig, dispersal = dispersal)
          species.rasters <- replace(species.rasters, next.free.edgetable,species.rasters[[i]])
          species.rasters <- replace(species.rasters, next.free.edgetable+1,parapatric.speciation$species.rasters[[1]])
          newbr1 <- c(edgetable[i, 2], next.free.edgetable, stepsize, NA,
                      sum(species.rasters[[i]]@data@values, na.rm=T), time, position,
                      breadth, "parapatric", "X")
          newbr2 <- c(edgetable[i, 2], next.free.edgetable + 1, stepsize, NA,
                      sum(parapatric.speciation$species.rasters[[1]]@data@values, na.rm=T), time, parapatric.speciation$pos[[1]],
                      parapatric.speciation$breadth[[1]], "parapatric", "X")
          edgetable[c(next.free.edgetable,next.free.edgetable+1),] <- rbind(newbr1, newbr2)
        }
        # Dispersal speciation
        if(event == 4) {
          edgetable[i, 4] <- "dispersal"
          position <- as.numeric(edgetable[i,7])
          breadth <- as.numeric(edgetable[i,8])
          dispersal.speciation <- speciateDispersal(position, breadth, species.rasters[[i]], env, phylo.sig, dispersal)
          species.rasters <- replace(species.rasters, next.free.edgetable, species.rasters[[i]])
          species.rasters <- replace(species.rasters, next.free.edgetable+1, dispersal.speciation$species.rasters[[1]])
          newbr1 <- c(edgetable[i, 2], next.free.edgetable, stepsize, NA,
                      sum(species.rasters[[i]]@data@values, na.rm=T), time, position,
                      breadth, "dispersal", "X")
          newbr2 <- c(edgetable[i, 2], next.free.edgetable +1, stepsize, NA,
                      sum(dispersal.speciation$species.rasters[[1]]@data@values, na.rm=T), time, dispersal.speciation$pos[[1]],
                      dispersal.speciation$breadth[[1]], "dispersal", "X")
          edgetable[c(next.free.edgetable,next.free.edgetable+1),] <- rbind(newbr1, newbr2)
        }
        # No-Event (niche evolution)
        if(event == 5) {
          edgetable[i,3] <- as.numeric(edgetable[i,3]) + stepsize
          # we evolve the niche at each timestep
          niche <- nicheEvolution(position, breadth, pos.ev.rate = niche.ev.rate, breadth.ev.rate=breadth.ev.rate, species.raster =current.species, env=env)
          edgetable[i,7] <- niche[[2]]
          edgetable[i,8] <- niche[[3]]
          species.rasters[[i]] <- niche[[1]]
          edgetable[i, 5] <- sum(species.rasters[[i]]@data@values, na.rm=T)
        }
        # Extinction
        if (event == 6) {
          extinct[i] <- TRUE
          edgetable[i, 10] <- "extinct"
        }
      }#end of for loop

    # if all species are extinct reset to the starting values of the simulation
      if (all(edgetable[, 10][which(!edgetable[, 2] %in% edgetable[, 1])] =="extinct")) {
        initial.species <- seedSpecies(env, dispersal=dispersal)
        edgetable <- matrix(ncol=10, nrow=10000)
        edgetable[1,] <- c(0, 1, 0, NA, NA, 1, NA, NA, NA, "X")
        edgetable[1,5] <- sum(initial.species[[1]]@data@values, na.rm=T)
        edgetable[1,7] <- initial.species[[2]]
        edgetable[1,8] <- breadth
        species.rasters <- vector('list', 10000)
        species.rasters[[1]] <- initial.species[[1]]
        time <- 1
        stepsize <- stepsize
        tips <- 1
        extinct <- vector("logical", 10000)
        # Record mass extinction event. When 5 mass extinction events occur (extinction =5) in a row the simulation will sample new parameters and start again
        extinction <- extinction + 1
        extinct.number <- 0
        if(extinction == 5) {return(list(params, "too many extinctions with these parameters"))}
        next
      }
      extinct.species <- which(extinct==TRUE)
      extinct.number <- length(extinct.species)
      terminal.branches <- which(!edgetable[, 2] %in% edgetable[, 1])
      tips <- length(terminal.branches)
    }# end of while loop

  # species.rasters.2 are extant species rasters (extinct rasters are empty)
  species.rasters.2 <- species.rasters[terminal.branches[!terminal.branches %in% extinct.species]]
  # Build phylogeny
  phy.table <- buildPhyInSim(edgetable, species.rasters.2@layers, extinct, tips)

  # plot species ranges
  if (plot == TRUE ) {
    cols <- sample(rainbow(100,end=1, start=0, alpha=0.6 ), 100, replace = TRUE)
    par(mfrow = c(2, 1), mar=c(2,2,2,2))
    plot(env)
    for (i in 1:length(species.rasters.2)) {
      par(new=T)
      sp<-rasterToPolygons(species.rasters.2[[i]], dissolve=T)
      plot(sp, add=T, col=cols[i])
    }
    plot(phy.table[[2]]); axisPhylo()
  }

  species.rasters.2 <- stack(species.rasters.2)
  # record proportion of domain occupied by the whole clade
  clade.area <- sum(stackApply(species.rasters.2, indices = rep(1, length(species.rasters.2@layers)), fun=mean)@data@values, na.rm=T)
  #generate summary statistics (must have ENMTools loaded)
  if(generateSummaryStatistics == TRUE) {
    summaryStats <- try(generateSummaryStatistics(phy.table[[1]], species.rasters.2, phy.table[[2]], drop.fossil = T))
    ARC.clade <- createENMToolsclade(species.rasters.2, phy.table[[1]], phy.table[[2]], env, drop.fossil = T)
    ARC <- try(enmtools.aoc(clade = ARC.clade,  nreps = 0, overlap.source = "range"))
    if(class(summaryStats) == "try-error") {
      summary.vector=NA
    } else {
      if(class(ARC)=="try-error"){
        summary.vector = c(summaryStats$analysis, NA, NA)
        ARC <- NA
      } else {
        summary.vector = c(summaryStats$analysis, ARC$coefficients[[1,1]], ARC$coefficients[[1,2]])
      }
    }
    results <- list("df" = phy.table[[1]], "phy" = phy.table[[2]], "rasters" = species.rasters.2, "env" = env, "params"=params, "ARC"=ARC, "summaryStats"=summaryStats, "summaryVector"= summary.vector , "extinct"= length(which(extinct==TRUE)), "clade.area"= clade.area)
  } else {
  results <- list("df" = phy.table[[1]], "phy" = phy.table[[2]], "rasters" = species.rasters.2, "env" = env, "params"=params, "extinct"= length(which(extinct==TRUE)), "clade.area"=clade.area)
  }
# returns a list with the following elements 1) data frame with species information, 2) phylogeney, 3) environment (as was at end of simulation), 4) parameters used in that simulation
# 5) the Age-Range_correlation object crated by ENMTools, 6) a list of summary statistics (see generateSummaryStatistics function for details), 7) number of extinct lineages, 8) total area occupied by clade
  return(results)
}



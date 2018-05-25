### This version forked by Dan Rosauer to modify and prototype an alternative macroevolutionary simulation

Don't expect this code to work properly for some time
-----------------------------------------------------


DynamicRangeEvolutionAndDiversification (DREaD)
===============================================

### Alexander Skeels & Marcel Cardillo

Load scripts and packages
-------------------------

``` r
# set to directory containing source scripts
setwd()

scripts <- c("rangeDispersal.R", "nicheEvolution.R","speciateAllopatric.R","speciateSympatric.R",
             "speciateParapatric.R","speciateDispersal.R","seedSpecies.R","environmentalChange.R",
             "nicheRecenter.R","DREaD.R","generateSummaryStatistics.R", "helperFunctions.R")

lapply(scripts, source)

required.packages <- (c("raster","gstat","ape","phytools","geiger","phyloclim","ggplot2","gridExtra","moments",
                        "apTreeshape","parallel", "doSNOW", "rgeos", "data.table", "fossil", "ENMTools"))

lapply(required.packages, require, character.only=T)
```

Run Simulation
--------------

``` r
#sample parameters
tips <- runif(1, 5, 100)
D <- runif(1, 1, 10)
ENVa <- runif(1, 0.25, 2)
ENVf <- runif(1, 0.25, 2)
NEp <- runif(1, 0.005, 2)
NEb <- runif(1, 0.0025, 1)
PS <- runif(1, 0.25, 1)
m <- runif(1, 50, 250)

# Run model
simulation.1 <- DREaD(totaltips=tips, dispersal=D, amp=ENVa, freq=ENVf,
                             niche.ev.rate = NEp, breadth.ev.rate=NEb,
                             phylo.sig=PS, Me=m,  enviro.hetero =T, geo.mode = "dispersal",
                             enviro.mode = "sine")
```

Results
-------

``` r
# model parameters
print(simulation.1$params)
```

    ##   totaltips dispersal amp freq niche.ev.rate breath.ev.rate phylo.sig  Me
    ## 1         8         3 1.5  1.5             1            0.5       0.5 100
    ##    geo.mode slope enviro.mode enviro.hetero
    ## 1 dispersal     1        sine         FALSE

``` r
# measured summary statistics
print(simulation.1$summaryStats$analysis[1:30])
```

    ##   RO0 RO50 RO75 RO90 RO100    ROmean   asymmean bimod100 bimod90 bimod75
    ## 1 0.5    0    0    0     0 0.1036866 0.02441617        0       0       0
    ##   bimod50 ROkurt ROskew     TOmean      TOsd RSskew   RSmean      RSsd
    ## 1       0      1      0 -0.4493919 0.7786774 1.1646 0.261453 0.4049734
    ##        RDmean       RDsd  ROslope ROintercept asymslope asymintercept
    ## 1 0.007142495 0.01010101 16.52074 -0.06912442  3.483402   -0.01202109
    ##     RDslope RDintercept      Beta CI SI    gamma
    ## 1 -1.138038  0.01904665 -1.243582 11 29 2.321321

``` r
# Age-Range-Correlation
print(simulation.1$ARC$reps$rep.1)
```

    ## 
    ## Call:
    ## lm(formula = overlap ~ age, data = rep.df)
    ## 
    ## Coefficients:
    ## (Intercept)          age  
    ##      0.3156      -0.1608

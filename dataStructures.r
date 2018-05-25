
makeEdgeTable <- function(rowcount=10000,
                          dynamicSpeciation=FALSE) {

  # this function creates and returns an edgetable matrix it is placed in a
  # separate function for clarity and flexibility

  #edgeTable is the matrix that stores data about the characteristics of each
  #species or branch, it's environmental niche, ancestry and desendents

  # the dynamicSpeciation argument, when true, requests columns in edgeTable
  # required for dynamic speciation, such as a speciesID
  # to maintain compatibility - for now - none of the original columns are changed
  # so references by column number are the same, and unused columns are retained.

  #in the future this function may include arguments for
  # a) number of environmental dimensions
  # b) whether to structure for instantaneous or dynamic speciation
  # c) whether to return column numbers

  edgetable       <- matrix(data = NA, nrow = rowcount, ncol = 10)
  edgetable[1, ]  <- c(0, 1, 0, NA, NA, 1, NA, NA, NA, 1)

  # add column names to edgetable

  column.names <- c("parent",
                    "daughter",
                    "branch.length",
                    "split.mode",
                    "range.size",
                    "time.of.speciation",
                    "niche.position",
                    "niche.breadth",
                    "birth.mode",
                    "extant")             # 1=extant species, 0=extinct
  colnames(edgetable) <- column.names

  if (dynamicSpeciation) {
    speciesID <- rep(NA, rowcount)
    edgetable <- cbind(edgetable, speciesID)
    column.names <- c("speciesID")  # a single row for each species in edgtable links to 1 row per occupuied
                      # cell in demeTable

    colnames(edgetable)[ncol(edgetable)] <- column.names  # name the last column
  }

  return((edgetable))

}


makeDemeTable <- function(genomeDimensions=3, columnInfo=FALSE) {

  # this function creates and returns an empty demes matrix
  # it is placed in a separate function for clarity and flexibility

  #demesTable is the matrix that stores data about the location and
  #characteristics of each population or deme, and their changes over time as
  #the simulation progresses. it is linked to EdgeTable by the speciesID, and
  #to the region raster by cellID

  #in the future this function may include arguments for
  # a) number of environmental dimensions

  demestable <- data.frame(cellID=0,
                           speciesID=0,
                           x=0.0,
                           y=0.0,
                           amount=0.0,
                           niche1.position=0.0,
                           niche1.breadth=0.0,
                           niche2.position=0.0,
                           niche2.breadth=0.0,
                           gene.pos1=0.0,
                           gene.pos2=0.0)

  if (genomeDimensions > 2) {
    for (i in 3:genomeDimensions) {
      genecol.df <- data.frame(gene.pos=0.0)
      names(genecol.df) <- paste("gene.pos", i, sep="")
      demestable <- cbind(demestable, genecol.df)
    }
  }

  # returns vectors of column numbers for niche and genetic parameters for easy
  # reference to the columns
  if (columnInfo) {
    niche.position.columns <- 6:9
    gene.position.columns  <- 10:ncol(demestable)
    columns <- list(niche.position.columns=niche.position.columns,
                    gene.position.columns=gene.position.columns)
    return(columns)
  } else {
    return(demestable)
  }
}


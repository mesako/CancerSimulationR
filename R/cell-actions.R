
MarkDeadCells <- function(cell.states, current.map) {
  updated.map <- current.map
  dead.ind <- which(`dim<-`(grepl(cell.states, pattern = "F"),
                            dim(cell.states)), arr.ind = TRUE)
  if (nrow(dead.ind) > 0) {
    for (i in 1:nrow(dead.ind)) {
      this.coord <- paste0(dead.ind[i, ], collapse = ".")
      this.ind <- which(updated.map$data[, "id"] == this.coord)
      updated.map$data[this.ind, "value"] <- "None"
    }
    cell.states[dead.ind] <- NA
  }
  return(list(updated.map, cell.states))
}

RunBaselineCellActions <- function(current.map, cell.states, cell.divisions, max.divisions) {
  results <- MarkDeadCells(cell.states, current.map)
  current.map <- results[[1]]
  cell.states <- results[[2]]
  results <- RunCellDivision(current.map, cell.states, cell.divisions, max.divisions)
  current.map <- results[[1]]
  cell.states <- results[[2]]
  cell.divisions <- results[[3]]
  results <- CheckEnergyNeed(current.map, cell.states)
  current.map <- results[[1]]
  cell.states <- results[[2]]
  return(list(current.map, cell.states, cell.divisions))
}

CheckEnergyNeed <- function(current.map, cell.states) {
  updated.map <- current.map
  blood.ind <- current.map$data$id[current.map$data$value == "Blood"]
  blood.ind.num <- as.numeric(as.character(blood.ind))
  if (length(unique(floor(blood.ind.num))) > 1) {
    neighbors <- blood.ind.num - 0.1
    neighbors <- c(neighbors, (blood.ind.num + 0.1))
    neighbors <- c(neighbors, (blood.ind.num + 0.2))
    neighbors <- c(neighbors, (blood.ind.num - 0.2))
  } else {
    neighbors <- blood.ind.num - 1
    neighbors <- c(neighbors, (blood.ind.num + 1))
    neighbors <- c(neighbors, (blood.ind.num + 2))
    neighbors <- c(neighbors, (blood.ind.num - 2))
  }
  neighbors <- as.character(neighbors)
  neighbors[duplicated(neighbors)] <- paste(neighbors[duplicated(neighbors)], "0", sep = "")
  non.neighbors <- which(!(current.map$data$id %in% neighbors))
  all.cells <- which(current.map$data$value == "Cell")
  remove.candidates <- intersect(non.neighbors, all.cells)
  energy.ind1 <- which(`dim<-`(grepl(cell.states, pattern = mutation.encoding$metabolic),
                              dim(cell.states)), arr.ind = TRUE)
  energy.ind2 <- which(`dim<-`(grepl(cell.states, pattern = mutation.encoding$angio),
                              dim(cell.states)), arr.ind = TRUE)
  energy.ind <- rbind(energy.ind1, energy.ind2)
  remove.candidates <- current.map$data$id[remove.candidates]
  if (length(energy.ind) > 0) {
    energy.ind <- paste(energy.ind[seq(1, length(energy.ind), by = 2)],
                        energy.ind[seq(2, length(energy.ind), by = 2)], sep = ".")
    remove.candidates <- setdiff(remove.candidates, energy.ind)
  }
  
  if (length(remove.candidates) > 0) {
    num.to.remove <- round(length(remove.candidates) * 0.2)
    to.remove <- sample(remove.candidates, num.to.remove,  replace = FALSE)
    remove.ind <- which(updated.map$data[, "id"] %in% to.remove)
    updated.map$data[remove.ind, "value"] <- "None"
    for (i in 1:length(to.remove)) {
      this.ind <- to.remove[i]
      this.ind <- as.character(this.ind)
      this.ind <- unlist(strsplit(this.ind, split = "\\."))
      this.ind <- as.numeric(this.ind)
      cell.states[this.ind[1], this.ind[2]] <- NA
    }
  }
  return(list(updated.map, cell.states))
}

ChangeCellMutationRate <- function(cell.mut.rate, which.cell.mut) {
  changed.cell.mut.rate <- cell.mut.rate
  changed.cell.mut.rate[which.cell.mut] <- 2 * changed.cell.mut.rate[which.cell.mut]
  return(changed.cell.mut.rate)
}

FindNeighborCell <- function(current.map, cell.coord) {
  cell.coord <- as.character(cell.coord)
  cell.coord.digits <- unlist(strsplit(cell.coord, split = "\\."))
  first.digit <- cell.coord.digits[1]
  second.digit <- cell.coord.digits[2]
  first.digit <- as.numeric(first.digit)
  second.digit <- as.numeric(second.digit)
  possible.first.digit <- sort(c(first.digit + 1, first.digit, first.digit - 1))
  possible.second.digit <- sort(c(second.digit + 1, second.digit, second.digit - 1))
  possible.options <- c()
  for (x in possible.first.digit) {
    new.add <- paste(x, ".", possible.second.digit, sep = "")
    possible.options <- c(possible.options, new.add)
  }
  possible.options <- setdiff(possible.options, cell.coord)
  possible.options <- possible.options[possible.options %in% as.character(current.map$data$id)]
  return(possible.options)
}


ReturnEmptySpot <- function(current.map, possible.options) {
  find.ind <- as.character(current.map$data$id) %in% possible.options
  find.ind <- which(find.ind == TRUE)
  valid.ind <- which(current.map$data$value == "None")
  match.ind <- intersect(find.ind, valid.ind)
  match.ind <- current.map$data$id[match.ind]
  return(match.ind)
}

RunCellDivision <- function(current.map, cell.states, cell.divisions, max.divisions) {
  extract.dim <- as.numeric(current.map$data$id)
  extract.dim <- floor(max(extract.dim))
  temp1 <- matrix(current.map$data$id, nrow = extract.dim, byrow = TRUE)
  temp1 <- temp1[c(nrow(temp1):1), , drop = FALSE]
  temp2 <- matrix(current.map$data$value, nrow = extract.dim, byrow = TRUE)
  temp2 <- temp2[c(nrow(temp2):1), , drop = FALSE]
  updated.map <- current.map
  cellcycle <- which(`dim<-`(grepl(cell.states, pattern = mutation.encoding$cellcycle),
                             dim(cell.states)), arr.ind = TRUE)
  if (nrow(cellcycle) > 0) {
    for (i in 1:nrow(cellcycle)) {
      this.row <- cellcycle[i, ]
      temp <- cell.states[this.row[1], this.row[2]]
      if (cell.divisions[this.row[1], this.row[2]] < max.divisions | grepl(temp, pattern = mutation.encoding$telomerase)) {
        num.div <- 1
        if (grepl(temp, pattern = mutation.encoding$growthact)) {
          num.div <- 2
        }
      } else {
        num.div <- 0
      }
      cell.divisions[this.row[1], this.row[2]] <- cell.divisions[this.row[1], this.row[2]] + num.div
      cell.coord <- paste(this.row[1], this.row[2], sep = ".")
      results <- MakeNewCells(updated.map, cell.states, cell.divisions, cell.coord, num.div)
      updated.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
    }
  }
  return(list(updated.map, cell.states, cell.divisions))
}

MakeNewCells <- function(current.map, cell.states, cell.divisions, cell.coord, num.div) {
  updated.map <- current.map
  possible.options <- FindNeighborCell(current.map, cell.coord)
  possible.options <- ReturnEmptySpot(current.map, possible.options)
  if (length(possible.options) > num.div) {
    chosen.spot <- sample(possible.options, num.div, replace = FALSE)
  } else {
    chosen.spot <- possible.options
  }
  updated.map$data[updated.map$data$id %in% chosen.spot, "value"] <- "Cell"
  overwrite.state <- GetOriginalData(cell.states, cell.divisions, cell.coord)[[1]]
  overwrite.division <- GetOriginalData(cell.states, cell.divisions, cell.coord)[[2]]
  
  for (each in chosen.spot) {
    x <- as.numeric(unlist(strsplit(each, split = "\\."))[1])
    y <- as.numeric(unlist(strsplit(each, split = "\\."))[2])
    cell.states[x, y] <- overwrite.state
    cell.divisions[x, y] <- overwrite.division
  }
  return(list(updated.map, cell.states, cell.divisions))
}

GetOriginalData <- function(cell.states, cell.divisions, cell.coord) {
  x <- as.numeric(unlist(strsplit(as.character(cell.coord), split = "\\."))[1])
  y <- as.numeric(unlist(strsplit(as.character(cell.coord), split = "\\."))[2])
  this.state <- cell.states[x, y]
  this.division <- cell.divisions[x, y]
  return(list(this.state, this.division))
}

RunCellMove <- function(cell.states) {
  grep(cell.states, pattern = mutation.encoding$emt)
  return()
}

RollMutations <- function(map.dim, cell.mut.rate) {
  dice.roll <- runif(n = (map.dim[1] * map.dim[2]), min = 0, max = 1)
  match.ind <- which(dice.roll < as.vector(cell.mut.rate))
  return(match.ind)
}

GetMutations <- function(match.ind, mutation.set) {
  which.mut <- runif(n = length(match.ind), min = 0, max = 1)
  temp.prob <- cumsum(mutation.set)
  temp.prob <- c(0, temp.prob)
  temp.prob <- temp.prob[1:(length(temp.prob) - 1)]
  names(temp.prob) <- names(mutation.set)
  round.mut <- names(mutation.set)[findInterval(which.mut, temp.prob)]
  return(round.mut)
}

EncodeMutation <- function(mutation.encoding, new.mutations) {
  short.names <- c()
  for (i in 1:length(new.mutations)) {
    short.names <- c(short.names, mutation.encoding[[new.mutations[i]]])
  }
  return(short.names)
}

AddMutations <- function(cell.states, cell.mut.rate,
                         map.dim, mutation.encoding) {
  match.ind <- RollMutations(map.dim, cell.mut.rate)
  updated.cell.states <- cell.states
  if (length(match.ind) > 0) {
    new.mutations <- GetMutations(match.ind, mutation.set)
    new.mutations <- EncodeMutation(mutation.encoding, new.mutations)
    is.na.test <- !is.na(updated.cell.states[match.ind])
    is.mutated <- match.ind[is.na.test]
    not.mutated <- match.ind[!is.na.test]
    updated.cell.states[not.mutated] <- new.mutations[!is.na.test]
    updated.cell.states[is.mutated] <- paste(updated.cell.states[is.mutated],
                                             new.mutations[is.na.test], sep = "")
  }
  return(updated.cell.states)
}

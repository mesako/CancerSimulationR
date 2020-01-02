
SwitchPrePost <- function(cell.states, current.map, mutation.encoding) {
  metabolic.ind <- grep(cell.states, pattern = mutation.encoding$metabolic)
  num.cells <- sum(current.map$data$value == "Cell")
  if (length(metabolic.ind) > round(nrow(current.map$data) * 0.1)) {
    print("metabolic signal exceeded")
    patient.state <- "post"
  } else if (num.cells > round(nrow(current.map$data) * 0.25)) {
    print("cell number threshold exceeded")
    patient.state <- "post"
  } else {
    patient.state <- "pre"
  }
  return(patient.state)
}

KillCellsCycling <- function(cell.states, current.map, cell.divisions, cell.mut.rate, mutation.encoding) {
  updated.map <- current.map
  cycling.ind <- which(`dim<-`(grepl(cell.states, pattern = mutation.encoding$cellcycle),
                               dim(cell.states)), arr.ind = TRUE)
  if (nrow(cycling.ind) > 0) {
    kill.these <- runif(n = nrow(cycling.ind), min = 0, max = 1)
    kill.these <- kill.these < 0.33
    if (sum(kill.these) > 1) {
      kill.these <- cycling.ind[kill.these, ]
      for (i in 1:nrow(kill.these)) {
        this.coord <- paste0(kill.these[i, ], collapse = ".")
        this.ind <- which(updated.map$data[, "id"] == this.coord)
        updated.map$data[this.ind, "value"] <- "None"
      }
      cell.states[kill.these] <- NA
      cell.divisions[kill.these] <- NA
      cell.mut.rate[kill.these] <- 0
    }
  }
  # This drug disrupts the cytoskeleton, causing
  # cells to fail to divide or die if they are in the
  # middle of dividing. All cells that divided in the
  # last round must roll for chance of survival. No
  # new cell divisions are allowed during this round.

  return(list(updated.map, cell.states, cell.divisions, cell.mut.rate))
}

KillCellsDNA <- function(cell.states, current.map, cell.divisions, cell.mut.rate, mutation.encoding) {
  updated.map <- current.map
  repair.ind <- which(`dim<-`(!grepl(cell.states, pattern = mutation.encoding$dnarepair),
                              dim(cell.states)), arr.ind = TRUE)
  rearrange.names <- paste(repair.ind[, 1], repair.ind[, 2], sep = ".")
  only.cells <- current.map$data[current.map$data$value == "Cell", "id"]
  only.cells <- as.character(only.cells)
  repair.ind <- repair.ind[which(rearrange.names %in% only.cells), ]
  if (class(repair.ind) == "integer") {
    kill.flag <- runif(n = 1, min = 0, max = 1)
    kill.flag <- kill.flag < 0.5
    if (kill.flag) {
      this.coord <- paste0(repair.ind, collapse = ".")
      this.ind <- which(updated.map$data[, "id"] == this.coord)
      updated.map$data[this.ind, "value"] <- "None"
      cell.states[repair.ind[1], repair.ind[2]] <- NA
      cell.divisions[repair.ind[1], repair.ind[2]] <- NA
      cell.mut.rate[repair.ind[1], repair.ind[2]] <- 0
    }
  } else if (class(repair.ind) == "matrix") {
    if (nrow(repair.ind) > 0) {
      kill.these <- runif(n = nrow(repair.ind), min = 0, max = 1)
      kill.these <- kill.these < 0.5
      if (sum(kill.these) > 1) {
        kill.these <- repair.ind[kill.these, ]
        for (i in 1:nrow(kill.these)) {
          this.coord <- paste0(kill.these[i, ], collapse = ".")
          this.ind <- which(updated.map$data[, "id"] == this.coord)
          updated.map$data[this.ind, "value"] <- "None"
        }
        cell.states[kill.these] <- NA
        cell.divisions[kill.these] <- NA
        cell.mut.rate[kill.these] <- 0
      }
    }
  }
  # This therapy damages DNA, inducing death. All
  # cells without the DNA Repair Mutation must
  # roll for chance of survival. No new cell divisions
  # are allowed during this round. Any cells that
  # survive or have DNA Repair Mutation must roll
  # for a chance to mutate.

  return(list(updated.map, cell.states, cell.divisions, cell.mut.rate))
}

KillCellsSpace <- function(cell.states, current.map, cell.divisions, cell.mut.rate, mutation.encoding) {
  resist.ind <- which(`dim<-`(grepl(cell.states, pattern = mutation.encoding$deathinhibit),
                              dim(cell.states)), arr.ind = TRUE)
  resistant.cells <- paste(resist.ind[, 1], resist.ind[, 2], sep = ".")
  crowded.cells <- CheckSimilarNeighbors(current.map, search.type = "Cell")
  possible.options <- setdiff(crowded.cells, resistant.cells)
  kill.these <- runif(n = length(possible.options), min = 0, max = 1)
  kill.these <- kill.these < 0.33
  kill.these <- possible.options[kill.these]
  updated.map <- current.map
  remove.ind <- which(updated.map$data[, "id"] %in% kill.these)
  updated.map$data[remove.ind, "value"] <- "None"
  for (i in 1:length(kill.these)) {
    this.ind <- kill.these[i]
    this.ind <- as.character(this.ind)
    this.ind <- unlist(strsplit(this.ind, split = "\\."))
    this.ind <- as.numeric(this.ind)
    cell.states[this.ind[1], this.ind[2]] <- NA
    cell.divisions[this.ind[1], this.ind[2]] <- NA
    cell.mut.rate[this.ind[1], this.ind[2]] <- 0
  }
  # Normal cells in the surrounding lung will signal
  # to the cancerous cell to stop dividing or commit
  # apoptosis as they begin to compete for space.
  # All cells without the Death Inhibition Mutation
  # must roll for chance of survival.

  return(list(updated.map, cell.states, cell.divisions, cell.mut.rate))
}


KillCellsEnergy <- function(cell.states, current.map, cell.divisions, cell.mut.rate, mutation.encoding) {
  updated.map <- current.map
  angio.ind <- which(`dim<-`(grepl(cell.states, pattern = mutation.encoding$angio),
                             dim(cell.states)), arr.ind = TRUE)
  metabolic.ind <- which(`dim<-`(grepl(cell.states, pattern = mutation.encoding$metabolic),
                                 dim(cell.states)), arr.ind = TRUE)
  angio.cells <- paste(angio.ind[, 1], angio.ind[, 2], sep = ".")
  metabolic.cells <- paste(metabolic.ind[, 1], metabolic.ind[, 2], sep = ".")
  possible.options <- setdiff(angio.cells, metabolic.cells)
  if (length(possible.options) > 0) {
    kill.these <- runif(n = length(possible.options), min = 0, max = 1)
    kill.these <- kill.these < 0.5
    if (sum(kill.these) > 1) {
      kill.these <- possible.options[kill.these]
      for (i in 1:length(kill.these)) {
        this.ind <- which(updated.map$data[, "id"] == kill.these[i])
        updated.map$data[this.ind, "value"] <- "None"
        this.ind <- unlist(strsplit(as.character(kill.these[i]), split = "\\."))
        this.ind <- as.numeric(this.ind)
        cell.states[this.ind[1], this.ind[2]] <- NA
        cell.divisions[this.ind[1], this.ind[2]] <- NA
        cell.mut.rate[this.ind[1], this.ind[2]] <- 0
      }
    }
  }
  # Blood vessels that were created through the
  # Angiogenesis Mutation are destroyed. If this
  # results in cells being too far from a source of
  # food (blood vessel), these cells must roll for a
  # chance of survival unless they have a
  # Metabolic Mutation.

  return(list(updated.map, cell.states, cell.divisions, cell.mut.rate))
}

KillCellsGrowth <- function(cell.states, current.map, cell.divisions, cell.mut.rate, mutation.encoding) {
  updated.map <- current.map
  grow.ind <- which(`dim<-`(grepl(cell.states, pattern = mutation.encoding$growthact),
                            dim(cell.states)), arr.ind = TRUE)
  if (nrow(grow.ind) > 0) {
    kill.these <- runif(n = nrow(grow.ind), min = 0, max = 1)
    kill.these <- kill.these < 0.33
    if (sum(kill.these) > 1) {
      kill.these <- grow.ind[kill.these, ]
      for (i in 1:nrow(kill.these)) {
        this.coord <- paste0(kill.these[i, ], collapse = ".")
        this.ind <- which(updated.map$data[, "id"] == this.coord)
        updated.map$data[this.ind, "value"] <- "None"
      }
      cell.states[kill.these] <- NA
      cell.divisions[kill.these] <- NA
      cell.mut.rate[kill.these] <- 0
    }
  }
  # This inhibitor specifically binds to a signaling
  # protein that activates growth and blocks it from
  # keeping the cell alive. All cells with the Growth
  # Activation Mutation must roll for chance of survival

  return(list(updated.map, cell.states, cell.divisions, cell.mut.rate))
}

CheckSimilarNeighbors <- function(current.map, search.type = "Cell") {
  extract.dim <- as.numeric(current.map$data$id)
  extract.dim <- floor(max(extract.dim))
  map.types <- matrix(current.map$data$value, nrow = extract.dim, byrow = TRUE)
  map.types <- map.types[nrow(map.types):1, ]
  map.types <- as.data.frame(map.types)
  surrounded.tiles <- c()
  for (i in current.map$data$id) {
    this.ind <- unlist(strsplit(as.character(i), split = "\\."))
    this.ind <- as.numeric(this.ind)
    all.neighbors.row <- c(this.ind[1] - 1, this.ind[1],
                           this.ind[1] + 1)
    all.neighbors.col <- c(this.ind[2] - 1, this.ind[2],
                           this.ind[2] + 1)
    coord.pairs <- expand.grid(all.neighbors.row, all.neighbors.col)
    coord.pairs <- paste(coord.pairs[, 1], coord.pairs[, 2], sep = ".")
    coord.pairs <- intersect(coord.pairs, as.character(current.map$data$id))
    if (length(coord.pairs) < 5) {
      next
    } else {
      all.neighbors <- current.map$data[as.character(current.map$data$id) %in% coord.pairs, "value"]
      all.neighbors <- as.character(all.neighbors)
      if (sum(all.neighbors == search.type) == length(coord.pairs)) {
        surrounded.tiles <- c(surrounded.tiles, coord.pairs)
      }
    }
  }
  surrounded.tiles <- unique(surrounded.tiles)
  return(surrounded.tiles)
}

KillCellsSurgery <- function(cell.states, current.map, cell.divisions, cell.mut.rate) {
  updated.map <- current.map
  to.remove <- CheckSimilarNeighbors(current.map, search.type = "Cell")
  remove.ind <- which(updated.map$data[, "id"] %in% to.remove)
  updated.map$data[remove.ind, "value"] <- "None"
  for (i in 1:length(to.remove)) {
    this.ind <- to.remove[i]
    this.ind <- as.character(this.ind)
    this.ind <- unlist(strsplit(this.ind, split = "\\."))
    this.ind <- as.numeric(this.ind)
    cell.states[this.ind[1], this.ind[2]] <- NA
    cell.divisions[this.ind[1], this.ind[2]] <- NA
    cell.mut.rate[this.ind[1], this.ind[2]] <- 0
  }
  # Surgery will remove cancerous cells from the
  # body (removing them from the game). If there is
  # a density of cells in a certain area (above some
  # threshold), all cells in that area are killed.

  return(list(updated.map, cell.states, cell.divisions, cell.mut.rate))
}

GetPatientAction <- function(patient.state, patient.pre.action.set,
                             patient.post.action.set) {
  if (patient.state == "pre") {
    temp.prob <- cumsum(patient.pre.action.set)
    temp.prob <- c(0, temp.prob)
    temp.prob <- temp.prob[1:(length(temp.prob) - 1)]
    names(temp.prob) <- names(patient.pre.action.set)
    round.action <- names(patient.pre.action.set)[findInterval(runif(1, 0, 1), temp.prob)]
  } else if (patient.state == "post") {
    temp.prob <- cumsum(patient.post.action.set)
    temp.prob <- c(0, temp.prob)
    temp.prob <- temp.prob[1:(length(temp.prob) - 1)]
    names(temp.prob) <- names(patient.post.action.set)
    round.action <- names(patient.post.action.set)[findInterval(runif(1, 0, 1), temp.prob)]
  } else {
    stop("!!!")
  }
  return(round.action)
}

CheckPatientDeath <- function(cell.states, mutation.encoding) {
  required.mut <- c("emt", "telomerase", "growthact",
                    "deathinhibit", "cellcycle", "metastasis")
  first.mut <- mutation.encoding[[required.mut[1]]]
  cell.meet.cond <- grep(cell.states, pattern = first.mut)
  if (length(cell.meet.cond) > 0) {
    for (i in 2:length(required.mut)) {
      mut.encode <- mutation.encoding[[required.mut[i]]]
      meets.cond <- grep(cell.states, pattern = mut.encode)
      if (length(meets.cond) < 1) {
        cell.meet.cond <- NULL
        break
      } else {
        cell.meet.cond <- intersect(cell.meet.cond, meets.cond)
      }
      if (length(cell.meet.cond) < 1) {
        cell.meet.cond <- NULL
        break
      }
    }
  }
  # This cell must have the EMT Mutation,
  # Telomerase Mutation, Growth Activation
  # Mutation, Death Inhibition Mutation, and be
  # surrounded by cells. If this cell meets these
  # criteria, the cell is able to leave the organ and
  # spread to other organs, killing this patient.

  return(cell.meet.cond)
}

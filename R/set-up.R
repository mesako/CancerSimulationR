
GenerateMap <- function(map.dim) {
  if (sum(map.dim < 7)) {
    stop("Map dimensions are too small!")
  } else if (map.dim[1] != map.dim[2]) {
    stop("Map must be symmetrical/square!")
  }
  ids <- c()
  for (i in map.dim[2]:1) {
    for (j in 1:map.dim[1]) {
      ids <- c(ids, paste(i, ".", j, sep = ""))
    }
  }
  ids <- as.factor(ids)
  # print(ids)
  options <- factor(c("Cell", "None", "Blood"))
  value <- rep(options[2], times = length(ids))
  values <- data.frame(id = ids, value)
  id <- rep(ids, each = 4)
  
  CreateBloodVessel <- function(values) {
    fixed.values <- values
    to.pick <- fixed.values$id
    to.pick <- as.character(to.pick)
    to.pick <- unlist(strsplit(to.pick, split = "\\."))
    to.pick <- sort(as.numeric(unique(to.pick)))
    chosen.range <- to.pick[3:(length(to.pick) - 3)]
    chosen.coord <- sample(chosen.range, size = 1)
    which.direction <- sample(c("vertical", "horizontal"), size = 1)
    if (which.direction == "vertical") {
      replace.these <- paste(to.pick, ".", chosen.coord, sep = "")
    } else {
      replace.these <- paste(chosen.coord, ".", to.pick, sep = "")
    }
    fixed.values[fixed.values$id %in% replace.these, ]$value <- "Blood"
    return(fixed.values)
  }
  
  CreateCancerCells <- function(values) {
    blood.ind <- which(values$value == "Blood")
    blood.ind <- values$id[blood.ind]
    original.options <- as.character(blood.ind)
    first.digit <- unlist(strsplit(original.options, split = "\\."))
    first.digit <- first.digit[seq(1, length(first.digit), by = 2)]
    first.digit <- unique(first.digit)
    second.digit <- unlist(strsplit(original.options, split = "\\."))
    second.digit <- second.digit[seq(2, length(second.digit), by = 2)]
    first.digit <- as.numeric(first.digit)
    second.digit <- as.numeric(second.digit)
    possible.first.digit <- sort(c(first.digit + 1, first.digit - 1))
    possible.second.digit <- sort(c(second.digit + 1, second.digit - 1))
    possible.options <- c()
    for (x in possible.first.digit) {
      new.add <- paste(x, ".", possible.second.digit, sep = "")
      possible.options <- c(possible.options, new.add)
    }
    possible.options <- setdiff(possible.options, original.options)
    possible.options <- possible.options[possible.options %in% as.character(values$id)]
    change.these <- sample(possible.options, size = round(length(possible.options) * 0.5))
    fixed.values <- values
    fixed.values[fixed.values$id %in% change.these, ]$value <- "Cell"
    return(fixed.values)
  }
  
  values <- CreateBloodVessel(values)
  values <- CreateCancerCells(values)
  
  GenerateCoord <- function(id, map.dim) {
    base.x <- c(1, 0, 0, 1)
    x <- base.x
    if (map.dim[1] > 1) {
      for (i in 1:(map.dim[1] - 1)) {
        new.x <- c(base.x + i)
        x <- c(x, new.x)
      }
    }
    base.y <- c(0, 0, 1, 1)
    y <- rep(base.y, times = map.dim[1])
    if (map.dim[2] > 1) {
      for (j in 1:(map.dim[2] - 1)) {
        new.y <- c(base.y + j)
        y <- c(y, rep(new.y, times = map.dim[1]))
      }
    }
    x <- rep(x, times = map.dim[2])
    positions <- data.frame(id, x, y)
    return(positions)
  }
  positions <- GenerateCoord(id, map.dim)
  label.x <- positions$x[seq(1, length(positions$x), by = 4)]
  label.y <- positions$y[seq(1, length(positions$y), by = 4)]
  label.x <- label.x - 0.5
  label.y <- label.y + 0.5
  my.plot <- ggplot(values, aes(fill = value, color = "black")) +
    geom_map(aes(map_id = id), map = positions) +
    expand_limits(positions) + scale_color_manual(values = c("black"))
  my.plot <- my.plot + theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_line(colour = "black"),
          panel.grid.major.y = element_line(colour = "black")) + 
    scale_x_continuous(limits = c(0, map.dim[1]), breaks = seq(0, map.dim[1], 1)) +
    scale_y_continuous(limits = c(0, map.dim[2]), breaks = seq(0, map.dim[2], 1))
  my.plot <- my.plot + scale_fill_manual(values = c("magenta", "yellow", "cyan"))
  my.plot <- my.plot + guides(color = FALSE)
  my.plot <- my.plot + geom_text(aes(label = id, x = label.x, y = label.y))
  return(my.plot)
}

GenerateDivisionMatrix <- function(map.dim, starting.map) {
  cell.divisions <- matrix(0, nrow = map.dim[1],
                           ncol = map.dim[2])
  save.ind <- which(starting.map$data$value != "Cell")
  cell.divisions[save.ind] <- NA
  cell.divisions <- t(cell.divisions)
  cell.divisions <- cell.divisions[c(nrow(cell.divisions):1), , drop = FALSE] 
  return(cell.divisions)
}

GenerateMutationMatrix <- function(map.dim, starting.map, mutation.rate) {
  cell.mut.rate <- matrix(mutation.rate, nrow = map.dim[1],
                          ncol = map.dim[2])
  save.ind <- which(starting.map$data$value != "Cell")
  cell.mut.rate[save.ind] <- 0
  cell.mut.rate <- t(cell.mut.rate)
  cell.mut.rate <- cell.mut.rate[c(nrow(cell.mut.rate):1), , drop = FALSE] 
  return(cell.mut.rate)
}

EstablishMutations <- function() {
  mutation.set <- c("noeffect", "cellcycle", "growthact", "deathinhibit",
                    "metabolic", "emt", "angio", "telomerase",
                    "fatal", "dnarepair", "metastasis")
  mutation.set.prob <- c(0.2, 0.1, 0.15, 0.1, 0.05, 0.1, 0.05,
                         0.1, 0.05, 0.05, 0.05)
  names(mutation.set.prob) <- mutation.set
  mutation.set <- mutation.set.prob
  return(mutation.set)
}

EstablishPatientPreActions <- function() {
  patient.pre.action.set <- c("noaction", "smoke", "asbestos", "radexpose",
                              "doctor", "immune", "growthinhibit",
                              "deathinduct")
  patient.pre.action.prob <- c(0.5, 0.05, 0.05, 0.05, 0.1, 0.1, 0.05, 0.1)
  names(patient.pre.action.prob) <- patient.pre.action.set
  patient.pre.action.set <- patient.pre.action.prob
  return(patient.pre.action.set)
}

EstablishPatientPostActions <- function() {
  patient.post.action.set <- c("noaction", "chemo1", "chemo2",
                               "radtherapy", "targeted", "angblock",
                               "surgery")
  patient.post.action.prob <- c(0.3, 0.15, 0.1, 0.1, 0.1, 0.05, 0.2)
  names(patient.post.action.prob) <- patient.post.action.set
  patient.post.action.set <- patient.post.action.prob
  return(patient.post.action.set)
}

MakeStartingCellStateMatrix <- function(map.dim) {
  starting.cell.states <- matrix(NA, nrow = map.dim[1],
                                 ncol = map.dim[2])
  return(starting.cell.states)
}

MakeMutationEncoding <- function(mutation.set) {
  abbrev <- substr(toupper(names(mutation.set)), 1, 1)
  fix.abbrev <- which(duplicated(abbrev))
  for (to.fix in fix.abbrev) {
    this.one <- names(mutation.set)[to.fix]
    avail.letters <- unlist(strsplit(toupper(this.one), split = ""))
    avail.letters <- setdiff(avail.letters, abbrev)
    if (length(avail.letters) > 0) {
      abbrev[to.fix] <- avail.letters[1]
    } else {
      abbrev[to.fix] <- setdiff(LETTERS, abbrev)[1]
    }
  }
  mutation.encoding <- list()
  for (i in 1:length(mutation.set)) {
    mutation.encoding[[i]] <- abbrev[i]
  }
  names(mutation.encoding) <- names(mutation.set)
  return(mutation.encoding)
}

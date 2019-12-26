
DeterminePatientAction <- function(current.map, patient.state, num.surgery,
                                   patient.pre.action.set, patient.post.action.set) {
  if (num.surgery < 2 & sum(current.map$data$value == "Cell") > round(nrow(current.map$data) * 0.66)) {
    if (patient.state == "post") {
      patient.action <- "surgery"
    } else {
      patient.action <- "doctor"
    }
  } else {
    patient.action <- GetPatientAction(patient.state, patient.pre.action.set,
                                       patient.post.action.set)
  }
  while (patient.action == "chemo1" | patient.action == "chemo2" | patient.action == "radtherapy" &
         sum(current.map$data$value == "Cell") < round(nrow(current.map$data) * 0.1)) {
    patient.action <- GetPatientAction(patient.state, patient.pre.action.set,
                                       patient.post.action.set)
  }
  while (patient.action == "surgery" & sum(current.map$data$value == "Cell") < round(nrow(current.map$data) * 0.25)) {
    patient.action <- GetPatientAction(patient.state, patient.pre.action.set,
                                       patient.post.action.set)
  }
  while (num.surgery >= 2 & patient.action == "surgery") {
    patient.action <- GetPatientAction(patient.state, patient.pre.action.set,
                                       patient.post.action.set)
  }
  # print(patient.action)
  return(patient.action)
}

RunPatientSimulation <- function(num.rounds, starting.map, mutation.set,
                                 patient.pre.action.set, patient.post.action.set,
                                 mutation.encoding, cell.mut.rate,
                                 cell.divisions, cell.states, max.divisions) {
  num.surgery <- 0
  patient.state <- "pre"
  current.map <- starting.map
  round.cell.num <- c()
  round.average.mut.num <- c()
  tumor.burden.count <- 0
  for (i in 1:num.rounds) {
    patient.action <- DeterminePatientAction(current.map, patient.state, num.surgery,
                                             patient.pre.action.set, patient.post.action.set)
    if (patient.action == "noaction") {
      cell.states <- AddMutations(cell.states, cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- RunBaselineCellActions(current.map, cell.states, cell.divisions, max.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
    } else if (patient.action == "doctor") {
      cell.states <- AddMutations(cell.states, cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- RunBaselineCellActions(current.map, cell.states, cell.divisions, max.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      patient.state <- SwitchPrePost(cell.states, current.map, mutation.encoding)
    } else if (patient.action == "smoke" | patient.action == "asbestos" | patient.action == "radexpose") {
      temp.cell.mut.rate <- cell.mut.rate * 20
      cell.states <- AddMutations(cell.states, temp.cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- RunBaselineCellActions(current.map, cell.states, cell.divisions, max.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
    } else if (patient.action == "deathinduct") {
      cell.states <- AddMutations(cell.states, cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- RunBaselineCellActions(current.map, cell.states, cell.divisions, max.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      results <- KillCellsSpace(cell.states, current.map, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
    } else if (patient.action == "immune" | patient.action == "growthinhibit" | patient.action == "targeted") {
      cell.states <- AddMutations(cell.states, cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- RunBaselineCellActions(current.map, cell.states, cell.divisions, max.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      results <- KillCellsGrowth(cell.states, current.map, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
    } else if (patient.action == "chemo1") {
      cell.states <- AddMutations(cell.states, cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- MarkDeadCells(cell.states, current.map, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      results <- CheckEnergyNeed(current.map, cell.states, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      results <- KillCellsCycling(cell.states, current.map, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
    } else if (patient.action == "chemo2" | patient.action == "radtherapy") {
      temp.cell.mut.rate <- cell.mut.rate * 10
      cell.states <- AddMutations(cell.states, temp.cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- MarkDeadCells(cell.states, current.map, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      results <- CheckEnergyNeed(current.map, cell.states, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      results <- KillCellsDNA(cell.states, current.map, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
    } else if (patient.action == "angblock") {
      cell.states <- AddMutations(cell.states, cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- RunBaselineCellActions(current.map, cell.states, cell.divisions, max.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      results <- KillCellsEnergy(cell.states, current.map, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
    } else if (patient.action == "surgery") {
      cell.states <- AddMutations(cell.states, cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- RunBaselineCellActions(current.map, cell.states, cell.divisions, max.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      results <- KillCellsSurgery(cell.states, current.map, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      num.surgery <- num.surgery + 1
    }
    cell.meet.cond <- CheckPatientDeath(cell.states, mutation.encoding)
    if (length(cell.meet.cond) > round(length(cell.states) * 0.1) && !is.null(cell.meet.cond)) {
      cat("round", i, "\n")
      print("patient has died")
      # has mutations: E, T, D, G, C, S
      # EMT, telomerase, death inhibition,
      # growth activation, cell cycle, metastasis
      print("aggressive tumor has metastasized")
      print(cell.states[cell.meet.cond])
      break
    }
    if (sum(current.map$data$value == "Cell") > round(nrow(current.map$data) * 0.80)) {
      tumor.burden.count <- tumor.burden.count + 1
      if (tumor.burden.count >= 3) {
        cat("round", i, "\n")
        print("patient has died")
        print("too much tumor burden")
        break
      }
    }
    round.cell.num <- c(round.cell.num, GetCellNumbers(current.map))
    round.average.mut.num <- c(round.average.mut.num,
                               GetAverageMutationNum(cell.states, GetCellNumbers(current.map)))
  }
  return(list(current.map, cell.states, cell.divisions, cell.mut.rate,
              round.cell.num, round.average.mut.num))
}


RunMultiSimulation <- function(num.rounds, starting.map, mutation.set,
                               patient.pre.action.set, patient.post.action.set,
                               mutation.encoding, cell.mut.rate,
                               cell.divisions, cell.states, max.divisions) {
  
  # TO BE FILLED IN
}


RunSimulation <- function(num.rounds, starting.map, mutation.set,
                          patient.pre.action.set, patient.post.action.set,
                          mutation.encoding, cell.mut.rate,
                          cell.divisions, cell.states, max.divisions) {
  num.surgery <- 0
  patient.state <- "pre"
  current.map <- starting.map
  tumor.burden.count <- 0
  for (i in 1:num.rounds) {
    patient.action <- DeterminePatientAction(current.map, patient.state, num.surgery,
                                             patient.pre.action.set, patient.post.action.set)
    if (patient.action == "noaction") {
      cell.states <- AddMutations(cell.states, cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- RunBaselineCellActions(current.map, cell.states, cell.divisions, max.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
    } else if (patient.action == "doctor") {
      cell.states <- AddMutations(cell.states, cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- RunBaselineCellActions(current.map, cell.states, cell.divisions, max.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      patient.state <- SwitchPrePost(cell.states, current.map, mutation.encoding)
    } else if (patient.action == "smoke" | patient.action == "asbestos" | patient.action == "radexpose") {
      temp.cell.mut.rate <- cell.mut.rate * 20
      cell.states <- AddMutations(cell.states, temp.cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- RunBaselineCellActions(current.map, cell.states, cell.divisions, max.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
    } else if (patient.action == "deathinduct") {
      cell.states <- AddMutations(cell.states, cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- RunBaselineCellActions(current.map, cell.states, cell.divisions, max.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      results <- KillCellsSpace(cell.states, current.map, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
    } else if (patient.action == "immune" | patient.action == "growthinhibit" | patient.action == "targeted") {
      cell.states <- AddMutations(cell.states, cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- RunBaselineCellActions(current.map, cell.states, cell.divisions, max.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      results <- KillCellsGrowth(cell.states, current.map, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
    } else if (patient.action == "chemo1") {
      cell.states <- AddMutations(cell.states, cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- MarkDeadCells(cell.states, current.map, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      results <- CheckEnergyNeed(current.map, cell.states, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      results <- KillCellsCycling(cell.states, current.map, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
    } else if (patient.action == "chemo2" | patient.action == "radtherapy") {
      temp.cell.mut.rate <- cell.mut.rate * 10
      cell.states <- AddMutations(cell.states, temp.cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- MarkDeadCells(cell.states, current.map, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      results <- CheckEnergyNeed(current.map, cell.states, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      results <- KillCellsDNA(cell.states, current.map, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
    } else if (patient.action == "angblock") {
      cell.states <- AddMutations(cell.states, cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- RunBaselineCellActions(current.map, cell.states, cell.divisions, max.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      results <- KillCellsEnergy(cell.states, current.map, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
    } else if (patient.action == "surgery") {
      cell.states <- AddMutations(cell.states, cell.mut.rate,
                                  map.dim, mutation.encoding)
      results <- RunBaselineCellActions(current.map, cell.states, cell.divisions, max.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      results <- KillCellsSurgery(cell.states, current.map, cell.divisions, cell.mut.rate)
      current.map <- results[[1]]
      cell.states <- results[[2]]
      cell.divisions <- results[[3]]
      cell.mut.rate <- results[[4]]
      num.surgery <- num.surgery + 1
    }
    cell.meet.cond <- CheckPatientDeath(cell.states, mutation.encoding)
    if (length(cell.meet.cond) > round(length(cell.states) * 0.1) && !is.null(cell.meet.cond)) {
      cat("round", i, "\n")
      print("patient has died")
      # has mutations: E, T, D, G, C, S
      # EMT, telomerase, death inhibition,
      # growth activation, cell cycle, metastasis
      print("aggressive tumor has metastasized")
      print(cell.states[cell.meet.cond])
      break
    }
    if (sum(current.map$data$value == "Cell") > round(nrow(current.map$data) * 0.80)) {
      tumor.burden.count <- tumor.burden.count + 1
      if (tumor.burden.count >= 3) {
        cat("round", i, "\n")
        print("patient has died")
        print("too much tumor burden")
        break
      }
    }
  }
  return(list(current.map, cell.states, cell.divisions, cell.mut.rate))
}
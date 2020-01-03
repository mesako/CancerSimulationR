
GetCellNumbers <- function(current.map) {
  total.cells <- sum(current.map$data$value == "Cell")
  return(total.cells)
}

GetCellMutNum <- function(cell.states) {
  if (sum(!is.na(cell.states)) < 1) {
    num.mut.per.cell <- 0
  } else {
    num.mut.per.cell <- sapply(cell.states[!is.na(cell.states)], nchar)
  }
  return(num.mut.per.cell)
}

GetAverageMutationNum <- function(cell.states, total.cells) {
  num.mut <- sum(GetCellMutNum(cell.states))
  average.mut.num <- num.mut / total.cells
  average.mut.num <- round(average.mut.num)
  return(average.mut.num)
}

PlotAverageMutNum <- function(all.average.mut.num) {
  this.data <- cbind(1:length(all.average.mut.num), all.average.mut.num)
  this.data <- as.data.frame(this.data)
  colnames(this.data) <- c("Round", "AvgMutNum")
  mut.plot <- ggplot(data = this.data, mapping = aes(x = Round, y = AvgMutNum))
  mut.plot <- mut.plot + geom_point() + geom_line()
  if (max(this.data$Round) > 20) {
    mut.plot <- mut.plot + theme(legend.position = "none")
  }
  return(mut.plot)
}

PlotTotalCellNum <- function(all.total.cell.num) {
  this.data <- cbind(1:length(all.total.cell.num), all.total.cell.num)
  this.data <- as.data.frame(this.data)
  colnames(this.data) <- c("Round", "CellNum")
  cell.plot <- ggplot(data = this.data, mapping = aes(x = Round, y = CellNum))
  cell.plot <- cell.plot + geom_point() + geom_line()
  if (max(this.data$Round) > 20) {
    cell.plot <- cell.plot + theme(legend.position = "none")
  }
  return(cell.plot)
}

PlotMultiPatientSummary <- function(patient.summary.statements) {
  num.alive <- sum(grepl("alive", patient.summary.statements))
  num.dead.metastatic <- sum(grepl("metastastatic", patient.summary.statements))
  num.dead.burden <- sum(grepl("tumor burden", patient.summary.statements))
  patient.data <- data.frame("Alive" = num.alive,
                             "Metastatic" = num.dead.metastatic,
                             "TumorBurden" = num.dead.burden)
  patient.data <- t(patient.data)
  patient.data <- cbind(rownames(patient.data), patient.data[, 1])
  patient.data <- as.data.frame(patient.data)
  patient.data[, 2] <- as.numeric(as.character(patient.data[, 2]))
  rownames(patient.data) <- NULL
  colnames(patient.data) <- c("Status", "NumPatients")
  patient.plot <- ggplot(patient.data, aes(x = Status, y = NumPatients, fill = Status)) +
    geom_bar(stat = "identity") + theme(legend.position = "none")
  return(patient.plot)
}

PlotMultiSurvivalCurve <- function(patient.summary.statements,
                                   num.rounds) {
  total.num <- length(patient.summary.statements)
  when.death <- grep("round", patient.summary.statements)
  if (length(when.death) < 1) {
    round.num <- sprintf("%02d", 1:num.rounds)
    round.data <- round.num
    round.data <- as.data.frame(round.data)
    round.data <- cbind(round.data, rep(total.num, times = nrow(round.data)))
    colnames(round.data) <- c("RoundNum", "NumAlive")
    patient.plot <- ggplot(round.data, aes(x = RoundNum, y = NumAlive))
    patient.plot <- patient.plot + geom_point() + geom_line(group = 1)
    patient.plot <- patient.plot + ylim(0, total.num)
  } else {
    when.death <- patient.summary.statements[when.death]
    when.death <- regmatches(when.death, gregexpr("[[:digit:]]+", when.death))
    when.death <- as.numeric(unlist(when.death))
    when.death <- table(when.death)
    num.dead <- cumsum(when.death)
    num.alive <- total.num - num.dead
    round.num <- sprintf("%02d", 1:num.rounds)
    round.data <- round.num
    round.data <- as.data.frame(round.data)
    round.data <- cbind(round.data, rep(total.num, times = nrow(round.data)))
    colnames(round.data) <- c("RoundNum", "NumAlive")
    change.num.ind <- as.numeric(names(num.alive))
    change.num.ind <- c(change.num.ind, num.rounds + 1)
    for (x in 1:(length(change.num.ind) - 1)) {
      new.num <- unname(num.alive[as.character(change.num.ind[x])])
      round.data[change.num.ind[x]:(change.num.ind[x + 1] - 1), 2] <- new.num
    }
    patient.plot <- ggplot(round.data, aes(x = RoundNum, y = NumAlive))
    patient.plot <- patient.plot + geom_point() + geom_line(group = 1)
    patient.plot <- patient.plot + ylim(0, total.num)
    if (num.rounds > 20) {
      patient.plot <- patient.plot + theme(axis.text.x = element_blank())
    }
  }
  return(patient.plot)
}

MultiPlotAverageMutNum <- function(all.average.mut.nums) {
  patient.nums <- 1:nrow(all.average.mut.nums)
  patient.nums <- sprintf("%02d", patient.nums)
  this.data <- cbind(paste0("Patient", patient.nums), all.average.mut.nums)
  this.data <- as.data.frame(this.data)
  colnames(this.data)[1] <- c("Patient")
  round.num <- sprintf("%02d", 1:ncol(all.average.mut.nums))
  colnames(this.data)[2:ncol(this.data)] <- round.num
  this.data.m <- melt(this.data, id.vars = "Patient")
  this.data.m$value <- as.numeric(as.character(this.data.m$value))
  mut.plot <- ggplot(data = this.data.m, mapping = aes(x = variable, y = value,
                                                       group = Patient, color = Patient))
  mut.plot <- mut.plot + geom_point() + geom_line()
  mut.plot <- mut.plot + xlab("Round Number") + ylab("Average Number of Mutations")
  if (ncol(all.average.mut.nums) > 20) {
    mut.plot <- mut.plot + theme(axis.text.x = element_blank())
  }
  return(mut.plot)
}

MultiPlotTotalCellNum <- function(all.total.cell.nums) {
  patient.nums <- 1:nrow(all.total.cell.nums)
  patient.nums <- sprintf("%02d", patient.nums)
  this.data <- cbind(paste0("Patient", patient.nums), all.total.cell.nums)
  this.data <- as.data.frame(this.data)
  colnames(this.data)[1] <- c("Patient")
  round.num <- sprintf("%02d", 1:ncol(all.total.cell.nums))
  colnames(this.data)[2:ncol(this.data)] <- round.num
  this.data.m <- melt(this.data, id.vars = "Patient")
  this.data.m$value <- as.numeric(as.character(this.data.m$value))
  cell.plot <- ggplot(data = this.data.m, mapping = aes(x = variable, y = value,
                                                       group = Patient, color = Patient))
  cell.plot <- cell.plot + geom_point() + geom_line()
  cell.plot <- cell.plot + xlab("Round Number") + ylab("Total Number of Cells")
  if (ncol(all.total.cell.nums) > 20) {
    cell.plot <- cell.plot + theme(axis.text.x = element_blank())
  }
  return(cell.plot)
}

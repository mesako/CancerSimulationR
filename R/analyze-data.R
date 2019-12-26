
GetCellNumbers <- function(current.map) {
  total.cells <- sum(current.map$data$value == "Cell")
  return(total.cells)
}

GetCellMutNum <- function(cell.states) {
  # print(cell.states)
  # print(sum(!is.na(cell.states)))
  if (sum(!is.na(cell.states)) < 1) {
    num.mut.per.cell <- 0
  } else {
    num.mut.per.cell <- sapply(cell.states[!is.na(cell.states)], nchar)
  }
  # print(cell.states[!is.na(cell.states)])
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
  return(mut.plot)
}

PlotTotalCellNum <- function(all.total.cell.num) {
  this.data <- cbind(1:length(all.total.cell.num), all.total.cell.num)
  this.data <- as.data.frame(this.data)
  colnames(this.data) <- c("Round", "CellNum")
  cell.plot <- ggplot(data = this.data, mapping = aes(x = Round, y = CellNum))
  cell.plot <- cell.plot + geom_point() + geom_line()
  return(cell.plot)
}


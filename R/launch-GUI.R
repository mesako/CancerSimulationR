#' GUI for launching shiny Simulation Application
#' @examples
#' \dontrun{LaunchGUI()}
#' @export
LaunchGUI <- function() {
  # requireNamespace("shiny")
  # requireNamespace("???")
  # requireNamespace("???")
  runApp(appDir = file.path(system.file(package = "CancerSimulationR"), "app"))
}

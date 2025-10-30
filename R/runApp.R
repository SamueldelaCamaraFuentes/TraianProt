#' @title Launch the TraianProt Shiny Application
#' @description Launches the Shiny application for the TraianProt package.
#'
#' @details This function finds the Shiny app directory within the installed
#' package and runs it.
#'
#' @return Nothing. This function is called for its side effect of
#'         launching the Shiny application.
#'
#' @export
#' @importFrom shiny runApp
#'
#' @examples
#' \dontrun{
#'   # This is the correct way to show an app example
#'   runTraianProt()
#' }
runTraianProt <- function() {
  appDir <- system.file(package = "TraianProt", "app.R")
  if (appDir == "") {
    stop(
      "Could not find the app.R file. ",
      "Try re-installing 'TraianProt'."
    )
  }
  shiny::runApp(appDir, display.mode = "normal")
}

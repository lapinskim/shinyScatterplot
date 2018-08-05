#' Run Shiny scatterplot application
#'
#' Launch Shiny application that displays the interactive scatterplot.
#'
#' Recreation of Figure 2. from the publication \emph{Cytoplasmic
#' polyadenylation-mediated translational control of maternal mRNAs directs
#' maternal to zygotic transition}, Winata CL, Development, 11 December 2017
#' \url{http://dev.biologists.org/content/early/2017/12/08/dev.159566}
#'
#' Measurements of transcription (∆TS) and translation (∆TL) rates were obtained
#' as a fraction of total and polysome-bound expression values at particular
#' stages compared with the egg stage as baseline. Translational regulation is
#' defined as a non-linear relationship  between  transcription and translation
#' rates.
#'
#' \deqn{transcription rate (∆TS)=\frac{total(stage x)}{total(egg)}}{transcription rate (∆TS) = total(stage x) / total(egg)}
#' \deqn{translation rate(∆TL)=\frac{polysome(stage x)}{polysome(egg)}}{translation rate(∆TL) = polysome(stage x) / polysome(egg)}
#'
#' @examples
#' ## Run in interactive session only
#' if (interactive()) {
#'     runPlot()
#' }
#'
#' @seealso Winata CL, Development, 11 December 2017
#' \url{http://dev.biologists.org/content/early/2017/12/08/dev.159566}
#'
#' @export
runPlot <- function() {
    appDir <- system.file("application", package = "shinyScatterplot")

    shiny::runApp(appDir, display.mode = "normal")
}

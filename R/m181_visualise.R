#' Plots an xts-object and highlight preassigned FALSE data using static
#' ggplot2 style.
#'
#' @title Plot xts object and highlight FALSE data
#' @param x An \code{xts} object to be plotted.
#' @param y An \code{xts} object with TRUE and FALSE values. 
#' @param plot_only_false_flags logical. If TRUE, only FALSE flags are plotted.
#' @param title The title of the plot.
#' @rdname visualise_gg
#' @export
#' @seealso \code{\link[xts]{xts}}, \code{\link[ggplot2]{ggplot}}.
visualise_gg <- function(x, y, title="flag", plot_only_false_flags=TRUE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 needed for this function to work. Please install it.",
         call. = FALSE)
  }
    
  # checks if x and y are really xts objects and stops if not
  stopifnot(xts::is.xts(x),xts::is.xts(y))
  
  # print only false flags?
  if (plot_only_false_flags) y <- y[!y]
  
  # merge both xts to prepare conversion to data frame
  xts <- merge(x,y)  
  # conversion to data frame 
  df_of_xts <- data.frame(Index = zoo::index(xts), zoo::coredata(xts))
  # rename col
  colnames(df_of_xts) <- c("Index","Value","Flag")
  # ggplot part
  p <- ggplot2::ggplot(data = df_of_xts,
                       mapping = ggplot2::aes(x = Index,
                                              y = Value)) 
  p <- p + ggplot2::geom_line()
  p <- p + ggplot2::geom_point(data = df_of_xts,
                               mapping = ggplot2::aes(x = Index,
                                                      y = Value, 
                                                      colour = ifelse(Flag == 0,
                                                                      'red',
                                                                      'black'),
                                                      size = ifelse(Flag == 0,
                                                                    0.2, 0.1))) 
  p <- p + ggplot2::scale_color_identity()
  p <- p + ggplot2::theme(legend.position = "none")
  p <- p + ggplot2::ggtitle(title)
  #return ggplot object
  return(p)
}

#' Plots an xts-object and highlight preassigned FALSE data using 
#' interactive dygraphs style.
#'
#' @title Plot xts object and highlight FALSE data
#' @param x An \code{xts} object to be plotted.
#' @param y An \code{xts} object with TRUE and FALSE values. 
#' @rdname visualise_dy
#' @importFrom magrittr "%>%"
#' @export
#' @seealso \code{\link[xts]{xts}}, \code{\link[dygraphs]{dygraph}}.
visualise_dy <- function(x, y) {

  if (!requireNamespace("dygraphs", quietly = TRUE)) {
    stop("dygraphs needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  # select only indices of x, where checks are flagged FALSE
  y <- x[zoo::index(y[!y])]
  
  # merge objects to have one xts object which can be handled by dygraph
  data <- merge(x,y)
  
  graph <- dygraphs::dygraph(data) %>%
    dygraphs::dyOptions(useDataTimezone = TRUE) %>% 
    dygraphs::dySeries(colnames(data[,1]),
                       color = "black") %>% 
    dygraphs::dySeries(colnames(data[,2]),
                       label = "False Flags",
                       color = "red",
                       drawPoints = TRUE,
                       strokeWidth = 0,
                       pointSize = 5) %>% 
    dygraphs::dyRangeSelector()
  
  #return dygraph object
  return(graph) 
  
}

#' Plots an xts-object and highlight preassigned FALSE data
#'
#' @title Plot xts object and highlight FALSE data
#' @param x An \code{xts} object to be plotted.
#' @param y An \code{xts} object with TRUE and FALSE values. 
#' @param a An \code{xts} object to be plotted.
#' @param b An \code{xts} object to be plotted.
#' @param c An \code{xts} object to be plotted.
#' @rdname visualise_dy3
#' @importFrom magrittr "%>%"
#' @seealso \code{\link[xts]{xts}}, \code{\link[dygraphs]{dygraph}}.
visualise_dy3 <- function(x, y, a=NA, b=NA, c=NA) {
  
  message("This method has still prototype status!")
  
  if (!requireNamespace("dygraphs", quietly = TRUE)) {
    stop("dygraphs needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  # select only indices of x, where checks are flagged FALSE
  y <- x[zoo::index(y[!y])]
  
  # merge objects to have one xts object which can be handled by dygraph
  data <- merge(x, y, a, b, c)
  
  graph <- dygraphs::dygraph(data) %>% 
    dygraphs::dySeries(colnames(data[,1]),
                       color = "black") %>% 
    dygraphs::dySeries(colnames(data[,2]),
                       label = "False Flags",
                       color = "red",
                       drawPoints = TRUE,
                       strokeWidth = 0,
                       pointSize = 5) %>%
    dygraphs::dySeries(colnames(data[,3]),
                       label = "a",
                       color = "blue",
                       drawPoints = FALSE,
                       strokeWidth = 1) %>%
    dygraphs::dySeries(colnames(data[,4]),
                       label = "b",
                       color = "green",
                       drawPoints = FALSE,
                       strokeWidth = 1) %>%
    dygraphs::dySeries(colnames(data[,5]),
                       label = "c",
                       color = "green",
                       drawPoints = FALSE,
                       strokeWidth = 1) %>%
    dygraphs::dyRangeSelector(.)
  
  #return dygraph object
  return(graph) 
  
}
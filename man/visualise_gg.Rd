% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m181_visualise.R
\name{visualise_gg}
\alias{visualise_gg}
\title{Plot xts object and highlight FALSE data}
\usage{
visualise_gg(x, y, title = "flag", plot_only_false_flags = TRUE)
}
\arguments{
\item{x}{An \code{xts} object to be plotted.}

\item{y}{An \code{xts} object with TRUE and FALSE values.}

\item{title}{The title of the plot.}

\item{plot_only_false_flags}{logical. If TRUE, only FALSE flags are plotted.}
}
\description{
Plots an xts-object and highlight preassigned FALSE data using static
ggplot2 style.
}
\seealso{
\code{\link[xts]{xts}}, \code{\link[ggplot2]{ggplot}}.
}


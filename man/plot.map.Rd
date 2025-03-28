% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot.map}
\alias{plot.map}
\title{Nice scatter plots for georeferenced data}
\usage{
plot.map(
  x,
  y = NULL,
  xlim = NULL,
  ylim = NULL,
  domain = NULL,
  center = 0,
  hires = FALSE,
  land.col = "darkolivegreen4",
  sea.col = "aliceblue",
  boundaries.col = "black",
  grid.col = "white",
  grid = TRUE,
  grid.lwd = 0.5,
  cex = 0.5,
  pch = 19,
  main = NULL,
  add = FALSE,
  axes = TRUE,
  land = TRUE,
  border = !axes,
  asp = NA,
  axs = "i",
  xaxs = axs,
  yaxs = axs,
  cex.axis = 0.75,
  interior = FALSE,
  fill = TRUE,
  countries = FALSE,
  nx = NULL,
  ny = nx,
  ...
)
}
\arguments{
\item{x}{vector of longitude values, or a data.frame containing variables named 'lon'
and 'lat'}

\item{y}{vector of latitude values.}

\item{xlim}{}

\item{ylim}{}

\item{center}{}

\item{hires}{}

\item{land.col}{}

\item{sea.col}{}

\item{boundaries.col}{}

\item{grid.col}{}

\item{grid}{}

\item{cex}{}

\item{pch}{}

\item{main}{}

\item{add}{}

\item{axes}{}

\item{land}{}

\item{border}{}

\item{asp}{}

\item{axs}{}

\item{xaxs}{}

\item{yaxs}{}

\item{cex.axis}{}

\item{interior}{}

\item{fill}{}

\item{countries}{}

\item{nx}{}

\item{ny}{}

\item{...}{}
}
\description{
This function makes and scatter plot given the longitude and latitude
of the data points.
}
\examples{
load(mackerel)
plot.map(mackerel)
plot.map(domain="peru")
}
\author{
Ricardo Oliveros-Ramos
}

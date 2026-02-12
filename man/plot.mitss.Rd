\name{plot.mitss}
\alias{plot.mitss}
\title{Plot method for mitss objects}
\usage{
\method{plot}{mitss}(x, type = c("estimates", "propensity", "balance"), ...)
}
\arguments{
\item{x}{A mitss object}

\item{type}{Type of plot: "estimates" for distribution of MI estimates,
  "propensity" for propensity score distribution, or "balance" for covariate balance}

\item{...}{Additional arguments passed to plot functions}
}
\description{
Creates diagnostic plots for MITSS results.
}

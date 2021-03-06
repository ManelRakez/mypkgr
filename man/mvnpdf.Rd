% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvnpdfRAW.R
\name{mvnpdf}
\alias{mvnpdf}
\title{Density for multivariate normal distribution}
\usage{
mvnpdf(x, mean = rep(0, nrow(x)), varcovM = diag(nrow(x)), Log = TRUE)
}
\arguments{
\item{x}{a matrix with n columns (the observations) and p rows}

\item{mean}{a vetor of means}

\item{varcovM}{a variance-covariance matrix}

\item{Log}{a logical parameter, with default value to TRUE}
}
\value{
a list containing the matrix x, and a vector of length n of the multivariate normal distribution density
}
\description{
Description
}

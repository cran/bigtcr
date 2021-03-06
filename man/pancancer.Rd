% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bigtcr_package.R
\name{pancancer}
\alias{pancancer}
\title{Example Pancreatic Cancer Dataset}
\format{A dataframe with 3 variables:
\describe{
  \item{obs.y}{Observed time to failure events or censoring in months}
  \item{min.type}{Type of failure events
\describe{
\item{0}{Censored}
\item{1}{death with metastasis limited to lung only}
\item{2}{death with metastasis that involves sites other than lung (e.g. liver)}
\item{3}{death without disease recurrence}
}
}
\item{v}{Time to recurrence. NA if no recurrence observed}
}}
\description{
Simulated data used in \pkg{bigtcr} examples.
}
\details{
Data simulated based on the patients who had surgical resection of pancreatic
adenocarcinomas and had postoperative follow-up at the Johns Hopkins Hospital
between 1998 and 2007.
}


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CVEP_IPW.R
\name{CVEP_IPW}
\alias{CVEP_IPW}
\title{Fit a linear mixed effects model with inverse probability weighting}
\usage{
CVEP_IPW(data, control = list(wyear = 3, wcycle = 3), yes.unwgt = F)
}
\arguments{
\item{data}{a data frame with the following variable names: Variety, Year, Loc, Rep, Cycle, Checks, TT_Trait,t DM_Trait}

\item{control}{a list of control parameters including \code{wyear} which specifies the window size in years of the missing entries for regular varieties and \code{wcycle} which specifies the window size in cycles of the missing entries for control varieties}

\item{yes.unwgt}{a logical variable indicating whether to fit an unweighted linear mixed effects model}
}
\value{
a list of data used to fit the linear mixed effects model, the fitted model and summarized results
}
\description{
Fit a linear mixed effects model with inverse probability weighting
}

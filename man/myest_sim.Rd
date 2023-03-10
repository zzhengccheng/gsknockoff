% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/myest_sim.R
\name{myest_sim}
\alias{myest_sim}
\title{This function implement our proposed Knockoff method for testing union null}
\usage{
myest_sim(
  Xlist,
  Ylist,
  familylist,
  Kmethod,
  Zmethod,
  Wmethod,
  q,
  offset,
  Zfun = NULL,
  Wfun = NULL
)
}
\arguments{
\item{Xlist}{list of length K for predictor matrices from each experiment}

\item{Ylist}{list of length K for outcomes from each experiment}

\item{familylist}{vector of length K for type of outcomes from each experiment (Gaussian, Binomial or Cox)}

\item{Kmethod}{vecttor of length K for knockoff construction methods for each experiment (fixed or second order)}

\item{Zmethod}{vecttor of length K for Z statistics calculation methods for each experiment. Current implemented method is glmnet}

\item{Wmethod}{Current implemented method for general K is Diff, method="Other" can be used for other user defined function}

\item{q}{FDR level for control}

\item{offset}{0 for Knockoff and 1 for Knockoff+}

\item{Zfun}{User defined function on how to construct Z statistics from a list of length K for Z and a list of length K for Ztilde (required if select Zmethod="Other")}

\item{Wfun}{User defined function on how to construct W statistics from a list of length K for Z and a list of length K for Ztilde (required if select Wmethod="Other")}
}
\value{
A list of selected variables
}
\description{
This function implement our proposed Knockoff method for testing union null
}

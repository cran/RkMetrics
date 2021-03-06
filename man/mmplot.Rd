% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RkMetrics_functions.R
\name{mmplot}
\alias{mmplot}
\title{A Plotting Function}
\usage{
mmplot(n, x, young, old)
}
\arguments{
\item{n}{the length of the vector Defaults to TRUE.}

\item{x}{the vector arguement.}

\item{young}{the age at which the accident hump begins. Must be entered}

\item{old}{age at which, either mortality experience between males and females converge, or rapid acceleration of mortality. This is typically over 80 years.}
}
\description{
Produces a plot of the difference between the area-under-the-curve for the mortality data and the extended mortality boundaries
}
\examples{
#Examples
m1 <- Mortality$D.Male[which(Mortality$Year == 2008)]
m2 <- Mortality$E.Male[which(Mortality$Year == 2008)]
male.1 <- m1/m2
male.2 <- log(male.1[!is.na(male.1)])
lplot(1:length(male.2),male.2)


mmplot(1:length(male.2),male.2,young=17,old=80)

}
\keyword{mmplot}

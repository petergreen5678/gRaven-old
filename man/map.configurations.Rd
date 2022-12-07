\name{map.configurations}
\alias{map.configurations}
\title{Get belief in a gRaven domain}
\description{Get belief (marginal probability) in a gRaven domain}
\usage{map.configurations(dom, nodes, pmin)}
\arguments{
  \item{dom}{name of gRaven domain}
  \item{nodes}{character vector of names of nodes}
  \item{pmin}{a single numeric value between 0 and 1 specifying the minimum probability for the most probable configurations}
}
\details{Emulates function of the same name in the \code{RHugin} package by calls to \code{gRain} functions}
\value{a data.frame with one column for each node in \code{nodes}. Each row contains a most probable configuration. The final column of the data.frame (\code{Prob}) gives the probability of the configuration.}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}
\examples{
demo('chest',package="gRaven",echo=FALSE)
compile(chest)
set.finding(chest,'asia','yes')
set.finding(chest,'dysp','no')
map.configurations(chest,c('lung','bronc'),0.0001)
}

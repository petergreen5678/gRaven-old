\name{get.normalization.constant}
\alias{get.normalization.constant}
\title{Get the normalisation constant of a gRaven domain}
\description{Get the normalisation constant of a gRaven domain}
\usage{
get.normalization.constant(dom, log = FALSE)
}
\arguments{
  \item{dom}{
name of a gRaven domain
}
  \item{log}{
logical, should constant be returned on log scale?
}
}
\value{numeric, the value of the normalisation constant (or its logarithm)}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}
\examples{
yn <- c("yes","no")
chest<-hugin.domain()

add.node(chest,"asia",yn)
add.node(chest,"tub",yn)
add.node(chest,"smoke",yn)
add.node(chest,"lung",yn)
add.node(chest,"bronc",yn)
add.node(chest,"either",yn)
add.node(chest,"xray",yn)
add.node(chest,"dysp",yn)

add.edge(chest,"tub","asia")
add.edge(chest,"lung","smoke")
add.edge(chest,"bronc","smoke")
add.edge(chest,"either","lung")
add.edge(chest,"either","tub")
add.edge(chest,"xray","either")
add.edge(chest,"dysp","bronc")
add.edge(chest,"dysp","either")

set.table(chest,"asia",c(0.01,0.99))
set.table(chest,"tub",c(0.05,0.95,0.01,0.99))
set.table(chest,"smoke",c(0.5,0.5))
set.table(chest,"lung",c(0.1,0.9,0.01,0.99))
set.table(chest,"bronc",c(0.6,0.4,0.3,0.7))
set.table(chest,"either",c(1,0,1,0,1,0,0,1))
set.table(chest,"xray",c(0.98,0.02,0.05,0.95))
set.table(chest,"dysp",c(0.9,0.1,0.7,0.3,0.8,0.2,0.1,0.9))

compile(chest)
get.normalization.constant(chest)

set.finding(chest,"asia","yes")
set.finding(chest,"dysp","yes")
get.normalization.constant(chest)
}

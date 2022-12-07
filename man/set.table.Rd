\name{set.table}
\alias{set.table}
\title{Set CPT in a gRaven domain}
\description{Set CPT in a gRaven domain}
\usage{
set.table(dom, n, tab = 1, type = "cpt")
}
\arguments{
  \item{dom}{name of gRaven domain}
  \item{n}{name of node}
  \item{tab}{values of conditional probabilties}  
  \item{type}{only "cpt" is accepted}
}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}

\value{a NULL value is invisibly returned.}

\examples{
demo(chest,package="gRaven")
chest
}

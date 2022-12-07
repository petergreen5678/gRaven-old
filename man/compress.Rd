\name{compress}
\alias{compress}
\title{Dummy network compression routine}
\description{Dummy routine.}
\usage{compress(dom)}
\arguments{
  \item{dom}{name of gRaven domain}
}
\value{1}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}
\examples{
demo('chest',package="gRaven",echo=FALSE)
compile(chest)
compress(chest)
}

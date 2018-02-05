\name{get.smoothed.enrichment.mle}
\alias{get.smoothed.enrichment.mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate chromosome-wide profiles of smoothed enrichment estimate }
\description{
  Given signal and control tag positions, the method calculates log2
  signal to control enrichment esimates (maximum likelihood) for each
  chromosome, based on the smoothed tag density profile (see \link{get.smoothed.tag.density}).
}
\usage{

get.smoothed.enrichment.mle(signal.tags, 
  control.tags, 
  tag.shift = 146/2, 
  background.density.scaling = F, 
  pseudocount = 1,
  bg.weight = NULL, 
  rngl = NULL, 
  chrl = NULL, 
  ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{signal.tags}{ signal chromosome tag coordinate vectors (e.g. output
    of \code{\link{select.informative.tags}} }
  \item{control.tags}{ control (input) tags }
  \item{tag.shift}{ number of base pairs by which positive and negative
    tag coordinates should be shifted towards eachother (half of binding
    peak separation distance)}
  \item{background.density.scaling}{background.density.scaling}
  \item{pseudocount}{ pseudocount value to be added to tag density -
    defaults to 1 }
  \item{bg.weight}{ optional weight by which the background density
    should be multipled for scaling. If not supplied, the weight is
    calculated based on the ratio of the reduced ChIP to input dataset sizes. }
  \item{rngl}{rngl}
  \item{chrl}{chrl}
  \item{\dots}{ additional parameters should be the same as those passed
      to the \code{\link{get.smoothed.tag.density}}, such as for example
      bandwidth (default value 150) and step (default value 50).
      see appropriate reference for
      \code{\link{get.smoothed.tag.density}} for details.}
}
\value{
  A list of elements corresponding to chromosomes, with each element
  being an $x/$y data.frame giving the position and associated
  log2 signal/control enrichment estimate. 
}
\seealso{ \code{\link{writewig}} }
\examples{
\dontrun{
  # get smoothed enrichment estimate profile using 500bp bandwidth at
  # 50bp steps
  smoothed.M <- get.smoothed.enrichment.mle(chip.data,bandwidth=500,step=50);
  writewig(smoothed.M,"example.smoothedM.wig","Example smoothed log2 intensity ratio estimate");
  }
}

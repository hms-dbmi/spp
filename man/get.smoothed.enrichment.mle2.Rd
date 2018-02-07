\name{get.smoothed.enrichment.mle2}
\alias{get.smoothed.enrichment.mle2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate background input controlled chromosome-wide
       profiles of smoothed enrichment estimate }
\description{
  Given signal and control tag positions, the method calculates log2
  signal to control enrichment esimates (maximum likelihood) for each
  chromosome, based on the smoothed tag density profile
  (see \link{get.smoothed.tag.density}).
}
\usage{
get.smoothed.enrichment.mle2(signal.tags1, control.tags1,
                             signal.tags2, control.tags2,
                             tag.shift = 146/2, background.density.scaling = F,
                             pseudocount = 1, bg.weight1 = NULL,
                             bg.weight2 = NULL,
                             rngl = NULL, chrl = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{signal.tags1}{
Parameter
}
  \item{control.tags1}{
Parameter
}
  \item{signal.tags2}{
Parameter
}
  \item{control.tags2}{
Parameter
}
  \item{tag.shift}{
Parameter
}
  \item{background.density.scaling}{
Parameter
}
  \item{pseudocount}{
Parameter
}
  \item{bg.weight1}{
Parameter
}
  \item{bg.weight2}{
Parameter
}
  \item{rngl}{
Parameter
}
  \item{chrl}{
Parameter
}
  \item{\dots}{
Parameter
}
}
\value{
  A list of elements corresponding to chromosomes, with each element
  being an $x/$y data.frame giving the position and associated
  log2 signal/control enrichment estimate. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{get.smoothed.enrichment.mle}}
}
\examples{
\dontrun{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (signal.tags1, control.tags1, signal.tags2, control.tags2, 
    tag.shift = 146/2, background.density.scaling = F, pseudocount = 1, 
    bg.weight1 = NULL, bg.weight2 = NULL, rngl = NULL, chrl = NULL, 
    ...) 
{
    if (is.null(chrl)) {
        chrl <- intersect(names(signal.tags1), names(signal.tags2))
        names(chrl) <- chrl
    }
    if (is.null(rngl)) {
        rngl <- lapply(chrl, function(chr) range(c(range(abs(signal.tags1[[chr]] + 
            tag.shift)), range(abs(signal.tags2[[chr]] + tag.shift)))))
    }
    else {
        chrl <- names(rngl)
        names(chrl) <- chrl
    }
    ssd1 <- get.smoothed.tag.density(signal.tags1, rngl = rngl, 
        ..., scale.by.dataset.size = F)
    ssd2 <- get.smoothed.tag.density(signal.tags2, rngl = rngl, 
        ..., scale.by.dataset.size = F)
    csd1 <- get.smoothed.tag.density(control.tags1, rngl = rngl, 
        ..., scale.by.dataset.size = F)
    csd2 <- get.smoothed.tag.density(control.tags2, rngl = rngl, 
        ..., scale.by.dataset.size = F)
    if (is.null(bg.weight1)) {
        bg.weight1 <- dataset.density.ratio(signal.tags1, control.tags1, 
            background.density.scaling = background.density.scaling)
    }
    if (is.null(bg.weight2)) {
        bg.weight2 <- dataset.density.ratio(signal.tags2, control.tags2, 
            background.density.scaling = background.density.scaling)
    }
    cmle <- lapply(chrl, function(chr) {
        d <- ssd1[[chr]]
        d$y <- log2(ssd1[[chr]]$y + pseudocount * bg.weight1) - 
            log2(csd1[[chr]]$y + pseudocount) - log2(bg.weight1) - 
            log2(ssd2[[chr]]$y + pseudocount * bg.weight2) + 
            log2(csd2[[chr]]$y + pseudocount) + log2(bg.weight2)
        return(d)
    })
  }
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{get.conservative.fold.enrichment.profile2}
\alias{get.conservative.fold.enrichment.profile2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Return Conservative fold enrichment profile controlling for input
       and a single background scale  }
\description{
Returns a conservative upper/lower bound profile (log2) given signal tag list,
  background tag list and window scales controlling for input,
  and supporting only a single background scale.
  
  Novel version of get.conservative.fold.enrichment.profile()
  supporting a single background scale.
}
\usage{
get.conservative.fold.enrichment.profile2(ftl1, ftl2, btl1, btl2, fws,
           bws = 1 * fws, step = 50, tag.shift = 146/2, 
           alpha = 0.05, background.density.scaling = T,
           bg.weight1 = NULL, bg.weight2 = NULL, 
           posl = NULL, return.mle = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ftl1}{
Parameter
}
  \item{ftl2}{
Parameter
}
  \item{btl1}{
Parameter
}
  \item{btl2}{
Parameter
}
  \item{fws}{
Parameter
}
  \item{bws}{
Parameter
}
  \item{step}{
Parameter
}
  \item{tag.shift}{
Parameter
}
  \item{alpha}{
Parameter
}
  \item{background.density.scaling}{
Parameter
}
  \item{bg.weight1}{
Parameter
}
  \item{bg.weight2}{
Parameter
}
  \item{posl}{
Parameter
}
  \item{return.mle}{
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
function (ftl1, ftl2, btl1, btl2, fws, bws = 1 * fws, step = 50, 
    tag.shift = 146/2, alpha = 0.05, background.density.scaling = T, 
    bg.weight1 = NULL, bg.weight2 = NULL, posl = NULL, return.mle = F) 
{
    ftl1 <- ftl1[unlist(lapply(ftl1, length)) > 2]
    chrl <- names(ftl1)
    names(chrl) <- chrl
    if (!is.null(posl)) {
        chrl <- chrl[chrl \%in\% names(posl)]
    }
    if (is.null(bg.weight1)) {
        bg.weight1 <- dataset.density.ratio(ftl1, btl1,
                      background.density.scaling = background.density.scaling)
    }
    if (is.null(bg.weight2)) {
        bg.weight2 <- dataset.density.ratio(ftl2, btl2,
                      background.density.scaling = background.density.scaling)
    }
    lapply(chrl, function(chr) {
        x <- binomial.proportion.ratio.bounds(abs(ftl1[[chr]] + 
            tag.shift), abs(btl1[[chr]] + tag.shift), abs(ftl2[[chr]] + 
            tag.shift), abs(btl2[[chr]] + tag.shift), fws = fws, 
            bws = bws, step = step, bg.weight1 = bg.weight1, 
            bg.weight2 = bg.weight2, alpha = alpha, pos = if (is.null(posl)) {
                NULL
            }
            else {
                posl[[chr]]
            })
        ps <- rep(0, length(x$mle))
        vi <- which(!is.na(x$lb) & x$lb > 0)
        ps[vi] <- x$lb[vi]
        vi <- which(!is.na(x$ub) & x$ub < 0)
        ps[vi] <- x$ub[vi]
        if (is.null(posl)) {
            if (return.mle) {
                return(data.frame(x = x$x, y = ps, mle = x$mle, 
                  lb = x$lb, ub = x$ub))
            }
            else {
                return(data.frame(x = x$x, y = ps))
            }
        }
        else {
            if (return.mle) {
                return(data.frame(x = posl[[chr]], y = ps, mle = x$mle, 
                  lb = x$lb, ub = x$ub))
            }
            else {
                return(data.frame(x = posl[[chr]], y = ps))
            }
        }
    })
  }
}
}

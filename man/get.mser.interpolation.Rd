\name{get.mser.interpolation}
\alias{get.mser.interpolation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ will do}
\description{
 A concise (1-5 lines) description of what the function does. 
}
\usage{
get.mser.interpolation(signal.data, control.data, target.fold.enrichment = 5, n.chains = 10, n.steps = 6, step.size = 1e+05, chains = NULL, test.agreement = 0.99, return.chains = F, enrichment.background.scales = c(1), plot = T, excluded.steps = c(seq(2, n.steps - 2)), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{signal.data}{descibe the code}
  \item{control.data}{ descibe the code }
  \item{target.fold.enrichment}{ descibe the code }
  \item{n.chains}{ descibe the code }
  \item{n.steps}{ tesing the code}
  \item{step.size}{tesing the code}
  \item{chains}{tesing the code}
  \item{test.agreement}{ tesing the code}
  \item{return.chains}{ tesing the code}
  \item{enrichment.background.scales}{ tesing the code }
  \item{plot}{tesing the code }
  \item{excluded.steps}{ tesing the code }
  \item{\dots}{ tesing the code}
}
\details{
not required
}
\value{
  
retuns on 0 or 1 for fail and success respectively}

\references{"http://compbio.med.harvard.edu/spp"}
\author{Peterk}
\note{none }

\examples{

## To run the with all the default settings.
x<-someexmaple
}


##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (signal.data, control.data, target.fold.enrichment = 5, 
    n.chains = 10, n.steps = 6, step.size = 1e+05, chains = NULL, 
    test.agreement = 0.99, return.chains = F, enrichment.background.scales = c(1), 
    plot = T, excluded.steps = c(seq(2, n.steps - 2)), ...) 
{
    msers <- get.mser(signal.data, control.data, n.chains = n.chains, 
        n.steps = n.steps, step.size = step.size, chains = chains, 
        test.agrement = test.agreement, return.chains = T, enrichment.background.scales = enrichment.background.scales, 
        excluded.steps = excluded.steps, ...)
    mser <- mser.chain.interpolation(chains = msers$chains, enrichment.background.scales = c(1), 
        test.agreement = test.agreement, return.lists = T)
    sr <- sum(unlist(lapply(signal.data, length)))/mser[[1]][[1]]$n[1]
    intpn <- lapply(mser, function(ms) {
        lmvo <- do.call(rbind, ms)
        lmvo$n <- lmvo$n * sr
        lmvo <- lmvo[lmvo$nd == lmvo$nd[1], ]
        lmvo <- na.omit(lmvo)
        if (any(lmvo$me == 1)) {
            return(list(prediction = NA))
        }
        lmvo$n <- log10(lmvo$n)
        lmvo$me <- log10(lmvo$me - 1)
        emvf <- lm(me ~ n, data = lmvo)
        tfe <- (log10(target.fold.enrichment - 1) - coef(emvf)[[1]])/coef(emvf)[[2]]
        tfen <- 10^tfe
        return(list(prediction = tfen, log10.fit = emvf))
    })
    if (return.chains) {
        return(list(interpolation = intpn, chains = msers$chains))
    }
    else {
        return(intpn)
    }
    return(msers)
    msers$mser
    unlist(lapply(msers$chains[[3]], function(d) sum(unlist(lapply(d$npl, 
        function(x) length(x$x))))))
    mser <- mser.chain.interpolation(chains = msers$chains, enrichment.background.scales = c(1), 
        test.agreement = test.agreement, return.lists = F)
    mser
    ms <- mser[[1]]
    names(mser[[1]])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{array}
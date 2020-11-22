\name{read.arachne.tags}
\alias{read.arachne.tags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read in Arachne tags }
\description{
Read in ARACHNE Tag file
}
\usage{
read.arachne.tags(filename, fix.chromosome.names = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
filename
}
  \item{fix.chromosome.names}{
do we fix chromosome names
}
}

\value{
A list like element
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Fill in later
}
\examples{
\dontrun{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (filename, fix.chromosome.names = F) 
{
    tl <- lapply(.Call("read_arachne_long", path.expand(filename)), 
        function(d) {
            xo <- order(abs(d$t))
            d$t <- d$t[xo]
            d$n <- d$n[xo]
            d$l <- d$l[xo]
            return(d)
        })
    if (fix.chromosome.names) {
        names(tl) <- gsub("\\.fa", "", names(tl))
    }
    return(list(tags = lapply(tl, function(d) d$t), quality = lapply(tl, 
        function(d) d$n), length = lapply(tl, function(d) d$l)))
  }
}
}

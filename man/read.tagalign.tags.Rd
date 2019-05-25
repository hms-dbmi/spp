\name{read.tagalign.tags}
\alias{read.tagalign.tags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read in tagalign tags }
\description{
Fill in later
}
\usage{
read.tagalign.tags(filename, fix.chromosome.names = T, fix.quality = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
Filename of tag file
}
  \item{fix.chromosome.names}{
chromosome names
}
  \item{fix.quality}{
fix quality
}
}
\details{
...
}
\value{
a list like structure
}
\references{
spp by kharchenko
}
\author{
Peter K.
}
\note{
Needs further editing
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Buh!
}
\examples{
\dontrun{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (filename, fix.chromosome.names = T, fix.quality = T) 
{
    tl <- lapply(.Call("read_tagalign", path.expand(filename)), 
        function(d) {
            xo <- order(abs(d$t))
            d$t <- d$t[xo]
            d$n <- d$n[xo]
            if (fix.quality) {
                if (min(d$n) < 0.5) {
                  d$n = ceiling(1000/4^d$n)
                }
                break.vals <- unique(sort(c(0, unique(d$n))))
                d$n <- length(break.vals) - 1 - cut(d$n, breaks = break.vals, 
                  labels = F)
            }
            return(d)
        })
    if (fix.chromosome.names) {
        names(tl) <- gsub("\\.fa", "", names(tl))
    }
    return(list(tags = lapply(tl, function(d) d$t), quality = lapply(tl, 
        function(d) d$n)))
  }
}
}


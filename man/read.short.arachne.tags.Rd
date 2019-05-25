\name{read.short.arachne.tags}
\alias{read.short.arachne.tags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read in ARACHNE short tags }
\description{
Read in short arachne reads
}
\usage{
read.short.arachne.tags(filename, fix.chromosome.names = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
filename 
}
  \item{fix.chromosome.names}{
Fix chromosome names
}
}
\details{
Not necessary
}
\value{
A list like structure
}
\references{
spp by Peter Kharchenko
}
\author{
Peter Kharchenko
}
\note{
No Notes
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Nothing to see here
}
\examples{
\dontrun{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (filename, fix.chromosome.names = F) 
{
    tl <- lapply(.Call("read_arachne", path.expand(filename)), 
        function(d) {
            xo <- order(abs(d$t))
            d$t <- d$t[xo]
            d$n <- d$n[xo]
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

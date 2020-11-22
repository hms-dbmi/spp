\name{read.helicos.tags}
\alias{read.helicos.tags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read in helicos tags }
\description{
Read in Helicos tags
}
\usage{
read.helicos.tags(filename, read.tag.names = F,
                  fix.chromosome.names = F, include.length.info = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
filename
}
  \item{read.tag.names}{
Read in tag names
}
  \item{fix.chromosome.names}{
Do we fix chromosome names
}
  \item{include.length.info}{
include length information
}
}
\value{
A list like structure
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
function (filename, read.tag.names = F, fix.chromosome.names = F, 
    include.length.info = T) 
{
    if (read.tag.names) {
        rtn <- as.integer(1)
    }
    else {
        rtn <- as.integer(0)
    }
    tl <- lapply(.Call("read_helicostabf", path.expand(filename), 
        rtn), function(d) {
        xo <- order(abs(d$t))
        d$t <- d$t[xo]
        d$n <- d$n[xo]
        d$l <- d$l[xo]
        if (read.tag.names) {
            d$s <- d$s[xo]
        }
        return(d)
    })
    if (fix.chromosome.names) {
        names(tl) <- gsub("\\.fa", "", names(tl))
    }
    if (read.tag.names) {
        return(list(tags = lapply(tl, function(d) d$t), quality = lapply(tl, 
            function(d) d$n), length = lapply(tl, function(d) d$l), 
            names = lapply(tl, function(d) d$s)))
    }
    else {
        return(list(tags = lapply(tl, function(d) d$t), quality = lapply(tl, 
            function(d) d$n), length = lapply(tl, function(d) d$l)))
    }
  }
}
}

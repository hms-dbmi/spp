\name{points_within}
\alias{points_within}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Find points within
}
\description{
points_within substitutes the deprecated function 
"points.within"
}
\usage{
points_within(x, fs, fe, return.list = F, return.unique = F,
              sorted = F, return.point.counts = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Parameter
}
  \item{fs}{
Parameter
}
  \item{fe}{
Parameter
}
  \item{return.list}{
Parameter
}
  \item{return.unique}{
Parameter
}
  \item{sorted}{
Parameter
}
  \item{return.point.counts}{
Parameter
}
  \item{\dots}{
Parameter
}
}
\value{
Parameter
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
\dontrun{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, fs, fe, return.list = F, return.unique = F, sorted = F, 
    return.point.counts = F, ...) 
{
    if (is.null(x) | length(x) < 1) {
        return(c())
    }
    if (!sorted) {
        ox <- rank(x, ties.method = "first")
        x <- sort(x)
    }
    se <- c(fs, fe)
    fi <- seq(1:length(fs))
    fi <- c(fi, -1 * fi)
    fi <- fi[order(se)]
    se <- sort(se)
    storage.mode(x) <- storage.mode(fi) <- storage.mode(se) <- "integer"
    if (return.unique) {
        iu <- 1
    }
    else {
        iu <- 0
    }
    if (return.list) {
        il <- 1
    }
    else {
        il <- 0
    }
    if (return.point.counts) {
        rpc <- 1
    }
    else {
        rpc <- 0
    }
    storage.mode(iu) <- storage.mode(il) <- storage.mode(rpc) <- "integer"
    result <- .Call("points_withinC", x, se, fi, il, iu, rpc)
    if (!sorted & !return.point.counts) {
        result <- result[ox]
    }
    return(result)
  }
}
}

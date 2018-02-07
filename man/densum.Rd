\name{densum}
\alias{densum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Do Something
}
\description{
Densum
}
\usage{
densum(vin, bw = 5, dw = 3, match.wt.f = NULL, return.x = T,
       from = min(vin), to = max(vin), step = 1, new.code = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vin}{
Parameter
}
  \item{bw}{
Parameter
}
  \item{dw}{
Parameter
}
  \item{match.wt.f}{
Parameter
}
  \item{return.x}{
Parameter
}
  \item{from}{
Parameter
}
  \item{to}{
Parameter
}
  \item{step}{
Parameter
}
  \item{new.code}{
Parameter
}
}
\value{
Some sum
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (vin, bw = 5, dw = 3, match.wt.f = NULL, return.x = T, 
    from = min(vin), to = max(vin), step = 1, new.code = T) 
{
    tc <- table(vin[vin >= from & vin <= to])
    pos <- as.numeric(names(tc))
    storage.mode(pos) <- "double"
    tc <- as.numeric(tc)
    storage.mode(tc) <- "double"
    n <- length(pos)
    if (!is.null(match.wt.f)) {
        tc <- tc * match.wt.f(pos)
    }
    rng <- c(from, to)
    if (rng[1] < 0) {
        stop("range extends into negative values")
    }
    if (range(pos)[1] < 0) {
        stop("position vector contains negative values")
    }
    storage.mode(n) <- storage.mode(rng) <- storage.mode(bw)
          <- storage.mode(dw) <- storage.mode(step) <- "integer"
    spos <- rng[1]
    storage.mode(spos) <- "double"
    dlength <- floor((rng[2] - rng[1])/step) + 1
    if (dlength < 1) {
        stop("zero data range")
    }
    if (new.code) {
        storage.mode(step) <- storage.mode(dlength)
              <- storage.mode(bw) <- storage.mode(dw) <- "integer"
        dout <- .Call("ccdensum", pos, tc, spos, bw, dw, dlength, 
            step)
    }
    else {
        dout <- numeric(dlength)
        storage.mode(dout) <- "double"
        storage.mode(dlength) <- "integer"
        .C("cdensum", n, pos, tc, spos, bw, dw, dlength, step, 
            dout)
    }
    if (return.x) {
        return(list(x = c(rng[1], rng[1] + step * (dlength - 
            1)), y = dout, step = step))
    }
    else {
        return(dout)
    }
  }
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

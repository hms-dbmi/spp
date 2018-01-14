\name{write.narrowpeak.binding}
\alias{write.narrowpeak.binding}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Write out determined binding peaks using narrowPeak format }
\description{
  Writes out determined binding positions into a narrowPeak file.
  The region will correspond to associated broad enrichment region, if
  such were added using add.broad.peak.regions method. Otherwise the
  region size will be determined using margin (which defaults to the
  window half size that was used to determine binding positions)
  Note: since v1.13, FDR is written out in -log10() scale.
}
\usage{
write.narrowpeak.binding(bd, fname, margin=bd$whs, npeaks)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bd}{ output of the \code{\link{find.binding.positions}} }
  \item{fname}{ file name }
  \item{margin}{ explicit value of the margin to be used if the broad
    region information is absent (defaults to peak detection window
    half-size}
  \item{npeaks}{ optionally, limit the output to the specified number of
    top peaks }
}

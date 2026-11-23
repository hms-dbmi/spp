# Helpers for the tag-reader tests. The C++ readers print progress with Rprintf,
# so every call goes through quiet().

quiet <- function(expr) {
	invisible(capture.output(res <- expr))
	res
}

# Write lines (always newline terminated) to a fresh temporary file.
write_tmp <- function(lines, fileext = "") {
	f <- tempfile(fileext = fileext)
	writeLines(lines, f)
	f
}

# Join fields with tabs, one line per row (rows are character vectors).
tabline <- function(...) paste(c(...), collapse = "\t")

# Strip dimnames/attributes that tapply() adds so that values compare plainly.
plain <- function(x) as.numeric(unlist(x, use.names = FALSE))

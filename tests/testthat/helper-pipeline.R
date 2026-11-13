# Helpers shared by the test-pipeline-*.R files.
# Tag lists are named lists (one element per chromosome) of position-sorted
# integer 5' positions; the sign is the strand (negative = minus strand).

# Evaluate expr with a fixed seed under the "Rejection" sample.kind, then put
# back the previous RNG kind and state (or its absence), so tests are neither
# sensitive to the session's RNGkind nor leak a seed into the global state.
pl.seeded <- function(seed, expr) {
	genv <- globalenv()
	old.kind <- RNGkind()
	had.seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
	if (had.seed) old.seed <- get(".Random.seed", envir = genv, inherits = FALSE)
	on.exit({
		suppressWarnings(RNGkind(old.kind[1], old.kind[2], old.kind[3]))
		if (had.seed) {
			assign(".Random.seed", old.seed, envir = genv)
		} else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
			rm(".Random.seed", envir = genv)
		}
	})
	suppressWarnings(RNGkind(sample.kind = "Rejection"))
	set.seed(seed)
	expr
}

# Run an expression, dropping printed progress dots, messages and the known
# harmless warnings of the C code and of tiny inputs; other warnings pass through.
pl.quiet <- function(expr, ignore = c("stack imbalance", "too few tags for shuffling", "no non-missing arguments", "adjusted down wig")) {
	res <- NULL
	utils::capture.output(res <- suppressMessages(withCallingHandlers(expr, warning = function(w) {
		known <- vapply(ignore, function(p) grepl(p, conditionMessage(w), fixed = TRUE), logical(1))
		if (any(known)) invokeRestart("muffleWarning")
	})))
	res
}

# A ChIP-seq-like pile around binding site c: n plus-strand tags centred on
# c-d and n minus-strand tags centred on c+d, offsets spread evenly over
# -spread..spread (deterministic and symmetric).
pl.pile <- function(c, n = 30, d = 50, spread = 5) {
	o <- as.integer(round(seq(-spread, spread, length.out = n)))
	list(pos = as.integer(c - d + o), neg = as.integer(c + d - o))
}

# signed, position-sorted tag vector made of one pile per centre
pl.tags <- function(centers, n = 30, ...) {
	pos <- integer(0)
	neg <- integer(0)
	for (c in centers) {
		p <- pl.pile(c, n, ...)
		pos <- c(pos, p$pos)
		neg <- c(neg, p$neg)
	}
	tg <- c(pos, -neg)
	tg[order(abs(tg))]
}

# uniform background tags with random strand
pl.bg <- function(n, L = 12000, seed = 1) {
	pl.seeded(seed, {
		x <- sample(L, n)
		s <- sample(c(-1L, 1L), n, replace = TRUE)
		y <- as.integer(x * s)
		y[order(abs(y))]
	})
}

# the standard data set: piles at 3000 and 8000 on chr1, 5000 on chr2, and an
# empty chrE
pl.signal <- function() {
	list(chr1 = pl.tags(c(3000, 8000)), chr2 = pl.tags(5000), chrE = integer(0))
}

# Control: uniform background, with the tags within `gap` bp of a pile centre
# removed so that background subtraction does not move the called peak positions.
pl.control <- function(gap = 160) {
	clear <- function(x, centers) x[vapply(abs(x), function(p) all(abs(p - centers) > gap), logical(1))]
	list(chr1 = clear(pl.bg(200, seed = 1), c(3000, 8000)), chr2 = clear(pl.bg(100, seed = 2), 5000), chrE = integer(0))
}

pl.ntags <- function(tl) sum(vapply(tl, length, numeric(1)))

# Brute-force Gaussian kernel sum exp(-d^2/(2 bw^2)) over tags (shifted, folded
# to positive) inside [from, to].  densum truncates the kernel near 3 bw (and
# unevenly on the two sides, see BUGS-pipeline.txt), so it agrees with this up
# to the dropped tail; see pl.close.
pl.density <- function(tags, x, bw, shift, from, to) {
	t0 <- abs(tags + shift)
	t0 <- t0[t0 >= from & t0 <= to]
	vapply(x, function(p) sum(exp(-(p - t0)^2 / (2 * bw^2))), numeric(1))
}

# TRUE when two vectors differ by at most tol everywhere
pl.close <- function(a, b, tol) {
	length(a) == length(b) && all(abs(a - b) <= tol)
}

# number of values in a window centred at p of width w (inclusive ends)
pl.count <- function(t, p, w) {
	vapply(p, function(q) sum(t >= q - w / 2 & t <= q + w / 2), numeric(1))
}

# number of values within +-h of p (inclusive)
pl.count.within <- function(t, p, h) {
	vapply(p, function(q) sum(t >= q - h & t <= q + h), numeric(1))
}

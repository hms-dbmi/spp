# Shared helpers and pure-R reference implementations for the kernel tests.

spp_internal <- function(name) {
	get(name, envir = asNamespace("spp"))
}

# The wtd/lwcc native routines leave one PROTECT unbalanced when
# return.peaks is FALSE, which makes R print "stack imbalance" warnings.
# Those are listed in BUGS-kernels.txt and are not what the tests check.
quiet <- function(expr) {
	suppressWarnings(expr)
}

# Brute-force Gaussian kernel sum used by densum().  Contributions are
# summed over every output position, with no window cutoff.
ref_density <- function(v, bw, from, to, step, wt = NULL) {
	v <- v[v >= from & v <= to]
	npos <- floor((to - from) / step) + 1
	grid <- from + step * (seq_len(npos) - 1)
	out <- numeric(npos)
	for (p in unique(v)) {
		ic <- sum(v == p)
		w <- ic
		if (!is.null(wt)) w <- ic * wt(p)
		out <- out + w * exp(-0.5 * ((grid - p) / bw)^2)
	}
	out
}

# Windowed Gaussian sum as written in ccdensum: each unique position p with
# multiplicity ic adds ic * exp(-0.5 * ((j * step - in) / bw)^2) to output j for
# j in [max(0, floor((in - whs) / step)), min(floor((in + whs) / step), npos - 1) - 1],
# where in = p - from and whs = dw * bw * ic.  This mirrors the C window rule
# (upper end exclusive); the final output element is never written, so it is 0.
ref_density_win <- function(v, bw, dw, from, to, step) {
	v <- v[v >= from & v <= to]
	npos <- floor((to - from) / step) + 1
	out <- numeric(npos)
	for (p in unique(v)) {
		ic <- sum(v == p)
		inn <- p - from
		whs <- dw * bw * ic
		ws <- max(0, floor((inn - whs) / step))
		we <- min(floor((inn + whs) / step), npos - 1)
		for (j in seq_len(max(0, we - ws)) + ws - 1) {
			out[j + 1] <- out[j + 1] + ic * exp(-0.5 * ((j * step - inn) / bw)^2)
		}
	}
	out
}

# NOTE on ref_wtd / ref_lwcc below: they mirror the arithmetic of the C code
# (background windows are bg.whs + 1 wide, scoring starts at whs + 1 for wtd and
# at os for lwcc, the lwcc means divide by 2 * os + 1 although 2 * (os - is + 1)
# positions are sampled, round.up applies to different count terms in wtd and
# lwcc, masked pairs are skipped together).  They are consistency checks of the
# native code against a plain-R transcription of its formulas, not independent
# definitions.  Independent hand-derived cases are in the test files; no
# hand derivation was possible for the lwcc mean divisor, since which divisor is
# intended is not documented.

# Window sums over a pair of strand histograms (0-based index i, windows
# [i - whs, i - 1] on the left and [i, i + whs - 1] on the right).
ref_wtd <- function(hx, hy, whs, tag.weight = 1, direct.count = FALSE,
		round.up = FALSE, bg = NULL, bg.whs = whs, bg.weight = 1,
		masked = FALSE) {
	n <- length(hx)
	out <- numeric(n)
	for (i in seq(whs + 1, n - whs - 2)) {
		pn1 <- pn2 <- nn1 <- nn2 <- 0
		for (k in 0:(whs - 1)) {
			a <- hx[i - k]; b <- hy[i - k]       # position i-k-1 (0-based)
			c <- hx[i + k + 1]; d <- hy[i + k + 1] # position i+k
			if (masked && (a == -1 || b == -1 || c == -1 || d == -1)) next
			pn1 <- pn1 + a; nn1 <- nn1 + b
			pn2 <- pn2 + c; nn2 <- nn2 + d
		}
		spn1 <- pn1 * tag.weight; snn1 <- nn1 * tag.weight
		spn2 <- pn2 * tag.weight; snn2 <- nn2 * tag.weight
		if (round.up) {
			if (pn1 > 0 && spn1 < 1) spn1 <- 1
			if (nn2 > 0 && snn2 < 1) snn2 <- 1
		}
		if (!is.null(bg)) {
			lo <- max(0, i - bg.whs); hi <- min(n - 1, i + bg.whs)
			bgp1 <- sum(bg$x[(lo:i) + 1]); bgn1 <- sum(bg$y[(lo:i) + 1])
			bgp2 <- sum(bg$x[(i:hi) + 1]); bgn2 <- sum(bg$y[(i:hi) + 1])
		}
		if (direct.count) {
			val <- spn1 + snn2
			if (round.up && val < 1) val <- 1
			if (!is.null(bg)) val <- val - (bgp1 + bgn2) * bg.weight
		} else if (!is.null(bg)) {
			spn1 <- spn1 - bgp1 * bg.weight; snn1 <- snn1 - bgn1 * bg.weight
			spn2 <- spn2 - bgp2 * bg.weight; snn2 <- snn2 - bgn2 * bg.weight
			spn2 <- max(spn2, 0); snn1 <- max(snn1, 0)
			if (spn1 > 0 && snn2 > 0) {
				val <- 2 * sqrt(spn1 * snn2) - (spn2 + snn1 + 1)
			} else {
				val <- -(spn2 + snn1 + 1)
			}
		} else {
			val <- 2 * sqrt(spn1 * snn2) - (spn2 + snn1 + tag.weight)
		}
		out[i + 1] <- val
	}
	out
}

# Strand histogram over positions rx[1]..rx[2]
ref_hist <- function(v, rx) {
	tabulate(v - rx[1] + 1, diff(rx) + 1)
}

# Weighted lwcc score over (0-based) positions os .. n-os-2
ref_lwcc <- function(hx, hy, os, is, tag.weight = 1, round.up = FALSE,
		bg = NULL, bg.whs = os, bg.weight = 1, masked = FALSE) {
	n <- length(hx)
	rs <- 2 * os + 1
	out <- numeric(n)
	for (i in seq(os, n - os - 2)) {
		ks <- 0:(os - is)
		xp1 <- hx[i - os + ks + 1]; xp2 <- hx[i + os - ks + 1]
		xn1 <- hy[i + os - ks + 1]; xn2 <- hy[i - os + ks + 1]
		ok1 <- rep(TRUE, length(ks)); ok2 <- ok1
		if (masked) {
			ok1 <- xp1 != -1 & xn1 != -1
			ok2 <- xp2 != -1 & xn2 != -1
		}
		pn1 <- sum(xp1[ok1]); nn1 <- sum(xn1[ok1])
		pn2 <- sum(xp2[ok2]); nn2 <- sum(xn2[ok2])
		mp <- (pn1 + pn2) / rs
		mn <- (nn1 + nn2) / rs
		val <- -1e3
		if (mp > 0 && mn > 0) {
			a <- c(xp1[ok1], xp2[ok2]) - mp
			b <- c(xn1[ok1], xn2[ok2]) - mn
			num <- sum(a * b); varp <- sum(a * a); varn <- sum(b * b)
			spn1 <- pn1 * tag.weight; snn1 <- nn1 * tag.weight
			spn2 <- pn2 * tag.weight; snn2 <- nn2 * tag.weight
			if (round.up) {
				if (pn1 > 0 && spn1 < 1) spn1 <- 1
				if (nn1 > 0 && snn1 < 1) snn1 <- 1
			}
			if (!is.null(bg)) {
				lo <- max(0, i - bg.whs); hi <- min(n - 1, i + bg.whs)
				bgp1 <- sum(bg$x[(lo:i) + 1]); bgn1 <- sum(bg$y[(lo:i) + 1])
				bgp2 <- sum(bg$x[(i:hi) + 1]); bgn2 <- sum(bg$y[(i:hi) + 1])
				spn1 <- spn1 - bgp1 * bg.weight
				snn1 <- snn1 - bgn2 * bg.weight
				spn2 <- spn2 - bgp2 * bg.weight
				snn2 <- snn2 - bgn1 * bg.weight
				spn2 <- max(spn2, 0); snn2 <- max(snn2, 0)
				if (spn1 > 0 && snn1 > 0) {
					tagw <- 2 * sqrt(spn1 * snn1) - (spn2 + snn2 + 1)
				} else {
					tagw <- -(spn2 + snn2 + 1)
				}
			} else {
				tagw <- 2 * sqrt(spn1 * snn1) - (spn2 + snn2)
			}
			if (tagw < 0) {
				val <- 0
			} else if (num == 0) {
				val <- tagw
			} else {
				val <- num / sqrt(varp * varn) * sqrt(tagw) + tagw
			}
		}
		out[i + 1] <- val
	}
	out
}

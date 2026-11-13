# get.binding.characteristics(): strand cross-correlation (tag.scc, in R),
# peak, and window half-size.

tag.scc <- spp_internal("tag.scc")

# Naive strand cross-correlation of position histograms.  For a shift s it is
# the Pearson-like sum over all positions p in [lo, hi] of
# (pos[p] - mean) * (neg[p + s] - mean), normalised by the strand variances,
# with the negative-strand histogram taken as zero outside [lo, hi].
naive_scc <- function(pos, neg, shifts) {
	lo <- min(pos, neg); hi <- max(pos, neg); n <- hi - lo + 1
	hp <- tabulate(pos - lo + 1, n); hn <- tabulate(neg - lo + 1, n)
	mp <- sum(hp) / n; mn <- sum(hn) / n
	ss <- sqrt(sum((hp - mp)^2) * sum((hn - mn)^2))
	vapply(shifts, function(s) {
		q <- seq_len(n) + s
		hq <- ifelse(q >= 1 & q <= n, hn[pmin(pmax(q, 1), n)], 0)
		sum((hp - mp) * (hq - mn)) / ss
	}, numeric(1))
}

# Tags where every minus-strand tag lies exactly `shift` bp after a plus-strand
# tag.  A lone plus-strand tag at position 1 fixes the left edge of the range so
# that no minus-strand tag falls closer than 200 bp to it.
make_tags <- function(n, shift, seed) {
	set.seed(seed)
	p <- sample(200:20000, n)
	list(pos = c(1, p), neg = p + shift)
}
signed <- function(t) c(t$pos, -t$neg)

test_that("tag.scc agrees with a naive cross-correlation", {
	t <- make_tags(800, 90, 1)
	r <- tag.scc(signed(t), srange = c(40, 180), bin = 1)
	expect_equal(as.numeric(names(r)), 40:180)
	expect_equal(as.numeric(r), naive_scc(t$pos, t$neg, 40:180), tolerance = 1e-10)
})

test_that("the cross-correlation peak is at the fragment shift", {
	t <- make_tags(1500, 100, 5)
	d <- list(tags = list(chr1 = signed(t)))
	r <- get.binding.characteristics(d, srange = c(50, 200), bin = 1,
		remove.tag.anomalies = FALSE)
	expect_named(r, c("cross.correlation", "peak", "whs"))
	expect_equal(r$peak$x, 100)
	expect_equal(r$cross.correlation$y[r$cross.correlation$x == 100], r$peak$y)
	# a one-to-one shifted copy correlates almost perfectly (the single lone
	# plus-strand tag lowers it slightly)
	expect_gt(r$peak$y, 0.99)
	expect_lte(r$peak$y, 1 + 1e-10)
})

test_that("the profile equals the naive cross-correlation at every shift", {
	t <- make_tags(1500, 100, 5)
	d <- list(tags = list(chr1 = signed(t)))
	r <- get.binding.characteristics(d, srange = c(50, 200), bin = 1,
		remove.tag.anomalies = FALSE)
	cc <- r$cross.correlation
	expect_equal(cc$x, 50:200)
	expect_equal(cc$y, naive_scc(t$pos, t$neg, cc$x), tolerance = 1e-10)
})

test_that("srange sets the range of lags reported", {
	t <- make_tags(1000, 70, 6)
	d <- list(tags = list(chr1 = signed(t)))
	r <- get.binding.characteristics(d, srange = c(30, 150), bin = 1,
		remove.tag.anomalies = FALSE)
	expect_equal(range(r$cross.correlation$x), c(30, 150))
	expect_equal(r$peak$x, 70)
})

# whs is defined relative to the last lag of the profile (its baseline), not to
# zero; this test restates that definition (it pins the implementation).
test_that("whs is the widest lag with correlation >= last lag + (peak - last lag) / 3", {
	t <- make_tags(1500, 100, 5)
	d <- list(tags = list(chr1 = signed(t)))
	r <- get.binding.characteristics(d, srange = c(50, 200), bin = 1,
		remove.tag.anomalies = FALSE)
	cc <- r$cross.correlation
	last <- cc$y[length(cc$y)]
	th <- (r$peak$y - last) / 3 + last
	expect_equal(r$whs, max(cc$x[cc$y >= th]))
})

test_that("bin = 5 evaluates lags on a 5 bp grid and finds the same peak", {
	# the profile values themselves are not compared with naive_scc on binned
	# positions: tag.scc scales the strand means by bin (suspected bug 10 in
	# BUGS-kernels.txt), which offsets the whole profile.
	t <- make_tags(1500, 100, 5)
	d <- list(tags = list(chr1 = signed(t)))
	r <- get.binding.characteristics(d, srange = c(50, 200), bin = 5,
		remove.tag.anomalies = FALSE)
	expect_equal(r$cross.correlation$x, seq(50, 200, by = 5))
	expect_equal(r$peak$x, 100)
	expect_equal(which.max(r$cross.correlation$y), which(r$cross.correlation$x == 100))
})

test_that("chromosome profiles are averaged weighted by tag counts", {
	t1 <- make_tags(1500, 100, 5)
	t2 <- make_tags(500, 100, 8)
	d <- list(tags = list(chr1 = signed(t1), chr2 = signed(t2)))
	r <- get.binding.characteristics(d, srange = c(50, 200), bin = 1,
		remove.tag.anomalies = FALSE)
	n1 <- length(signed(t1)); n2 <- length(signed(t2))
	shifts <- r$cross.correlation$x
	want <- (n1 * naive_scc(t1$pos, t1$neg, shifts) +
		n2 * naive_scc(t2$pos, t2$neg, shifts)) / (n1 + n2)
	expect_equal(r$cross.correlation$y, want, tolerance = 1e-10)
	expect_equal(r$peak$x, 100)
})

test_that("removing tag anomalies leaves the peak of clean data unchanged", {
	t <- make_tags(1500, 100, 5)
	d <- list(tags = list(chr1 = signed(t)))
	r <- get.binding.characteristics(d, srange = c(50, 200), bin = 1)
	expect_equal(r$peak$x, 100)
	expect_gt(r$peak$y, 0.9)
})

test_that("piled-up positions are removed when remove.tag.anomalies is TRUE", {
	t <- make_tags(1500, 100, 5)
	# 8 + tags at 30000 and 8 - tags at 30060: a pair 60 bp apart inside srange
	# that is not caught by tag.scc's own cap (10 x the mean count)
	tags <- c(signed(t), rep(30000, 8), rep(-30060, 8))
	d <- list(tags = list(chr1 = tags))
	on <- get.binding.characteristics(d, srange = c(50, 200), bin = 1)
	off <- get.binding.characteristics(d, srange = c(50, 200), bin = 1,
		remove.tag.anomalies = FALSE)
	# with the pile removed the profile is that of the clean tags
	expect_equal(on$peak$y, naive_scc(t$pos, t$neg, 100), tolerance = 1e-10)
	# with the pile kept the correlation at the true shift is lower
	expect_lt(off$peak$y, on$peak$y - 0.01)
})

test_that("remove.tag.anomalies drops positions far above the typical count", {
	remove.tag.anomalies <- spp_internal("remove.tag.anomalies")
	t <- make_tags(1500, 100, 5)
	tags <- c(signed(t), rep(30000, 8), rep(-30060, 8))
	r <- remove.tag.anomalies(list(tags = list(chr1 = tags)))
	expect_equal(length(r$tags$chr1), length(tags) - 16)
	expect_false(any(abs(r$tags$chr1) >= 30000))
})

test_that("quality scores add the bin acceptance table unless accept.all.tags", {
	t <- make_tags(1000, 100, 5)
	tags <- signed(t)
	d <- list(tags = list(chr1 = tags), quality = list(chr1 = rep(1, length(tags))))
	r <- get.binding.characteristics(d, srange = c(50, 200), bin = 1,
		remove.tag.anomalies = FALSE)
	expect_named(r, c("cross.correlation", "peak", "whs", "quality.bin.acceptance"))
	expect_equal(r$peak$x, 100)
	r2 <- get.binding.characteristics(d, srange = c(50, 200), bin = 1,
		remove.tag.anomalies = FALSE, accept.all.tags = TRUE)
	expect_named(r2, c("cross.correlation", "peak", "whs"))
	expect_equal(r2$peak$x, 100)
})

test_that("only the lowest quality bin is used for the profile", {
	t <- make_tags(1000, 100, 5)
	tags <- signed(t)
	# decoy tags in quality bin 2 would destroy the profile if they were used
	set.seed(77)
	decoy <- c(sample(1000:15000, 600), -sample(1000:15000, 600))
	d <- list(tags = list(chr1 = c(tags, decoy)),
		quality = list(chr1 = c(rep(1, length(tags)), rep(2, length(decoy)))))
	r <- get.binding.characteristics(d, srange = c(50, 200), bin = 1,
		remove.tag.anomalies = FALSE, min.tag.count = 1e6)
	expect_equal(r$cross.correlation$y,
		naive_scc(t$pos, t$neg, 50:200), tolerance = 1e-10)
	expect_named(r$quality.bin.acceptance, c("informative.bins", "quality.cc"))
})

test_that("empty tag lists give an empty profile and no peak", {
	d <- list(tags = list(chr1 = numeric(0)))
	r <- suppressWarnings(get.binding.characteristics(d, srange = c(1, 5),
		bin = 1, remove.tag.anomalies = FALSE))
	expect_equal(nrow(r$cross.correlation), 0)
	expect_length(r$peak$x, 0)
	expect_length(r$peak$y, 0)
})

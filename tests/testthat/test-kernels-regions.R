# get.relative.coordinates() (get_relative_coordinates), regionset.intersection.c()
# (region_intersection), tag.enrichment.clusters() and
# find.significantly.enriched.regions() (find_poisson_enrichment_clusters).
# All in peaks.cpp.

get.relative.coordinates <- spp_internal("get.relative.coordinates")
regionset.intersection.c <- spp_internal("regionset.intersection.c")
tag.enrichment.clusters <- spp_internal("tag.enrichment.clusters")
find.significantly.enriched.regions <- spp_internal("find.significantly.enriched.regions")

# ---- get.relative.coordinates ----

# Tests below use a single centre or a single data point on purpose: with
# several centres and several data points the native code can lose matches
# (suspected bug 3 in BUGS-kernels.txt), so such inputs are not asserted.
#
# brute-force pairs (relative coordinate, centre index), ordered by data point
# and then by |centre|
ref_relative <- function(d, pos, size) {
	d <- sort(d)
	op <- order(abs(pos))
	xs <- c(); is <- c()
	for (k in seq_along(d)) {
		for (j in op) {
			pd <- d[k] - abs(pos[j])
			if (abs(pd) <= size) {
				xs <- c(xs, sign(pos[j]) * pd)
				is <- c(is, j)
			}
		}
	}
	list(x = as.integer(xs), i = as.integer(is))
}

test_that("get.relative.coordinates gives the hand-computed coordinates", {
	# centre 8 (+): data 5 -> -3, data 10 -> 2; centre -12 (-): data 10 -> -(10-12) = 2
	r <- get.relative.coordinates(c(5, 10, 20), c(8, -12), 5)
	expect_named(r, c("x", "i"))
	expect_equal(r$x, c(-3, 2, 2))
	expect_equal(r$i, c(1, 1, 2))
})

test_that("get.relative.coordinates refers to the original centre order", {
	r <- get.relative.coordinates(c(5, 10, 20), c(-12, 8), 5)
	expect_equal(r$x, c(-3, 2, 2))
	expect_equal(r$i, c(2, 2, 1))
})

test_that("get.relative.coordinates sorts unsorted data", {
	a <- get.relative.coordinates(c(20, 5, 10), c(8, -12), 5)
	b <- get.relative.coordinates(c(5, 10, 20), c(8, -12), 5)
	expect_equal(a, b)
})

test_that("the size window includes distance == size and excludes size + 1", {
	expect_equal(get.relative.coordinates(13, 8, 5)$x, 5)
	expect_equal(get.relative.coordinates(3, 8, 5)$x, -5)
	expect_length(get.relative.coordinates(14, 8, 5)$x, 0)
	expect_length(get.relative.coordinates(2, 8, 5)$x, 0)
})

test_that("the sign of a centre flips the coordinate", {
	expect_equal(get.relative.coordinates(10, 7, 5)$x, 3)
	expect_equal(get.relative.coordinates(10, -7, 5)$x, -3)
})

test_that("get.relative.coordinates matches a brute force for a single centre", {
	set.seed(3)
	for (it in 1:50) {
		d <- sort(sample(1:200, sample(1:30, 1)))
		pos <- sample(c(1:200, -(1:200)), 1)
		size <- sample(1:30, 1)
		got <- get.relative.coordinates(d, pos, size)
		want <- ref_relative(d, pos, size)
		expect_equal(got$x, want$x)
		expect_equal(got$i, want$i)
	}
})

test_that("get.relative.coordinates matches a brute force for a single data point", {
	set.seed(4)
	for (it in 1:50) {
		d <- sample(1:200, 1)
		pos <- sample(c(1:200, -(1:200)), sample(1:10, 1))
		size <- sample(1:30, 1)
		got <- get.relative.coordinates(d, pos, size)
		want <- ref_relative(d, pos, size)
		expect_equal(got$x, want$x)
		expect_equal(got$i, want$i)
	}
})

test_that("get.relative.coordinates with no data points returns empty vectors", {
	r <- get.relative.coordinates(integer(0), c(8, -12), 5)
	expect_length(r$x, 0)
	expect_length(r$i, 0)
})

# ---- regionset.intersection.c ----

test_that("intersection of two region sets is computed by hand", {
	a <- list(s = c(1, 10, 30), e = c(5, 20, 40))
	b <- list(s = c(3, 15), e = c(12, 35))
	r <- regionset.intersection.c(list(a, b))
	expect_named(r, c("s", "e"))
	expect_equal(r$s, c(3, 10, 15, 30))
	expect_equal(r$e, c(5, 12, 20, 35))
})

test_that("intersection of three region sets", {
	r <- regionset.intersection.c(list(list(s = 1, e = 10), list(s = 5, e = 20),
		list(s = 8, e = 12)))
	expect_equal(r$s, 8)
	expect_equal(r$e, 10)
})

test_that("union of region sets merges overlapping regions", {
	r <- regionset.intersection.c(list(list(s = 1, e = 10), list(s = 5, e = 20),
		list(s = 8, e = 12)), do.union = TRUE)
	expect_equal(r$s, 1)
	expect_equal(r$e, 20)
	r2 <- regionset.intersection.c(list(list(s = 1, e = 5), list(s = 10, e = 20)),
		do.union = TRUE)
	expect_equal(r2$s, c(1, 10))
	expect_equal(r2$e, c(5, 20))
})

test_that("disjoint region sets have an empty intersection", {
	r <- regionset.intersection.c(list(list(s = 1, e = 5), list(s = 10, e = 20)))
	expect_equal(nrow(r), 0)
	expect_named(r, c("s", "e"))
})

test_that("a single region set is returned unchanged", {
	a <- list(s = c(1, 10, 30), e = c(5, 20, 40))
	r <- regionset.intersection.c(list(a))
	expect_equal(r$s, a$s)
	expect_equal(r$e, a$e)
})

test_that("the intersection does not depend on the order of the sets", {
	a <- list(s = c(1, 10, 30), e = c(5, 20, 40))
	b <- list(s = c(3, 15), e = c(12, 35))
	r1 <- regionset.intersection.c(list(a, b))
	r2 <- regionset.intersection.c(list(b, a))
	expect_equal(r1, r2)
})

test_that("region values are reported as the min (max.val = -1) or max (1) overlapping value", {
	a <- list(s = c(1, 10, 30), e = c(5, 20, 40), v = c(1, 2, 3))
	b <- list(s = c(3, 15), e = c(12, 35), v = c(5, 1))
	rmin <- regionset.intersection.c(list(a, b), max.val = -1)
	expect_named(rmin, c("s", "e", "v"))
	expect_equal(rmin$s, c(3, 10, 15, 30))
	expect_equal(rmin$e, c(5, 12, 20, 35))
	expect_equal(rmin$v, c(1, 2, 1, 1))
	rmax <- regionset.intersection.c(list(a, b), max.val = 1)
	expect_equal(rmax$s, rmin$s)
	expect_equal(rmax$e, rmin$e)
	expect_equal(rmax$v, c(5, 5, 2, 3))
	# the default is the minimum
	expect_equal(regionset.intersection.c(list(a, b)), rmin)
})

test_that("values of a union are the max over the merged regions", {
	a <- list(s = c(1, 10, 30), e = c(5, 20, 40), v = c(1, 2, 3))
	b <- list(s = c(3, 15), e = c(12, 35), v = c(5, 1))
	r <- regionset.intersection.c(list(a, b), max.val = 1, do.union = TRUE)
	expect_equal(r$s, 1)
	expect_equal(r$e, 40)
	expect_equal(r$v, 5)
})

# ---- tag.enrichment.clusters ----

tec <- function(signal, background, ...) {
	tag.enrichment.clusters(signal, background, thr = 0, tag.shift = 0, ...)
}

test_that("empty signal and background give no clusters", {
	r <- suppressWarnings(tag.enrichment.clusters(numeric(0), numeric(0)))
	expect_named(r, c("s", "e"))
	expect_length(r$s, 0)
	expect_length(r$e, 0)
})

test_that("an enrichment ratio no window can reach gives no clusters", {
	r <- tec(c(100, 110, 500, 505, 510), numeric(0), wsize = 20, mcs = 0,
		bg.weight = 1000)
	expect_length(r$s, 0)
	expect_length(r$e, 0)
})

test_that("a signal pair gives one cluster inside the span of qualifying windows", {
	# thr = 0, no background: a window qualifies when it holds a signal tag, i.e.
	# when its start lies in [80, 110]; window centres are start + 10.
	r <- tec(c(100, 110), numeric(0), wsize = 20, mcs = 0, bg.weight = 2)
	expect_length(r$s, 1)
	expect_gte(r$s, 90)
	expect_lte(r$e, 120)
	expect_lte(r$s, r$e)
})

test_that("a regular run of signal tags gives one cluster spanning the tags", {
	# tags every 5 bp from 1000 to 1100, windows of 20 bp, every window holds tags
	r <- tec(seq(1000, 1100, by = 5), numeric(0), wsize = 20, mcs = 0, bg.weight = 2)
	expect_length(r$s, 1)
	expect_gte(r$s, 990)
	expect_lte(r$s, 1020)
	expect_gte(r$e, 1090)
	expect_lte(r$e, 1110)
})

test_that("mcs drops clusters shorter than the minimum size", {
	sig <- c(100, 110, 500, 505, 510, 900, 905, 910, 915, 920, 925, 930)
	base <- tec(sig, numeric(0), wsize = 20, mcs = 0, bg.weight = 2)
	expect_gt(length(base$s), 1)
	for (mcs in c(1, 3, 5, 8, 20)) {
		r <- tec(sig, numeric(0), wsize = 20, mcs = mcs, bg.weight = 2)
		keep <- base$e - base$s >= mcs
		expect_equal(r$s, base$s[keep])
		expect_equal(r$e, base$e[keep])
	}
})

test_that("tag.shift moves all clusters by the shift", {
	sig <- c(100, 110, 500, 505, 510)
	bg <- c(105, 300, 502)
	a <- tag.enrichment.clusters(sig, bg, wsize = 20, thr = 0, mcs = 0,
		bg.weight = 2, tag.shift = 0)
	b <- tag.enrichment.clusters(sig, bg, wsize = 20, thr = 0, mcs = 0,
		bg.weight = 2, tag.shift = 10)
	expect_gt(length(a$s), 0)
	expect_equal(b$s, a$s + 10)
	expect_equal(b$e, a$e + 10)
})

test_that("negative-strand tags are folded onto their absolute position", {
	a <- tec(c(100, 110, 500, 505, 510), numeric(0), wsize = 20, mcs = 0, bg.weight = 2)
	b <- tec(c(-100, -110, 500, 505, 510), numeric(0), wsize = 20, mcs = 0, bg.weight = 2)
	expect_equal(a, b)
})

test_that("either = TRUE accepts windows that meet only the tag-count criterion", {
	sig <- c(100, 110, 500, 505, 510)
	strict <- tec(sig, numeric(0), wsize = 20, mcs = 0, bg.weight = 1000)
	loose <- tec(sig, numeric(0), wsize = 20, mcs = 0, bg.weight = 1000, either = TRUE)
	expect_length(strict$s, 0)
	expect_gt(length(loose$s), 0)
})

test_that("an unreachable minimum tag count removes all clusters", {
	sig <- c(100, 110, 500, 505, 510)
	r <- tec(sig, numeric(0), wsize = 20, mcs = 0, bg.weight = 2,
		min.tag.count.z = 1, min.tag.count.ratio = 1e6, tag.av.den = 1)
	expect_length(r$s, 0)
	# with a negligible expected count the threshold is zero and nothing changes
	same <- tec(sig, numeric(0), wsize = 20, mcs = 0, bg.weight = 2,
		min.tag.count.z = 1, min.tag.count.ratio = 1e-6, tag.av.den = 1)
	base <- tec(sig, numeric(0), wsize = 20, mcs = 0, bg.weight = 2)
	expect_equal(same, base)
})

# ---- find.significantly.enriched.regions ----

test_that("find.significantly.enriched.regions pads clusters by half the masking window", {
	sig <- list(chr1 = c(100, 110))
	ctl <- list(chr1 = numeric(0))
	base <- tag.enrichment.clusters(sig$chr1, ctl$chr1, wsize = 20, thr = 0, mcs = 0,
		bg.weight = 2, tag.shift = 0)
	expect_gt(length(base$s), 0)
	r <- find.significantly.enriched.regions(sig, ctl, window.size = 20, z.thr = 0,
		bg.weight = 2, tag.shift = 0)
	expect_named(r, "chr1")
	expect_equal(r$chr1$s, base$s - 10)
	expect_equal(r$chr1$e, base$e + 10)
	r2 <- find.significantly.enriched.regions(sig, ctl, window.size = 20, z.thr = 0,
		bg.weight = 2, tag.shift = 0, masking.window.size = 40)
	expect_equal(r2$chr1$s, base$s - 20)
	expect_equal(r2$chr1$e, base$e + 20)
})

test_that("find.significantly.enriched.regions handles several chromosomes and multiplier", {
	sig <- list(chr1 = c(100, 110), chr2 = c(300, 310, 700))
	ctl <- list(chr1 = numeric(0), chr2 = numeric(0))
	r <- find.significantly.enriched.regions(sig, ctl, window.size = 20, z.thr = 0,
		bg.weight = 2, tag.shift = 0)
	expect_named(r, c("chr1", "chr2"))
	for (chr in names(sig)) {
		base <- tag.enrichment.clusters(sig[[chr]], ctl[[chr]], wsize = 20, thr = 0,
			mcs = 0, bg.weight = 2, tag.shift = 0)
		expect_equal(r[[chr]]$s, base$s - 10)
		expect_equal(r[[chr]]$e, base$e + 10)
	}
	expect_gt(length(r$chr1$s), 0)
	expect_gt(length(r$chr2$s), 0)
	hi <- find.significantly.enriched.regions(sig, ctl, window.size = 20, z.thr = 0,
		bg.weight = 2, multiplier = 1000, tag.shift = 0)
	expect_length(hi$chr1$s, 0)
	expect_length(hi$chr2$s, 0)
})

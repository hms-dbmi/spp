# get.smoothed.tag.density, get.smoothed.enrichment.mle, get.smoothed.enrichment.mle2

# random signed tags whose shifted positions (shift 10 or 73) are all distinct,
# so that no position carries more than one tag
pl.dens.tags <- function(n, seed, lo = 100, hi = 2000) {
	repeat {
		a <- pl.seeded(seed, as.integer(sample(lo:hi, n) * sample(c(-1, 1), n, replace = TRUE)))
		ok <- all(vapply(c(10, 73), function(s) !anyDuplicated(abs(a + s)), logical(1)))
		if (ok) return(a[order(abs(a))])
		seed <- seed + 1000
	}
}

pl.dens.data <- function() {
	list(a = pl.dens.tags(60, 21), b = c(500L, -700L, 900L, 1200L))
}

pl.dens.ctl <- function() {
	list(a = pl.dens.tags(30, 22), b = c(600L, -800L, 1000L, 1500L))
}

# smoothed density with the defaults get.smoothed.enrichment.mle uses (shift 73)
pl.sd <- function(tags, rngl, bw = 30, step = 1) {
	get.smoothed.tag.density(tags, rngl = rngl, bandwidth = bw, step = step)
}

test_that("a single tag gives a Gaussian bump at its shifted position", {
	d <- get.smoothed.tag.density(list(a = 1000L), bandwidth = 30, tag.shift = 0, step = 10, rngl = list(a = c(900, 1100)))
	expect_equal(d$a$x, seq(900, 1100, by = 10))
	y <- setNames(d$a$y, d$a$x)
	# values at 0, 1 and 2 bandwidths (the exact kernel; a 3% change of bandwidth moves exp(-2) by 6%)
	expect_equal(unname(y["1000"]), 1)
	expect_equal(unname(y[c("970", "1030")]), rep(exp(-0.5), 2))
	expect_equal(unname(y[c("940", "1060")]), rep(exp(-2), 2))
	expect_lte(max(y), 1 + 1e-12)
	expect_true(all(y >= 0))
})

test_that("tag.shift moves plus-strand tags up and minus-strand tags down", {
	rng <- list(a = c(500, 1300))
	plus <- get.smoothed.tag.density(list(a = 900L), bandwidth = 30, tag.shift = 100, step = 10, rngl = rng)
	minus <- get.smoothed.tag.density(list(a = -900L), bandwidth = 30, tag.shift = 100, step = 10, rngl = rng)
	expect_equal(plus$a$x[which.max(plus$a$y)], 1000)
	expect_equal(minus$a$x[which.max(minus$a$y)], 800)
	expect_equal(max(plus$a$y), 1)
})

test_that("the default range and step follow the data and the bandwidth", {
	s <- pl.dens.data()
	d <- get.smoothed.tag.density(s, bandwidth = 30, tag.shift = 10)
	expect_equal(names(d), c("a", "b"))
	for (chr in names(s)) {
		t0 <- abs(s[[chr]] + 10)
		expect_equal(d[[chr]]$x[1], min(t0))
		expect_true(all(diff(d[[chr]]$x) == 10))
		expect_lte(max(d[[chr]]$x), max(t0))
		expect_gt(max(d[[chr]]$x), max(t0) - 10)
	}
	d2 <- get.smoothed.tag.density(s, bandwidth = 30, tag.shift = 10, step = 7)
	expect_true(all(diff(d2$a$x) == 7))
})

test_that("the density matches a brute-force Gaussian kernel sum", {
	s <- pl.dens.data()
	d <- get.smoothed.tag.density(s, bandwidth = 30, tag.shift = 10, step = 1)
	for (chr in names(s)) {
		x <- d[[chr]]$x
		rng <- range(abs(s[[chr]] + 10))
		ref <- pl.density(s[[chr]], x, 30, 10, rng[1], rng[2])
		inner <- x >= rng[1] + 90 & x <= rng[2] - 90
		expect_gt(sum(inner), 5)
		# densum drops the kernel tail beyond ~3 bandwidths: up to ~0.05 per point on this sparse data
		expect_true(pl.close(d[[chr]]$y[inner], ref[inner], 0.06))
	}
	# a 10% wrong bandwidth would be far outside that tolerance
	x <- d$a$x
	rng <- range(abs(s$a + 10))
	inner <- x >= rng[1] + 90 & x <= rng[2] - 90
	wrong <- pl.density(s$a, x, 33, 10, rng[1], rng[2])
	expect_false(pl.close(d$a$y[inner], wrong[inner], 0.06))
})

test_that("rngl restricts the chromosomes and the range", {
	s <- pl.dens.data()
	d <- get.smoothed.tag.density(s, bandwidth = 30, tag.shift = 10, step = 10, rngl = list(a = c(500, 1000)))
	expect_equal(names(d), "a")
	expect_equal(d$a$x, seq(500, 1000, by = 10))
	d1 <- get.smoothed.tag.density(s, bandwidth = 30, tag.shift = 10, step = 1, rngl = list(a = c(500, 1000)))
	ref <- pl.density(s$a, d1$a$x, 30, 10, 500, 1000)
	inner <- d1$a$x >= 590 & d1$a$x <= 910
	expect_true(pl.close(d1$a$y[inner], ref[inner], 0.06))
})

test_that("control density is subtracted with the tag-count ratio as weight", {
	s <- pl.dens.data()
	cc <- pl.dens.ctl()
	rngl <- list(a = c(300, 1800))
	sig <- get.smoothed.tag.density(s, bandwidth = 30, tag.shift = 10, step = 10, rngl = rngl)
	ctl <- get.smoothed.tag.density(cc, bandwidth = 30, tag.shift = 10, step = 10, rngl = rngl)
	bgw <- sum(vapply(s, length, numeric(1))) / sum(vapply(cc, length, numeric(1)))
	d <- get.smoothed.tag.density(s, cc, bandwidth = 30, tag.shift = 10, step = 10, rngl = rngl, background.density.scaling = FALSE)
	expect_equal(d$a$y, sig$a$y - bgw * ctl$a$y)
	d2 <- get.smoothed.tag.density(s, cc, bandwidth = 30, tag.shift = 10, step = 10, rngl = rngl, bg.weight = 0.25)
	expect_equal(d2$a$y, sig$a$y - 0.25 * ctl$a$y)
	d0 <- get.smoothed.tag.density(s, cc, bandwidth = 30, tag.shift = 10, step = 10, rngl = rngl, bg.weight = 0)
	expect_equal(d0$a$y, sig$a$y)
})

test_that("identical signal and control with weight 1 subtract to zero", {
	s <- pl.dens.data()
	d <- get.smoothed.tag.density(s, s, bandwidth = 30, tag.shift = 10, step = 10, background.density.scaling = FALSE)
	expect_true(all(abs(d$a$y) < 1e-12))
	expect_true(all(abs(d$b$y) < 1e-12))
})

test_that("scale.by.dataset.size scales to counts per million", {
	s <- pl.dens.data()
	plain <- get.smoothed.tag.density(s, bandwidth = 30, tag.shift = 10, step = 10, background.density.scaling = FALSE)
	scaled <- get.smoothed.tag.density(s, bandwidth = 30, tag.shift = 10, step = 10, background.density.scaling = FALSE, scale.by.dataset.size = TRUE)
	expect_equal(scaled$a$y, plain$a$y * 1e6 / 64)
	expect_equal(scaled$b$y, plain$b$y * 1e6 / 64)
})

test_that("a chromosome absent from the control is left unsubtracted", {
	s <- pl.dens.data()
	cc <- pl.dens.ctl()["a"]
	d <- get.smoothed.tag.density(s, cc, bandwidth = 30, tag.shift = 10, step = 10, bg.weight = 1)
	plain <- get.smoothed.tag.density(s, bandwidth = 30, tag.shift = 10, step = 10)
	expect_equal(d$b$y, plain$b$y)
})

test_that("get.smoothed.enrichment.mle is the log2 pseudocounted ratio of the two densities", {
	s <- pl.dens.data()
	cc <- pl.dens.ctl()
	rngl <- list(a = c(300, 1800), b = c(520, 1200))
	m <- get.smoothed.enrichment.mle(s, cc, bandwidth = 30, step = 1, rngl = rngl)
	expect_equal(names(m), c("a", "b"))
	ssd <- pl.sd(s, rngl)
	csd <- pl.sd(cc, rngl)
	bgw <- 64 / 34
	for (chr in c("a", "b")) {
		expect_equal(m[[chr]]$x, ssd[[chr]]$x)
		expect_equal(m[[chr]]$y, log2(ssd[[chr]]$y + bgw) - log2(csd[[chr]]$y + 1) - log2(bgw))
	}
	# bandwidth reaches both densities: a 10% different bandwidth changes the result clearly
	wide <- get.smoothed.enrichment.mle(s, cc, bandwidth = 33, step = 1, rngl = rngl)
	expect_gt(max(abs(wide$a$y - m$a$y)), 0.05)
	m2 <- get.smoothed.enrichment.mle(s, cc, bandwidth = 30, step = 1, rngl = rngl, bg.weight = 2, pseudocount = 0.5)
	expect_equal(m2$a$y, log2(ssd$a$y + 0.5 * 2) - log2(csd$a$y + 0.5) - log2(2))
})

test_that("identical signal and control give an enrichment of 1 (log2 0)", {
	s <- pl.dens.data()
	m <- get.smoothed.enrichment.mle(s, s, bandwidth = 30, step = 10)
	expect_equal(names(m), c("a", "b"))
	for (chr in names(m)) expect_true(all(abs(m[[chr]]$y) < 1e-12))
})

test_that("enrichment is highest over a signal pile and flat elsewhere", {
	flat <- as.integer(seq(5, 2000, by = 5))
	pile <- as.integer(round(seq(1000, 1100, length.out = 150)))
	sig <- list(a = sort(c(flat, pile)))
	ctl <- list(a = as.integer(seq(3, 2000, by = 5)))
	m <- get.smoothed.enrichment.mle(sig, ctl, bandwidth = 30, step = 10)
	best <- m$a$x[which.max(m$a$y)]
	# plus-strand tags are moved up by the default 73 bp shift
	expect_gte(best, 1073)
	expect_lte(best, 1173)
	expect_gt(max(m$a$y), 1)
	# away from the pile both samples are flat with 1 tag per 5 bp:
	# density D = bw * sqrt(2 pi) / 5, so log2((D + w) / (D + 1) / w) with w = n_sig / n_ctl
	w <- length(sig$a) / length(ctl$a)
	D <- 30 * sqrt(2 * pi) / 5
	expected <- log2((D + w) / (D + 1) / w)
	far <- m$a$x > 300 & m$a$x < 800
	expect_true(sum(far) > 20)
	expect_true(pl.close(m$a$y[far], rep(expected, sum(far)), 0.1))
	expect_equal(m$a$x[1], min(abs(sig$a + 73), abs(ctl$a + 73)))
})

test_that("chromosomes missing from the control are not reported", {
	s <- pl.dens.data()
	cc <- pl.dens.ctl()["a"]
	m <- get.smoothed.enrichment.mle(s, cc, bandwidth = 30, step = 10)
	expect_equal(names(m), "a")
	m2 <- get.smoothed.enrichment.mle(s, pl.dens.ctl(), bandwidth = 30, step = 10, chrl = c(a = "a"))
	expect_equal(names(m2), "a")
})

# pure-R reference of the experiment-corrected log2 ratio, as a function of the
# four smoothed densities and the two signal/control weights
pl.mle2.ref <- function(ss1, cs1, ss2, cs2, b1, b2, pc = 1) {
	log2(ss1 + pc * b1) - log2(cs1 + pc) - log2(b1) - (log2(ss2 + pc * b2) - log2(cs2 + pc) - log2(b2))
}

test_that("get.smoothed.enrichment.mle2 follows the two-experiment formula", {
	s1 <- pl.dens.data()
	c1 <- pl.dens.ctl()
	s2 <- list(a = pl.dens.tags(50, 31), b = c(450L, 750L, 950L, 1250L, 1300L))
	c2 <- list(a = pl.dens.tags(40, 32), b = c(550L, -750L, 1100L, 1400L, 1450L))
	rngl <- list(a = c(300, 1800), b = c(520, 1200))
	m <- get.smoothed.enrichment.mle2(s1, c1, s2, c2, bandwidth = 30, step = 1, rngl = rngl)
	expect_equal(names(m), c("a", "b"))
	b1 <- 64 / 34
	b2 <- 55 / 45
	for (chr in c("a", "b")) {
		d <- list(ss1 = pl.sd(s1, rngl)[[chr]], cs1 = pl.sd(c1, rngl)[[chr]], ss2 = pl.sd(s2, rngl)[[chr]], cs2 = pl.sd(c2, rngl)[[chr]])
		expect_equal(m[[chr]]$x, d$ss1$x)
		expect_equal(m[[chr]]$y, pl.mle2.ref(d$ss1$y, d$cs1$y, d$ss2$y, d$cs2$y, b1, b2))
	}
	# the experiments are not interchangeable: swapping them negates the ratio
	sw <- get.smoothed.enrichment.mle2(s2, c2, s1, c1, bandwidth = 30, step = 1, rngl = rngl)
	expect_equal(sw$a$y, -m$a$y)
	expect_gt(max(abs(m$a$y)), 0.05)
	# explicit weights and pseudocount
	mw <- get.smoothed.enrichment.mle2(s1, c1, s2, c2, bandwidth = 30, step = 1, rngl = rngl, bg.weight1 = 1, bg.weight2 = 2, pseudocount = 0.5)
	d <- list(ss1 = pl.sd(s1, rngl)$a$y, cs1 = pl.sd(c1, rngl)$a$y, ss2 = pl.sd(s2, rngl)$a$y, cs2 = pl.sd(c2, rngl)$a$y)
	expect_equal(mw$a$y, pl.mle2.ref(d$ss1, d$cs1, d$ss2, d$cs2, 1, 2, 0.5))
})

test_that("get.smoothed.enrichment.mle2 of an experiment against itself is zero", {
	s <- pl.dens.data()
	cc <- pl.dens.ctl()
	m <- get.smoothed.enrichment.mle2(s, cc, s, cc, bandwidth = 30, step = 10)
	expect_equal(names(m), c("a", "b"))
	for (chr in names(m)) expect_true(all(abs(m[[chr]]$y) < 1e-12))
	m2 <- get.smoothed.enrichment.mle2(s, cc, s, cc, bandwidth = 30, step = 10, bg.weight1 = 0.5, bg.weight2 = 0.5, pseudocount = 2)
	for (chr in names(m2)) expect_true(all(abs(m2[[chr]]$y) < 1e-12))
})

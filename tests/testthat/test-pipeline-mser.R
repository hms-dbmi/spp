# get.mser and get.mser.interpolation.
# A saturation "chain" is a list of find.binding.positions-like results, one per
# subsampling step, each with $n (tags used) and $npl (peaks with $x and $enr).
# The minimal saturated enrichment (MSER) between a step and the next one is the
# lowest enrichment such that the peaks above it overlap (within 50 bp) with
# peaks of the next step at the test.agreement rate.

pl.step <- function(n, x, enr) {
	list(n = n, npl = list(chr1 = data.frame(x = x, enr = enr)))
}

# Hand-derived chain.  Step 0 has peaks (100: 8), (1000: 2); step 1 has (100: 6),
# (5000: 3); step 2 has (100: 4), (9000: 2.5).
# 0 -> 1: peak 100 agrees, peak 1000 does not.  Sorted by min(ref, other) enrichment
#   the agreement is 1/1 at 6 and 1/2 below, so the MSER at agreement .99 is 6
#   (the floor max(min enr of both steps, 1) is 3).
# 1 -> 2: peaks (100: 6, other 4), (5000: 3, none) give 4 in the same way.
pl.chain.a <- function() {
	list("1" = list(
		"0" = pl.step(1000, c(100, 1000), c(8, 2)),
		"1" = pl.step(800, c(100, 5000), c(6, 3)),
		"2" = pl.step(600, c(100, 9000), c(4, 2.5))
	))
}

test_that("get.mser on supplied chains matches the hand-derived values", {
	m <- get.mser(NULL, NULL, chains = pl.chain.a(), n.steps = 3)
	expect_equal(names(m), "1")
	expect_equal(m[["1"]]$n, c(1000, 800))
	expect_equal(m[["1"]]$me, c(6, 4))
})

test_that("get.mser with a single step returns the mser vector", {
	chain <- list("1" = list("0" = pl.step(1000, c(100, 1000), c(8, 2)), "1" = pl.step(800, c(100, 5000), c(6, 3))))
	m <- get.mser(NULL, NULL, chains = chain, n.steps = 1)
	expect_equal(unname(m), 6)
})

test_that("test.agreement lowers or raises the MSER", {
	# peaks (100: 10), (200: 8), (300: 5), (400: 3); the next step reproduces
	# 100 (9) and 200 (7.5) and adds the unrelated peak 500 (4).
	chain <- list("1" = list(
		"0" = pl.step(1000, c(100, 200, 300, 400), c(10, 8, 5, 3)),
		"1" = pl.step(800, c(110, 210, 500), c(9, 7.5, 4))
	))
	# at 0.99 only the top two peaks (agreement 2/2) qualify: MSER 7.5
	strict <- get.mser(NULL, NULL, chains = chain, n.steps = 1, test.agreement = 0.99)
	expect_equal(unname(strict), 7.5)
	# at 0.5 everything qualifies (2/4), so the MSER falls to 3, below the floor
	# max(min enr) = 4, and is reported as 1 (saturated)
	loose <- get.mser(NULL, NULL, chains = chain, n.steps = 1, test.agreement = 0.5)
	expect_equal(unname(loose), 1)
})

test_that("get.mser averages over chains", {
	ref <- pl.step(1000, c(100, 200, 300, 400), c(10, 8, 5, 3))
	# chain 1 agrees on 100, 200 (MSER 7.5, floor 4); chain 2 also on 300 (MSER 4.5, floor 4)
	chain1 <- list("0" = ref, "1" = pl.step(800, c(110, 210, 500), c(9, 7.5, 4)))
	chain2 <- list("0" = ref, "1" = pl.step(800, c(110, 210, 310, 500), c(9, 7.5, 4.5, 4)))
	m <- get.mser(NULL, NULL, chains = list("1" = chain1, "2" = chain2), n.steps = 2)
	expect_equal(m[["1"]]$n, 1000)
	expect_equal(m[["1"]]$me, mean(c(7.5, 4.5)))
	r <- get.mser(NULL, NULL, chains = list("1" = chain1, "2" = chain2), n.steps = 2, return.chains = TRUE)
	expect_equal(names(r), c("mser", "chains"))
	expect_equal(length(r$chains), 2)
})

test_that("a fully reproduced peak set is reported as saturated (MSER 1)", {
	same <- pl.step(1000, c(100, 200), c(8, 5))
	chain <- list("1" = list("0" = same, "1" = pl.step(800, c(102, 198), c(8, 5))))
	m <- get.mser(NULL, NULL, chains = chain, n.steps = 1)
	expect_equal(unname(m), 1)
})

test_that("get.mser.interpolation fits log10(mser - 1) against log10(n)", {
	sig <- list(chr1 = 1:1000)
	res <- get.mser.interpolation(sig, NULL, chains = pl.chain.a(), n.steps = 3, excluded.steps = c(), target.fold.enrichment = 5)
	expect_equal(names(res), "1")
	# two points: (n = 1000, mser 6) and (n = 800, mser 4); the line in log10 space
	# is solved for mser = 5 (target fold enrichment)
	x <- log10(c(1000, 800))
	y <- log10(c(6, 4) - 1)
	slope <- diff(y) / diff(x)
	tfe <- x[1] + (log10(5 - 1) - y[1]) / slope
	expect_equal(res[["1"]]$prediction, 10^tfe)
	expect_equal(unname(coef(res[["1"]]$log10.fit)), c(y[1] - slope * x[1], slope))
	# tag counts are rescaled by the ratio of the signal size to the first step
	sig2 <- list(chr1 = 1:2000)
	res2 <- get.mser.interpolation(sig2, NULL, chains = pl.chain.a(), n.steps = 3, excluded.steps = c(), target.fold.enrichment = 5)
	expect_equal(res2[["1"]]$prediction, 2 * res[["1"]]$prediction)
})

test_that("get.mser.interpolation predicts NA when a step is fully saturated", {
	sat <- pl.step(1000, c(100, 200), c(8, 5))
	chain <- list("1" = list("0" = sat, "1" = pl.step(800, c(102, 198), c(8, 5)), "2" = pl.step(600, c(101, 199), c(8, 5))))
	# cause: every peak is reproduced at every step, so each MSER is 1 (saturated) ...
	m <- get.mser(NULL, NULL, chains = chain, n.steps = 3)
	expect_equal(m[["1"]]$me, c(1, 1))
	# ... and a log10(mser - 1) fit is impossible, so the prediction is NA without a fit
	res <- get.mser.interpolation(list(chr1 = 1:1000), NULL, chains = chain, n.steps = 3, excluded.steps = c())
	expect_true(is.na(res[["1"]]$prediction))
	expect_null(res[["1"]]$log10.fit)
})

test_that("get.mser.interpolation can return the chains", {
	res <- get.mser.interpolation(list(chr1 = 1:1000), NULL, chains = pl.chain.a(), n.steps = 3, excluded.steps = c(), return.chains = TRUE)
	expect_equal(names(res), c("interpolation", "chains"))
	expect_equal(length(res$chains), 1)
})

test_that("get.mser runs end to end on subsampled synthetic data", {
	cs <- seq(2000, 20000, by = 3000)
	sig <- list(chr1 = pl.tags(cs, n = 30, spread = 10), chr2 = pl.tags(c(4000, 9000), n = 30, spread = 10))
	# background kept 160 bp away from every pile so that it does not move the peaks
	clear <- function(x, centers) x[vapply(abs(x), function(p) all(abs(p - centers) > 160), logical(1))]
	ctl <- list(chr1 = clear(pl.bg(300, L = 22000, seed = 1), cs), chr2 = clear(pl.bg(200, L = 12000, seed = 2), c(4000, 9000)))
	centers <- list(chr1 = cs, chr2 = c(4000, 9000))
	# debug, min.ntags, excluded.steps and test.chromosomes are named so that the chain
	# index get.mser passes positionally cannot land on one of them (see BUGS-pipeline.txt);
	# tec.filter is off so that background masking cannot hide a pile
	res <- pl.seeded(1, pl.quiet(get.mser(sig, ctl, n.chains = 2, step.size = 200, n.steps = 1, whs = 100, e.value = 1e6, return.chains = TRUE, debug = FALSE, min.ntags = 10, excluded.steps = c(), test.chromosomes = NULL, tec.filter = FALSE)))
	expect_equal(length(res$chains), 2)
	# step 0 uses all 540 tags; step 1 keeps floor(n * (540 - 200) / 540) tags per chromosome
	n1 <- floor(420 * 340 / 540) + floor(120 * 340 / 540)
	expect_equal(n1, 339)
	for (ch in res$chains) {
		expect_equal(unname(vapply(ch, function(d) d$n, numeric(1))), c(540, n1))
		# every pile is called at step 0 and every called peak has a partner (within 50 bp)
		# at step 1, so the peaks are saturated and the MSER is 1
		for (chr in c("chr1", "chr2")) {
			x0 <- ch[["0"]]$npl[[chr]]$x
			x1 <- ch[["1"]]$npl[[chr]]$x
			expect_equal(length(x0), length(centers[[chr]]))
			expect_true(all(abs(sort(x0) - centers[[chr]]) <= 5))
			expect_true(all(vapply(x0, function(p) any(abs(p - x1) <= 50), logical(1))))
		}
	}
	expect_equal(unname(res$mser), 1)
})

# wtd() / lwcc() and their wrappers tag.wtd() / tag.lwcc() (spp_wtd and
# spp_lwcc in wdl.cpp).
#
# wtd() is exercised through its masked code path (a mask position outside the
# scanned range makes ignore.masking FALSE).  The unmasked fast path
# disagrees with the masked path; see BUGS-kernels.txt.  A mask that does not
# fall in the scanned range changes nothing else.

wtd <- spp_internal("wtd")
lwcc <- spp_internal("lwcc")
tag.wtd <- spp_internal("tag.wtd")
tag.lwcc <- spp_internal("tag.lwcc")

far <- -5000 # mask position outside every scanned range

wtd_m <- function(x, y, s, e, ...) {
	quiet(wtd(x, y, s, e, mask.x = far, mask.y = far, ...))
}
lwcc_b <- function(x, y, s, e, ...) {
	quiet(lwcc(x, y, s, e, bg.weight = 1, ...))
}

# random tag sets shared by the reference comparisons
set.seed(2)
rx_ <- sample(100:300, 60, replace = TRUE)
ry_ <- sample(100:300, 60, replace = TRUE)
rbx <- sample(100:300, 80, replace = TRUE)
rby <- sample(100:300, 80, replace = TRUE)

test_that("wtd returns the hand-computed scores for one tag per strand", {
	# + tag at 10, - tag at 11, whs = 3, scanned 5..15 -> positions 2..18
	r <- wtd_m(10L, 11L, 5, 15, whs = 3, return.peaks = FALSE)
	expect_equal(r$x, c(2, 18))
	expect_equal(r$y, c(rep(0, 4), -1, -1, -2, -2, -2, 1, -2, -2, -2, rep(0, 4)))
})

test_that("wtd round.up and tag.weight follow hand-computed scores", {
	# score at position 11 (index 10): + tag at 10 in the left window, - tag at 11
	# in the right window, tag.weight 0.5.
	# without round.up: 2 * sqrt(0.5 * 0.5) - 0.5 = 0.5
	# with round.up the fractional counts 0.5 are raised to 1: 2 * 1 - 0.5 = 1.5
	a <- wtd_m(10L, 11L, 5, 15, whs = 3, return.peaks = FALSE, tag.weight = 0.5)
	b <- wtd_m(10L, 11L, 5, 15, whs = 3, return.peaks = FALSE, tag.weight = 0.5,
		round.up = TRUE)
	expect_equal(a$y[10], 0.5)
	expect_equal(b$y[10], 1.5)
})

test_that("wtd background subtraction and the bg.whs rescale follow hand-computed scores", {
	# two + and two - background tags at the same positions as the signal tags.
	# bg.weight 0.25: signal 1 - 2 * 0.25 = 0.5 on each strand, score
	# 2 * sqrt(0.5 * 0.5) - 1 = 0; bg.weight 0.5: the + signal 1 - 2 * 0.5 = 0 is not
	# positive, score = -1.  bg.whs = 6 with bg.weight 0.5 is rescaled by
	# whs / bg.whs = 0.5 to 0.25 and gives 0 again.
	bgx <- c(10L, 10L); bgy <- c(11L, 11L)
	g <- function(...) wtd_m(10L, 11L, 5, 15, whs = 3, return.peaks = FALSE,
		bg.x = bgx, bg.y = bgy, ...)$y[10]
	expect_equal(g(bg.weight = 0.25), 0)
	expect_equal(g(bg.weight = 0.5), -1)
	expect_equal(g(bg.weight = 0.5, bg.whs = 6), 0)
})

test_that("wtd masking skips every window pair that touches a masked position", {
	r <- quiet(wtd(10L, 11L, 5, 15, whs = 3, return.peaks = FALSE,
		mask.x = 11, mask.y = far))
	expect_equal(r$y, c(rep(0, 4), -1, -1, -2, -2, -2, -1, -1, -1, -1, rep(0, 4)))
})

test_that("a mask on a position that holds a tag is ignored", {
	a <- quiet(wtd(10L, 11L, 5, 15, whs = 3, return.peaks = FALSE,
		mask.x = 10, mask.y = 11))
	b <- wtd_m(10L, 11L, 5, 15, whs = 3, return.peaks = FALSE)
	expect_equal(a$y, b$y)
})

test_that("wtd matches the window-sum reference without background", {
	s <- 100; e <- 300; whs <- 10
	rx <- c(s - whs, e + whs)
	hx <- ref_hist(rx_, rx); hy <- ref_hist(ry_, rx)
	r <- wtd_m(rx_, ry_, s, e, whs = whs, return.peaks = FALSE)
	expect_equal(r$x, rx)
	expect_length(r$y, diff(rx) + 1)
	expect_equal(r$y, ref_wtd(hx, hy, whs), tolerance = 1e-10)
})

test_that("wtd tag.weight, direct.count and round.up follow the formulas", {
	s <- 100; e <- 300; whs <- 10
	rx <- c(s - whs, e + whs)
	hx <- ref_hist(rx_, rx); hy <- ref_hist(ry_, rx)
	expect_equal(
		wtd_m(rx_, ry_, s, e, whs = whs, return.peaks = FALSE, tag.weight = 0.5)$y,
		ref_wtd(hx, hy, whs, tag.weight = 0.5), tolerance = 1e-10)
	expect_equal(
		wtd_m(rx_, ry_, s, e, whs = whs, return.peaks = FALSE, direct.count = TRUE)$y,
		ref_wtd(hx, hy, whs, direct.count = TRUE), tolerance = 1e-10)
	expect_equal(
		wtd_m(rx_, ry_, s, e, whs = whs, return.peaks = FALSE, tag.weight = 0.1,
			round.up = TRUE)$y,
		ref_wtd(hx, hy, whs, tag.weight = 0.1, round.up = TRUE), tolerance = 1e-10)
	expect_equal(
		wtd_m(rx_, ry_, s, e, whs = whs, return.peaks = FALSE, direct.count = TRUE,
			tag.weight = 0.02, round.up = TRUE)$y,
		ref_wtd(hx, hy, whs, tag.weight = 0.02, direct.count = TRUE, round.up = TRUE),
		tolerance = 1e-10)
})

test_that("wtd background subtraction follows the formulas", {
	s <- 100; e <- 300; whs <- 10
	rx <- c(s - whs, e + whs)
	hx <- ref_hist(rx_, rx); hy <- ref_hist(ry_, rx)
	bg <- list(x = ref_hist(rbx, rx), y = ref_hist(rby, rx))
	expect_equal(
		wtd_m(rx_, ry_, s, e, whs = whs, return.peaks = FALSE, bg.x = rbx,
			bg.y = rby, bg.weight = 0.3)$y,
		ref_wtd(hx, hy, whs, bg = bg, bg.weight = 0.3), tolerance = 1e-10)
	# a larger background window rescales bg.weight by whs / bg.whs
	expect_equal(
		wtd_m(rx_, ry_, s, e, whs = whs, return.peaks = FALSE, bg.x = rbx,
			bg.y = rby, bg.weight = 0.3, bg.whs = 20)$y,
		ref_wtd(hx, hy, whs, bg = bg, bg.whs = 20, bg.weight = 0.15),
		tolerance = 1e-10)
	expect_equal(
		wtd_m(rx_, ry_, s, e, whs = whs, return.peaks = FALSE, bg.x = rbx,
			bg.y = rby, bg.weight = 0.3, direct.count = TRUE, tag.weight = 2)$y,
		ref_wtd(hx, hy, whs, tag.weight = 2, direct.count = TRUE, bg = bg,
			bg.weight = 0.6), tolerance = 1e-10)
})

test_that("wtd masking matches the reference with several masked positions", {
	s <- 100; e <- 300; whs <- 10
	rx <- c(s - whs, e + whs)
	mx <- c(120, 150, 151); my <- c(130, 200)
	hx <- ref_hist(rx_, rx); hy <- ref_hist(ry_, rx)
	hx[setdiff(mx, rx_) - rx[1] + 1] <- -1
	hy[setdiff(my, ry_) - rx[1] + 1] <- -1
	r <- quiet(wtd(rx_, ry_, s, e, whs = whs, return.peaks = FALSE,
		mask.x = mx, mask.y = my))
	expect_equal(r$y, ref_wtd(hx, hy, whs, masked = TRUE), tolerance = 1e-10)
})

test_that("wtd peaks are the plateau centres of two strand-separated clusters", {
	# + tags 100..110 (40 of them), - tags 125..135; a plateau of maximal score
	# for window centres 116..120 (left window holds all + tags, right window
	# holds all - tags) with score 2 * sqrt(40 * 40) - 1.
	x <- c(rep(100:110, length.out = 40), rep(400:410, length.out = 30))
	y <- c(rep(125:135, length.out = 40), rep(425:435, length.out = 30))
	p <- wtd_m(x, y, 50, 500, whs = 20, min.thr = 1, min.dist = 30)
	expect_equal(p$x, c(118, 418))
	expect_equal(p$y, c(79, 59))
	expect_named(p, c("x", "y"))
})

test_that("wtd peak of an even-width plateau is its middle", {
	# - tags 125..134 instead of 125..135: the maximal plateau is centres 115..120
	x <- rep(100:110, length.out = 40)
	y <- rep(125:134, length.out = 40)
	full <- wtd_m(x, y, 50, 300, whs = 20, return.peaks = FALSE)
	top <- which(full$y == max(full$y)) + full$x[1] - 1
	expect_equal(top, 115:120)
	p <- wtd_m(x, y, 50, 300, whs = 20, min.thr = 1, min.dist = 30)
	expect_equal(p$y, max(full$y))
	expect_equal(max(full$y), 79)
	# either of the two middle positions is acceptable
	expect_true(p$x %in% c(117, 118))
})

test_that("wtd min.thr and min.dist filter peaks", {
	x <- c(rep(100:110, length.out = 40), rep(400:410, length.out = 30))
	y <- c(rep(125:135, length.out = 40), rep(425:435, length.out = 30))
	none <- wtd_m(x, y, 50, 500, whs = 20, min.thr = 1e6, min.dist = 30)
	expect_equal(nrow(none), 0)
	only_best <- wtd_m(x, y, 50, 500, whs = 20, min.thr = 1, min.dist = 1000)
	expect_equal(only_best$x, 118)
	expect_equal(only_best$y, 79)
	# a threshold between the two scores keeps only the stronger peak
	thr <- wtd_m(x, y, 50, 500, whs = 20, min.thr = 70, min.dist = 30)
	expect_equal(thr$x, 118)
})

test_that("wtd step bins coordinates and returns them on the original scale", {
	x <- c(rep(100:110, length.out = 40), rep(400:410, length.out = 30))
	y <- c(rep(125:135, length.out = 40), rep(425:435, length.out = 30))
	p <- wtd_m(x, y, 50, 500, whs = 20, min.thr = 1, min.dist = 30, step = 5)
	# in bins of 5 the + tags fall in bins 20..22 and the - tags in 25..27, so the
	# only window centre that holds all of them is bin 24 (= 120)
	expect_equal(p$x, c(120, 420))
	expect_equal(p$y, c(79, 59))
	full <- wtd_m(x, y, 50, 500, whs = 20, return.peaks = FALSE, step = 5)
	# range: s / step - whs / step = 10 - 4 and e / step + whs / step = 100 + 4, in bp
	expect_equal(full$x, c(30, 520))
})

test_that("tag.wtd splits signed tags by strand and trims to [s, e]", {
	ctv <- c(10, 10, -11, 12, -30, 40)
	a <- quiet(tag.wtd(ctv, 5, 15, return.peaks = FALSE, whs = 3,
		mask.ctv = numeric(0)))
	b <- wtd_m(c(10L, 10L, 12L), 11L, 5, 15, whs = 3, return.peaks = FALSE)
	expect_equal(a$y, b$y)
	expect_equal(a$x, b$x)
})

test_that("tag.wtd returns an empty data frame when a strand has no tags", {
	r <- tag.wtd(c(10, 12), 5, 15)
	expect_equal(nrow(r), 0)
	r2 <- tag.wtd(c(-10, -12), 5, 15)
	expect_equal(nrow(r2), 0)
	r3 <- tag.wtd(c(100, -100), 5, 15)
	expect_equal(nrow(r3), 0)
})

test_that("tag.wtd with one strand empty and return.peaks = FALSE gives zeros", {
	r <- tag.wtd(c(10, 12), 5, 15, return.peaks = FALSE)
	expect_equal(r$x, c(10, 12))
	expect_equal(r$y, c(0, 0, 0))
})

test_that("tag.wtd passes bg.ctv and mask.ctv through", {
	set.seed(21)
	pos <- sample(100:200, 40, replace = TRUE)
	neg <- sample(100:200, 40, replace = TRUE)
	bgp <- sample(100:200, 30, replace = TRUE)
	bgn <- sample(100:200, 30, replace = TRUE)
	a <- quiet(tag.wtd(c(pos, -neg), 100, 200, return.peaks = FALSE, whs = 8,
		bg.ctv = c(bgp, -bgn), bg.weight = 0.4, mask.ctv = numeric(0)))
	b <- quiet(wtd(pos, neg, 100, 200, return.peaks = FALSE, whs = 8,
		bg.x = bgp, bg.y = bgn, bg.weight = 0.4, mask.x = far, mask.y = far))
	expect_equal(a$y, b$y)
})

# ---- lwcc ----

test_that("lwcc returns the hand-computed score for a lone tag pair", {
	# + tag at 8, - tag at 12, outer 2, inner 2: at position 10 the correlation
	# is 1 and the tag weight is 2 * sqrt(1 * 1) = 2, so score = sqrt(2) + 2;
	# every other scored position has no tags on one strand (-1000).
	r <- lwcc_b(8L, 12L, 5, 15, whs = 2, isize = 2, return.peaks = FALSE)
	expect_equal(r$x, c(3, 17))
	want <- c(0, 0, rep(-1000, 5), sqrt(2) + 2, rep(-1000, 4), 0, 0, 0)
	expect_equal(r$y, want, tolerance = 1e-10)
})

test_that("lwcc matches the correlation reference", {
	s <- 100; e <- 300
	rx <- c(s - 20, e + 20)
	hx <- ref_hist(rx_, rx); hy <- ref_hist(ry_, rx)
	r <- lwcc_b(rx_, ry_, s, e, whs = 20, isize = 5, return.peaks = FALSE)
	expect_equal(r$x, rx)
	expect_equal(r$y, ref_lwcc(hx, hy, 20, 5), tolerance = 1e-10)
	r2 <- lwcc_b(rx_, ry_, s, e, whs = 20, isize = 5, return.peaks = FALSE,
		tag.weight = 0.4)
	expect_equal(r2$y, ref_lwcc(hx, hy, 20, 5, tag.weight = 0.4), tolerance = 1e-10)
	r3 <- lwcc_b(rx_, ry_, s, e, whs = 20, isize = 5, return.peaks = FALSE,
		tag.weight = 0.02, round.up = TRUE)
	expect_equal(r3$y, ref_lwcc(hx, hy, 20, 5, tag.weight = 0.02, round.up = TRUE),
		tolerance = 1e-10)
})

test_that("lwcc background subtraction follows the formulas", {
	s <- 100; e <- 300
	rx <- c(s - 20, e + 20)
	hx <- ref_hist(rx_, rx); hy <- ref_hist(ry_, rx)
	bg <- list(x = ref_hist(rbx, rx), y = ref_hist(rby, rx))
	# bg.weight is scaled by (whs - isize) / bg.whs
	r <- quiet(lwcc(rx_, ry_, s, e, whs = 20, isize = 5, return.peaks = FALSE,
		bg.x = rbx, bg.y = rby, bg.weight = 0.3))
	expect_equal(r$y, ref_lwcc(hx, hy, 20, 5, bg = bg, bg.whs = 20,
		bg.weight = 0.3 * 15 / 20), tolerance = 1e-10)
	r2 <- quiet(lwcc(rx_, ry_, s, e, whs = 20, isize = 5, return.peaks = FALSE,
		bg.x = rbx, bg.y = rby, bg.weight = 0.3, bg.whs = 30))
	expect_equal(r2$y, ref_lwcc(hx, hy, 20, 5, bg = bg, bg.whs = 30,
		bg.weight = 0.3 * 15 / 30), tolerance = 1e-10)
})

test_that("lwcc masking matches the reference", {
	s <- 100; e <- 300
	rx <- c(s - 20, e + 20)
	mx <- c(120, 150, 151); my <- c(130, 200)
	hx <- ref_hist(rx_, rx); hy <- ref_hist(ry_, rx)
	hx[setdiff(mx, rx_) - rx[1] + 1] <- -1
	hy[setdiff(my, ry_) - rx[1] + 1] <- -1
	r <- lwcc_b(rx_, ry_, s, e, whs = 20, isize = 5, return.peaks = FALSE,
		mask.x = mx, mask.y = my)
	expect_equal(r$y, ref_lwcc(hx, hy, 20, 5, masked = TRUE), tolerance = 1e-10)
})

test_that("lwcc peaks sit at the centre between + and - tag clusters", {
	x <- c(rep(100:110, length.out = 40), rep(400:410, length.out = 30))
	y <- c(rep(125:135, length.out = 40), rep(425:435, length.out = 30))
	full <- lwcc_b(x, y, 50, 500, whs = 40, isize = 5, return.peaks = FALSE)
	p <- lwcc_b(x, y, 50, 500, whs = 40, isize = 5, min.thr = 1, min.dist = 30)
	# midpoint between the + cluster (centre 105) and the - cluster (centre 130)
	expect_length(p$x, 2)
	expect_true(all(abs(p$x - c(117.5, 417.5)) <= 1.5))
	# each peak score is the score vector at that position
	expect_equal(p$y, full$y[p$x - full$x[1] + 1], tolerance = 1e-12)
	# and is the highest score within 30 bp of the peak
	for (k in 1:2) {
		win <- (p$x[k] - 30):(p$x[k] + 30) - full$x[1] + 1
		expect_equal(max(full$y[win]), p$y[k], tolerance = 1e-12)
	}
	expect_equal(nrow(lwcc_b(x, y, 50, 500, whs = 40, isize = 5, min.thr = 1e6,
		min.dist = 30)), 0)
	best <- lwcc_b(x, y, 50, 500, whs = 40, isize = 5, min.thr = 1, min.dist = 1000)
	expect_length(best$x, 1)
	expect_true(abs(best$x - 117.5) <= 1.5)
})

test_that("lwcc step bins coordinates", {
	x <- c(rep(100:110, length.out = 40), rep(400:410, length.out = 30))
	y <- c(rep(125:135, length.out = 40), rep(425:435, length.out = 30))
	p <- lwcc_b(x, y, 50, 500, whs = 40, isize = 5, min.thr = 1, min.dist = 30,
		step = 5)
	full <- lwcc_b(x, y, 50, 500, whs = 40, isize = 5, return.peaks = FALSE, step = 5)
	# range: (s - whs) / step = 10 - 8 and (e + whs) / step = 100 + 8, in bp
	expect_equal(full$x, c(10, 540))
	expect_length(full$y, 107)
	expect_length(p$x, 2)
	expect_true(all(p$x %% 5 == 0))
	# bins of 5 bp: + tags in bins 20..22, - tags in bins 25..27, mid-point 23.5
	expect_true(all(abs(p$x - c(117.5, 417.5)) <= 5))
	expect_equal(p$y, full$y[(p$x - full$x[1]) / 5 + 1], tolerance = 1e-12)
})

test_that("tag.lwcc splits signed tags by strand and trims to [s, e]", {
	ctv <- c(8, 8, -12, 30, -40)
	a <- quiet(tag.lwcc(ctv, 5, 15, return.peaks = FALSE, whs = 2, isize = 2,
		bg.weight = 1))
	b <- lwcc_b(c(8L, 8L), 12L, 5, 15, whs = 2, isize = 2, return.peaks = FALSE)
	expect_equal(a$y, b$y)
	expect_equal(a$x, b$x)
})

test_that("tag.lwcc returns an empty data frame when a strand has no tags", {
	expect_equal(nrow(tag.lwcc(c(8, 9), 5, 15)), 0)
	expect_equal(nrow(tag.lwcc(c(-8, -9), 5, 15)), 0)
	r <- tag.lwcc(c(8, 9), 5, 15, return.peaks = FALSE)
	expect_equal(r$x, c(8, 9))
	expect_equal(r$y, c(0, 0))
})

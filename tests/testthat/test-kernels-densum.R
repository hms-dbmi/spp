# densum(): Gaussian-kernel density of tag positions (ccdensum in cdensum.c).
# The last element of the output is never written by the native code (see
# BUGS-kernels.txt), so comparisons against the brute-force reference leave it out.

densum <- spp_internal("densum")

# Comparison against the brute-force reference over all but the last element.
expect_density_equal <- function(vin, bw, from, to, step = 1, ...) {
	got <- densum(vin, bw = bw, dw = 10L, from = from, to = to, step = step,
		return.x = FALSE, ...)
	want <- ref_density(vin, bw, from, to, step)
	n <- length(want)
	expect_length(got, n)
	expect_equal(got[-n], want[-n], tolerance = 1e-10)
}

test_that("a single tag gives the hand-computed Gaussian profile", {
	# tag at 10, window 10..14, bw = 2: y[j] = exp(-0.5 * (j / 2)^2)
	y <- densum(10, bw = 2, dw = 3, from = 10, to = 14, return.x = FALSE)
	expect_length(y, 5)
	expect_equal(y[1:4], c(1, exp(-1 / 8), exp(-0.5), exp(-9 / 8)),
		tolerance = 1e-10)
})

test_that("return.x gives the x range, the step and the profile", {
	r <- densum(c(10, 20), bw = 3, from = 0, to = 30, return.x = TRUE)
	expect_named(r, c("x", "y", "step"))
	expect_equal(r$x, c(0, 30))
	expect_equal(r$step, 1)
	expect_length(r$y, 31)
	expect_equal(r$y, densum(c(10, 20), bw = 3, from = 0, to = 30, return.x = FALSE))
})

test_that("density matches a brute-force sum for several bandwidths", {
	vin <- c(12, 12, 15, 30, 31, 31, 31, 55, 70)
	for (bw in c(1, 2, 5, 8)) {
		expect_density_equal(vin, bw, from = 0, to = 90)
	}
})

test_that("duplicated positions are counted as multiplicities", {
	vin <- c(40, 40, 40, 60)
	expect_density_equal(vin, 3, from = 20, to = 80)
	y <- densum(vin, bw = 3, dw = 10L, from = 20, to = 80, return.x = FALSE)
	# at the triple position the three tags dominate the single tag 20 bp away
	expect_equal(y[21], 3 + exp(-0.5 * (20 / 3)^2), tolerance = 1e-10)
})

test_that("from and to restrict which tags contribute", {
	vin <- c(5, 10, 50, 90)
	expect_density_equal(vin, 2, from = 8, to = 60)
	# the tag at 90 and the tag at 5 are outside the range
	y <- densum(vin, bw = 2, dw = 10L, from = 8, to = 60, return.x = FALSE)
	expect_equal(y, densum(c(10, 50), bw = 2, dw = 10L, from = 8, to = 60,
		return.x = FALSE))
})

test_that("step subsamples the output grid", {
	vin <- c(10, 20, 27)
	for (step in c(2L, 5L, 7L)) {
		expect_density_equal(vin, 4, from = 0, to = 40, step = step)
		r <- densum(vin, bw = 4, dw = 10L, from = 0, to = 40, step = step)
		n <- floor(40 / step) + 1
		expect_length(r$y, n)
		expect_equal(r$x, c(0, step * (n - 1)))
		expect_equal(r$step, step)
	}
})

test_that("integer and numeric inputs give identical output", {
	vin <- c(10L, 12L, 12L, 25L)
	a <- densum(vin, bw = 3L, dw = 4L, from = 0L, to = 40L, step = 2L)
	b <- densum(as.numeric(vin), bw = 3, dw = 4, from = 0, to = 40, step = 2)
	expect_identical(a, b)
})

test_that("match.wt.f weights each position's count", {
	vin <- c(20, 20, 50)
	wf <- function(p) ifelse(p == 20, 0.5, 2)
	got <- densum(vin, bw = 3, dw = 10L, from = 0, to = 70, match.wt.f = wf,
		return.x = FALSE)
	want <- ref_density(vin, 3, 0, 70, 1, wt = wf)
	n <- length(want)
	expect_equal(got[-n], want[-n], tolerance = 1e-10)
})

test_that("the dw window cuts off contributions beyond dw * bw", {
	# whole-vector comparison against the windowed reference (last element is
	# never written, see BUGS-kernels.txt)
	for (cfg in list(
		list(v = 20, bw = 2, dw = 1, from = 0, to = 40, step = 1),
		list(v = c(10, 10, 10, 30), bw = 2, dw = 1, from = 0, to = 50, step = 1),
		list(v = c(10, 10, 31), bw = 3, dw = 2, from = 0, to = 60, step = 3),
		list(v = c(25, 25), bw = 2, dw = 2, from = 5, to = 60, step = 2))) {
		got <- do.call(densum, c(list(cfg$v, bw = cfg$bw, dw = cfg$dw, from = cfg$from,
			to = cfg$to, step = cfg$step, return.x = FALSE)))
		want <- ref_density_win(cfg$v, cfg$bw, cfg$dw, cfg$from, cfg$to, cfg$step)
		n <- length(want)
		expect_length(got, n)
		expect_equal(got[-n], want[-n], tolerance = 1e-10)
	}
})

test_that("a tag with multiplicity ic widens its window to dw * bw * ic", {
	one <- densum(20, bw = 2, dw = 1, from = 0, to = 60, return.x = FALSE)
	three <- densum(c(20, 20, 20), bw = 2, dw = 1, from = 0, to = 60, return.x = FALSE)
	# single tag: window half size 2; three tags: half size 6
	expect_true(all(one[1:17] == 0))
	expect_true(any(three[15:17] > 0))
	expect_true(all(three[1:13] == 0))
	expect_equal(three[21], 3)
})

test_that("the last output element is non-negative", {
	y <- densum(c(10, 20), bw = 5, dw = 3, from = 0, to = 30, return.x = FALSE)
	expect_gte(y[length(y)], 0)
})

test_that("empty ranges and bad ranges raise errors", {
	expect_error(densum(c(1, 2), from = -1), "negative")
	expect_error(suppressWarnings(densum(c(1, 2), from = 5, to = 3)),
		"zero data range")
	expect_error(densum(c(1, 2), new.code = FALSE), "new.code")
	expect_error(suppressWarnings(densum(numeric(0))), "missing value|TRUE/FALSE")
})

test_that("tags outside the range are ignored, including negative ones", {
	y <- densum(c(-3, 2), bw = 5, dw = 3, from = 0, to = 5, return.x = FALSE)
	expect_equal(y[1:5], densum(2, bw = 5, dw = 3, from = 0, to = 5,
		return.x = FALSE)[1:5])
})

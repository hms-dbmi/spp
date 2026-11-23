# points_withinFunction() / points.within() (points_within in peaks.cpp).
# Membership convention of the native code: a point x lies in fragment
# (fs, fe) when fs < x <= fe (start exclusive, end inclusive); the start
# convention is pinned in the boundary test below.

pw <- spp::points_withinFunction

# brute-force reference: sorted fragment indices containing each point
ref_membership <- function(x, fs, fe) {
	lapply(x, function(p) which(fs < p & p <= fe))
}

test_that("default output is the containing fragment, or -1 when there is none", {
	r <- pw(c(1, 5, 10), c(2, 8), c(6, 12))
	expect_equal(r, c(-1, 1, 2))
})

test_that("results come back in the order of the input points", {
	x <- c(10, 5, 1, 7, 11)
	fs <- c(2, 8); fe <- c(6, 12)
	r <- pw(x, fs, fe)
	expect_equal(r, c(2, 1, -1, -1, 2))
	expect_equal(pw(sort(x), fs, fe, sorted = TRUE), c(-1, 1, -1, 2, 2))
})

test_that("points in several fragments report the lowest index by default", {
	r <- pw(c(5, 20), c(1, 2), c(10, 10))
	expect_equal(r, c(1, -1))
})

test_that("return.unique drops points that fall in more than one fragment", {
	r <- pw(c(5, 5, 12), c(1, 2, 11), c(10, 10, 15), return.unique = TRUE)
	expect_equal(r, c(-1, -1, 3))
})

test_that("return.list gives all fragment indices per point", {
	r <- pw(c(5, 1, 9), c(1, 2), c(10, 8), return.list = TRUE)
	expect_type(r, "list")
	expect_length(r, 3)
	expect_equal(r[[1]], c(1, 2))
	expect_null(r[[2]])
	expect_equal(r[[3]], 1)
})

test_that("return.list with return.unique keeps only single-fragment points", {
	r <- pw(c(5, 1, 9), c(1, 2), c(10, 8), return.list = TRUE, return.unique = TRUE)
	expect_null(r[[1]])
	expect_null(r[[2]])
	expect_equal(r[[3]], 1)
})

test_that("return.point.counts counts points per fragment", {
	r <- pw(c(5, 5, 5, 30), c(1, 2), c(10, 10), return.point.counts = TRUE)
	expect_equal(r, c(3, 3))
	r2 <- pw(c(3, 9, 11), c(1, 8, 20), c(5, 12, 21), return.point.counts = TRUE)
	expect_equal(r2, c(1, 2, 0))
})

test_that("start is exclusive and end is inclusive (pins the current output)", {
	r <- pw(c(2, 6, 7), c(2, 6), c(6, 9))
	# x = 2: not in (2, 6]; x = 6: in (2, 6] only; x = 7: in (6, 9]
	expect_equal(r, c(-1, 1, 2))
})

# The brute-force reference uses the start-exclusive convention fs < x <= fe, so
# this test pins that convention (see the boundary test above).
test_that("results agree with a start-exclusive brute force on random fragments and points", {
	set.seed(11)
	for (it in 1:20) {
		n <- sample(1:30, 1); m <- sample(1:8, 1)
		x <- sample(0:100, n, replace = TRUE)
		fs <- sample(0:90, m, replace = TRUE)
		fe <- fs + sample(1:30, m, replace = TRUE)
		want <- ref_membership(x, fs, fe)
		lst <- pw(x, fs, fe, return.list = TRUE)
		for (k in seq_len(n)) {
			if (length(want[[k]]) == 0) {
				expect_null(lst[[k]])
			} else {
				expect_equal(lst[[k]], want[[k]])
			}
		}
		first <- vapply(want, function(w) if (length(w)) w[1] else -1L, numeric(1))
		expect_equal(pw(x, fs, fe), first)
		uniq <- vapply(want, function(w) if (length(w) == 1) w else -1L, numeric(1))
		expect_equal(pw(x, fs, fe, return.unique = TRUE), uniq)
		counts <- vapply(seq_len(m), function(f) sum(vapply(want,
			function(w) f %in% w, logical(1))), numeric(1))
		expect_equal(pw(x, fs, fe, return.point.counts = TRUE), counts)
	}
})

test_that("integer and numeric inputs agree", {
	expect_identical(pw(c(3L, 1L), 0L, 10L), pw(c(3, 1), 0, 10))
})

test_that("empty point vectors return an empty result", {
	expect_length(pw(c(), 1, 2), 0)
	expect_length(pw(integer(0), 1, 2), 0)
})

test_that("a single point and a single fragment work", {
	expect_equal(pw(5, 1, 10), 1)
	expect_equal(pw(5, 6, 10), -1)
	expect_equal(pw(11, 1, 10), -1)
})

test_that("points beyond all fragment ends and before all starts get -1", {
	expect_equal(pw(c(-5, 0, 100), c(1, 20), c(10, 30)), c(-1, -1, -1))
})

test_that("the deprecated points.within warns and forwards its arguments", {
	expect_warning(r <- spp::points.within(c(1, 5), 0, 3), "deprecated")
	expect_equal(r, c(1, -1))
	suppressWarnings(r2 <- spp::points.within(c(1, 5), c(0, 4), c(3, 8),
		return.point.counts = TRUE))
	expect_equal(r2, c(1, 1))
})

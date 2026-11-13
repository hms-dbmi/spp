# window.tag.count() and window.tag.count.around() (cwindow_n_tags and
# cwindow_n_tags_around in cdensum.c).

window.tag.count <- spp_internal("window.tag.count")
window.tag.count.around <- spp_internal("window.tag.count.around")

# Windows are [L, L + size] (both ends inclusive) with L = from + step * i -
# floor(size / 2).
ref_window_count <- function(tags, size, step, from, nsteps) {
	vapply(seq_len(nsteps) - 1, function(i) {
		L <- from + step * i - floor(size / 2)
		sum(tags >= L & tags <= L + size)
	}, numeric(1))
}

test_that("window.tag.count matches the hand-computed counts", {
	# tags 1 3 4 8 9 15, size 4: from = 1 + 2 = 3, to = 15 - 2 = 13
	# window i spans [1 + i, 5 + i] inclusive
	r <- window.tag.count(c(1, 3, 4, 8, 9, 15), 4)
	expect_named(r, c("x", "y", "step"))
	expect_equal(r$x, c(3, 12))
	expect_equal(r$step, 1)
	expect_equal(r$y, c(3, 2, 2, 2, 2, 2, 2, 2, 1, 0))
})

test_that("window.tag.count is insensitive to input order and storage type", {
	tags <- c(1, 3, 4, 8, 9, 15)
	a <- window.tag.count(tags, 4)
	expect_identical(window.tag.count(rev(tags), 4), a)
	expect_identical(window.tag.count(as.integer(tags), 4L), a)
})

test_that("window.tag.count matches a brute-force count for sizes and steps", {
	tags <- c(3, 3, 7, 10, 11, 11, 11, 25, 26, 40, 41, 44, 58, 59, 60, 90)
	for (size in c(1, 4, 5, 10, 21)) {
		for (step in c(1, 2, 3, 7)) {
			from <- min(tags) + floor(size / 2)
			to <- max(tags) - floor(size / 2)
			nsteps <- ceiling((to - from) / step)
			r <- window.tag.count(tags, size, window.step = step)
			expect_equal(r$y, ref_window_count(tags, size, step, from, nsteps))
			expect_equal(r$x, c(from, from + (nsteps - 1) * step))
			expect_equal(r$step, step)
		}
	}
})

test_that("from and to select the scanned range", {
	tags <- c(1, 3, 4, 8, 9, 15, 20, 21, 22)
	r <- window.tag.count(tags, 4, from = 8, to = 20)
	expect_equal(r$x, c(8, 19))
	expect_equal(r$y, ref_window_count(tags, 4, 1, 8, 12))
})

test_that("return.x = FALSE returns only the counts", {
	tags <- c(1, 3, 4, 8, 9, 15)
	y <- window.tag.count(tags, 4, return.x = FALSE)
	expect_equal(y, window.tag.count(tags, 4)$y)
	expect_type(y, "integer")
})

test_that("tv replaces the tag vector that is counted", {
	tags <- c(1, 3, 4, 8, 9, 15)
	full <- window.tag.count(tags, 5, window.step = 3)
	same <- window.tag.count(tags, 5, window.step = 3, tv = sort(tags))
	expect_equal(same, full)
	# a deliberately trimmed tv (tag at 1 dropped) lowers only the first window,
	# which spans [1, 6] for from = 3
	trimmed <- window.tag.count(tags, 5, window.step = 3, tv = c(3, 4, 8, 9, 15))
	expect_equal(trimmed$y[1], full$y[1] - 1)
	expect_equal(trimmed$y[-1], full$y[-1])
})

test_that("odd window sizes use floor(size / 2) and an inclusive window (pins the current output)", {
	# size 5 -> window i spans [1 + i, 6 + i] (size + 1 positions, both ends inclusive)
	r <- window.tag.count(c(1, 3, 4, 8, 9, 15), 5)
	expect_equal(r$x, c(3, 12))
	expect_equal(r$y, c(3, 2, 3, 3, 2, 2, 2, 2, 1, 1))
})

test_that("a zero-length scan range is an error", {
	expect_error(window.tag.count(c(5, 5, 5), 4), "zero data range")
	expect_error(window.tag.count(c(1, 3, 4, 8, 9, 15), 4, from = 3, to = 3),
		"zero data range")
})

test_that("window.tag.count.around counts tags within +/- floor(size / 2)", {
	tags <- c(1, 3, 4, 8, 9, 15)
	r <- window.tag.count.around(tags, 4, c(3, 9))
	expect_equal(r$x, c(3, 9))
	# [1, 5] holds 1 3 4; [7, 11] holds 8 9
	expect_equal(r$y, c(3, 2))
})

test_that("window.tag.count.around matches a brute force and keeps query order", {
	tags <- c(2, 2, 5, 9, 9, 9, 14, 30, 31, 31, 60)
	pos <- c(31, 2, 60, 10, 0, 100, 14, 9)
	for (size in c(0, 1, 6, 10, 21)) {
		whs <- floor(size / 2)
		want <- vapply(pos, function(p) sum(tags >= p - whs & tags <= p + whs),
			numeric(1))
		r <- window.tag.count.around(tags, size, pos)
		expect_equal(r$x, pos)
		expect_equal(r$y, want)
		expect_equal(window.tag.count.around(tags, size, pos, return.x = FALSE),
			want)
	}
})

test_that("window.tag.count.around handles sorted = TRUE and a precomputed table", {
	tags <- c(1, 3, 4, 8, 9, 15)
	r1 <- window.tag.count.around(tags, 4, c(3, 9), sorted = TRUE)
	expect_equal(r1$x, c(3, 9))
	expect_equal(r1$y, c(3, 2))
	expect_equal(window.tag.count.around(tags, 4, c(3, 9), sorted = TRUE,
		return.x = FALSE), c(3, 2))
	r2 <- window.tag.count.around(tags, 4, c(3, 9), tc = table(tags))
	expect_equal(r2$y, c(3, 2))
})

test_that("window.tag.count.around handles duplicated query positions", {
	tags <- c(1, 3, 4, 8, 9, 15)
	r <- window.tag.count.around(tags, 4, c(9, 3, 9))
	expect_equal(r$y, c(2, 3, 2))
})

test_that("window.tag.count.around handles empty tags and empty queries", {
	r <- window.tag.count.around(numeric(0), 4, c(1, 2))
	expect_equal(r$y, c(0, 0))
	r2 <- window.tag.count.around(c(1, 3, 4), 4, numeric(0))
	expect_equal(nrow(r2), 0)
})

test_that("tags at the ends of the tag vector are counted", {
	tags <- c(10, 11, 12)
	# size 2 -> half size 1
	expect_equal(window.tag.count.around(tags, 2, c(10, 12), return.x = FALSE),
		c(2, 2))
	expect_equal(window.tag.count.around(tags, 2, c(9, 13), return.x = FALSE),
		c(1, 1))
	expect_equal(window.tag.count.around(tags, 2, c(8, 14), return.x = FALSE),
		c(0, 0))
})

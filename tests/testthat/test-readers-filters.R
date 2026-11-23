# Tests of remove.local.tag.anomalies and select.informative.tags.

test_that("remove.local.tag.anomalies returns short or empty vectors unchanged", {
	r <- remove.local.tag.anomalies(list(a = integer(0), b = 5L))
	expect_equal(r$a, integer(0))
	expect_equal(r$b, 5L)
})

# Local window: window.size = 200, so +/- 100 around each position. Background
# is 100 single tags at every other position 1000..1198; a hot position at 1101
# (odd, so it does not collide) carries `count` tags. The window around 1101
# holds the 99 background tags at 1002..1198, so the local density is
# (99 + count - count + 1) / 200 = 0.5. With z.threshold 3 the tail probability
# is p = pnorm(3, lower.tail = FALSE). The hot position is removed when its
# count exceeds qpois(p, 0.5 * 10, lower.tail = FALSE) = 13 and is capped to
# 8 when it exceeds qpois(p, 0.5 * 4, lower.tail = FALSE) = 7.
local_tags <- function(count) {
	bg <- seq(1000, 1198, by = 2)
	sort(c(bg, rep(1101, count)))
}

test_that("a position far above local density is removed entirely", {
	# 30 tags > 13 (eliminate threshold): all hot tags dropped, background kept
	tags <- local_tags(30)
	r <- remove.local.tag.anomalies(list(chr1 = tags))
	expect_equal(r$chr1, seq(1000, 1198, by = 2))
})

test_that("a position moderately above local density is capped", {
	# 10 tags: above the cap threshold (7) but not above the elimination one (13);
	# the count is reduced to 7 + 1 = 8
	tags <- local_tags(10)
	r <- remove.local.tag.anomalies(list(chr1 = tags))
	expect_equal(sum(r$chr1 == 1101), 8L)
	expect_equal(r$chr1[r$chr1 != 1101], seq(1000, 1198, by = 2))
})

test_that("positions at the local density are left alone", {
	tags <- local_tags(3)
	r <- remove.local.tag.anomalies(list(chr1 = tags))
	expect_equal(r$chr1, tags)
})

test_that("remove.local.tag.anomalies works per chromosome and keeps strands sorted by abs", {
	pos <- local_tags(1)
	neg <- -pos
	r <- remove.local.tag.anomalies(list(a = pos, b = c(pos, neg)))
	expect_named(r, c("a", "b"))
	expect_equal(r$a, pos)
	expect_equal(length(r$b), length(pos) * 2)
	expect_false(is.unsorted(abs(r$b)))
})

test_that("explicit default arguments give the same exact filtering", {
	# the defaults (200, 10, 4, 3) equal the constants used internally
	tags <- local_tags(30)
	r <- remove.local.tag.anomalies(list(chr1 = tags), window.size = 200,
		eliminate.fold = 10, cap.fold = 4, z.threshold = 3)
	expect_equal(r$chr1, seq(1000, 1198, by = 2))
})

test_that("an anomaly on one strand only is removed and the other strand is untouched", {
	bg <- seq(1000, 1198, by = 2)
	tags <- c(bg, -bg, rep(1101, 30))
	r <- remove.local.tag.anomalies(list(chr1 = tags))
	# negative coordinates form their own windows; ties in abs() keep the
	# negative tag first because positions are tabulated in increasing order
	expect_equal(r$chr1, as.numeric(rbind(-bg, bg)))
})

test_that("select.informative.tags returns all tags without binding characteristics", {
	data <- list(tags = list(chr1 = c(10L, -20L, 30L)), quality = list(chr1 = c(0L, 1L, 0L)))
	expect_equal(select.informative.tags(data), data$tags)
})

test_that("select.informative.tags keeps only tags in informative quality bins", {
	data <- list(
		tags = list(chr1 = c(10L, -20L, 30L), chr2 = c(5L, 6L)),
		quality = list(chr1 = c(0L, 1L, 2L), chr2 = c(1L, 0L))
	)
	bc <- list(quality.bin.acceptance = list(informative.bins = c("0" = TRUE, "1" = FALSE, "2" = TRUE)))
	r <- select.informative.tags(data, bc)
	expect_equal(r$chr1, c(10L, 30L))
	expect_equal(r$chr2, 6L)
})

test_that("select.informative.tags without quality.bin.acceptance accepts all tags", {
	data <- list(tags = list(chr1 = c(10L, 20L)), quality = list(chr1 = c(0L, 1L)))
	out <- capture.output(r <- select.informative.tags(data, list(whs = 100)))
	expect_equal(r, data$tags)
	expect_match(out, "accepting all tags")
})

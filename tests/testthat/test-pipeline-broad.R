# add.broad.peak.regions and get.broad.enrichment.clusters
# Data are fully deterministic: signal and control both have a tag every ~100 bp
# (alternating strands, the control offset by 37 bp), and the signal of chr1 has
# 500 extra tags evenly spread over 8000-9500.  With the default 73 bp shift
# (plus-strand tags move up, minus-strand tags down) the block occupies 7927-9573.

pl.broad.tags <- function(offset, block) {
	x <- round(seq(50, 19950, length.out = 200)) + offset
	if (!is.null(block)) x <- c(x, round(seq(block[1], block[2], length.out = block[3])))
	x <- sort(as.integer(x))
	y <- as.integer(x * rep(c(1L, -1L), length.out = length(x)))
	y[order(abs(y))]
}

pl.broad.data <- function() {
	list(
		sig = list(chr1 = pl.broad.tags(0, c(8000, 9500, 500)), chr2 = pl.broad.tags(0, NULL)),
		ctl = list(chr1 = pl.broad.tags(37, NULL), chr2 = pl.broad.tags(37, NULL))
	)
}

pl.block <- c(7927, 9573)

test_that("get.broad.enrichment.clusters finds the enriched block", {
	d <- pl.broad.data()
	bg <- get.broad.enrichment.clusters(d$sig, d$ctl, window.size = 500, z.thr = 2)
	expect_equal(names(bg), c("chr1", "chr2"))
	# exactly one cluster overlaps the block and it contains it
	hit <- bg$chr1$s <= pl.block[2] & bg$chr1$e >= pl.block[1]
	expect_equal(sum(hit), 1)
	expect_lte(bg$chr1$s[hit], pl.block[1])
	expect_gte(bg$chr1$e[hit], pl.block[2])
	expect_gt(bg$chr1$rv[hit], 1)
	# no enrichment on the chromosome with the same density in signal and control
	expect_equal(length(bg$chr2$s), 0)
	expect_equal(length(bg$chr2$e), 0)
})

test_that("clusters extend past the block by a fraction of the window size", {
	d <- pl.broad.data()
	for (w in c(200, 500, 1000)) {
		bg <- get.broad.enrichment.clusters(d$sig, d$ctl, window.size = w, z.thr = 3)
		hit <- bg$chr1$s <= pl.block[2] & bg$chr1$e >= pl.block[1]
		expect_equal(sum(hit), 1)
		left <- pl.block[1] - bg$chr1$s[hit]
		right <- bg$chr1$e[hit] - pl.block[2]
		# the cluster is built from windows of w tags positions and padded by w/2 on
		# each side, trimmed where the window is not yet significant
		expect_true(left >= 0.25 * w && left <= w)
		expect_true(right >= 0.25 * w && right <= w)
	}
})

test_that("a larger window size widens the clusters", {
	d <- pl.broad.data()
	width <- function(w) {
		bg <- get.broad.enrichment.clusters(d$sig, d$ctl, window.size = w, z.thr = 3)
		hit <- bg$chr1$s <= pl.block[2] & bg$chr1$e >= pl.block[1]
		bg$chr1$e[hit] - bg$chr1$s[hit]
	}
	w2 <- width(200)
	w5 <- width(500)
	w10 <- width(1000)
	expect_lt(w2, w5)
	expect_lt(w5, w10)
	# every extra bp of window size adds about 1.8 bp to the cluster (both ends)
	expect_true((w10 - w5) / 500 > 1.4 && (w10 - w5) / 500 < 2.2)
	expect_true((w5 - w2) / 300 > 1.4 && (w5 - w2) / 300 < 2.2)
})

test_that("cluster scores are the log2 of the pseudocounted signal/control ratio", {
	d <- pl.broad.data()
	shift <- 146 / 2
	st <- abs(d$sig$chr1 + shift)
	ct <- abs(d$ctl$chr1 + shift)
	for (bds in c(FALSE, TRUE)) {
		bg <- get.broad.enrichment.clusters(d$sig, d$ctl, window.size = 500, z.thr = 2, background.density.scaling = bds)
		s <- bg$chr1$s
		e <- bg$chr1$e
		# points_within treats region ends unevenly, so make sure no tag sits on one
		expect_false(any(st %in% c(s, e)) || any(ct %in% c(s, e)))
		# weight: tag-count ratio, or the ratio of tags outside enriched regions with scaling
		bgw <- if (bds) {
			spp:::dataset.density.ratio(d$sig, d$ctl, background.density.scaling = TRUE)
		} else {
			pl.ntags(d$sig) / pl.ntags(d$ctl)
		}
		expect_gte(length(s), 1)
		for (i in seq_along(s)) {
			sc <- sum(st > s[i] & st < e[i])
			cc <- sum(ct > s[i] & ct < e[i])
			expect_equal(bg$chr1$rv[i], log2((sc + 1) / (cc + 1) / bgw))
		}
	}
})

test_that("a very high z threshold leaves no clusters", {
	d <- pl.broad.data()
	bg <- get.broad.enrichment.clusters(d$sig, d$ctl, window.size = 500, z.thr = 1000)
	expect_equal(length(bg$chr1$s), 0)
	expect_equal(length(bg$chr2$s), 0)
})

test_that("identical signal and control have no broad clusters", {
	d <- pl.broad.data()
	bg <- get.broad.enrichment.clusters(d$ctl, d$ctl, window.size = 500, z.thr = 2)
	expect_equal(length(bg$chr1$s), 0)
	expect_equal(length(bg$chr2$s), 0)
})

test_that("add.broad.peak.regions annotates peaks inside enriched regions", {
	d <- pl.broad.data()
	bp <- list(npl = list(
		chr1 = data.frame(x = c(1000, 8700, 9100, 15000), y = 1:4),
		chr2 = data.frame(x = numeric(0), y = numeric(0))
	))
	res <- add.broad.peak.regions(d$sig, d$ctl, bp, window.size = 500, z.thr = 2)
	expect_equal(names(res$npl), c("chr1", "chr2"))
	npl <- res$npl$chr1
	expect_equal(npl$x, bp$npl$chr1$x)
	expect_equal(npl$y, 1:4)
	expect_true(all(c("rs", "re") %in% colnames(npl)))
	# the two peaks inside the block share one region that covers the block
	expect_equal(npl$rs[2], npl$rs[3])
	expect_equal(npl$re[2], npl$re[3])
	expect_lte(npl$rs[2], pl.block[1])
	expect_gte(npl$re[2], pl.block[2])
	# peaks in quiet areas are not annotated
	expect_true(all(is.na(npl$rs[c(1, 4)]) & is.na(npl$re[c(1, 4)])))
	# the empty chromosome is returned untouched
	expect_equal(nrow(res$npl$chr2), 0)
})

test_that("add.broad.peak.regions equals a brute-force region lookup, including at the margins", {
	d <- pl.broad.data()
	se <- spp:::find.significantly.enriched.regions(d$sig, d$ctl, window.size = 500, z.thr = 2, poisson.z = 0, poisson.ratio = 0, either = FALSE)$chr1
	expect_gte(length(se$s), 1)
	# peaks just inside and just outside both ends of every region, and far away
	xs <- sort(unique(c(1000, 15000, se$s - 5, se$s + 5, se$e - 5, se$e + 5)))
	bp <- list(npl = list(chr1 = data.frame(x = xs, y = seq_along(xs))))
	res <- add.broad.peak.regions(d$sig, d$ctl, bp, window.size = 500, z.thr = 2)$npl$chr1
	for (i in seq_along(xs)) {
		inside <- which(se$s <= xs[i] & se$e >= xs[i])
		if (length(inside) == 0) {
			expect_true(is.na(res$rs[i]))
			expect_true(is.na(res$re[i]))
		} else {
			expect_equal(res$rs[i], min(se$s[inside]))
			expect_equal(res$re[i], max(se$e[inside]))
		}
	}
	expect_true(is.na(res$rs[res$x == min(se$s) - 5]))
	expect_false(is.na(res$rs[res$x == min(se$s) + 5]))
	expect_true(is.na(res$re[res$x == max(se$e) + 5]))
	expect_false(is.na(res$re[res$x == max(se$e) - 5]))
})

test_that("broad region boundaries flow into the narrowPeak output", {
	d <- pl.broad.data()
	bp <- list(
		npl = list(chr1 = data.frame(x = c(1017, 8713, 9117), y = c(10, 30, 20), evalue = 0.01, fdr = c(0.01, 0.001, 0.01), enr = 2, enr.mle = 3)),
		whs = 100
	)
	bd <- add.broad.peak.regions(d$sig, d$ctl, bp, window.size = 500, z.thr = 2)
	f <- tempfile()
	write.narrowpeak.binding(bd, f)
	np <- read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
	expect_equal(nrow(np), 3)
	expect_equal(np$V7, c(30, 20, 10))
	# the peak in the quiet area uses the default margin (whs = 100)
	expect_equal(np$V2[3], 1017 - 100)
	expect_equal(np$V3[3], 1017 + 100)
	# the two block peaks use the region ends and the offset of the peak inside it
	expect_equal(np$V2[1], bd$npl$chr1$rs[2])
	expect_equal(np$V3[1], bd$npl$chr1$re[2])
	expect_equal(np$V10[1], 8713 - bd$npl$chr1$rs[2])
	expect_equal(np$V10[2], 9117 - bd$npl$chr1$rs[3])
})

test_that("broad clusters can be written as broadPeak rows", {
	d <- pl.broad.data()
	# a second block on chr2 so that the file has several rows
	d$sig$chr2 <- pl.broad.tags(0, c(3000, 4000, 400))
	bg <- get.broad.enrichment.clusters(d$sig, d$ctl, window.size = 500, z.thr = 2)
	expect_gte(length(bg$chr1$s) + length(bg$chr2$s), 2)
	f <- tempfile()
	write.broadpeak.info(bg, f)
	out <- read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
	expect_equal(ncol(out), 9)
	expect_equal(nrow(out), length(bg$chr1$s) + length(bg$chr2$s))
	expect_true(!is.unsorted(rev(out$V7)))
	expect_equal(sort(out$V7, decreasing = TRUE), sort(c(bg$chr1$rv, bg$chr2$rv), decreasing = TRUE))
	expect_equal(sort(out$V2), sort(c(bg$chr1$s, bg$chr2$s)))
})

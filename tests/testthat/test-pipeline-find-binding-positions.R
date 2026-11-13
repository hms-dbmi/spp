# find.binding.positions on tiny synthetic pile-ups.  Each pile is a block of
# plus-strand tags 50 bp upstream and minus-strand tags 50 bp downstream of a
# binding site, so the true peak position is the site itself.

pl.fbp <- function(signal, ..., seed = 1) {
	pl.seeded(seed, pl.quiet(find.binding.positions(signal, debug = FALSE, ...)))
}

pl.peak.x <- function(res, chr) {
	d <- res$npl[[chr]]
	if (is.null(d)) numeric(0) else sort(d$x)
}

pl.n.peaks <- function(res) length(unlist(lapply(res$npl, function(d) d$x)))

test_that("a symmetric pile gives a peak at its centre (no control)", {
	res <- pl.fbp(pl.signal(), e.value = 1e6, whs = 100)
	expect_equal(length(pl.peak.x(res, "chr1")), 2)
	expect_equal(length(pl.peak.x(res, "chr2")), 1)
	expect_true(all(abs(pl.peak.x(res, "chr1") - c(3000, 8000)) <= 2))
	expect_true(abs(pl.peak.x(res, "chr2") - 5000) <= 2)
})

test_that("empty chromosomes yield no peaks and no error", {
	res <- pl.fbp(pl.signal(), e.value = 1e6, whs = 100)
	expect_equal(NROW(res$npl[["chrE"]]), 0)
	only.empty <- list(chrE = integer(0), chrF = integer(0))
	res0 <- pl.fbp(c(pl.signal()["chr1"], only.empty), e.value = 1e6, whs = 100)
	expect_equal(NROW(res0$npl[["chrE"]]), 0)
	expect_equal(NROW(res0$npl[["chrF"]]), 0)
})

test_that("result bookkeeping fields record the inputs", {
	sig <- pl.signal()
	ctl <- pl.control()
	res <- pl.fbp(sig, control.data = ctl, e.value = 1, whs = 100, background.density.scaling = FALSE)
	expect_equal(res$whs, 100)
	expect_equal(res$f, 1)
	expect_equal(res$n, 60 * 3)
	expect_equal(res$n.bg, pl.ntags(ctl))
	expect_equal(res$thr$type, "Evalue")
	expect_equal(res$thr$e.value, 1)
	expect_equal(res$enrichment.whs, 100)
	expect_equal(res$enrichment.z, 2)
	expect_equal(res$enrichment.bg.weight, 180 / pl.ntags(ctl))
})

test_that("with a control, peaks are still at the pile centres", {
	res <- pl.fbp(pl.signal(), control.data = pl.control(), e.value = 1, whs = 100, background.density.scaling = FALSE)
	expect_equal(length(pl.peak.x(res, "chr1")), 2)
	expect_true(all(abs(pl.peak.x(res, "chr1") - c(3000, 8000)) <= 5))
	expect_equal(length(pl.peak.x(res, "chr2")), 1)
	expect_true(abs(pl.peak.x(res, "chr2") - 5000) <= 5)
	expect_equal(NROW(res$npl[["chrE"]]), 0)
})

test_that("background.density.scaling TRUE and FALSE both recover the piles", {
	for (bds in c(TRUE, FALSE)) {
		res <- pl.fbp(pl.signal(), control.data = pl.control(), e.value = 1, whs = 100, background.density.scaling = bds)
		expect_equal(length(pl.peak.x(res, "chr1")), 2)
		expect_true(all(abs(pl.peak.x(res, "chr1") - c(3000, 8000)) <= 5))
		expect_true(abs(pl.peak.x(res, "chr2") - 5000) <= 5)
	}
})

test_that("tag.wtd and tag.lwcc both locate the pile centres", {
	for (m in list(spp:::tag.wtd, spp:::tag.lwcc)) {
		res <- pl.fbp(pl.signal(), control.data = pl.control(), fdr = 1, whs = 100, method = m)
		expect_equal(length(pl.peak.x(res, "chr1")), 2)
		expect_true(all(abs(pl.peak.x(res, "chr1") - c(3000, 8000)) <= 5))
		expect_true(abs(pl.peak.x(res, "chr2") - 5000) <= 5)
		expect_equal(res$thr$type, "FDR")
		expect_equal(res$thr$fdr, 1)
	}
})

test_that("whs smaller than the half distance between strand piles finds nothing", {
	# plus and minus piles are 100 bp apart, a window half size of 10 never sees both
	res <- pl.fbp(pl.signal(), control.data = pl.control(), e.value = 1, whs = 10, background.density.scaling = FALSE)
	expect_equal(pl.n.peaks(res), 0)
	res <- pl.fbp(pl.signal(), control.data = pl.control(), e.value = 1, whs = 100, background.density.scaling = FALSE)
	expect_gt(pl.n.peaks(res), 0)
})

test_that("an fdr nothing can satisfy gives no peaks, a satisfiable one gives them", {
	ctl <- pl.control()
	none <- pl.fbp(pl.signal(), control.data = ctl, fdr = 1e-9, whs = 100)
	expect_equal(none$thr$type, "FDR")
	expect_equal(pl.n.peaks(none), 0)
	some <- pl.fbp(pl.signal(), control.data = ctl, fdr = 1, whs = 100)
	expect_equal(pl.n.peaks(some), 3)
})

test_that("a looser e.value keeps a superset of the peaks of a stricter one", {
	sig <- pl.signal()
	# unequal pile sizes so that scores differ between peaks
	sig$chr1 <- c(pl.tags(c(3000, 8000, 11000), n = 30), pl.tags(14000, n = 12))
	sig$chr1 <- sig$chr1[order(abs(sig$chr1))]
	ctl <- pl.control()
	strict <- pl.fbp(sig, control.data = ctl, e.value = 0.2, whs = 100, background.density.scaling = FALSE)
	loose <- pl.fbp(sig, control.data = ctl, e.value = 1e6, whs = 100, background.density.scaling = FALSE)
	for (chr in c("chr1", "chr2")) {
		expect_true(all(pl.peak.x(strict, chr) %in% pl.peak.x(loose, chr)))
	}
	expect_gte(pl.n.peaks(loose), pl.n.peaks(strict))
	expect_equal(pl.n.peaks(loose), 5)
	expect_gte(strict$thr$root, loose$thr$root)
	for (d in loose$npl) expect_true(all(d$y >= loose$thr$root))
	for (d in strict$npl) expect_true(all(d$y >= strict$thr$root))
})

test_that("peak statistics columns are consistent with brute-force counts", {
	sig <- pl.signal()
	ctl <- pl.control()
	res <- pl.fbp(sig, control.data = ctl, e.value = 1, whs = 100, tag.count.whs = 100, enrichment.background.scales = c(1, 5), background.density.scaling = FALSE)
	bgw <- 180 / pl.ntags(ctl)
	for (chr in c("chr1", "chr2")) {
		d <- res$npl[[chr]]
		# points_within handles tags exactly on a window end unevenly, so make sure none is there
		expect_false(any(outer(abs(c(sig[[chr]], ctl[[chr]])), d$x, "-") %in% c(-100, 100, -500, 500)))
		nt <- pl.count.within(abs(sig[[chr]]), d$x, 100)
		nb <- pl.count.within(abs(ctl[[chr]]), d$x, 100)
		nb5 <- pl.count.within(abs(ctl[[chr]]), d$x, 500)
		expect_equal(d$nt, nt)
		expect_equal(d$nbgt, nb)
		expect_equal(d$nbgt.5, nb5)
		expect_equal(d$enr.mle, (nt + 0.5) / (nb + 0.5) / bgw)
		expect_equal(d$enr.mle.5, (nt + 0.5) / (nb5 + 0.5) / (bgw / 5))
		expect_true(all(d$enr <= d$enr.mle & d$enr.mle <= d$enr.ub))
		expect_true(all(d$enr.5 <= d$enr.mle.5 & d$enr.mle.5 <= d$enr.ub.5))
		expect_true(all(d$evalue > 0 & d$fdr > 0))
	}
	expect_gt(sum(unlist(lapply(res$npl, function(d) d$nbgt.5))), 0)
})

test_that("tag.count.whs sets the window of the tag counts", {
	sig <- pl.signal()
	ctl <- pl.control()
	res <- pl.fbp(sig, control.data = ctl, e.value = 1, whs = 100, tag.count.whs = 70, background.density.scaling = FALSE)
	expect_equal(res$enrichment.whs, 70)
	for (chr in c("chr1", "chr2")) {
		d <- res$npl[[chr]]
		expect_false(any(outer(abs(c(sig[[chr]], ctl[[chr]])), d$x, "-") %in% c(-70, 70)))
		expect_equal(d$nt, pl.count.within(abs(sig[[chr]]), d$x, 70))
		expect_equal(d$nbgt, pl.count.within(abs(ctl[[chr]]), d$x, 70))
	}
})

test_that("without a control the enrichment is estimated from the tag count", {
	sig <- pl.signal()
	res <- pl.fbp(sig, e.value = 1e6, whs = 100, tag.count.whs = 100)
	lambda <- pl.ntags(sig) * 2 * 100 / 3e9
	d <- res$npl$chr1
	expect_false(any(outer(abs(sig$chr1), d$x, "-") %in% c(-100, 100)))
	expect_equal(d$nt, pl.count.within(abs(sig$chr1), d$x, 100))
	expect_equal(d$enr.mle, d$nt / lambda)
	expect_true(all(d$enr < d$enr.mle & d$enr.mle < d$enr.ub))
})

test_that("min.dist suppresses weaker neighbours closer than the distance", {
	sig <- list(chr1 = c(pl.tags(3000, n = 40), pl.tags(3400, n = 20)))
	sig$chr1 <- sig$chr1[order(abs(sig$chr1))]
	ctl <- list(chr1 = pl.bg(200, seed = 3))
	near <- pl.fbp(sig, control.data = ctl, e.value = 1e6, whs = 100, min.dist = 100, background.density.scaling = FALSE)
	far <- pl.fbp(sig, control.data = ctl, e.value = 1e6, whs = 100, min.dist = 1000, background.density.scaling = FALSE)
	expect_equal(length(pl.peak.x(near, "chr1")), 2)
	expect_equal(length(pl.peak.x(far, "chr1")), 1)
	expect_true(abs(pl.peak.x(far, "chr1") - 3000) <= 50)
	expect_true(all(diff(pl.peak.x(near, "chr1")) >= 100))
})

test_that("topN returns only the top-scoring peaks and records the threshold", {
	sig <- pl.signal()
	# distinct pile sizes give distinct scores (ties at the threshold are all dropped)
	sig$chr1 <- c(pl.tags(3000, n = 30), pl.tags(8000, n = 26), pl.tags(11000, n = 14))
	sig$chr1 <- sig$chr1[order(abs(sig$chr1))]
	sig$chr2 <- pl.tags(5000, n = 22)
	ctl <- pl.control()
	all.peaks <- pl.fbp(sig, control.data = ctl, topN = 100, whs = 100, background.density.scaling = FALSE)
	top <- pl.fbp(sig, control.data = ctl, topN = 2, whs = 100, background.density.scaling = FALSE)
	ay <- unname(sort(unlist(lapply(all.peaks$npl, function(d) d$y)), decreasing = TRUE))
	ty <- unname(sort(unlist(lapply(top$npl, function(d) d$y)), decreasing = TRUE))
	expect_equal(length(ay), 4)
	expect_equal(top$thr$type, "topN")
	expect_equal(top$thr$topN, 2)
	expect_equal(unname(top$thr$root), ay[2])
	# what is returned is a non-empty prefix of the ranking with at most N peaks
	expect_gte(length(ty), 1)
	expect_lte(length(ty), 2)
	expect_equal(ty, ay[seq_along(ty)])
	expect_lte(top$thr$root, min(ty))
})

test_that("only chromosomes present in both signal and control are called", {
	sig <- pl.signal()
	ctl <- pl.control()
	ctl$chr2 <- NULL
	res <- pl.fbp(sig, control.data = ctl, e.value = 1, whs = 100, background.density.scaling = FALSE)
	expect_equal(length(pl.peak.x(res, "chr2")), 0)
	expect_equal(length(pl.peak.x(res, "chr1")), 2)
	expect_equal(res$n, length(sig$chr1))
	expect_equal(res$n.bg, length(ctl$chr1))
})

test_that("f < 1 subsamples each chromosome by floor(length * f)", {
	res <- pl.fbp(pl.signal(), control.data = pl.control(), e.value = 1e6, whs = 100, f = 0.5, background.density.scaling = FALSE)
	expect_equal(res$f, 0.5)
	expect_equal(res$n, 60 + 30 + 0)
})

test_that("select.informative.tags feeds find.binding.positions", {
	pile <- pl.signal()
	noise <- list(chr1 = pl.bg(150, seed = 11), chr2 = pl.bg(60, seed = 12), chrE = integer(0))
	# the pile is placed first, noise after it; sort jointly and label by origin
	tags <- list()
	quality <- list()
	for (chr in names(pile)) {
		x <- c(pile[[chr]], noise[[chr]])
		q <- c(rep("1", length(pile[[chr]])), rep("2", length(noise[[chr]])))
		o <- order(abs(x), q)
		tags[[chr]] <- x[o]
		quality[[chr]] <- q[o]
	}
	data <- list(tags = tags, quality = quality)
	bc <- list(quality.bin.acceptance = list(informative.bins = c("1" = TRUE, "2" = FALSE)))
	kept <- select.informative.tags(data, bc)
	for (chr in names(tags)) {
		expect_equal(kept[[chr]], tags[[chr]][quality[[chr]] == "1"])
		expect_equal(sort(kept[[chr]]), sort(pile[[chr]]))
	}
	# with no characteristics, or none carrying quality info, all tags are kept
	expect_equal(pl.quiet(select.informative.tags(data, NULL)), tags)
	expect_equal(pl.quiet(select.informative.tags(data, list(peak = 1))), tags)

	via.select <- pl.fbp(kept, control.data = pl.control(), e.value = 1, whs = 100, background.density.scaling = FALSE)
	expect_equal(length(pl.peak.x(via.select, "chr1")), 2)
	expect_true(all(abs(pl.peak.x(via.select, "chr1") - c(3000, 8000)) <= 5))
	expect_true(abs(pl.peak.x(via.select, "chr2") - 5000) <= 5)
})

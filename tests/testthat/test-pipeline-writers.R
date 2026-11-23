# writewig, output.binding.results, write.narrowpeak.binding, write.broadpeak.info

pl.peaks <- function() {
	list(
		npl = list(
			chr1 = data.frame(x = c(1007, 5013, 9021), y = c(50, 80, 65), evalue = c(0.1, 0.01, 0.05), fdr = c(0.1, 0.001, 0.01), enr = c(3, 6, 4), enr.mle = c(4, 8, 5)),
			chr2 = data.frame(x = 2011, y = 70, evalue = 0.02, fdr = 0.01, enr = 5, enr.mle = 6)
		),
		thr = list(root = 50, type = "Evalue", e.value = 1),
		whs = 100
	)
}

# ---- writewig -------------------------------------------------------------

test_that("writewig output parses back to the input values", {
	dat <- list(
		chr1 = data.frame(x = c(100, 200, 300, 400), y = c(1.5, NA, -2.25, 0)),
		chr2 = data.frame(x = c(50, 150), y = c(3, 4))
	)
	f <- tempfile(fileext = ".wig")
	out <- writewig(dat, f, "myfeat")
	expect_equal(out, f)
	lines <- readLines(f)
	expect_equal(length(lines), 1 + 3 + 2)
	expect_true(grepl("^track type=wiggle_0 name=\"Bed Format\" description=\"myfeat\"", lines[1]))
	tab <- read.table(f, skip = 1, header = FALSE, stringsAsFactors = FALSE)
	expect_equal(tab$V1, c("chr1", "chr1", "chr1", "chr2", "chr2"))
	# the NA position is skipped
	expect_equal(tab$V2, c(100, 300, 400, 50, 150))
	expect_equal(tab$V4, c(1.5, -2.25, 0, 3, 4))
	# probe length 35 when positions are far apart
	expect_equal(tab$V3 - tab$V2, rep(35, 5))
})

test_that("writewig shortens the probe to fit the closest spacing", {
	dat <- list(chr1 = data.frame(x = c(100, 120, 400), y = c(1, 2, 3)))
	f <- tempfile()
	pl.quiet(writewig(dat, f, "a"))
	tab <- read.table(f, skip = 1, stringsAsFactors = FALSE)
	# min spacing is 20, so probes are 19 long
	expect_equal(tab$V3 - tab$V2, rep(19, 3))
	expect_equal(tab$V2, c(100, 120, 400))
})

test_that("writewig writes the header once and the description is the feature", {
	dat <- list(a = data.frame(x = c(10, 100), y = c(1, 2)), b = data.frame(x = c(10, 100), y = c(3, 4)), c = data.frame(x = c(10, 100), y = c(5, 6)))
	f <- tempfile()
	writewig(dat, f, "feature one")
	lines <- readLines(f)
	expect_equal(sum(grepl("^track", lines)), 1)
	expect_true(grepl("description=\"feature one\"", lines[1], fixed = TRUE))
	expect_equal(length(lines), 7)
})

test_that("writewig overwrites an existing file", {
	dat <- list(a = data.frame(x = c(10, 100), y = c(1, 2)))
	f <- tempfile()
	writewig(dat, f, "a")
	writewig(dat, f, "a")
	expect_equal(length(readLines(f)), 3)
})

# ---- output.binding.results ----------------------------------------------------

test_that("output.binding.results writes the documented table", {
	bd <- pl.peaks()
	f <- tempfile()
	output.binding.results(bd, f)
	lines <- readLines(f)
	expect_equal(lines[1], "chr\tpos\tscore\tEvalue\tFDR\tenrichment.lb\tenrichment.mle")
	tab <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
	expect_equal(nrow(tab), 4)
	expect_equal(tab$chr, c("chr1", "chr1", "chr1", "chr2"))
	expect_equal(tab$pos, c(1007, 5013, 9021, 2011))
	expect_equal(tab$score, c(50, 80, 65, 70))
	expect_equal(tab$Evalue, c(0.1, 0.01, 0.05, 0.02))
	expect_equal(tab$FDR, c(0.1, 0.001, 0.01, 0.01))
	expect_equal(tab$enrichment.lb, c(3, 6, 4, 5))
	expect_equal(tab$enrichment.mle, c(4, 8, 5, 6))
})

test_that("output.binding.results skips chromosomes with no rows", {
	bd <- pl.peaks()
	bd$npl$chr3 <- bd$npl$chr1[0, ]
	f <- tempfile()
	output.binding.results(bd, f)
	tab <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
	expect_equal(unique(tab$chr), c("chr1", "chr2"))
	expect_equal(nrow(tab), 4)
	# nothing at all: header only
	bd0 <- list(npl = list(chr1 = bd$npl$chr1[0, ]), thr = bd$thr)
	output.binding.results(bd0, f)
	expect_equal(readLines(f), "chr\tpos\tscore\tEvalue\tFDR\tenrichment.lb\tenrichment.mle")
})

test_that("output.binding.results overwrites an existing file", {
	bd <- pl.peaks()
	f <- tempfile()
	output.binding.results(bd, f)
	output.binding.results(bd, f)
	expect_equal(length(readLines(f)), 5)
})

test_that("output.binding.results round-trips a real find.binding.positions result", {
	res <- pl.seeded(1, pl.quiet(find.binding.positions(pl.signal(), debug = FALSE, e.value = 1e6, whs = 100)))
	f <- tempfile()
	output.binding.results(res, f)
	tab <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
	expect_equal(nrow(tab), 3)
	expect_equal(sort(tab$pos[tab$chr == "chr1"]), sort(res$npl$chr1$x))
	expect_equal(tab$score, unname(unlist(lapply(res$npl, function(d) d$y))))
	expect_equal(tab$enrichment.mle, unname(unlist(lapply(res$npl, function(d) d$enr.mle))), tolerance = 1e-6)
})

# ---- write.narrowpeak.binding ------------------------------------------------------

pl.read.np <- function(f) read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

test_that("narrowPeak rows have the documented ten columns", {
	bd <- pl.peaks()
	f <- tempfile()
	write.narrowpeak.binding(bd, f)
	np <- pl.read.np(f)
	expect_equal(ncol(np), 10)
	expect_equal(nrow(np), 4)
	# sorted by signal value (column 7), strongest first
	expect_equal(np$V7, c(80, 70, 65, 50))
	expect_equal(np$V1, c("chr1", "chr2", "chr1", "chr1"))
	expect_equal(np$V4, rep(".", 4))
	expect_equal(np$V5, rep(0, 4))
	expect_equal(np$V6, rep(".", 4))
	expect_equal(np$V8, rep(-1, 4))
	# the default margin is the window half size recorded in the result (100)
	x <- c(5013, 2011, 9021, 1007)
	expect_equal(np$V2, x - 100)
	expect_equal(np$V3, x + 100)
	# qValue is -log10(fdr); the peak offset is the position relative to the start
	expect_equal(np$V9, -log10(c(0.001, 0.01, 0.01, 0.1)))
	expect_equal(np$V10, rep(100, 4))
})

test_that("narrowPeak margin, npeaks and the missing-whs default", {
	bd <- pl.peaks()
	f <- tempfile()
	write.narrowpeak.binding(bd, f, margin = 10, npeaks = 2)
	np <- pl.read.np(f)
	expect_equal(nrow(np), 2)
	expect_equal(np$V7, c(80, 70))
	expect_equal(np$V3 - np$V2, c(20, 20))
	expect_equal(np$V10, c(10, 10))
	# npeaks larger than the number of peaks keeps all of them
	write.narrowpeak.binding(bd, f, npeaks = 100)
	expect_equal(nrow(pl.read.np(f)), 4)
	# no whs recorded: the margin falls back to 50
	bd$whs <- NULL
	write.narrowpeak.binding(bd, f)
	np <- pl.read.np(f)
	expect_equal(np$V3 - np$V2, rep(100, 4))
	expect_equal(np$V10, rep(50, 4))
})

test_that("narrowPeak uses broad-region boundaries where present", {
	bd <- pl.peaks()
	bd$npl$chr1$rs <- c(907, NA, 8821)
	bd$npl$chr1$re <- c(1307, NA, 9521)
	bd$npl$chr2$rs <- NA
	bd$npl$chr2$re <- NA
	f <- tempfile()
	write.narrowpeak.binding(bd, f)
	np <- pl.read.np(f)
	# order: 5013 (80), 2011 (70), 9021 (65), 1007 (50)
	expect_equal(np$V2, c(4913, 1911, 8821, 907))
	expect_equal(np$V3, c(5113, 2111, 9521, 1307))
	expect_equal(np$V10, c(100, 100, 200, 100))
})

# ---- write.broadpeak.info -------------------------------------------------------

test_that("broadPeak rows hold region bounds and the log2 ratio, strongest first", {
	bp <- list(
		chr1 = list(s = c(1111, 5123), e = c(2222, 6234), rv = c(1.5, 3.25)),
		chr2 = list(s = 7111, e = 8222, rv = 2.5),
		chr3 = list(s = numeric(0), e = numeric(0), rv = numeric(0))
	)
	f <- tempfile()
	write.broadpeak.info(bp, f)
	bpf <- read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
	expect_equal(ncol(bpf), 9)
	expect_equal(nrow(bpf), 3)
	expect_equal(bpf$V1, c("chr1", "chr2", "chr1"))
	expect_equal(bpf$V2, c(5123, 7111, 1111))
	expect_equal(bpf$V3, c(6234, 8222, 2222))
	expect_equal(bpf$V4, rep(".", 3))
	expect_equal(bpf$V5, rep(0, 3))
	expect_equal(bpf$V6, rep(".", 3))
	expect_equal(bpf$V7, c(3.25, 2.5, 1.5))
	expect_equal(bpf$V8, rep(-1, 3))
	expect_equal(bpf$V9, rep(-1, 3))
})

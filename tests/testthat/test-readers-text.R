# Tests of the text-format tag readers: bowtie, tagalign, arachne (short and
# long), helicos, maq text map, meland. Expected values are worked out by hand.
# For the arachne, helicos, meland and maq column layouts and the mismatch and
# length rules the C++ source (src/bed2vector.cpp) is the specification, so
# those expectations are not independent of it.
# Conventions: positive tag = forward strand at its 5' coordinate, negative tag
# = reverse strand, stored as -(pos + len - 1) where the format gives a start.
# The R wrappers return the tags of each chromosome sorted by abs(position).

# ---------------------------------------------------------------- bowtie
# columns: name, strand, chr, pos, seq, qual, count, mismatch list
bowtie_lines <- function() {
	c(
		tabline("r1", "+", "chr1.fa", "100", "ACGTACGTAC", "IIIIIIIIII", "0", ""),
		tabline("r2", "-", "chr1.fa", "50", "ACGTACGTAC", "IIIIIIIIII", "3", "3:A>T,5:C>G"),
		tabline("r3", "+", "chr2", "7", "ACGT", "IIII", "0", "1:A>T")
	)
}

test_that("read.bowtie.tags parses positions, strands and mismatch counts", {
	f <- write_tmp(bowtie_lines())
	r <- quiet(read.bowtie.tags(f))
	expect_named(r, c("tags", "quality"))
	expect_equal(sort(names(r$tags)), c("chr1.fa", "chr2"))
	# r2 is on the minus strand: -(50 + 10 - 1) = -59; sorted by abs value
	expect_equal(r$tags[["chr1.fa"]], c(-59L, 100L))
	# mismatches: r2 has two comma-separated entries, r1 an empty field
	expect_equal(r$quality[["chr1.fa"]], c(2L, 0L))
	expect_equal(r$tags[["chr2"]], 7L)
	expect_equal(r$quality[["chr2"]], 1L)
})

test_that("read.bowtie.tags returns tag names in sorted order when requested", {
	f <- write_tmp(bowtie_lines())
	r <- quiet(read.bowtie.tags(f, read.tag.names = TRUE))
	expect_named(r, c("tags", "quality", "names"))
	expect_equal(r$names[["chr1.fa"]], c("r2", "r1"))
	expect_equal(r$names[["chr2"]], "r3")
})

test_that("read.bowtie.tags fix.chromosome.names strips .fa", {
	f <- write_tmp(bowtie_lines())
	expect_equal(sort(names(quiet(read.bowtie.tags(f))$tags)), c("chr1.fa", "chr2"))
	r <- quiet(read.bowtie.tags(f, fix.chromosome.names = TRUE))
	expect_equal(sort(names(r$tags)), c("chr1", "chr2"))
	expect_equal(r$tags$chr1, c(-59L, 100L))
})

test_that("read.bowtie.tags handles empty and missing files", {
	r <- quiet(read.bowtie.tags(write_tmp(character(0))))
	expect_named(r, c("tags", "quality"))
	expect_length(r$tags, 0)
	expect_length(r$quality, 0)
	r <- quiet(read.bowtie.tags(file.path(tempdir(), "no-such-file.bowtie")))
	expect_length(r$tags, 0)
})

test_that("read.bowtie.tags keeps every line of a larger file", {
	n <- 100
	pos <- seq_len(n) * 3
	lines <- vapply(seq_len(n), function(i) {
		tabline(paste0("r", i), "+", paste0("chr", 1 + i %% 3), pos[i], "ACGT", "IIII", "0", "")
	}, "")
	r <- quiet(read.bowtie.tags(write_tmp(lines)))
	expect_equal(sort(names(r$tags)), c("chr1", "chr2", "chr3"))
	expect_equal(sum(lengths(r$tags)), n)
	# chr1 gets i with i %% 3 == 0, positions 3 * i
	expect_equal(r$tags$chr1, as.integer(pos[seq_len(n) %% 3 == 0]))
	for (chr in names(r$tags)) expect_false(is.unsorted(r$tags[[chr]]))
})

# ---------------------------------------------------------------- tagalign
# columns: chr, start, end, sequence, score, strand
tagalign_lines <- function(scores = c(0, 1, 2, 0)) {
	c(
		tabline("chr1.fa", 100, 135, "N", scores[1], "+"),
		tabline("chr1.fa", 200, 235, "N", scores[2], "-"),
		tabline("chr1.fa", 50, 85, "N", scores[3], "+"),
		tabline("chr2.fa", 10, 45, "N", scores[4], "-")
	)
}

test_that("read.tagalign.tags uses start for + and -end for - strand tags", {
	f <- write_tmp(tagalign_lines())
	r <- quiet(read.tagalign.tags(f, fix.quality = FALSE))
	expect_named(r, c("tags", "quality"))
	expect_equal(names(r$tags), c("chr1", "chr2"))
	expect_equal(r$tags$chr1, c(50L, 100L, -235L))
	expect_equal(r$tags$chr2, -45L)
	# raw scores follow the tags through the sort
	expect_equal(r$quality$chr1, c(2L, 0L, 1L))
	expect_equal(r$quality$chr2, 0L)
})

test_that("read.tagalign.tags fix.chromosome.names = FALSE keeps .fa", {
	f <- write_tmp(tagalign_lines())
	r <- quiet(read.tagalign.tags(f, fix.chromosome.names = FALSE, fix.quality = FALSE))
	expect_equal(names(r$tags), c("chr1.fa", "chr2.fa"))
})

test_that("read.tagalign.tags fix.quality converts mismatch-like scores into ordered bins", {
	# Per chromosome: scores containing 0 are mapped to ceiling(1000 / 4^score) =
	# 1000, 250, 63 for scores 0, 1, 2; the three distinct values become bins
	# numbered from the highest value down, so score 0 -> bin 0, 1 -> 1, 2 -> 2.
	f <- write_tmp(tagalign_lines())
	r <- quiet(read.tagalign.tags(f))
	expect_equal(r$tags$chr1, c(50L, 100L, -235L))
	expect_equal(as.integer(r$quality$chr1), c(2L, 0L, 1L))
	# a chromosome with a single distinct value gets bin 0
	expect_equal(as.integer(r$quality$chr2), 0L)
})

test_that("read.tagalign.tags fix.quality handles mapq-like scores (no zeros)", {
	# scores 40, 10, 40: breaks (0, 10, 40], bins 2, 1, 2 -> 2 - bin = 0, 1, 0
	f <- write_tmp(c(
		tabline("chr1", 100, 135, "N", 40, "+"),
		tabline("chr1", 200, 235, "N", 10, "+"),
		tabline("chr1", 300, 335, "N", 40, "+")
	))
	r <- quiet(read.tagalign.tags(f))
	expect_equal(r$tags$chr1, c(100L, 200L, 300L))
	expect_equal(as.integer(r$quality$chr1), c(0L, 1L, 0L))
})

test_that("read.tagalign.tags accepts space separated fields", {
	f <- write_tmp(c("chr1 100 135 N 0 +", "chr1 200 235 N 0 -"))
	r <- quiet(read.tagalign.tags(f, fix.quality = FALSE))
	expect_equal(r$tags$chr1, c(100L, -235L))
})

test_that("read.tagalign.tags reads files with more than 30 lines completely", {
	n <- 45
	lines <- vapply(seq_len(n), function(i) tabline("chr1", i * 10, i * 10 + 35, "N", 0, "+"), "")
	r <- quiet(read.tagalign.tags(write_tmp(lines), fix.quality = FALSE))
	expect_equal(r$tags$chr1, as.integer(seq_len(n) * 10))
})

test_that("read.tagalign.tags returns empty lists for an empty file", {
	r <- quiet(read.tagalign.tags(write_tmp(character(0))))
	expect_named(r, c("tags", "quality"))
	expect_length(r$tags, 0)
	expect_length(r$quality, 0)
})

# ---------------------------------------------------------------- short arachne
# columns: chr, position, optional mismatch count. There is no strand column,
# the position is used as given.
test_that("read.short.arachne.tags reads position and optional mismatches", {
	f <- write_tmp(c("chr1.fa 300 1", "chr1.fa 100", "chr2.fa -40 0"))
	r <- quiet(read.short.arachne.tags(f))
	expect_named(r, c("tags", "quality"))
	# default fix.chromosome.names = FALSE
	expect_equal(names(r$tags), c("chr1.fa", "chr2.fa"))
	expect_equal(r$tags[["chr1.fa"]], c(100L, 300L))
	expect_equal(r$quality[["chr1.fa"]], c(0L, 1L))
	expect_equal(r$tags[["chr2.fa"]], -40L)
	r2 <- quiet(read.short.arachne.tags(f, fix.chromosome.names = TRUE))
	expect_equal(names(r2$tags), c("chr1", "chr2"))
})

test_that("read.short.arachne.tags accepts tab separators and empty input", {
	r <- quiet(read.short.arachne.tags(write_tmp(c("chrA\t5\t2\n"))))
	expect_equal(r$tags$chrA, 5L)
	expect_equal(r$quality$chrA, 2L)
	r <- quiet(read.short.arachne.tags(write_tmp(character(0))))
	expect_named(r, c("tags", "quality"))
	expect_length(r$tags, 0)
})

test_that("read.short.arachne.tags reads more than 30 lines", {
	n <- 40
	r <- quiet(read.short.arachne.tags(write_tmp(paste("chr1", seq_len(n) * 7, seq_len(n) %% 4))))
	expect_equal(r$tags$chr1, as.integer(seq_len(n) * 7))
	expect_equal(r$quality$chr1, as.integer(seq_len(n) %% 4))
})

# Files without a final newline: only values that are the same whether or not
# the last character of the last line is kept are asserted.
test_that("tag positions survive a missing final newline", {
	f <- tempfile()
	cat("chr1.fa 100 5\nchr1.fa 300 7", file = f)
	expect_equal(quiet(read.short.arachne.tags(f))$tags[["chr1.fa"]], c(100L, 300L))
	f <- tempfile()
	cat(tabline("r1", "+", "chr1", "10", "ACGT", "IIII", "0", ""), "\n",
		tabline("r2", "+", "chr1", "20", "ACGT", "IIII", "0", "1:A>T,2:C>G"), sep = "", file = f)
	expect_equal(quiet(read.bowtie.tags(f))$tags$chr1, c(10L, 20L))
})

# ---------------------------------------------------------------- arachne (long)
# QUERY lines: QUERY, 4 skipped fields, strand (1 = reverse), chr, start, end,
# 1 skipped field, number of blocks, then (gaps, length, mismatches) per block.
# length = sum of block lengths; mismatches = sum(|gaps| + mm) + blocks - 1.
arachne_lines <- function() {
	c(
		"HEADER should be ignored",
		"QUERY a b c d 0 chr1.fa 100 150 x 1 0 50 2",
		"QUERY a b c d 1 chr1.fa 300 360 x 2 0 30 1 -1 30 0",
		"QUERY a b c d 0 chr2.fa 10 30 x 1 0 20 0"
	)
}

test_that("read.arachne.tags computes position, mismatches and length", {
	r <- quiet(read.arachne.tags(write_tmp(arachne_lines())))
	expect_named(r, c("tags", "quality", "length"))
	expect_equal(names(r$tags), c("chr1.fa", "chr2.fa"))
	# the reverse-strand read is placed at -end
	expect_equal(r$tags[["chr1.fa"]], c(100L, -360L))
	# first read: |0| + 2 + (1 - 1) = 2; second: (0 + 1) + (1 + 0) + (2 - 1) = 3
	expect_equal(r$quality[["chr1.fa"]], c(2L, 3L))
	expect_equal(r$length[["chr1.fa"]], c(50L, 60L))
	expect_equal(r$tags[["chr2.fa"]], 10L)
	expect_equal(r$quality[["chr2.fa"]], 0L)
	expect_equal(r$length[["chr2.fa"]], 20L)
})

test_that("read.arachne.tags fix.chromosome.names and empty input", {
	r <- quiet(read.arachne.tags(write_tmp(arachne_lines()), fix.chromosome.names = TRUE))
	expect_equal(names(r$tags), c("chr1", "chr2"))
	r <- quiet(read.arachne.tags(write_tmp(c("HEADER only"))))
	expect_named(r, c("tags", "quality", "length"))
	expect_length(r$tags, 0)
	expect_length(r$length, 0)
})

test_that("read.arachne.tags reads more than 30 QUERY lines", {
	n <- 40
	lines <- vapply(seq_len(n), function(i) sprintf("QUERY a b c d 0 chr1 %d %d x 1 0 25 1", i * 10, i * 10 + 25), "")
	r <- quiet(read.arachne.tags(write_tmp(lines)))
	expect_equal(r$tags$chr1, as.integer(seq_len(n) * 10))
	expect_equal(r$length$chr1, rep(25L, n))
	expect_equal(r$quality$chr1, rep(1L, n))
})

# ---------------------------------------------------------------- helicos
# columns: chr, name, start, end, tstart, tend, 2 skipped, ndel, nins, nsub, strand
# length = tend - tstart; mismatches = ndel + nins + nsub; minus strand at -end.
helicos_lines <- function() {
	c(
		"# a comment line",
		tabline("Reference_ID", "Read_ID", "Start", "End", "TStart", "TEnd", "a", "b", "Del", "Ins", "Sub", "Strand"),
		tabline("chr1.fa", "r1", 100, 150, 1, 51, "x", "x", 0, 1, 1, "+"),
		tabline("chr1.fa", "r2", 300, 340, 1, 41, "x", "x", 1, 0, 0, "-"),
		tabline("chr2.fa", "r3", 5, 25, 3, 23, "x", "x", 0, 0, 0, "+")
	)
}

test_that("read.helicos.tags parses coordinates, quality and length", {
	r <- quiet(read.helicos.tags(write_tmp(helicos_lines())))
	expect_named(r, c("tags", "quality", "length"))
	expect_equal(names(r$tags), c("chr1.fa", "chr2.fa"))
	expect_equal(r$tags[["chr1.fa"]], c(100L, -340L))
	expect_equal(r$quality[["chr1.fa"]], c(2L, 1L))
	expect_equal(r$length[["chr1.fa"]], c(50L, 40L))
	expect_equal(r$tags[["chr2.fa"]], 5L)
	expect_equal(r$quality[["chr2.fa"]], 0L)
	expect_equal(r$length[["chr2.fa"]], 20L)
})

test_that("read.helicos.tags tag names and chromosome name options", {
	f <- write_tmp(helicos_lines())
	r <- quiet(read.helicos.tags(f, read.tag.names = TRUE, fix.chromosome.names = TRUE))
	expect_named(r, c("tags", "quality", "length", "names"))
	expect_equal(names(r$tags), c("chr1", "chr2"))
	expect_equal(r$names$chr1, c("r1", "r2"))
	expect_equal(r$names$chr2, "r3")
})

test_that("read.helicos.tags returns nothing for header-only or empty files", {
	r <- quiet(read.helicos.tags(write_tmp(helicos_lines()[1:2])))
	expect_length(r$tags, 0)
	r <- quiet(read.helicos.tags(write_tmp(character(0))))
	expect_named(r, c("tags", "quality", "length"))
	expect_length(r$tags, 0)
})

test_that("read.helicos.tags reads more than 30 lines", {
	n <- 40
	lines <- vapply(seq_len(n), function(i) tabline("chr1", paste0("r", i), i * 10, i * 10 + 20, 1, 21, "x", "x", 0, 0, i %% 3, "+"), "")
	r <- quiet(read.helicos.tags(write_tmp(lines)))
	expect_equal(r$tags$chr1, as.integer(seq_len(n) * 10))
	expect_equal(r$quality$chr1, as.integer(seq_len(n) %% 3))
	expect_equal(r$length$chr1, rep(20L, n))
})

# ---------------------------------------------------------------- maq text map
# columns (0-based): 0 name, 1 chr, 2 pos, 3 strand, 4-8 skipped, 9 mismatches,
# 10-12 skipped, 13 read length.
maq_line <- function(name, chr, pos, strand, nm, len) {
	tabline(name, chr, pos, strand, "a", "b", "c", "d", "e", nm, "f", "g", "h", len)
}

test_that("read.maqmap.tags reads position, strand, mismatches and read length", {
	f <- write_tmp(c(
		maq_line("r1", "chr1.fa", 100, "+", 1, 36),
		maq_line("r2", "chr1.fa", 50, "-", 0, 36),
		maq_line("r3", "chr2.fa", 7, "+", 3, 20)
	))
	r <- quiet(read.maqmap.tags(f))
	expect_named(r, c("tags", "quality"))
	# default fix.chromosome.names = TRUE; -(50 + 36 - 1) = -85
	expect_equal(names(r$tags), c("chr1", "chr2"))
	expect_equal(r$tags$chr1, c(-85L, 100L))
	expect_equal(r$quality$chr1, c(0L, 1L))
	expect_equal(r$tags$chr2, 7L)
	expect_equal(r$quality$chr2, 3L)
	r <- quiet(read.maqmap.tags(f, read.tag.names = TRUE, fix.chromosome.names = FALSE))
	expect_equal(names(r$tags), c("chr1.fa", "chr2.fa"))
	expect_equal(r$names[["chr1.fa"]], c("r2", "r1"))
})

test_that("read.maqmap.tags returns empty lists for an empty file", {
	r <- quiet(read.maqmap.tags(write_tmp(character(0))))
	expect_named(r, c("tags", "quality"))
	expect_length(r$tags, 0)
	expect_length(r$quality, 0)
})

test_that("read.maqmap.tags reads more than 30 lines", {
	n <- 40
	lines <- vapply(seq_len(n), function(i) maq_line(paste0("r", i), "chr1", i * 10, "+", i %% 4, 36), "")
	r <- quiet(read.maqmap.tags(write_tmp(lines)))
	expect_equal(r$tags$chr1, as.integer(seq_len(n) * 10))
	expect_equal(r$quality$chr1, as.integer(seq_len(n) %% 4))
})

# ---------------------------------------------------------------- meland
# columns: 0 name, 1 skipped, 2 U<mismatches>, 3-5 skipped, 6 match length,
# 7 chromosome (first three characters are dropped), 8 position. Lines whose
# third field does not start with U are ignored. Quality is
# (max length over all tags) - length + mismatches / 10.
meland_lines <- function() {
	c(
		"r1 x U0 a b c 30 chr1.fa 100",
		"r2 x U1 a b c 25 chr1.fa 50",
		"r3 x U2 a b c 30 chr2.fa 10",
		"r4 x NM a b c 30 chr1.fa 5"
	)
}

test_that("read.meland.tags combines length and mismatches into the quality score", {
	r <- quiet(read.meland.tags(write_tmp(meland_lines())))
	expect_named(r, c("tags", "quality"))
	# 'chr' prefix dropped, then .fa removed by the default fix.chromosome.names
	expect_equal(names(r$tags), c("1", "2"))
	expect_equal(r$tags[["1"]], c(50L, 100L))
	expect_equal(r$tags[["2"]], 10L)
	# max length = 30: r2 -> 5 + 0.1, r1 -> 0 + 0, r3 -> 0 + 0.2
	expect_equal(as.numeric(r$quality[["1"]]), c(5.1, 0))
	expect_equal(as.numeric(r$quality[["2"]]), 0.2)
})

test_that("read.meland.tags options: names and fix.chromosome.names", {
	r <- quiet(read.meland.tags(write_tmp(meland_lines()), read.tag.names = TRUE,
		fix.chromosome.names = FALSE))
	expect_named(r, c("tags", "quality", "names"))
	expect_equal(names(r$tags), c("1.fa", "2.fa"))
	expect_equal(r$names[["1.fa"]], c("r2", "r1"))
})

test_that("read.meland.tags tab separated input is accepted", {
	f <- write_tmp(gsub(" ", "\t", meland_lines()))
	r <- quiet(read.meland.tags(f))
	expect_equal(r$tags[["1"]], c(50L, 100L))
})

test_that("read.meland.tags reads more than 30 lines", {
	n <- 40
	lines <- sprintf("r%d x U%d a b c 30 chr1.fa %d", seq_len(n), seq_len(n) %% 3, seq_len(n) * 10)
	r <- quiet(read.meland.tags(write_tmp(lines)))
	expect_equal(r$tags[["1"]], as.integer(seq_len(n) * 10))
	expect_equal(as.numeric(r$quality[["1"]]), (seq_len(n) %% 3) / 10)
})

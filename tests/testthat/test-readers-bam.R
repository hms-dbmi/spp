# Tests of read.bam.tags (needs Rsamtools). BAM files are built from SAM text
# with Rsamtools::asBam.
# Single-end: forward reads are placed at pos, reverse reads at a negative
# coordinate. Paired-end (any read has flag bit 0x1): only forward-strand reads
# are kept, tags = pos and flen = abs(isize).

sam_file <- function(reads, refs = c(chr1 = 1000, chr2 = 1000)) {
	hdr <- c("@HD\tVN:1.0\tSO:unsorted",
		paste0("@SQ\tSN:", names(refs), "\tLN:", refs))
	sam <- tempfile(fileext = ".sam")
	writeLines(c(hdr, vapply(reads, function(r) paste(r, collapse = "\t"), "")), sam)
	bam <- Rsamtools::asBam(sam, destination = tempfile(), indexDestination = FALSE,
		overwrite = TRUE)
	unname(bam)
}

# qname flag rname pos mapq cigar rnext pnext tlen seq qual
sam_read <- function(name, flag, chr, pos, mapq, len, tlen = 0, mate = "*", mpos = 0) {
	c(name, flag, chr, pos, mapq, paste0(len, "M"), mate, mpos, tlen,
		strrep("A", len), strrep("I", len))
}

test_that("read.bam.tags reads single-end forward reads and mapq", {
	skip_if_not_installed("Rsamtools")
	bam <- sam_file(list(
		sam_read("a", 0, "chr1", 300, 30, 20),
		sam_read("b", 0, "chr1", 100, 10, 20),
		sam_read("c", 0, "chr2", 50, 60, 20)
	))
	r <- quiet(read.bam.tags(bam))
	expect_named(r, c("tags", "quality"))
	expect_equal(names(r$tags), c("chr1", "chr2"))
	# tags follow the order of the file (no sorting in read.bam.tags)
	expect_equal(plain(r$tags$chr1), c(300, 100))
	expect_equal(plain(r$quality$chr1), c(30, 10))
	expect_equal(plain(r$tags$chr2), 50)
	expect_equal(plain(r$quality$chr2), 60)
})

test_that("read.bam.tags places reverse-strand reads at negative coordinates", {
	skip_if_not_installed("Rsamtools")
	bam <- sam_file(list(
		sam_read("a", 0, "chr1", 100, 30, 20),
		sam_read("b", 16, "chr1", 200, 30, 20)
	))
	r <- quiet(read.bam.tags(bam))
	t1 <- plain(r$tags$chr1)
	expect_length(t1, 2)
	expect_equal(t1[1], 100)
	# only the sign of the reverse read is checked, its exact coordinate is
	# convention dependent (see BUGS-readers.txt, item 5a)
	expect_lt(t1[2], 0)
})

test_that("read.bam.tags drops unmapped reads and returns names when asked", {
	skip_if_not_installed("Rsamtools")
	bam <- sam_file(list(
		sam_read("a", 0, "chr1", 100, 30, 20),
		sam_read("u", 4, "*", 0, 0, 20),
		sam_read("b", 0, "chr2", 10, 40, 20),
		sam_read("c", 0, "chr1", 150, 40, 20),
		sam_read("d", 0, "chr2", 60, 40, 20)
	))
	r <- quiet(read.bam.tags(bam, read.tag.names = TRUE))
	expect_named(r, c("tags", "quality", "names"))
	expect_equal(names(r$tags), c("chr1", "chr2"))
	# unmapped read "u" is excluded; names follow the file order per chromosome
	expect_equal(as.character(r$names$chr1), c("a", "c"))
	expect_equal(as.character(r$names$chr2), c("b", "d"))
	expect_equal(plain(r$tags$chr2), c(10, 60))
})

test_that("read.bam.tags keeps forward reads and fragment lengths for paired-end data", {
	skip_if_not_installed("Rsamtools")
	# flags 99 / 147: first read forward (pair flag 0x1 set), mate reverse
	bam <- sam_file(list(
		sam_read("p1", 99, "chr1", 100, 30, 20, tlen = 150, mate = "=", mpos = 230),
		sam_read("p1", 147, "chr1", 230, 30, 20, tlen = -150, mate = "=", mpos = 100),
		sam_read("p2", 99, "chr1", 400, 30, 20, tlen = 200, mate = "=", mpos = 580)
	))
	r <- quiet(read.bam.tags(bam))
	expect_true(all(c("tags", "flen") %in% names(r)))
	expect_equal(plain(r$tags$chr1), c(100, 400))
	expect_equal(plain(r$flen$chr1), c(150, 200))
})

test_that("read.bam.tags returns empty elements for a BAM without mapped reads", {
	skip_if_not_installed("Rsamtools")
	bam <- sam_file(list(sam_read("u", 4, "*", 0, 0, 20)))
	r <- quiet(read.bam.tags(bam))
	expect_named(r, c("tags", "quality"))
	expect_length(r$tags, 0)
	expect_length(r$quality, 0)
})

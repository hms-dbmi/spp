# Tests of read.bin.maqmap.tags on hand-built binary maq map files (gzip).
# Layout used by src/maqmap.h with MAQ_LONGREADS (set in src/Makevars):
#   header: int32 format (-1), int32 n_ref, per reference (int32 len incl NUL,
#   name bytes with NUL), uint64 n_mapped_reads
#   record (184 bytes): 128 byte seq, then size, map_qual, info1, info2, c[2],
#   flag, alt_qual (8 bytes), uint32 seqid, uint32 pos, int32 dist, char name[36].
# Dependencies: src/Makevars passes -DMAQ_LONGREADS (MAX_READLEN 128), giving
# 184-byte records with 4-byte ints; the files are written little-endian, so
# the tests are skipped on big-endian platforms.
# pos = 2 * (0-based position) + strand bit; mismatches = info1 & 0xf.

maq_record <- function(name, seqid, pos0, reverse, size, info1) {
	nm <- raw(36)
	nb <- charToRaw(name)
	nm[seq_along(nb)] <- nb
	c(
		raw(128),
		as.raw(c(size, 0, info1, 0, 0, 0, 0, 0)),
		writeBin(as.integer(seqid), raw(), size = 4, endian = "little"),
		writeBin(as.integer(2 * pos0 + reverse), raw(), size = 4, endian = "little"),
		writeBin(0L, raw(), size = 4, endian = "little"),
		nm
	)
}

write_maq <- function(refs, records) {
	f <- tempfile(fileext = ".map")
	con <- gzfile(f, "wb")
	on.exit(close(con))
	int <- function(x) writeBin(as.integer(x), raw(), size = 4, endian = "little")
	writeBin(int(-1), con)
	writeBin(int(length(refs)), con)
	for (r in refs) {
		writeBin(int(nchar(r) + 1), con)
		writeBin(c(charToRaw(r), as.raw(0)), con)
	}
	writeBin(int(c(length(records), 0)), con)
	for (rec in records) writeBin(rec, con)
	f
}

maq_example <- function() {
	write_maq(c("chr1.fa", "chr2.fa", "chr3.fa"), list(
		maq_record("rA", 0, 99, 0, 36, 0x21),
		maq_record("rB", 0, 49, 1, 36, 0x00),
		maq_record("rC", 1, 6, 0, 20, 0x03)
	))
}

test_that("read.bin.maqmap.tags decodes position, strand and mismatches", {
	skip_if(.Platform$endian != "little")
	r <- quiet(read.bin.maqmap.tags(maq_example()))
	expect_named(r, c("tags", "quality"))
	# default fix.chromosome.names = TRUE; references without reads are absent
	expect_equal(names(r$tags), c("chr1", "chr2"))
	# rA: 0-based 99 -> 100; rB reverse: -(50 + 36 - 1) = -85
	expect_equal(r$tags$chr1, c(-85L, 100L))
	# info1 & 0xf: rB 0, rA 1 (the high nibble of 0x21 is ignored)
	expect_equal(r$quality$chr1, c(0L, 1L))
	expect_equal(r$tags$chr2, 7L)
	expect_equal(r$quality$chr2, 3L)
})

test_that("read.bin.maqmap.tags tag names and chromosome name options", {
	skip_if(.Platform$endian != "little")
	r <- quiet(read.bin.maqmap.tags(maq_example(), read.tag.names = TRUE,
		fix.chromosome.names = FALSE))
	expect_named(r, c("tags", "quality", "names"))
	expect_equal(names(r$tags), c("chr1.fa", "chr2.fa"))
	expect_equal(r$names[["chr1.fa"]], c("rB", "rA"))
	expect_equal(r$names[["chr2.fa"]], "rC")
})

test_that("read.bin.maqmap.tags returns empty lists for a map without reads", {
	skip_if(.Platform$endian != "little")
	f <- write_maq(c("chr1.fa"), list())
	r <- quiet(read.bin.maqmap.tags(f))
	expect_named(r, c("tags", "quality"))
	expect_length(r$tags, 0)
	expect_length(r$quality, 0)
})

test_that("read.bin.maqmap.tags keeps reads on one chromosome in sorted order", {
	skip_if(.Platform$endian != "little")
	recs <- lapply(c(30, 10, 20), function(p) maq_record("r", 0, p - 1, 0, 25, 0))
	r <- quiet(read.bin.maqmap.tags(write_maq("c1", recs)))
	expect_equal(r$tags$c1, c(10L, 20L, 30L))
})

test_that("read.bin.maqmap.tags reads more than 30 records", {
	skip_if(.Platform$endian != "little")
	n <- 45
	recs <- lapply(seq_len(n), function(i) maq_record("r", i %% 2, i * 10 - 1, 0, 25, i %% 5))
	r <- quiet(read.bin.maqmap.tags(write_maq(c("c0", "c1"), recs)))
	i1 <- seq_len(n)[seq_len(n) %% 2 == 1]
	i0 <- seq_len(n)[seq_len(n) %% 2 == 0]
	expect_equal(r$tags$c1, as.integer(i1 * 10))
	expect_equal(r$tags$c0, as.integer(i0 * 10))
	expect_equal(r$quality$c1, as.integer(i1 %% 5))
})

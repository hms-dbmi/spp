# get.conservative.fold.enrichment.profile and get.conservative.fold.enrichment.profile2
# The profile keeps the lower confidence bound where it is above 1 (log2 > 0),
# the upper bound where it is below 1 (log2 < 0) and 0 elsewhere.

# signal with an extra pile over 1000-1100 on a uniform control (and the reverse)
pl.prof.data <- function() {
	flat <- as.integer(seq(5, 2000, by = 5))
	pile <- as.integer(round(seq(1000, 1100, length.out = 150)))
	sg <- c(flat, pile)
	sg <- sg[order(sg)]
	sg <- as.integer(sg * rep(c(1L, -1L), length.out = length(sg)))
	sg <- sg[order(abs(sg))]
	list(sig = list(a = sg), ctl = list(a = as.integer(seq(3, 2000, by = 5))), flat = flat)
}

pl.prof.ref.mle <- function(sg, ct, x, shift, fws, bws, bgw) {
	ft <- abs(sg + shift)
	bt <- abs(ct + shift)
	log2((pl.count(ft, x, fws) + 0.5) / (pl.count(bt, x, bws) + 0.5) * bws / fws / bgw)
}

test_that("grid profile mle matches a brute-force window count", {
	d <- pl.prof.data()
	r <- get.conservative.fold.enrichment.profile(d$sig, d$ctl, fws = 100, bwsl = 500, step = 50, tag.shift = 10, background.density.scaling = FALSE, return.mle = TRUE)
	expect_equal(names(r), "a")
	expect_equal(colnames(r$a), c("x", "y", "mle", "lb", "ub"))
	expect_true(all(diff(r$a$x) == 50))
	bgw <- length(d$sig$a) / length(d$ctl$a)
	ref <- pl.prof.ref.mle(d$sig$a, d$ctl$a, r$a$x, 10, 100, 500, bgw)
	expect_equal(r$a$mle, ref)
})

test_that("profile values follow the lower/upper-bound composition rule", {
	d <- pl.prof.data()
	r <- get.conservative.fold.enrichment.profile(d$sig, d$ctl, fws = 100, bwsl = 500, step = 50, tag.shift = 10, background.density.scaling = FALSE, return.mle = TRUE)$a
	expect_true(all(r$lb <= r$mle & r$mle <= r$ub))
	up <- r$lb > 0
	down <- r$ub < 0
	expect_true(any(up))
	expect_false(any(up & down))
	expect_equal(r$y[up], r$lb[up])
	expect_equal(r$y[down], r$ub[down])
	expect_true(all(r$y[!up & !down] == 0))
	# only the pile (1000-1100, shifted by 10) is significantly enriched
	expect_true(all(r$x[up] >= 950 & r$x[up] <= 1150))
	expect_true(any(r$x[up] >= 1000 & r$x[up] <= 1100))
	# without return.mle only x and y come back, with the same values
	r2 <- get.conservative.fold.enrichment.profile(d$sig, d$ctl, fws = 100, bwsl = 500, step = 50, tag.shift = 10, background.density.scaling = FALSE)$a
	expect_equal(colnames(r2), c("x", "y"))
	expect_equal(r2$y, r$y)
})

test_that("a depleted window gets the (negative) upper bound", {
	flat <- as.integer(seq(5, 2000, by = 5))
	hole <- list(a = flat[flat < 1000 | flat > 1100])
	ctl <- list(a = as.integer(seq(3, 2000, by = 2)))
	r <- get.conservative.fold.enrichment.profile(hole, ctl, fws = 100, bwsl = 500, step = 50, tag.shift = 0, background.density.scaling = FALSE, return.mle = TRUE)$a
	down <- r$ub < 0
	expect_true(any(down))
	expect_true(all(r$x[down] >= 950 & r$x[down] <= 1150))
	expect_true(all(r$y[down] < 0))
	expect_equal(r$y[down], r$ub[down])
	expect_true(all(r$lb[down] < 0))
})

test_that("a dataset against itself with one scale is exactly 1-fold", {
	x <- as.integer(round(seq(10, 2990, length.out = 300)) * rep(c(1L, -1L), 150))
	s <- list(a = x[order(abs(x))])
	r <- get.conservative.fold.enrichment.profile(s, s, fws = 100, bwsl = 100, step = 50, background.density.scaling = FALSE, return.mle = TRUE)$a
	expect_true(all(r$mle == 0))
	expect_true(all(r$y == 0))
	expect_true(all(r$lb <= 0 & r$ub >= 0))
})

test_that("explicit positions replace the grid", {
	d <- pl.prof.data()
	g <- get.conservative.fold.enrichment.profile(d$sig, d$ctl, fws = 100, bwsl = 500, step = 50, tag.shift = 10, background.density.scaling = FALSE, return.mle = TRUE)$a
	# do not assume where the grid starts: take the grid points nearest to the positions of interest
	near <- function(p) g$x[which.min(abs(g$x - p))]
	pp <- c(near(500), near(1050), near(1500))
	r <- get.conservative.fold.enrichment.profile(d$sig, d$ctl, fws = 100, bwsl = 500, tag.shift = 10, background.density.scaling = FALSE, return.mle = TRUE, posl = list(a = pp))
	expect_equal(names(r), "a")
	expect_equal(r$a$x, pp)
	bgw <- length(d$sig$a) / length(d$ctl$a)
	expect_equal(r$a$mle, pl.prof.ref.mle(d$sig$a, d$ctl$a, pp, 10, 100, 500, bgw))
	expect_gt(r$a$y[2], 0)
	expect_equal(r$a$y[c(1, 3)], c(0, 0))
	# the same positions on the grid and as explicit positions agree
	expect_equal(g$mle[match(pp, g$x)], r$a$mle)
})

test_that("posl restricts the chromosomes reported", {
	d <- pl.prof.data()
	d$sig$z <- d$sig$a
	d$ctl$z <- d$ctl$a
	r <- get.conservative.fold.enrichment.profile(d$sig, d$ctl, fws = 100, bwsl = 500, tag.shift = 10, background.density.scaling = FALSE, posl = list(z = c(600, 1050)))
	expect_equal(names(r), "z")
	expect_equal(r$z$x, c(600, 1050))
})

test_that("chromosomes with at most two tags are dropped", {
	d <- pl.prof.data()
	d$sig$tiny <- c(5L, -50L)
	d$ctl$tiny <- c(3L, 40L, 90L)
	r <- get.conservative.fold.enrichment.profile(d$sig, d$ctl, fws = 100, bwsl = 500, step = 50, tag.shift = 10, background.density.scaling = FALSE)
	expect_equal(names(r), "a")
	d$sig$three <- c(500L, 600L, -700L)
	d$ctl$three <- as.integer(seq(2, 2000, by = 10))
	r <- get.conservative.fold.enrichment.profile(d$sig, d$ctl, fws = 100, bwsl = 500, step = 50, tag.shift = 10, background.density.scaling = FALSE)
	expect_equal(names(r), c("a", "three"))
})

test_that("a chromosome missing from the control is profiled against an empty background", {
	d <- pl.prof.data()
	d$sig$z <- as.integer(seq(7, 1997, by = 20))
	r <- get.conservative.fold.enrichment.profile(d$sig, d$ctl, fws = 100, bwsl = 500, step = 50, tag.shift = 10, background.density.scaling = FALSE, return.mle = TRUE)
	expect_equal(names(r), c("a", "z"))
	# no control tags: the background count is the 0.5 pseudocount, and the weight is
	# the ratio of all signal tags to all control tags
	bgw <- pl.ntags(d$sig) / pl.ntags(d$ctl)
	expect_equal(r$z$mle, pl.prof.ref.mle(d$sig$z, numeric(0), r$z$x, 10, 100, 500, bgw))
})

test_that("bg.weight rescales the mle by a constant log2 offset", {
	d <- pl.prof.data()
	pp <- c(500, 1050, 1500)
	base <- get.conservative.fold.enrichment.profile(d$sig, d$ctl, fws = 100, bwsl = 500, tag.shift = 10, background.density.scaling = FALSE, return.mle = TRUE, posl = list(a = pp))$a
	w <- 2
	scaled <- get.conservative.fold.enrichment.profile(d$sig, d$ctl, fws = 100, bwsl = 500, tag.shift = 10, background.density.scaling = FALSE, return.mle = TRUE, posl = list(a = pp), bg.weight = w)$a
	bgw <- length(d$sig$a) / length(d$ctl$a)
	expect_equal(scaled$mle - base$mle, rep(log2(bgw / w), 3))
})

test_that("alpha controls the width of the interval", {
	d <- pl.prof.data()
	pp <- c(500, 1050, 1500)
	run <- function(alpha) get.conservative.fold.enrichment.profile(d$sig, d$ctl, fws = 100, bwsl = 500, tag.shift = 10, background.density.scaling = FALSE, return.mle = TRUE, posl = list(a = pp), alpha = alpha)$a
	a1 <- run(0.5)
	a2 <- run(0.05)
	a3 <- run(0.01)
	expect_equal(a1$mle, a3$mle)
	expect_true(all(a1$ub - a1$lb < a2$ub - a2$lb) && all(a2$ub - a2$lb < a3$ub - a3$lb))
	expect_true(all(a1$lb >= a3$lb) && all(a1$ub <= a3$ub))
	# on the log scale the interval half-width grows roughly like the normal quantile
	ratio <- (a3$ub - a3$lb) / (a2$ub - a2$lb)
	expect_true(all(abs(ratio / (qnorm(1 - 0.01 / 2) / qnorm(1 - 0.05 / 2)) - 1) < 0.25))
})

test_that("quick.calculation approximates the exact F-based bounds", {
	d <- pl.prof.data()
	pp <- c(500, 1050, 1500)
	quick <- get.conservative.fold.enrichment.profile(d$sig, d$ctl, fws = 100, bwsl = 500, tag.shift = 10, background.density.scaling = FALSE, return.mle = TRUE, posl = list(a = pp), quick.calculation = TRUE)$a
	exact <- get.conservative.fold.enrichment.profile(d$sig, d$ctl, fws = 100, bwsl = 500, tag.shift = 10, background.density.scaling = FALSE, return.mle = TRUE, posl = list(a = pp), quick.calculation = FALSE)$a
	expect_equal(quick$mle, exact$mle)
	expect_true(pl.close(quick$lb, exact$lb, 0.05))
	expect_true(pl.close(quick$ub, exact$ub, 0.05))
})

test_that("several background scales give the most conservative bounds", {
	d <- pl.prof.data()
	pp <- c(500, 1050, 1500)
	args <- list(d$sig, d$ctl, fws = 100, tag.shift = 10, background.density.scaling = FALSE, return.mle = TRUE, posl = list(a = pp))
	s1 <- do.call(get.conservative.fold.enrichment.profile, c(args, list(bwsl = 100)))$a
	s5 <- do.call(get.conservative.fold.enrichment.profile, c(args, list(bwsl = 500)))$a
	both <- do.call(get.conservative.fold.enrichment.profile, c(args, list(bwsl = c(100, 500))))$a
	expect_equal(both$mle, pmin(s1$mle, s5$mle))
	expect_equal(both$lb, pmin(s1$lb, s5$lb))
	expect_equal(both$ub, pmax(s1$ub, s5$ub))
 
})

test_that("the most informative scale is the one with the highest background density", {
	# 40 control tags sit in the 100 bp window around 500 and none elsewhere, so the
	# 100 bp scale has background density 40.5/100 against 40.5/500 for the 500 bp scale
	sig <- list(a = 470:529)
	ctl <- list(a = 460:499)
	args <- list(sig, ctl, fws = 100, bwsl = c(100, 500), tag.shift = 0, bg.weight = 1, return.mle = TRUE, posl = list(a = 500))
	inf <- do.call(get.conservative.fold.enrichment.profile, c(args, list(use.most.informative.scale = TRUE)))$a
	# 60 signal tags in the window; mle = (60 + .5) / (40 + .5) at the 100 bp scale (a wrong 500 bp pick would be 5x)
	expect_equal(inf$mle, log2(60.5 / 40.5))
	zal <- qnorm(0.975)
	fw <- 60.5
	nb <- 40.5
	rte <- fw + nb - 0.25 * zal^2
	dn <- nb - 0.25 * zal^2
	lb <- ((sqrt(fw * nb) - 0.5 * zal * sqrt(rte)) / dn)^2
	expect_equal(inf$lb, log2(lb))
	# the default (conservative over scales) picks the same lower bound here
	cons <- do.call(get.conservative.fold.enrichment.profile, args)$a
	expect_equal(cons$mle, inf$mle)
})

test_that("lb and ub equal the closed-form bounds on a fixed small-count case", {
	# 60 signal and 5 control tags in the window around 500, weight 1, one scale
	sig <- list(a = 461:520)
	ctl <- list(a = c(470L, 480L, 490L, 500L, 510L))
	fw <- 60.5
	nb <- 5.5
	zal <- qnorm(0.975)
	args <- list(sig, ctl, fws = 100, bwsl = 100, tag.shift = 0, bg.weight = 1, return.mle = TRUE, posl = list(a = 500))
	# Wilson-type quick bounds: (sqrt(fw nb) -+ z/2 sqrt(fw + nb - z^2/4)) / (nb - z^2/4), squared
	rte <- fw + nb - 0.25 * zal^2
	dn <- nb - 0.25 * zal^2
	quick.lb <- ((sqrt(fw * nb) - 0.5 * zal * sqrt(rte)) / dn)^2
	quick.ub <- ((sqrt(fw * nb) + 0.5 * zal * sqrt(rte)) / dn)^2
	quick <- do.call(get.conservative.fold.enrichment.profile, c(args, list(quick.calculation = TRUE)))$a
	expect_equal(quick$mle, log2(fw / nb))
	expect_equal(quick$lb, log2(quick.lb))
	expect_equal(quick$ub, log2(quick.ub))
	expect_equal(quick$y, log2(quick.lb))
	# exact path: Clopper-Pearson-type F quantiles of the Poisson ratio
	exact.lb <- fw / nb * qf(0.025, 2 * fw, 2 * nb)
	exact.ub <- fw / nb * qf(0.975, 2 * fw, 2 * nb)
	exact <- do.call(get.conservative.fold.enrichment.profile, c(args, list(quick.calculation = FALSE)))$a
	expect_equal(exact$lb, log2(exact.lb))
	expect_equal(exact$ub, log2(exact.ub))
	expect_equal(exact$y, log2(exact.lb))
})

# The reference below mirrors the estimator in binomial.proportion.ratio.bounds
# (a Beta(a = 1.25, b = 2.5) prior on each proportion and a delta-method variance
# for the log2 ratio); it checks the wiring and the window counts, not the
# statistical model itself.
test_that("profile2 matches the brute-force binomial proportion ratio", {
	d <- pl.prof.data()
	s2 <- list(a = d$flat)
	c2 <- list(a = as.integer(seq(2, 2000, by = 4)))
	shift <- 10
	fws <- 100
	r <- get.conservative.fold.enrichment.profile2(d$sig, s2, d$ctl, c2, fws = fws, bws = fws, step = 50, tag.shift = shift, background.density.scaling = FALSE, return.mle = TRUE)
	expect_equal(names(r), "a")
	expect_equal(colnames(r$a), c("x", "y", "mle", "lb", "ub"))
	ref.theta <- function(x) {
		bg1 <- length(d$sig$a) / length(d$ctl$a)
		bg2 <- length(s2$a) / length(c2$a)
		f1 <- pl.count(abs(d$sig$a + shift), x, fws)
		b1 <- pl.count(abs(d$ctl$a + shift), x, fws)
		f2 <- pl.count(abs(s2$a + shift), x, fws)
		b2 <- pl.count(abs(c2$a + shift), x, fws)
		a <- 1.25
		b <- 2.5
		lt <- log2((f1 + a - 1) / (f1 + b1 + a + b - 2)) - log2((f2 + a - 1) / (f2 + b2 + a + b - 2)) - log2(1 / (1 + 1 / bg1)) + log2(1 / (1 + 1 / bg2))
		vl <- 1 / ((f1 + a - 1) + (f1 + a - 1)^2 / (b1 + b - 1)) + 1 / ((f2 + a - 1) + (f2 + a - 1)^2 / (b2 + b - 1))
		list(mle = lt, v = vl)
	}
	ref <- ref.theta(r$a$x)
	expect_equal(r$a$mle, ref$mle)
	z <- qnorm(0.975)
	expect_equal(r$a$lb, ref$mle - z * sqrt(ref$v))
	expect_equal(r$a$ub, ref$mle + z * sqrt(ref$v))
	# composition rule
	up <- r$a$lb > 0
	down <- r$a$ub < 0
	expect_equal(r$a$y[up], r$a$lb[up])
	expect_equal(r$a$y[down], r$a$ub[down])
	expect_true(all(r$a$y[!up & !down] == 0))
	# the pile is enriched in experiment 1 relative to experiment 2
	expect_true(any(up))
	expect_true(all(r$a$x[up] >= 950 & r$a$x[up] <= 1150))
	# explicit positions
	pp <- c(500, 1050)
	rp <- get.conservative.fold.enrichment.profile2(d$sig, s2, d$ctl, c2, fws = fws, bws = fws, tag.shift = shift, background.density.scaling = FALSE, return.mle = TRUE, posl = list(a = pp))
	expect_equal(rp$a$x, pp)
	expect_equal(rp$a$mle, ref.theta(pp)$mle)
	expect_gt(rp$a$y[2], 0)
	expect_equal(rp$a$y[1], 0)
})

test_that("profile2 of an experiment against itself is 1-fold", {
	d <- pl.prof.data()
	r <- get.conservative.fold.enrichment.profile2(d$sig, d$sig, d$ctl, d$ctl, fws = 100, bws = 100, step = 50, tag.shift = 10, background.density.scaling = FALSE, return.mle = TRUE)$a
	expect_true(all(r$mle == 0))
	expect_true(all(r$y == 0))
	expect_true(all(r$lb < 0 & r$ub > 0))
	# swapping the experiments negates the log ratio
	s2 <- list(a = as.integer(round(seq(3, 2000, length.out = 300))))
	c2 <- list(a = as.integer(round(seq(1, 2000, length.out = 250))))
	ab <- get.conservative.fold.enrichment.profile2(d$sig, s2, d$ctl, c2, fws = 100, bws = 100, step = 50, tag.shift = 10, background.density.scaling = FALSE, return.mle = TRUE)$a
	ba <- get.conservative.fold.enrichment.profile2(s2, d$sig, c2, d$ctl, fws = 100, bws = 100, step = 50, tag.shift = 10, background.density.scaling = FALSE, return.mle = TRUE, posl = list(a = ab$x))$a
	expect_equal(ba$mle, -ab$mle)
	expect_equal(ba$lb, -ab$ub)
})

test_that("profile2 drops chromosomes with at most two tags", {
	d <- pl.prof.data()
	s2 <- list(a = d$flat, tiny = c(5L, 9L))
	d$sig$tiny <- c(5L, 9L)
	c2 <- list(a = as.integer(seq(2, 2000, by = 4)), tiny = c(2L, 4L))
	d$ctl$tiny <- c(2L, 4L)
	r <- get.conservative.fold.enrichment.profile2(d$sig, s2, d$ctl, c2, fws = 100, bws = 100, step = 50, tag.shift = 10, background.density.scaling = FALSE)
	expect_equal(names(r), "a")
})

test_that("profile2 interval width scales exactly with the normal quantile", {
	d <- pl.prof.data()
	s2 <- list(a = d$flat)
	c2 <- list(a = as.integer(seq(2, 2000, by = 4)))
	run <- function(alpha) get.conservative.fold.enrichment.profile2(d$sig, s2, d$ctl, c2, fws = 100, bws = 100, step = 50, tag.shift = 10, background.density.scaling = FALSE, return.mle = TRUE, alpha = alpha)$a
	a1 <- run(0.5)
	a2 <- run(0.05)
	a3 <- run(0.01)
	expect_equal((a1$ub - a1$lb) / qnorm(1 - 0.5 / 2), (a2$ub - a2$lb) / qnorm(1 - 0.05 / 2))
	expect_equal((a3$ub - a3$lb) / qnorm(1 - 0.01 / 2), (a2$ub - a2$lb) / qnorm(1 - 0.05 / 2))
})

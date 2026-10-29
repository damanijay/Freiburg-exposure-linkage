# ============================================================
# 12b_coverage.R - does the 95% CI contain the planted truth 95% of the time?
# ============================================================
library(dplyr); library(sandwich); library(lmtest); library(ggplot2)
set.seed(2026)

final <- read.csv("data_cache/exposure_child_final.csv")
bz    <- as.numeric(scale(final$built))          # same children, same built-up share
truth <- -0.5 / sqrt(1.25)                        # planted slope after rescaling

sim <- t(replicate(500, {
  ses <- as.numeric(scale(-0.5 * bz + rnorm(length(bz))))     # new noise each time
  d   <- data.frame(ses, bz, hex_id = final$hex_id)
  m   <- lm(ses ~ bz, data = d)
  ci_n <- confint(m)["bz", ]
  ci_c <- coefci(m, vcov. = vcovCL(m, cluster = ~hex_id))["bz", ]
  c(est      = unname(coef(m)["bz"]),
    in_naive = unname(truth >= ci_n[1] & truth <= ci_n[2]),
    in_clust = unname(truth >= ci_c[1] & truth <= ci_c[2]))
}))

coverage <- data.frame(
  truth              = round(truth, 3),
  mean_estimate      = round(mean(sim[, "est"]), 3),
  coverage_naive     = round(100 * mean(sim[, "in_naive"]), 1),
  coverage_clust     = round(100 * mean(sim[, "in_clust"]), 1),
  pct_as_far_as_ours = round(100 * mean(sim[, "est"] >= -0.363), 1))
print(coverage)
write.csv(coverage, "outputs/tables/T5b_coverage.csv", row.names = FALSE)

p15 <- ggplot(data.frame(est = sim[, "est"]), aes(est)) +
  geom_histogram(bins = 40, fill = "grey70") +
  geom_vline(xintercept = truth, linewidth = 1) +
  geom_vline(xintercept = -0.363, linetype = 2, colour = "firebrick") +
  labs(title = "500 repeated samples: estimates of the planted SES slope",
       subtitle = "Black = planted truth; red dashed = the one estimate from our cohort",
       x = "Estimated slope", y = "Samples") + theme_minimal()
print(p15); ggsave("outputs/figures/F15_coverage.png", p15, width = 7, height = 4, dpi = 200)


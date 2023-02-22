library(readr)
library(carbondate)

output_file = "output/kerr_predictive_density.csv"
data_file = "data/kerr.csv"
intcal_file = "data/intcal20.14c"

input_data <- read.csv(data_file)

intcal20 <- read.table(intcal_file, sep = ",", header = FALSE, skip = 11)
intcal20 = data.frame(calendar_age = intcal20[, 1], c14_age = intcal20[, 2], c14_sig = intcal20[, 3])

predictive_density <- read_csv(output_file, show_col_types = FALSE)

##############################################################################
# Initialise plotting parameters

SPD_colour <- grDevices::grey(0.1, alpha = 0.3)
calibration_curve_colour <- "blue"
calibration_curve_bg <- grDevices::rgb(0, 0, 1, .3)
output_colour <- "purple"

calendar_age_sequence <- predictive_density$calendar_age

##############################################################################
# Calculate plot scaling

.ScaleLimit <- function(lim, limscal) {
  lim <- lim + c(1, -1) * diff(lim) * (1 - limscal)
  return(lim)
}

xlim <- .ScaleLimit(rev(range(calendar_age_sequence)), 1)
ylim_calibration <- .ScaleLimit(
  range(input_data$c14_ages) + c(-2, 2) * stats::quantile(input_data$c14_sig, 0.9), 1)
ylim_density = c(0, 3 * max(predictive_density$mean))

##############################################################################
# Plot calibration curve

graphics::plot.default(
  intcal20$calendar_age,
  intcal20$c14_age,
  col = calibration_curve_colour,
  ylim = ylim_calibration,
  xlim = xlim,
  xlab = "Calendar Age (cal yr BP)",
  ylab = expression(paste(""^14, "C", " age (", ""^14, "C yr BP)")),
  type = "l",
  main = "Predictive density")
# multiplier for the confidence interval if you have a standard deviation
zquant <- - stats::qnorm((1 - 0.95) / 2)
intcal20$ub <- intcal20$c14_age + zquant * intcal20$c14_sig
intcal20$lb <- intcal20$c14_age - zquant * intcal20$c14_sig

graphics::lines(
  intcal20$calendar_age,
  intcal20$ub,
  lty = 2,
  col = calibration_curve_colour)
graphics::lines(
  intcal20$calendar_age,
  intcal20$lb,
  lty = 2,
  col = calibration_curve_colour)
graphics::polygon(
  c(rev(intcal20$calendar_age), intcal20$calendar_age),
  c(rev(intcal20$lb), intcal20$ub),
  col = calibration_curve_bg,
  border = NA)
graphics::rug(input_data$c14_ages, side = 2)

##############################################################################
# Plot density
graphics::par(new = TRUE)
graphics::plot.default(
  c(),
  c(),
  type = "n",
  ylim = ylim_density,
  xlim = xlim,
  axes = FALSE,
  xlab = NA,
  ylab = NA)

graphics::lines(
  predictive_density$calendar_age,
  predictive_density$mean,
  col = output_colour)
graphics::lines(
  predictive_density$calendar_age,
  predictive_density$ci_lower,
  col = output_colour,
  lty = 2)
graphics::lines(
  predictive_density$calendar_age,
  predictive_density$ci_upper,
  col = output_colour,
  lty = 2)

graphics::par(new = FALSE)

wo = WalkerBivarDirichlet(kerr$c14_ages, kerr$c14_sig, intcal20, 1e5, 10)
pd = FindPredictiveCalendarAgeDensity(wo, predictive_density$calendar_age, n_posterior_samples = 5000)
lines(pd$calendar_age, pd$density_mean, col="red", lty=3)
legend("topright", c("C++ executable", "R package"), col = c("purple", "red"), lty= c(1, 3))

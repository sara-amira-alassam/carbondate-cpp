library(readr)
library(carbondate)

ident = 4
output_file = paste("output/kerr_posterior_density_", ident, ".csv", sep="")

posterior_density <- read_csv(output_file, show_col_types = FALSE)

wo = WalkerBivarDirichlet(kerr$c14_ages, kerr$c14_sig, intcal20, 1e5, 10)
PlotCalendarAgeDensityIndividualSample(ident+1, wo, resolution = 5)
lines(posterior_density$calendar_age, posterior_density$probability, col="purple", lty=1)

legend("topright", c("C++ executable", "R package"), col = c("purple", "grey"), lty=c(1, -1), pch=c(NA, 15))

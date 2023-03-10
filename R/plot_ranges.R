plot_range = function(range, colour, height) {
  for (i in 1:length(range)) {
    lines(x = c(range[[i]][1], range[[i]][1]), y = c(0, height), col = colour)
    lines(x = c(range[[i]][1], range[[i]][2]), y = c(height, height), col = colour)
    lines(x = c(range[[i]][2], range[[i]][2]), y = c(0, height), col = colour)
  }
}

# Finds the area under a curve (defined by a set of calendar ages and probabilities) if we ignore
# all values below the cut-off
find_AUC = function(cal_age, prob, cut_off) {
  resolution = cal_age[2] - cal_age[1]
  ranges = list()
  areas = c()
  ind = 1
  for (i in 1:(length(cal_age) - 1)) {
    y1 = prob[i]
    y2 = prob[i+1]
    if (y1 <= cut_off && y2 >= cut_off) {
      dx = (cut_off - y1) / (y2 - y1) * resolution
      x_intercept_1 = cal_age[i] + dx
      areas[ind] = (cut_off + y2) * (resolution - dx) / 2
    } else if (y1 >= cut_off && y2 <= cut_off) {
      dx = (cut_off - y1) / (y2 - y1) * resolution
      x_intercept_2 = cal_age[i] + dx
      areas[ind] = areas[ind] + (cut_off + y1) * dx / 2
      ranges[[ind]] = c(x_intercept_1, x_intercept_2)
      ind = ind + 1
    } else if (y1 >= cut_off && y2 >= cut_off) {
      areas[ind] = areas[ind] + (y1 + y2) * resolution / 2
    }
  }
  return(list(area = sum(areas), ranges = ranges, areas = areas))
}

find_ranges = function(cal_age, prob, range_prob, print_ranges = FALSE) {
  # Use bisection method to find closest point
  auc = list(area = 0,ranges = c())
  a = 0
  b = max(prob)
  while (abs(auc$area - range_prob) > 1e-4) {
    p = (a + b) / 2
    auc = find_AUC(cal_age, prob, p)
    if (auc$area < range_prob) {
      b = p
    } else {
      a = p
    }
  }
  if (print_ranges) {
    print(paste(range_prob*100, "% probability"))
    for (i in 1:length(auc$ranges)) {
      print(paste("   ", round(auc$ranges[[i]][1]), "AD (", round(auc$areas[i] * 100, 1) ,"%) ", round(auc$ranges[[i]][2]),"AD", sep=""))
    }
  }
  return(auc$ranges)
}

## This is the data from point 5  - R-date 1270, 21
oxcal_prob = 0.012238 * c(0.000000, 0.000000, 0.000000, 0.000003, 0.000153, 0.003486, 0.038156, 0.141632, 0.329171, 0.567360, 0.770778, 0.923278, 0.950834, 0.869290, 1.000000, 0.763806, 0.800847, 0.801948, 0.955358, 0.988562, 0.979533, 0.995043, 0.626443, 0.425606, 0.273592, 0.679178, 0.756552, 0.506157, 0.219727, 0.000124, 0.003249, 0.069311, 0.183418, 0.106422, 0.088254, 0.103411, 0.146762, 0.101277, 0.042832, 0.013269, 0.009174, 0.004870, 0.004965, 0.008045, 0.015439, 0.022541, 0.021342, 0.018205, 0.010626, 0.002332, 0.000138, 0.000002, 0.000000, 0.000000, 0.000000)
oxcal_age = seq(635.5, by = 5, length.out = length(oxcal_prob))
oxcal_range_68.3 = list(c(685, 744))
oxcal_range_95.4 = list(c(670, 777.5), c(793, 799), c(813.5, 817))
oxcal_range_99.7 = list(c(661.5, 780.5), c(786, 833.5), c(852.5, 875.5))

## First compare calculating ranges with the same dataset - so for the Oxcal probability result,
## compare the OxCal rages and the ranges we calculate

plot(
  oxcal_age,
  oxcal_prob,
  col = "blue",
  xlab = "Calendar Age AD",
  ylab = "Probability",
  main = "Point 5 - R_Date(1270, 21)",
  type = "o")
plot_range(oxcal_range_68.3, "darkgreen", 2e-4)
plot_range(find_ranges(oxcal_age, oxcal_prob, 0.683, TRUE), "green", 3e-4)

plot_range(oxcal_range_95.4, "darkred", 6e-4)
plot_range(find_ranges(oxcal_age, oxcal_prob, 0.954, TRUE), "red", 7e-4)

plot_range(oxcal_range_99.7, "purple4", 9e-4)
plot_range(find_ranges(oxcal_age, oxcal_prob, 0.997, TRUE), "purple", 10e-4)

legend(
  "topright",
  c("OxCal density", "68.3% range OxCal", "68.3% range CD", "95.4% range OxCal", "95.4% range CD", "99.7% range OxCal", "99.7% range CD"),
  col = c("blue", "darkgreen", "green", "darkred", "red", "purple4", "purple"),
  lty = c(1, 1, 1, 1, 1, 1))


## Now compare calculating ranges with different datasets - compare range calculated with our
## datasets against the data from OxCal
plot(
  oxcal_age,
  oxcal_prob,
  col = "blue",
  xlab = "Calendar Age AD",
  ylab = "Probability",
  main = "Point 5 - R_Date(1270, 21)",
  type = "o")

# These results with resolution 1 give quite different ranges
carbondate_prob = 0.012462 * c(0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000002, 0.000006, 0.000013, 0.000032, 0.000066, 0.000133, 0.000306, 0.000578, 0.001063, 0.001905, 0.004136, 0.006878, 0.011154, 0.017640, 0.027207, 0.037038, 0.049714, 0.072007, 0.093518, 0.111798, 0.142613, 0.166373, 0.207076, 0.237648, 0.288603, 0.325826, 0.365446, 0.407204, 0.476696, 0.520500, 0.564773, 0.586918, 0.630869, 0.673867, 0.715289, 0.754505, 0.773090, 0.820681, 0.853614, 0.883696, 0.911309, 0.923068, 0.933381, 0.960445, 0.951091, 0.940154, 0.913766, 0.881819, 0.863958, 0.844954, 0.824899, 0.844954, 0.881819, 0.927691, 0.968166, 0.982046, 0.951091, 0.881819, 0.803888, 0.736112, 0.712275, 0.736112, 0.759394, 0.782020, 0.782020, 0.803888, 0.782020, 0.782020, 0.759394, 0.759394, 0.759394, 0.790038, 0.835085, 0.876198, 0.912564, 0.943437, 0.956607, 0.978063, 0.978063, 0.978063, 0.968170, 0.968170, 0.956607, 0.956607, 0.956607, 0.956607, 0.956607, 0.978063, 0.992634, 1.000000, 0.978063, 0.928730, 0.835085, 0.766306, 0.665874, 0.587768, 0.535880, 0.510263, 0.484974, 0.460087, 0.411788, 0.388496, 0.343878, 0.282430, 0.245415, 0.228135, 0.263516, 0.322633, 0.435671, 0.561744, 0.691590, 0.790038, 0.835085, 0.856186, 0.835085, 0.766306, 0.665874, 0.561744, 0.460087, 0.411788, 0.388496, 0.460087, 0.721094, 1.000000, 0.382111, 0.013829, 0.000307, 0.000032, 0.000018, 0.000067, 0.000200, 0.000475, 0.000917, 0.001260, 0.002002, 0.003129, 0.004809, 0.009485, 0.020163, 0.040105, 0.074631, 0.119020, 0.167235, 0.196049, 0.196049, 0.196049, 0.167235, 0.154030, 0.129942, 0.108814, 0.106212, 0.106212, 0.097084, 0.097084, 0.095438, 0.095438, 0.087246, 0.087246, 0.087246, 0.087246, 0.095438, 0.113617, 0.134337, 0.145700, 0.149467, 0.149467, 0.149467, 0.145700, 0.134337, 0.113617, 0.111798, 0.094274, 0.085849, 0.065792, 0.054667, 0.040920, 0.033471, 0.027207, 0.019705, 0.015767, 0.014785, 0.013221, 0.010523, 0.011219, 0.010024, 0.008943, 0.007968, 0.007968, 0.005793, 0.005793, 0.005793, 0.004707, 0.005350, 0.005350, 0.005350, 0.004967, 0.005658, 0.005658, 0.006435, 0.007307, 0.008283, 0.009375, 0.009375, 0.012538, 0.014071, 0.015767, 0.017640, 0.017640, 0.019705, 0.021976, 0.021976, 0.024471, 0.024471, 0.021295, 0.021295, 0.021295, 0.021295, 0.021295, 0.021295, 0.019034, 0.019034, 0.016984, 0.015130, 0.015767, 0.012538, 0.011154, 0.008786, 0.006878, 0.004707, 0.003627, 0.002425, 0.001599, 0.001039, 0.000666, 0.000360, 0.000190, 0.000115, 0.000039, 0.000018, 0.000009, 0.000002, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.00000);
carbondate_age = seq(645.5, by = 1, length.out = length(carbondate_prob))

#carbondate_prob = 0.012172 * c(0.000000, 0.000024, 0.000816, 0.013721, 0.074546, 0.213418, 0.429097, 0.649421, 0.836540, 0.958104, 0.910060, 0.910649, 0.891694, 0.772295, 0.795820, 0.854495, 0.989828, 0.984081, 1.000000, 0.854655, 0.528055, 0.342350, 0.370952, 0.820650, 0.586784, 0.604390, 0.002918, 0.000994, 0.015907, 0.154177, 0.154808, 0.102793, 0.090997, 0.130747, 0.141810, 0.084436, 0.028065, 0.012238, 0.007466, 0.005436, 0.006148, 0.010983, 0.018986, 0.023241, 0.021339, 0.016268, 0.007198, 0.001247, 0.000076, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000001, 0.000000, 0.000000, 0.000000)
#carbondate_age = seq(647.5, by = 5, length.out = length(carbondate_prob))

lines(carbondate_age, carbondate_prob, col="lightblue", type = "o")

plot_range(oxcal_range_68.3, "darkgreen", 2e-4)
plot_range(find_ranges(carbondate_age, carbondate_prob, 0.683, TRUE), "green", 3e-4)

plot_range(oxcal_range_95.4, "darkred", 6e-4)
plot_range(find_ranges(carbondate_age, carbondate_prob, 0.954, TRUE), "red", 7e-4)

plot_range(oxcal_range_99.7, "purple4", 9e-4)
plot_range(find_ranges(carbondate_age, carbondate_prob, 0.997), "purple", 10e-4)

legend(
  "topright",
  c("OxCal density", "Carbondate density", "68.3% range OxCal", "68.3% range CD", "95.4% range OxCal", "95.4% range CD", "99.7% range OxCal", "99.7% range CD"),
  col = c("blue", "lightblue", "darkgreen", "green", "darkred", "red", "purple4", "purple"),
  lty = c(1, 1, 1, 1, 1, 1))



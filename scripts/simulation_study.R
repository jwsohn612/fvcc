library(fvcc)

ID <- SimulationData$ID
W <- SimulationData[c('W1', 'W2')]
X <- SimulationData[c('X1', 'X2')]
Z <- SimulationData[c('Z1', 'Z2')]
t <- SimulationData$t
Y <- SimulationData$Y

output <- fvcc(
  ID = ID,
  W = W,
  X = X,
  Z = Z,
  t = t,
  Y = Y,
  thinning_unit = 5, 
  num_knots = 30,
  K_max = 10,
  num_iters = 10000
)

# If burn = NULL, it discard a half of samples as default.
par(mfrow = c(1,2))
plot(output$fixed_effect, position = 1, burn = NULL)
plot(output$fixed_effect, position = 2, burn = NULL)

par(mfrow = c(1,3))
plot(output$random_effect, row_position = 1, col_position = 1)
plot(output$random_effect, row_position = 1, col_position = 2)
plot(output$random_effect, row_position = 2, col_position = 2)

# check the estimated cluster labels for all subjects and then put cluster_number to draw varying coefficients
table(plot(output$cluster))
par(mfcol=c(2,3))
plot(output$varying_coefficient, cluster_number = 1, variable_number = 1)
plot(output$varying_coefficient, cluster_number = 1, variable_number = 2)
plot(output$varying_coefficient, cluster_number = 6, variable_number = 1)
plot(output$varying_coefficient, cluster_number = 6, variable_number = 2)
plot(output$varying_coefficient, cluster_number = 3, variable_number = 1)
plot(output$varying_coefficient, cluster_number = 3, variable_number = 2)

par(mfcol=c(2,3))
plot(output$latent_location, cluster_number = 1, variable_number = 1,time_range = output$time_range,knot_position = output$knot_position)
plot(output$latent_location, cluster_number = 1, variable_number = 2,time_range = output$time_range,knot_position = output$knot_position)
plot(output$latent_location, cluster_number = 6, variable_number = 1,time_range = output$time_range,knot_position = output$knot_position)
plot(output$latent_location, cluster_number = 6, variable_number = 2,time_range = output$time_range,knot_position = output$knot_position)
plot(output$latent_location, cluster_number = 3, variable_number = 1,time_range = output$time_range,knot_position = output$knot_position)
plot(output$latent_location, cluster_number = 3, variable_number = 2,time_range = output$time_range,knot_position = output$knot_position)


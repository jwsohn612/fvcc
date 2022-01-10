library(cvarpyp)

ID <- SimulationData$ID
W <- SimulationData[c('W1', 'W2')]
X <- SimulationData[c('X1', 'X2')]
Z <- SimulationData[c('Z1', 'Z2')]
t <- SimulationData$t
Y <- SimulationData$Y

output <- cvarpyp(
  ID = ID,
  W = W,
  X = X,
  Z = Z,
  t = t,
  Y = Y,
  num_knots = 15,
  K_max = 10,
  num_iters = 10000
)

plot(output$fixed_effect, position = 1)
plot(output$random_effect,
     row_position = 1,
     col_position = 1)
plot(output$cluster)
plot(
  output$varying_coefficient,
  cluster_number = 5,
  variable_number = 1
)
plot(
  output$latent_location,
  cluster_number = 5,
  variable_number = 1,
  time_range = output$time_range,
  knot_position = output$knot_position
)
plot(output$nu)
plot(output$lambda)

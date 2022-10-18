
library(fvcc)
library(rlist)

# Place your own directory to save model results
save_dir = "/Users/jinwonsohn/Desktop/packages/temp/"

sbj1 <- 400 
sbj2 <- 400
sbj3 <- 400

seeds <- sample(1:1000, size = 150, replace = FALSE)

for(seed in seeds){
  
  data <- get_simulation_data(sbj1, sbj2, sbj3, seed = seed)
  ID <- data$ID
  W <- data[c('W1','W2')]
  X <- data[c('X1','X2')]
  Z <- data[c('Z1','Z2')]
  t <- data$t
  Y <- data$Y
  tryCatch({
    output <- fvcc(ID = ID, W = W, X = X, Z = Z, t = t, Y = Y, thinning_unit = 5,  num_knots = 30, K_max = 10, num_iters = 10000)
    rlist::list.save(output, paste0(save_dir,"model_output_",seed,".rdata"))
  }, error = function(e){
    print(paste(seed, 'confront some computational errors. It needs to be re-run.'))
  }
  )
}
library(fvcc)

data("RealData")

# Data Preprocessing
dataset <- data %>%
  filter(ABITUR == 1) %>%
  filter(AGE <= 53)

# For a varying intercept
dataset$int00 <- rep(1, nrow(dataset))

# Reorder subjects id 
temp_tab <- table(dataset$id)
dataset$id_changed <- unlist(sapply(1:length(unique(dataset$id)),function(x) rep(x,temp_tab[x])))

ID <- dataset$id_changed
W <- dataset[c('int00','HHKIDS')]
X <- dataset[c("NEWHSAT",'HANDPER','MARRIED')]
Z <- dataset[c('int00')]
t <- dataset$AGE
Y <- dataset$WORKING

init_Psi1 <- matrix(c(2),nrow = 1)
init_beta1 <- c(-0.1, 0.1, 2.0)

output1 <- fvcc(ID = ID, W = W, X = X, Z = Z, t = t, Y = Y,
                thinning_unit = 50,  num_knots = 12, nu = 1,  
                Psi = init_Psi1, beta = init_beta1, K_max = 15, num_iters = 150000)

init_Psi2 <- matrix(c(16),nrow = 1)
init_beta2 <- c(0.1, -0.1, -2.0)

output2 <-  fvcc(ID = ID, W = W, X = X, Z = Z, t = t, Y = Y,
                 thinning_unit = 50,  num_knots = 12, nu = 1,  
                 Psi = init_Psi2, beta = init_beta2, K_max = 15, num_iters = 150000)

init_Psi3 <- matrix(c(10),nrow = 1)
init_beta3 <- c(-0.05, 0.05, 0.10)

output3 <-  fvcc(ID = ID, W = W, X = X, Z = Z, t = t, Y = Y,
                 thinning_unit = 50,  num_knots = 12, nu = 1,  
                 Psi = init_Psi3, beta = init_beta3, K_max = 15, num_iters = 150000)


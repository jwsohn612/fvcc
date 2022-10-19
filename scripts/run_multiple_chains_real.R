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


# ==== Calculate Gelman-Rubin Diagnostics ==== # 
get_gelman_score <- function(output_list){
  temp <- purrr::map(output_list, function(output){
    f1 <- purrr::map_dbl(output$fixed_effect, ~ .x[1])
    f2 <- purrr::map_dbl(output$fixed_effect, ~ .x[2])
    f3 <- purrr::map_dbl(output$fixed_effect, ~ .x[3])
    psi <- unlist(output$random_effect)
    
    mcmc(cbind(f1,f2,f3,psi))
  })
  comb.chain <- mcmc.list(temp)
  return(comb.chain)
}

mcmc_output <- get_gelman_score(list(output1, output2, output3))
gelman.diag(mcmc_output)

routput1 <- fvcc::relabel_chain(output1)
routput2 <- fvcc::relabel_chain(output2)
routput3 <- fvcc::relabel_chain(output3)

var_on_knots1 <- fvcc::make_vc_outputs_on_knots(routput1)
var_on_knots2 <- fvcc::make_vc_outputs_on_knots(routput2)
var_on_knots3 <- fvcc::make_vc_outputs_on_knots(routput3)

calculate_pointwise_gelman <- function(var_on_knots1, var_on_knots2, var_on_knots3){
  size_knots <- ncol(var_on_knots1[[1]][[1]])
  pt_wise_gelman_list <- list()
  
  for(k in 1:2){
    pt_wise_gel_cls <- list()
    
    for(p in 1:2){
      coef1 <- var_on_knots1[[k]][[p]]
      coef2 <- var_on_knots2[[k]][[p]]
      coef3 <- var_on_knots3[[k]][[p]]
      
      pt_wise_gelman <- purrr::map_dbl(1:size_knots, function(t){
        param1 <- mcmc(coef1[,t])
        param2 <- mcmc(coef2[,t])
        param3 <- mcmc(coef3[,t])
        comb.chain <- mcmc.list(param1, param2, param3)
        res <- gelman.diag(comb.chain)
        res$psrf[,1]  
      })
      pt_wise_gel_cls[[p]] <- pt_wise_gelman
    }
    pt_wise_gelman_list <- list.append(pt_wise_gelman_list, pt_wise_gel_cls)
  }
  return(pt_wise_gelman_list)
}

calculate_pointwise_gelman(var_on_knots1, var_on_knots2, var_on_knots3)


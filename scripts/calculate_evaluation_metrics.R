library(dplyr)
library(ggplot2)
library(cowplot)
library(fvcc)

# Place your own directory
save_dir <- "/Users/jinwonsohn/Desktop/packages/temp/"

file_list = list.files(save_dir)
replicates = length(file_list)

sbj1 <- 400
sbj2 <- 400
sbj3 <- 400

output <- fvcc::get_vc_results(save_dir, 400, 400, 400)

grid <- 100
timedomain <- ppoints(grid)
size_of_data <- length(save_dir)

# Get Coverage Ratio
g1_a1<-function(x) 2*exp(-200*(x-0.2)^2)+exp(-10*(x-0.6)^2)
g1_a2<-function(x) sin(2*pi*x^3)

g2_a1<-function(x) sin(8*(x-0.5))+1.5*exp(-(20^2)*(x-0.5)^2)
g2_a2<-function(x) 2*x + 0 

g3_a1<-function(x) -2*x + 0
g3_a2<-function(x) 0*x + 0

alpha<-list(c(g1_a1,g1_a2),c(g2_a1,g2_a2),c(g3_a1,g3_a2))

# === Calculate Coverages 
get_coverages <- function(varying_coefs_list, timedomain, size_of_data = 100){
  coverage <- list()
  for(k in 1:3){
    coverage_p <- list()
    for(p in 1:2){
      cov_matrix <- matrix(rep(0, length(timedomain)*size_of_data),nrow=size_of_data)
      for(i in 1:size_of_data){
        tmp <- alpha[[k]][[p]](timedomain)
        var <- varying_coefs_list[[i]][[k]][[1]][[p]]
        lower <- var[1,]
        upper <- var[3,]
        cov_matrix[i,] <- (lower <= tmp)*(tmp <= upper)
      }
      coverage_p[[p]] <- apply(cov_matrix,2,mean)
    }  
    coverage[[k]] <- coverage_p
  }
  return(coverage)
}

get_coverages(output$varying_coefs_list, timedomain, size_of_data)

# === Calculate Clutering Performance Measures
unify_factors <- function(classes, x){
  unex <- classes[!(classes %in% names(x))]
  if(length(unex) == 0){
    return(x)
  }else{
    vec <- rep(0,length(unex))
    names(vec) <- unex
    return(as.numeric(c(as.numeric(x),vec)))
  }
}

metrics <- function(cluster_list, sbj1, sbj2, sbj3){
  purrr::map_df(cluster_list, function(cls){
    
    tab <- table(c(rep(1,sbj1),rep(2,sbj2),rep(3,sbj3)), cls)
    
    row_classes <- rownames(tab)
    col_classes <- colnames(tab)
    classes <- unique(c(row_classes,col_classes))
    
    n <- sum(tab)
    di <- unify_factors(classes, diag(tab))
    rows <- unify_factors(classes, rowSums(tab))
    cols <- unify_factors(classes, colSums(tab))
    
    accus <- sum(di)/n
    pre <- di / cols
    re <- di / rows
    f1 <- 2 * pre * re / (pre + re)
    data.frame(accuracy = accus,
               precision1 = pre[1], recall1 = re[1], f1score1 = f1[1],
               precision2 = pre[2], recall2 = re[2], f1score2 = f1[2],
               precision3 = pre[3], recall3 = re[3], f1score3 = f1[3])
  })
}

metric <- metrics(output$cluster_list, sbj1, sbj2, sbj3)



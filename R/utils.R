#' @useDynLib cvarpyp
#' @importFrom Rcpp evalCpp
#' @import dplyr
#' @import purrr
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("cvarpyp", libpath)
}

tol <- 1e-200
solve <- function(x) base::solve(x, tol = tol)

preprocess_data <- function(ID, W, X, Z, t, Y) {
  id_mapper <- table(ID)
  ID <- map(1:length(id_mapper), ~ rep(.x, id_mapper[.x])) %>% unlist
  
  # Check X is null
  if (is.null(X)) {
    X <- data.frame(rep(0,length(ID)))
    q <- 1
  } else{
    q <- ncol(X)
  }
  # Change Name
  p <- ncol(W); r <- ncol(Z)
  data.names = map2(c(p,q,r),c('W','X','Z'), ~ paste0(.y, seq(1:.x)))
  
  temp <- map2(1:3, list(W,X,Z), function(x,y) {
    colnames(y) = data.names[[x]]
    y
  })
  
  temp <- bind_cols(temp) %>% mutate(
    ID = ID,
    t=t,
    Y=Y
  ) %>% arrange(ID, t)
  
  list(ID = temp %>% pull(ID),
       t = temp %>% pull(t),
       Y = temp %>% pull(Y),
       W = temp %>% select(all_of(data.names[[1]])),
       X = temp %>% select(all_of(data.names[[2]])),
       Z = temp %>% select(all_of(data.names[[3]])))
}

prod_v <- function(v,idx){
  rev_v <- prod((1-v)[1:idx-1])
  v[idx]*rev_v
}

get_DP_pi <- function(v, K){
  c(v[1], map_dbl(2:K, function(i) prod_v(v,i)))
}

get_sbj_clusters <- function(K, clusters){
  temp_cluster <- rep(0,K)
  names(temp_cluster) <- 1:K
  tab <- table(clusters)
  temp_cluster[names(tab)] <- tab %>% as.numeric
  temp_cluster
}


# ==== Unit Basis Data
B<- function(t, knot_position, sd){
  dist <- abs(t - knot_position)/sd
  basis <- dist^3
  return(c(1, t, basis))
}

make_sym <- function(A) (t(A)+A)/2
# ==== Generate Basis Data
get_basis_data <- function(tData, t, knot_position, sd){
  n_cols <- ncol(tData)
  UnitBasis <- lapply(t, function(x) B(x, knot_position, sd))
  out_data<-list()
  for(i in 1:n_cols){
    temp_data <- lapply(1:length(t), function(x) tData[x,i]*UnitBasis[[x]])
    out_data[[i]] <- do.call(rbind, temp_data)
  }
  return(out_data)
}

make_sym <- function(A) (t(A)+A)/2



flatten_gamma_with_fixC <- function(gamma_cls)  as.vector(t(gamma_cls))
# 
# # Function for location vectors
SameElements <- function (a,b){
  return(sum(a-b)>0)
}

# This does not consider a main effect term, but it is handled in other function.
devide_randomly <- function(current_gamma, num_knots){
  interact_term <- list(current_gamma[2])
  
  smoothing_term <- current_gamma[3:(num_knots+2)]
  
  blks <- sample(2:4, size = num_knots, replace = TRUE)
  cum_index <- cumsum(blks)
  end_index <- c(cum_index[cum_index<num_knots],num_knots)
  
  start_index <- c(1, end_index + 1)
  blocked_smoothing_term <- map(1:length(end_index), ~ smoothing_term[start_index[.x]:end_index[.x]])
  return(c(interact_term, blocked_smoothing_term))
}


propose_new_gamma_and_probs <- function(rand_seq, num_knots, a_pi, b_pi){
  
  which_blk <- sample(1:length(rand_seq), 1)
  current_blk <- rand_seq[[which_blk]]
  
  SCurrentBlk <- sum(current_blk) ; size_target <- length(current_blk)
  SRemainBlks <- sum(unlist(rand_seq[-which_blk]))
  
  blk_list <- blks_list[[size_target]]
  
  SCanditBlks <- lapply(blk_list,function(x) sum(x))
  
  denominator <- beta(SCurrentBlk+SRemainBlks+a_pi, (num_knots+1)-SCurrentBlk-SRemainBlks+b_pi)
  nominators <- vapply(SCanditBlks, function(x) beta(x+SRemainBlks+a_pi, (num_knots+1)-x-SRemainBlks+b_pi), numeric(1), USE.NAMES = FALSE)
  one_over_current_probs <- sum(nominators/denominator)
  probs_suggested <- nominators/denominator/one_over_current_probs
  candit_idx <- sample(1:length(blk_list),size=1,prob = probs_suggested)
  suggested_blk <- blk_list[[candit_idx]]
  with_prob <- probs_suggested[candit_idx]
  rand_seq[[which_blk]] <- suggested_blk
  changed <- SameElements(suggested_blk, current_blk)
  return(list(proposed_gamma=c(1,unlist(rand_seq)), changed=changed, with_prob_new=with_prob, with_prob_old=1/one_over_current_probs)) # 1 implies the intact varible itself
}

Gen_proposed_gamma<-function(current_gamma, num_knots, a_pi, b_pi){
  blocked_seq <- devide_randomly(current_gamma, num_knots)
  proposed_info <- propose_new_gamma_and_probs(blocked_seq, num_knots, a_pi, b_pi)
  proposed_gamma <- proposed_info$proposed_gamma ; changed <- proposed_info$changed
  proposing_prob <- proposed_info$with_prob_new ; remaining_prob <- proposed_info$with_prob_old
  return(proposed_gamma)
}

blk_1 <- list(0,1)
blk_2 <- list(c(0,0),c(0,1),c(1,0),c(1,1))
blk_3 <- list(c(0,0,0),c(0,0,1),c(0,1,0),c(1,0,0),c(1,1,0),c(1,0,1),c(0,1,1),c(1,1,1))
blk_4 <- list(c(0,0,0,0),c(0,0,0,1),c(0,0,1,0),c(0,1,0,0),c(1,0,0,0),c(0,0,1,1),c(0,1,0,1),c(1,0,0,1),c(1,0,1,0),c(0,1,1,0),c(1,1,0,0),        c(1,1,1,0),c(0,1,1,1),c(1,0,1,1),c(1,1,0,1),c(1,1,1,1))
blks_list <- list(blk_1, blk_2, blk_3, blk_4)

sample_gamma <- function(ID, K.max, p, active, num_knots, a_pi, b_pi, NSbj, gamma, tau, C_star_sbj, InvAdj, L_sbj, X_sbj, beta, Clusters) {
  for(k in 1:K.max){
    for(l in 1:p){
      # print(gamma[[k]][l,])
      if(k %in% active){
        temp_gamma <- gamma[[k]]
        current_gamma <- gamma[[k]][l,]
        proposed_gamma <- Gen_proposed_gamma(current_gamma, num_knots, a_pi, b_pi)
        temp_gamma[l,] <- proposed_gamma
        proposed_gamma_in <- flatten_gamma_with_fixC(temp_gamma)
        current_gamma_in <- flatten_gamma_with_fixC(gamma[[k]])
        
        changed <- sum(proposed_gamma_in!=current_gamma_in)
        
        if(changed==FALSE){
          gamma[[k]][l,] <- current_gamma
        }else{
          L_new <- KnotSelection(k, NSbj, tau[k], C_star_sbj, proposed_gamma_in, InvAdj, X_sbj, beta, L_sbj, ID-1, Clusters)[[1]]
          L_cur <- KnotSelection(k, NSbj, tau[k], C_star_sbj, current_gamma_in, InvAdj, X_sbj, beta, L_sbj, ID-1, Clusters)[[1]]
          accept_knotR <- exp(L_new - L_cur)
          if(runif(1)<accept_knotR) {
            gamma[[k]][l,] <- proposed_gamma 
          } else {
            gamma[[k]][l,] <- current_gamma
          }
        }
      }else{
        gamma[[k]][l,] <- c(1,rbinom(n=1+num_knots,size=1,prob=1/(num_knots+1)))
      }
    }
  }
  return(gamma)
}

sample_gamma_unit <- function(active, k, l, gamma, num_knots, a_pi, b_pi, NSbj, tau, C_star_sbj, InvAdj, L_sbj, ID, X_sbj, beta, Clusters){
  if(k %in% active){
    
    temp_gamma <- gamma[[k]]
    current_gamma <- gamma[[k]][l,]
    
    proposed_gamma <- Gen_proposed_gamma(current_gamma, num_knots, a_pi, b_pi)
    
    temp_gamma[l,] <- proposed_gamma
    proposed_gamma_in <- flatten_gamma_with_fixC(temp_gamma)
    current_gamma_in <- flatten_gamma_with_fixC(gamma[[k]])
    
    changed <- sum(proposed_gamma_in!=current_gamma_in)
    
    if(changed==FALSE){
      return(current_gamma)
    }else{
      L_new <- KnotSelection(k, NSbj, tau[k], C_star_sbj, proposed_gamma_in, InvAdj, X_sbj, beta, L_sbj, ID-1, Clusters)[[1]]
      L_cur <- KnotSelection(k, NSbj, tau[k], C_star_sbj, current_gamma_in, InvAdj, X_sbj, beta, L_sbj, ID-1, Clusters)[[1]]
      accept_knotR <- exp(L_new - L_cur)
      
      if (runif(1) < accept_knotR) {
        return(proposed_gamma)
      }else {
        return(current_gamma)
      }
    }
    
  }else{
    return(c(1,rbinom(n=1+num_knots,size=1,prob=1/(num_knots+1))))
  }
  
}


sample_nu <-function(current_nu, lambda, v, sd_of_nu, K.max, nu_param){
  trans_current <- log(current_nu+lambda)
  suggested <- rnorm(1, mean=trans_current, sd=sd_of_nu)
  proposed_nu <- exp(suggested)-lambda
  proposed <- map_dbl(1:(K.max-1), function(k) dbeta(v[k], 1-lambda, proposed_nu+k*lambda, log = T)) %>% reduce(`+`)  - nu_param*proposed_nu + suggested
  current <- map_dbl(1:(K.max-1), function(k) dbeta(v[k], 1-lambda, current_nu+k*lambda, log = T)) %>% reduce(`+`) - nu_param*current_nu  + trans_current
  prob <- exp(proposed-current)
  if(prob>=1){
    return(proposed_nu)
  }else{
    u<-runif(1)
    if(u<=prob){
      return(proposed_nu)
    }else{
      return(current_nu)}
  }
}

sample_lambda <- function(current_lambda, nu, v, K.max, sd_of_lambda){
  if(nu>0){
    trans_lambda <- log(current_lambda/(1-current_lambda))
    suggested <- rnorm(1, mean=trans_lambda, sd=sd_of_lambda)
    proposed_lambda <- 1/(1+exp(-suggested))
    log_jacobian_current <- -trans_lambda-2*log(1+exp(-trans_lambda))
    log_jacobian_proposed <- -suggested-2*log(1+exp(-suggested))
  }else{
    delta <- (current_lambda+nu)/(1+nu)
    trans_lambda <- log(delta/(1-delta))
    suggested <- rnorm(1, mean=trans_lambda, sd=sd_of_lambda)
    newly_delta <- 1/(1+exp(-suggested))
    proposed_lambda <- (1+nu)*newly_delta-nu
    
    log_jacobian_current <- -trans_lambda-2*log(1+exp(-trans_lambda))   # log(1+nu)?? ?????? ?????? ????.
    log_jacobian_proposed <- -suggested-2*log(1+exp(-suggested))
  }
  proposed <- map_dbl(1:(K.max-1), function(k) dbeta(v[k],1-proposed_lambda,nu + k*proposed_lambda, log = T)) %>% reduce(`+`) - proposed_lambda  + log_jacobian_proposed
  current <- map_dbl(1:(K.max-1), function(k) dbeta(v[k],1-current_lambda,nu + k*current_lambda, log = T)) %>% reduce(`+`) - current_lambda +log_jacobian_current
  
  prob <- exp(proposed-current)
  if(prob>=1){
    return(proposed_lambda)
  }else{
    u<-runif(1)
    if(u<=prob){
      return(proposed_lambda)
    }else{
      return(current_lambda)
    }
  }
}


for_theta <- function(num_k, active, XXi_k, Xi_k, Rk, tau) {
  
  if(num_k %in% active){
    cov_mat <- solve(make_sym( XXi_k[[num_k]] +  Rk[[num_k]]/tau[num_k] ))
    mean_vec <- cov_mat %*% Xi_k[[num_k]]
    return(as.vector(mvtnorm::rmvnorm(1, mean=as.vector(mean_vec), sigma=make_sym(cov_mat), method='eigen')))
    
  }else{
    cov_mat <- solve(Rk[[num_k]])
    return(as.vector(mvtnorm::rmvnorm(1, mean=rep(0,dim(cov_mat)[1]), sigma=tau[num_k]*make_sym(cov_mat), method='eigen')))
  }
}

sample_theta <- function(K.max, active, XXi_k, Xi_k, Rk, tau) {
  theta_k <- lapply(1:K.max, function(k) for_theta(k, active, XXi_k, Xi_k, Rk, tau))
  theta_k
}

get_trunc_variable <- function(mu,sd,a,b) qnorm(runif(length(mu),pnorm(a,mu,sd),pnorm(b,mu,sd)),mu,sd)

sample_L <- function(ID, C_star_sbj_k, theta_sbj_k, Z_sbj, X_sbj, bi, beta, lower.tn, upper.tn){
  mu.tn <-  map2(C_star_sbj_k, theta_sbj_k, ~ .x %*% .y) %>% unlist + 
    map2(Z_sbj, bi, ~ .x %*% t(.y)) %>% unlist + 
    map(X_sbj, function(x) x%*%beta) %>% unlist
  Li.vec <- get_trunc_variable(mu.tn, rep(1,length(mu.tn)), lower.tn, upper.tn)  
  L_sbj <- split(Li.vec, f=ID)
  L_sbj
}

sample_cluster <- function(out_mat, K.max) {
  out_mat <- do.call(cbind, map(1:K.max, ~ out_mat[[1]][[.x]] %>% unlist))
  obtain_max <- apply(out_mat, 1, function(x) max(x))
  out_mat <- exp(out_mat - obtain_max)
  Clusters <- apply(out_mat, 1, function(x) sample(1:K.max, 1, prob=x))  
  Clusters
}

sample_beta <- function(kappa, q, XSum, NSbj, X_sbj, L_sbj, C_star_sbj_k, theta_sbj_k, Z_sbj, bi) {
  Cov_beta <- solve(diag(q)/kappa+XSum)
  mu_beta <- map(1:NSbj, ~  t(X_sbj[[.x]])%*%(L_sbj[[.x]] - C_star_sbj_k[[.x]]%*%theta_sbj_k[[.x]] - Z_sbj[[.x]]%*%as.vector(bi[[.x]]) )) %>% reduce(`+`)
  beta <- as.vector(mvtnorm::rmvnorm(1,mean=Cov_beta%*%mu_beta,sigma=Cov_beta,method='eigen'))
  beta
}


get_theta_gamma_matrix <- function(p, unit_gamma, unit_theta){
  positions <- c(1,apply(unit_gamma, 1, sum))
  indexing <- map(1:p, function(i) {
    start <- positions[1:i] %>% sum
    end <- start + positions[i+1] - 1
    c(start,end)
  })
  
  do.call(rbind, map(1:p, function(i) {
    index <- indexing[[i]]
    temp_gamma <- unit_gamma[i,]
    temp_gamma[temp_gamma==1] <- unit_theta[index[1]:index[2]]
    temp_gamma
  }))
}

get_post_summary_vc <- function(t, time_range, knot_position, gamma, theta, grids=50, burn=NULL) {
  
  len <- length(gamma)
  burn <- ifelse(is.null(burn), round(len/2, 0), burn) 
  time_domain <- seq(from=time_range[1], to=time_range[2], l=grids)
  
  sd_ <- sd(t)
  K <- length(gamma[[1]])
  p <- dim(gamma[[1]][[1]])[[1]]
  
  
  basis_functions <- map(time_domain, ~ matrix(rep(B(.x, knot_position, sd_),p),nrow=p, byrow=T )) 
  
  # Time => MCMC => Cluster
  Bgamma <- map(1:length(time_domain), function(time) {
    map(gamma[(burn+1):len], function(gamma_)  map(gamma_, function(each_cls) {
      bgam_mat <- each_cls * basis_functions[[time]] # p times knots matrices
      bgam_mat
    }))
  })
  
  Tgamma <-  map((burn+1):len, function(MC) {
    map(1:K, function(k) get_theta_gamma_matrix(p, gamma[[MC]][[k]], theta[[MC]][[k]]) %>% t)  
  })
  
  vc_summary <-  map(1:K, function(k) {
    
    map(1:length(time_domain), function(time) {
      
      temp_mat<- do.call(rbind, 
                         map(1:length((burn+1):len), function(MC) {
                           diag(Bgamma[[time]][[MC]][[k]] %*% Tgamma[[MC]][[k]])
                         }))
      
      apply(temp_mat,2,function(x) quantile(x, c(0.025,0.5,0.975)))
      
    })
  })
  
  list(vc_summary = vc_summary,
       time_domain = time_domain)
}





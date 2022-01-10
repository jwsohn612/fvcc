#' @import dplyr
#' @import purrr
#' @importFrom truncnorm rtruncnorm
NULL


#' Clustering Variying Coefficients with a Pitman-Yor Process
#'
#' @param ID A vector of subjects' IDs.
#' @param W A part of design matrix for varying coefficients
#' @param X A part of design matrix for time-invariant fixed effects
#' @param Z A part of design matrix for random effects
#' @param t A vector of time observations
#' @param Y A vector of binary responses
#' @param num_knots The number of knots for the basis splines
#' @param K_max The upper bound of the number of clusters for the stick-breaking process
#' @param num_iters The number of iterations for the partially collapsed Gibbs sampler
#' @param thinning_unit Thinning unit of the pcg sampler
#' @param nu An initial value of the concentration parameter
#' @param lambda An initial value of the discount parameter
#' @param Psi An initial covariance matrix of the random effects
#' @param beta An initial vector of the fixed effects
#'
#' @examples
#' library(cvarpyp)
#'
#' #SimulationData is an example dataset.
#'
#' ID <- SimulationData$ID
#' W <- SimulationData[c('W1','W2')]
#' X <- SimulationData[c('X1','X2')]
#' Z <- SimulationData[c('Z1','Z2')]
#' t <- SimulationData$t
#' Y <- SimulationData$Y
#'
#' output <- fcvarpyp(ID = ID,
#'                W = W,
#'                X = X,
#'                Z = Z,
#'                t = t,
#'                Y = Y,
#'                num_knots = 15,
#'                K_max = 10,
#'                num_iters = 10000)
#'
#' plot(output$fixed_effect, position = 1)
#' plot(output$random_effect, row_position = 1, col_position = 1)
#' plot(output$cluster)
#' plot(output$varying_coefficient, cluster_number = 5, variable_number = 1)
#' plot(output$latent_location, cluster_number = 5, variable_number = 1, time_range = output$time_range, knot_position = output$knot_position)
#' plot(output$nu)
#' plot(output$lambda)
#'
#' @export
fcvarpyp <- function(ID,
                     W,
                     X = NULL,
                     Z,
                     t,
                     Y,
                     num_knots = 10,
                     K_max,
                     num_iters = 10000,
                     thinning_unit = 5,
                     nu = 1,
                     lambda = 0.5,
                     Psi = NULL,
                     beta = NULL,
                     kappa = 10,
                     stick_act = c(1, 10),
                     stick_nact = c(1, 10)) {
  temp <- preprocess_data(ID, W, X, Z, t, Y)
  ID <- temp$ID
  W <- temp$W
  X <- temp$X
  Z <- temp$Z
  t <- temp$t
  Y <- temp$Y
  rm(temp)
  
  num_obs_per_sbj <- table(ID)
  NSbj <- max(ID)
  one_obsr <- which(num_obs_per_sbj == 1)
  
  # translate to matrix values
  p <- ncol(W)
  r <- ncol(Z)
  q <- ncol(X)
  
  Z_sbj <- map(split(Z %>% data.matrix, f = ID), ~ matrix(.x, ncol = r))
  X_sbj <- map(split(X %>% data.matrix, f = ID), ~ matrix(.x, ncol = q))
  XSum <-
    map(1:NSbj, ~ t(X_sbj[[.x]]) %*% X_sbj[[.x]]) %>% reduce(`+`)
  
  
  lower.tn <- log(Y)
  upper.tn <- -log(1 - Y)
  
  sd_nu <- 0.5
  sd_lambda <- 0.5
  
  a_pi <- 1
  b_pi <- 1  # noninformative
  nu_param <- 1
  
  # rinvGamma
  g_k <- 1/2
  h_k <- NSbj/2
  
  knot_position <- quantile(t, ppoints(num_knots, a = 0))
  
  
  # Check Psi
  if (is.null(Psi)) {
    Psi <- diag(r)
  } else {
    if (is.matrix(Psi)) {
      if (((ncol(Psi) == nrow(Psi)) == r)) {
        Psi <- Psi
      } else {
        print(paste0(
          "Psi should be a positive definite matrix with",
          "(",
          r,
          ",",
          r,
          ")"
        ))
      }
    } else {
      print(paste0(
        "Psi should be a positive definite matrix with",
        "(",
        r,
        ",",
        r,
        ")"
      ))
    }
  }
  
  
  Psi_0 <- diag(r)
  
  if (is.null(beta)) {
    beta <- rnorm(q, sd = 1 / sqrt(q))
  }
  
  
  bi <- map(1:NSbj, ~ MASS::mvrnorm(1, mu = rep(0, r), Sigma = Psi))
  Li <-
    map_dbl(Y, function(sign)
      if (sign > 0) {
        rtruncnorm(1, a = 0)
      } else{
        rtruncnorm(1, b = 0)
      }) %>% as.matrix()
  L_sbj <- split(Li, f = ID)
  
  pu <- r
  D <- diag(r)
  
  g_vec <- rep(g_k, K_max)
  h_vec <- rep(h_k, K_max)
  tau <- MCMCpack::rinvgamma(n = K_max,
                             shape = g_k / 2,
                             scale = h_k / 2)
  gamma <- map(1:K_max, function(x) {
    temp <-
      matrix(rbinom(
        n = p * (num_knots + 2),
        size = 1,
        prob = 1 / (num_knots + 2)
      ), nrow = p)
    temp[, 1] <- 1
    temp
  })
  
  v <- c(map_dbl(1:(K_max - 1), ~ rbeta(1, 1, 1)), 1)
  Diri_p <- get_DP_pi(v = v, K = K_max)
  Clusters <- sample(1:K_max, size = NSbj, replace = T)
  
  # Cluster status
  active_ <- 1:K_max %in% Clusters
  active <- which(active_ == T)
  inactive <- which(active_ == F)
  
  theta_k <-
    map(1:K_max, function(k)
      MASS::mvrnorm(1, mu = rep(0, sum(gamma[[k]])), Sigma = diag(sum(gamma[[k]]))))
  
  W_star <-
    get_basis_data(
      tData = W %>% data.matrix,
      t = t,
      knot_position = knot_position,
      sd = sd(t)
    )
  C_star <- do.call(cbind, W_star)
  rm(W_star)
  
  C_star_sbj <-
    split(C_star %>% data.matrix, f = ID) %>% map(~ matrix(.x, ncol = dim(C_star)[2]))
  
  temp <-
    byproducts(
      num_obs_per_sbj,
      K_max,
      NSbj,
      C_star_sbj,
      gamma,
      X_sbj,
      beta,
      Z_sbj,
      L_sbj,
      Psi,
      Psi_0,
      ID - 1,
      Clusters,
      active,
      flatten_gamma = flatten_gamma_with_fixC
    )
  InvAdj <- temp[[1]]
  InvAdj_0 <- temp[[2]]
  XXi_k <- temp[[3]]
  Xi_k <- temp[[4]]
  
  temp <-
    MakeRk(K_max,
           tau,
           active,
           C_star_sbj,
           InvAdj,
           InvAdj_0,
           gamma,
           NSbj,
           flatten_gamma = flatten_gamma_with_fixC)
  Rk <- temp[[1]]
  
  updated_theta <- list()
  updated_beta <- list()
  updated_gamma <- list()
  updated_cluster <- list()
  updated_tau <- list()
  updated_Psi <- list()
  updated_nu <- c()
  updated_lambda <- c()
  
  update_indi <- 1
  
  pb <- progress::progress_bar$new(format = "  running the pcg [:bar] :percent eta: :eta",
                                   clear = FALSE,
                                   total = num_iters)
  
  
  for (iter in 1:num_iters) {
    pb$tick()
    
    # ================== activation_clusters,boolean ===================== #
    active_ <- 1:K_max %in% Clusters
    active <- which(active_ == T)
    inactive <- which(active_ == F)
    
    gamma <-
      sample_gamma(
        ID,
        K_max,
        p,
        active,
        num_knots,
        a_pi,
        b_pi,
        NSbj,
        gamma,
        tau,
        C_star_sbj,
        InvAdj,
        L_sbj,
        X_sbj,
        beta,
        Clusters
      )
    
    # Draw stick breaking beta
    num_of_ones_K_sbjs <- get_sbj_clusters(K_max, Clusters)
    v <-
      c(map_dbl(
        1:(K_max - 1),
        ~ rbeta(
          1,
          1 - lambda + num_of_ones_K_sbjs[.x],
          nu + .x * lambda + sum(num_of_ones_K_sbjs[(.x + 1):K_max])
        )
      ) , 1)
    Diri_p <- get_DP_pi(v, K_max)
    
    # Draw nu
    tryCatch({
      nu <- sample_nu(nu, lambda, v, sd_nu, K_max, nu_param)
    }, error = function(e) {
      nu <- tolerance::r2exp(n = 1,
                             rate = 1,
                             shift = -lambda)
      
    })
    
    # Draw lambda
    lambda <- tryCatch({
      lambda <- sample_lambda(lambda, nu, v, K_max, sd_lambda)
    }, error = function(e) {
      while (nu <= -lambda) {
        lambda <- runif(1)
      }
    })
    
    temp <-
      byproducts(
        num_obs_per_sbj,
        K_max,
        NSbj,
        C_star_sbj,
        gamma,
        X_sbj,
        beta,
        Z_sbj,
        L_sbj,
        Psi,
        Psi_0,
        ID - 1,
        Clusters,
        active,
        flatten_gamma = flatten_gamma_with_fixC
      )
    InvAdj <- temp[[1]]
    InvAdj_0 <- temp[[2]]
    XXi_k <- temp[[3]]
    Xi_k <- temp[[4]]
    
    temp <-
      MakeRk(K_max,
             tau,
             active,
             C_star_sbj,
             InvAdj,
             InvAdj_0,
             gamma,
             NSbj,
             flatten_gamma = flatten_gamma_with_fixC)
    Rk <- temp[[1]]
    
    C_star_sbj_k <-
      map(1:NSbj, ~ C_star_sbj[[.x]][, which(flatten_gamma_with_fixC(gamma[[Clusters[.x]]]) ==
                                               TRUE)])
    
    # Draw theta
    theta_k <- sample_theta(K_max, active, XXi_k, Xi_k, Rk, tau)
    theta_sbj_k <- map(1:NSbj, ~ theta_k[[Clusters[[.x]]]])
    
    # Draw tau
    num_of_ones_vars <- map_dbl(1:K_max, ~ sum(gamma[[.x]]))
    g_vec <-
      map_dbl(1:K_max, ~ ifelse((.x %in% active), g_k, g_k))
    h_vec <-
      map_dbl(1:K_max, ~ ifelse((.x %in% active), h_k, h_k))
    tau <-
      map_dbl(1:K_max,
              ~ MCMCpack::rinvgamma(
                1,
                shape = (num_of_ones_vars[.x] + g_vec[.x]) / 2,
                scale = h_vec[.x] / 2 + (theta_k[[.x]] %*% Rk[[.x]] %*% theta_k[[.x]]) /
                  2
              ))
    
    # Draw bi
    temp <-
      sample_randomeffects(
        1,
        Z_sbj,
        InvAdj,
        L_sbj,
        C_star_sbj,
        gamma,
        X_sbj,
        beta,
        Psi,
        theta_k,
        Clusters,
        r,
        flatten_gamma = flatten_gamma_with_fixC
      )
    bi <- temp[[1]]
    bi_mats <- temp[[2]]
    
    # Draw Psi
    proposed_Psi <- MCMCpack::riwish(pu + NSbj, D + bi_mats)
    probs_Psi <-
      MakeRkPsi(
        num_obs_per_sbj,
        K_max,
        Rk,
        tau,
        active,
        C_star_sbj,
        Z_sbj,
        proposed_Psi,
        gamma,
        theta_k,
        NSbj,
        flatten_gamma = flatten_gamma_with_fixC
      )
    if (runif(1) < sum(unlist(probs_Psi[[1]])))
      Psi <- proposed_Psi
    
    Psi_0 <- Psi
    
    # Draw Li
    L_sbj <-
      sample_L(ID,
               C_star_sbj_k,
               theta_sbj_k,
               Z_sbj,
               X_sbj,
               bi,
               beta,
               lower.tn,
               upper.tn)
    
    # Draw cluster
    out_mat <-
      ObtainOutput(
        NSbj,
        K_max,
        Diri_p ,
        L_sbj,
        C_star_sbj,
        X_sbj,
        beta,
        Z_sbj,
        bi,
        gamma,
        theta_k,
        flatten_gamma = flatten_gamma_with_fixC
      )
    Clusters <- sample_cluster(out_mat, K_max)
    
    # Draw beta
    if (!is.null(X)) {
      beta <-
        sample_beta(kappa,
                    q,
                    XSum,
                    NSbj,
                    X_sbj,
                    L_sbj,
                    C_star_sbj_k,
                    theta_sbj_k,
                    Z_sbj,
                    bi)
    } else {
      beta <- 0
    }
    
    # ======================= store variables ======================#
    if (iter %% thinning_unit == 0) {
      updated_gamma <- rlist::list.append(updated_gamma, gamma)
      updated_cluster <-
        rlist::list.append(updated_cluster, Clusters)
      updated_theta <- rlist::list.append(updated_theta, theta_k)
      updated_tau <- rlist::list.append(updated_tau, tau)
      updated_Psi <- rlist::list.append(updated_Psi, Psi)
      updated_lambda <- append(updated_lambda, lambda)
      updated_nu <- append(updated_nu, nu)
      
      if (!is.null(X))
        updated_beta <- rlist::list.append(updated_beta, beta)
    }
  }
  
  print("Generating model results..")
  varying_coefficient <- get_post_summary_vc(
    t = t,
    time_range = c(min(t), max(t)),
    knot_position = knot_position,
    gamma = updated_gamma,
    theta = updated_theta,
    grids = 100,
    burn = NULL
  )
  
  class(varying_coefficient) <- 'varying_coefficient'
  class(updated_gamma) <- 'latent_location'
  class(updated_cluster) <- 'cluster'
  class(updated_theta) <- 'varying_coefficient'
  class(updated_tau) <- 'dispersion'
  class(updated_Psi) <- 'random_effect'
  class(updated_beta) <- 'fixed_effect'
  class(updated_lambda) <- 'lambda'
  class(updated_nu) <- 'nu'
  
  list(
    latent_location = updated_gamma,
    cluster = updated_cluster,
    varying_coefficient = varying_coefficient,
    random_effect = updated_Psi,
    fixed_effect = updated_beta,
    dispersion = updated_tau,
    knot_position = knot_position,
    nu = updated_nu,
    lambda = updated_lambda,
    time_range = c(min(t), max(t)),
    p = p,
    q = ifelse(!is.null(X), dim(X)[2], 0),
    r = r
  )
}

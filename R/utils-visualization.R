#' @import dplyr
#' @import rlist
#' @import progress
#' @import purrr
#' @import label.switching
NULL

#' Draw the posterior results of random effects
#'
#' @export
plot.random_effect <-
  function(random_effect_samples,
           row_position,
           col_position,
           burn = NULL) {
    len <- length(random_effect_samples)
    if (is.null(burn))
      burn = round(len / 2)
    post_sample <-
      map_dbl(random_effect_samples, ~ .x[row_position, col_position])
    post_sample <- post_sample[(burn + 1):len]
    
    # par(mfrow=c(1,2))
    hist_ = hist(
      post_sample,
      breaks = 10,
      probability = T,
      main = "",
      xlab = "",
      ylab = "",
      col = gray(0:9 / 9)[8],
      border = gray(0:9 / 9)[8]
    )
    box()
    lines(hist_$breaks, c(hist_$density, 0), type = "s")
    abline(v = quantile(post_sample, c(0.025, 0.975)), lty = 2)
    abline(v = median(post_sample))
    mtext("Posterior density",
          side = 2,
          line = 2.3,
          cex = 0.9)
    
    # acf(post_sample,main="")
  }

#' Draw the posterior results of fixed effects
#'
#' @export
plot.fixed_effect <-
  function(fixed_effect_samples, position, burn = NULL) {
    len <- length(fixed_effect_samples)
    if (is.null(burn))
      burn = round(len / 2)
    post_sample <- map_dbl(fixed_effect_samples, ~ .x[position])
    post_sample <- post_sample[(burn + 1):len]
    
    # par(mfrow=c(1,2))
    hist_ = hist(
      post_sample,
      breaks = 10,
      probability = T,
      main = "",
      xlab = "",
      ylab = "",
      col = gray(0:9 / 9)[8],
      border = gray(0:9 / 9)[8]
    )
    box()
    lines(hist_$breaks, c(hist_$density, 0), type = "s")
    abline(v = quantile(post_sample, c(0.025, 0.975)), lty = 2)
    abline(v = median(post_sample))
    mtext("Posterior density",
          side = 2,
          line = 2.3,
          cex = 0.9)
    
    # acf(post_sample,main="")
  }


#' Draw the posterior results of varying coefficients
#'
#' @export
plot.varying_coefficient <-
  function(varying_coefficient_samples,
           cluster_number,
           variable_number) {
    vc_summary <- varying_coefficient_samples$vc_summary
    time_domain <- varying_coefficient_samples$time_domain
    
    vc_mat <-
      do.call(rbind, map(vc_summary[[cluster_number]], ~ .x[, variable_number]))
    L <- vc_mat[, 1]
    M <- vc_mat[, 2]
    U <- vc_mat[, 3]
    
    poly_range <- c(time_domain, rev(time_domain))
    poly_coef_UL <- c(L, rev(U))
    
    plot(
      NULL,
      type = "l",
      ylim = c(min(vc_mat), max(vc_mat)),
      xlim = c(min(time_domain), max(time_domain)),
      xlab = "",
      ylab = ""
    )
    polygon(poly_range,
            poly_coef_UL,
            col = gray(0:9 / 9)[8],
            border = F)
    lines(time_domain, M, lty = 1)
    
    # mtext(paste0('cluster',which_clusters,';varying',j),side=2,line=2.3,cex=0.9)
    # mtext("t",side=1,line=2.3,cex=0.9)
  }

# Time => Cluster => VC
# vc_object <- get_post_summary_vc(from=0, to=1, knots, gamma.sample, theta.sample, grids=50,  burn=NULL)
# summarise.varying.coefficient(vc_object, cluster_number = 1, variable_number = 3)

#' Draw the posterior results of latent location parameters
#'
#' @export
plot.latent_location <-
  function(latent_location_samples,
           cluster_number,
           variable_number,
           time_range,
           knot_position,
           burn = NULL) {
    len <- length(latent_location_samples)
    if (is.null(burn))
      burn = round(len / 2)
    
    location_mat <-
      do.call(rbind, map(latent_location_samples[(burn + 1):len], ~ .x[[cluster_number]][variable_number, ]))
    proba <- apply(location_mat, 2, mean)
    
    df <-
      data.frame(knot_position = c(time_range[1], knot_position),
                 freq = proba[2:(length(knot_position) + 2)])
    
    plot(
      df$knot_position[2:(length(knot_position) + 1)],
      df$freq[2:(length(knot_position) + 1)],
      type = 'h',
      ylim = c(0, 1),
      xlim = c(time_range[1], time_range[2]),
      lwd = 3,
      xlab = "",
      ylab = ""
    )
    points(0, df$freq[1], pch = 1)
    segments(0, 0, df$knot_position[1], df$freq[1], lty = 2)
    # mtext("klp",side=2,line=2.5,cex=0.9)
    # mtext("t",side=1,line=2.5,cex=0.9)
  }

# summarise.latent.location(time_range=c(0,1), gamma.sample, cluster_number=1, variable_number=2, burn=NULL)
#' Draw the posterior results of clusters
#'
#' @export
plot.cluster <- function(cluster_sample, burn = NULL) {
  len <- length(cluster_sample)
  if (is.null(burn))
    burn = round(len / 2)
  cluster_matrix <- do.call(rbind, cluster_sample)
  apply(cluster_matrix[(burn + 1):len, ], 2, function(x) {
    mode.value <- DescTools::Mode(x)
    mode.value[1]
  }) %>% unlist()
}

#' Draw the posterior histogram of the concentration parameter
#'
#' @export
plot.nu <- function(nu_sample, burn = NULL) {
  len <- length(nu_sample)
  if (is.null(burn))
    burn = round(len / 2)
  
  post_sample <- nu_sample[(burn + 1):len]
  
  # par(mfrow=c(1,2))
  hist_ = hist(
    post_sample,
    breaks = 10,
    probability = T,
    main = "",
    xlab = "",
    ylab = "",
    col = gray(0:9 / 9)[8],
    border = gray(0:9 / 9)[8]
  )
  box()
  lines(hist_$breaks, c(hist_$density, 0), type = "s")
  abline(v = quantile(post_sample, c(0.025, 0.975)), lty = 2)
  abline(v = median(post_sample))
  mtext("Posterior density",
        side = 2,
        line = 2.3,
        cex = 0.9)
  
  # acf(post_sample,main="")
}

#' Draw the posterior histogram of the discount parameter
#'
#' @export
plot.lambda <- function(lambda_sample, burn = NULL) {
  len <- length(lambda_sample)
  if (is.null(burn))
    burn = round(len / 2)
  
  post_sample <- lambda_sample[(burn + 1):len]
  
  hist_ = hist(
    post_sample,
    breaks = 10,
    probability = T,
    main = "",
    xlab = "",
    ylab = "",
    col = gray(0:9 / 9)[8],
    border = gray(0:9 / 9)[8]
  )
  box()
  lines(hist_$breaks, c(hist_$density, 0), type = "s")
  abline(v = quantile(post_sample, c(0.025, 0.975)), lty = 2)
  abline(v = median(post_sample))
  mtext("Posterior density",
        side = 2,
        line = 2.3,
        cex = 0.9)
}

init_prog_var <- function(n_iter){
  pb <- txtProgressBar(min = 0,      
                       max = n_iter, 
                       style = 3,    
                       width = 50,   
                       char = "=")   
}

#' @export
get_vc_results <- function(save_dir, sbj1, sbj2, sbj3){
  
  rdata_list <- list.files(save_dir)
  replicates <- length(rdata_list)
  
  pb <- init_prog_var(replicates)
  
  varying_coefs_list <- list()
  knot_selection_list <- list()
  cluster_list <- list()
  
  for(i in 1:length(rdata_list)){
    
    output <- rlist::list.load(paste0(save_dir, rdata_list[i]))
    saved_sequence_size <- length(output$cluster)
    
    p <- output$p
    K <- output$K_max
    
    vc <- output$varying_coefficient$vc_summary
    
    # Relabeling Cluster
    cluster_matrix <- do.call(rbind, output$cluster)
    cluster <- unlist(apply(cluster_matrix, 2,function(x) DescTools::Mode(x)[1]))
    
    es_label1 <- as.numeric(names(which.max(table(cluster[1:sbj1]))))
    es_label2 <- as.numeric(names(which.max(table(cluster[(sbj1+1):(sbj1+sbj2)]))))
    es_label3 <- as.numeric(names(which.max(table(cluster[(sbj1+sbj2+1):(sbj1+sbj2+sbj3)]))))
    
    trivial_cluster <- sort(unique(cluster[!(cluster == es_label1 | cluster == es_label2 | cluster == es_label3)]))
    
    if(is.null(trivial_cluster)){
      es_label <- unique(c(es_label1, es_label2, es_label3))
    }else{
      es_label <- unique(c(es_label1, es_label2, es_label3, trivial_cluster))
    }
    
    true_label <- 1:length(es_label)
    new_label <- sapply(cluster, function(x) true_label[es_label==x]) %>% as.vector()
    table(new_label, cluster)
    names(es_label) <- true_label
    
    cluster_list[[i]] <- new_label
    
    get_cluster_number <- function(x) es_label[names(es_label)==x]
    
    varying_coefs <- purrr::map(1:K, function(k){
      purrr::map(1:p, function(l){
        varing_coefs_each_covariate <- purrr::map(vc[[k]], ~ .x[,l])
        do.call(rbind, varing_coefs_each_covariate)
      })
    })
    
    
    # Summarize
    # Varying Coefficients
    varying_coefs_list[[i]] <- list(varying_coefs[get_cluster_number(1)],
                                    varying_coefs[get_cluster_number(2)],
                                    varying_coefs[get_cluster_number(3)])
    
    setTxtProgressBar(pb, i)
  }
  
  list(varying_coefs_list = varying_coefs_list, 
       cluster_list = cluster_list
  )
}

get_main_clusters <- function(cluster){
  
  Cluster_matrix = do.call(rbind, cluster)
  cclusters <- unlist(apply(Cluster_matrix,2,function(x) DescTools::Mode(x)[1]))
  
  active <- names(table(cclusters))
  num_sbjs_active <- table(cclusters)[active]
  active <- as.numeric(active)
  
  cclusters
}

#' @export
relabel_chain <- function(output){
  
  cluster <- output$cluster
  gamma <- output$gamma
  theta <- output$theta
  
  K <- length(theta[[1]])
  z <- do.call(rbind, cluster) 
  
  run <- label.switching::ecr.iterative.1(z = z, K = K)
  map <- run$permutations
  
  agamma <- purrr::map(1:nrow(map), function(i){
    gamma[[i]][map[i,]]
  })
  
  atheta <- purrr::map(1:nrow(map), function(i){
    theta[[i]][map[i,]]
  })
  
  acluster <- purrr::map(1:nrow(map), function(i){
    map[i,][cluster[[i]]]
  })
  
  return(list(gamma = agamma, 
              theta = atheta, 
              cluster = acluster,
              knots = output$knot_position,
              sdd = output$sd_scaler))
}

#' @export
make_vc_outputs_on_knots <- function(output){
  "
  All parameters should be burn in advance. 
  "
  gamma <- output$gamma
  theta <- output$theta
  cluster <- output$cluster
  
  K <- gamma[[1]] %>% length
  knots <- output$knots
  sdd <- output$sdd
  
  BasisFunctionsSets_onKnots <- purrr::map(knots, function(t) B(t, knots, sd = sdd))
  
  # In our example, the clustering order is consistent in accordance with the shape of varying coefficients.
  m_clusters <- get_main_clusters(cluster) 
  top2cls <- sort(table(m_clusters), decreasing = TRUE)[1:2]
  len <- gamma %>% length
  print(top2cls)
  
  outputs <- 
    purrr::map(names(top2cls) %>% as.numeric, function(k){
      each_K <- purrr::map(1:len, function(iter){
        along_time <- purrr::map(1:length(knots), function(t) {
          
          # Taking phis
          temp_gamma <- t(gamma[[iter]][[k]])
          temp_gamma[temp_gamma==1] <- theta[[iter]][[k]]
          
          # Taking Basis functions
          Btemp_gamma <- gamma[[iter]][[k]]
          Btemp_gamma[1,][Btemp_gamma[1, ]==1] <- BasisFunctionsSets_onKnots[[t]][gamma[[iter]][[k]][1,]==1]
          Btemp_gamma[2,][Btemp_gamma[2, ]==1] <- BasisFunctionsSets_onKnots[[t]][gamma[[iter]][[k]][2,]==1]
          
          diag(Btemp_gamma %*% temp_gamma)
        }) 
        do.call(rbind, along_time)
      }) 
      
      each_K <- list(
        do.call(rbind,purrr::map(each_K, function(mat) mat[,1])),
        do.call(rbind,purrr::map(each_K, function(mat) mat[,2]))
      )
    })
  
  graphics.off()
  par(mfrow=c(2,2),cex=0.5)
  drawing_varying_plots(outputs, cluster = 1, predictor = 1, cls_name = 1, timedomain = knots)
  drawing_varying_plots(outputs, cluster = 2, predictor = 1, cls_name = 2, timedomain = knots)
  drawing_varying_plots(outputs, cluster = 1, predictor = 2, cls_name = 1, timedomain = knots)
  drawing_varying_plots(outputs, cluster = 2, predictor = 2, cls_name = 2, timedomain = knots)
  return(outputs) 
}


drawing_varying_plots <- function(outputs, cluster, predictor, cls_name, y_range=NULL, timedomain){
  
  poly_range <- c(timedomain, rev(timedomain))
  varrr <- outputs[[cluster]][[predictor]] 
  
  med <- apply(varrr, 2, function(x) quantile(x, 0.5))
  lower <- apply(varrr, 2, function(x) quantile(x, 0.025))
  upper <- apply(varrr, 2, function(x) quantile(x, 0.975))
  
  poly_coef_UL <- c(lower, rev(upper))
  if(is.null(y_range)){
    plot(NULL,type="l",ylim=c(min(varrr),max(varrr)),xlim=range(timedomain), xlab="",ylab="")  
  }
  else{
    plot(NULL,type="l",xlim=range(timedomain), ylim=y_range, xlab="",ylab="")
  }
  polygon(poly_range, poly_coef_UL,col=gray(0:9/9)[8],border=F)
  lines(timedomain,med,lty=1)
  abline(h=0,lty=2)
  
  mtext(paste0('v',cls_name,predictor),side=2,line=2.8,cex=0.9)
  mtext("t",side=1,line=2.8,cex=0.9)
  
}



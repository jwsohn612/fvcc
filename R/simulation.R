#' Get simulation data
#' 
#' @param N The number of subjects
#' @param K The number of clusters
#' @param M The average number of observations
#' 
#' @export 
get_simulation_data <- function(N, K, M, seed=0){
  
  set.seed(seed)
  
  N <- N
  K <- 3 
  
  p <- 4 
  q <- 2 
  r <- 2
  
  ID <- 1:N
  ni = rpois(N, M)
  
  n = sum(ni)
  
  id=NULL
  for(i in 1:N) id=c(id, rep(i,ni[i]))
  dataset=data.frame(id,
                     rep(1,n),
                     matrix(rnorm(n*(p-1),0,1),n,p-1),
                     matrix(rnorm(n*(q),0,1),n,q),
                     rep(1,n),
                     matrix(rnorm(n*(r-1),0,1),n,r-1),
                     matrix(runif(n,0,1),n,1))
  names(dataset)=c("ID","var1","var2","var3","var4","fix1","fix2","re1","re2","t")
  # names(dataset)=c("ID","var1","var2","var3","fix1","fix2","re1","re2","t")

  sbj1 <- round(N/3, 0) 
  sbj2 <- round(N/3, 0)
  sbj3 <- round(N/3, 0)
  
  clusters <- c(rep(1,sbj1),rep(2,sbj2),rep(3,sbj3))
  dataset$G <- unlist(sapply(1:N, function(x) rep(clusters[x],ni[x])))
  
  g1_a1<-function(x) sin(2*pi*x)
  g1_a2<-function(x) 2*exp(-200*(x-0.2)^2)+exp(-10*(x-0.6)^2)
  g1_a3<-function(x) sin(2*pi*x^3)
  g1_a4<-function(x) -1+0*x
  
  g2_a1<-function(x) cos(2*pi*x)
  g2_a2<-function(x) 2*exp(-200*(x-0.2)^2)+exp(-10*(x-0.6)^2)
  g2_a3<-function(x) 2*x
  g2_a4<-function(x) 0*x
  
  g3_a1<-function(x) sin(8*(x-0.5))+1.5*exp(-(20^2)*(x-0.5)^2)
  g3_a2<-function(x) 2*exp(-200*(x-0.2)^2)+exp(-10*(x-0.6)^2)
  g3_a3<-function(x) 2+0*x
  g3_a4<-function(x) 1+0*x
  
  alpha<-list(
    list(g1_a1,g1_a2,g1_a3,g1_a4),
    list(g2_a1,g2_a2,g2_a3,g2_a4),
    list(g3_a1,g3_a2,g3_a3,g3_a4)
  )
  
  beta <- c(1, -1)
  D <- matrix(c(0.5 , 0.25 , 0.25 , 0.8), nrow=2)
  bi <- mvtnorm::rmvnorm(N, mean=c(0,0), sigma = D)
  
  sigma <- 1
  Li<-c() 
  Yi<-c()
  for(i in 1:n){
    the_cluster <- dataset$G[i]
    
    id_ <- id[i]
    Li[i] <- 
      dataset$var1[i]*alpha[[the_cluster]][[1]](dataset$t[i])+
      dataset$var2[i]*alpha[[the_cluster]][[2]](dataset$t[i])+
      dataset$var3[i]*alpha[[the_cluster]][[3]](dataset$t[i])+
      dataset$var4[i]*alpha[[the_cluster]][[4]](dataset$t[i])+
      dataset$fix1[i]*beta[1]+
      dataset$fix2[i]*beta[2]+
      dataset$re1[i]*bi[id_,1]+
      dataset$re2[i]*bi[id_,2] + 
      rnorm(1, mean=0, sd = sigma)
    
    if(Li[i]>=0){Yi[i]<-1}else{Yi[i]<-0}
  }
  dataset$Y <- Yi
  
  dataset
}

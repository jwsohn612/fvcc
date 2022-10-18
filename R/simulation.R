#' Get simulation data
#' 
#' @param sbj1 The number of subjects of the first group
#' @param sbj2 The number of subjects of the second group
#' @param sbj3 The number of subjects of the third group
#' @seed a random seed
#' 
get_simulation_data <- function(sbj1, sbj2, sbj3, seed=0){
  
  set.seed(seed)
  
  N <- sbj1 + sbj2 + sbj3
  K <- 3
  
  p <- 2
  q <- 2
  r <- 2
  
  ID <- 1:N
  ni = rpois(N, 10)
  
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
  names(dataset)=c("ID","var1","var2","fix1","fix2","re1","re2","t")
  # names(dataset)=c("ID","var1","var2","var3","fix1","fix2","re1","re2","t")

  clusters <- c(rep(1,sbj1),rep(2,sbj2),rep(3,sbj3))
  dataset$G <- unlist(sapply(1:N, function(x) rep(clusters[x],ni[x])))
  
  g1_a1<-function(x) 2*exp(-200*(x-0.2)^2)+exp(-10*(x-0.6)^2)
  g1_a2<-function(x) sin(2*pi*x^3)
  
  g2_a1<-function(x) sin(8*(x-0.5))+1.5*exp(-(20^2)*(x-0.5)^2)
  g2_a2<-function(x) 2*x + 0 
  
  g3_a1<-function(x) -2*x + 0
  g3_a2<-function(x) 0*x + 0
  
  alpha<-list(
    list(g1_a1,g1_a2),
    list(g2_a1,g2_a2),
    list(g3_a1,g3_a2)
  )
  
  beta <- c(1,-1)
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

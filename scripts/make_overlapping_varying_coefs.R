library(dplyr)
library(ggplot2)
library(cowplot)
library(fvcc)

# Place your directory
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

draw_overlapping_varcoefs <- function(varying_coefs_list, cls, pre, ylim=NULL, yname){
  len <- varying_coefs_list %>% length
  temp <- purrr::map_df(1:len, function(i){
    
    temp_df <- data.frame(varying_coefs_list[[i]][[cls]][[1]][[pre]])
    colnames(temp_df) <- c("L","M","U")
    
    temp_df$t <- timedomain
    temp_df$g <- rep(as.character(i),length(timedomain))
    
    temp_df
  })
  
  true_df <- data.frame(t=timedomain, y=alpha[[cls]][[pre]](timedomain))
  true_df$g <- 'true'
  
  temp_plt <- temp %>% 
    ggplot(aes(x=t, y=M, group=g)) + 
    theme_classic() + 
    geom_line(alpha=0.10) + ylab(yname) + 
    geom_line(data = true_df, mapping=aes(x=t, y=y, group=g), color='red')
  
  if(!is.null(ylim)){
    temp_plt <- temp_plt + ylim(ylim) 
  }
  return(temp_plt)
}

v11 <- draw_overlapping_varcoefs(output$varying_coefs_list, cls = 1, pre = 1, yname = 'v11')
v12 <- draw_overlapping_varcoefs(output$varying_coefs_list, cls = 1, pre = 2, yname = 'v12')
v21 <- draw_overlapping_varcoefs(output$varying_coefs_list, cls = 2, pre = 1, yname = 'v21')
v22 <- draw_overlapping_varcoefs(output$varying_coefs_list, cls = 2, pre = 2, yname = 'v22')
v31 <- draw_overlapping_varcoefs(output$varying_coefs_list, cls = 3, pre = 1, yname = 'v31')
v32 <- draw_overlapping_varcoefs(output$varying_coefs_list, cls = 3, pre = 2, yname = 'v32')

title_case1 <- ggdraw() + draw_label("         ", fontface='bold')
title_case0 <- ggdraw() + draw_label("", fontface='bold')
v1n2 <- plot_grid(v11, v12, nrow=2, align='v')
v1n2 <- plot_grid(title_case1, v1n2, ncol = 1, rel_heights=c(0.1, 1), align='v')
v2n2 <- plot_grid(v21, v22, nrow=2, align='v')
v2n2 <- plot_grid(title_case0, v2n2, ncol = 1, rel_heights=c(0.1, 1), align='v')
v3n2 <- plot_grid(v31, v32, nrow=2, align='v')
v3n2 <- plot_grid(title_case0, v3n2, ncol = 1, rel_heights=c(0.1, 1), align='v')

plot_grid(v1n2, v2n2, v3n2, nrow = 1, align = 'h')

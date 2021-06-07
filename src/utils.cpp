#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export()]]
List byproducts(arma::vec num_obs, int K, int NSbj, List C_star_sbj, List gamma, List X_sbj, arma::vec beta, List Z_sbj, List L_sbj, arma::mat Psi, arma::mat Psi0, arma::vec ID, arma::vec Clusters, arma::vec active, Function flatten_gamma){
  
  List IU(NSbj);
  List IU0(NSbj);
  List XXi_sbj(NSbj);
  List Xi_sbj(NSbj);
  
  for(int i=0 ; i<NSbj ; i++){
    arma::mat temp_gamma = as<arma::mat>(gamma[Clusters[i]-1]);
    NumericVector indZ = flatten_gamma(temp_gamma);
    arma::mat C_star_sbj_temp = as<arma::mat>(C_star_sbj[i]);
    arma::mat C_star_sbj_temp_sel = C_star_sbj_temp.cols(arma::find(as<arma::fvec>(indZ)>0.5));
    
    arma::mat Inv_middle_temp = inv_sympd(arma::eye(num_obs[i],num_obs[i])+as<arma::mat>(Z_sbj[i])*Psi*trans(as<arma::mat>(Z_sbj[i])));
    arma::mat Inv_middle_temp0 = inv_sympd(arma::eye(num_obs[i],num_obs[i])+as<arma::mat>(Z_sbj[i])*Psi0*trans(as<arma::mat>(Z_sbj[i])));
    
    IU[i] = Inv_middle_temp;
    IU0[i] = Inv_middle_temp0;
    
    arma::mat common_temp = trans(C_star_sbj_temp_sel) * Inv_middle_temp;
    arma::mat temp1 = common_temp * C_star_sbj_temp_sel;
    arma::mat temp2 = common_temp * (as<arma::vec>(L_sbj[i])-as<arma::mat>(X_sbj[i])*beta);
    
    XXi_sbj[i] = temp1;
    Xi_sbj[i] = temp2;
  }
  
  List XXi_sum(K);
  List Xi_sum(K);
  
  for(int k=0 ; k<K ; k++){
    if(std::find(active.begin(),active.end(),(k+1)) != active.end()){
      arma::uvec pos_k = find((Clusters-1)==k);
      
      int mat_size = as<arma::mat>(XXi_sbj[pos_k[0]]).n_cols;
      int pos_size = pos_k.size();
      
      arma::mat XXi_sum_temp = arma::zeros(mat_size, mat_size);
      arma::mat Xi_sum_temp = arma::zeros(mat_size, 1);
      for(int l=0 ; l < pos_size ; l++){
        XXi_sum_temp += as<arma::mat>(XXi_sbj[pos_k[l]]);
        Xi_sum_temp += as<arma::vec>(Xi_sbj[pos_k[l]]);
      }
      XXi_sum[k] = (trans(XXi_sum_temp)+XXi_sum_temp)/2;
      Xi_sum[k] = Xi_sum_temp;
    }else{
      XXi_sum[k] = 0;
      Xi_sum[k] = 0;
    }
    
  }
  
  return( List::create(IU, IU0, XXi_sum, Xi_sum) );
  
}

// [[Rcpp::export()]] 
List KnotSelection(int small_k, int NSbj, double tau, List C_star_sbj, NumericVector prsd_gamma, List Inv_middle, List X_sbj, arma::vec beta, List L_sbj, arma::vec ID, arma::vec Clusters){

  NumericVector indZ = prsd_gamma;
  List XXi_sbj(NSbj);
  List Xi_sbj(NSbj);
  
  arma::mat sumXUX_total = arma::zeros(sum(indZ),sum(indZ));
  arma::mat sumXUY_total = arma::zeros(sum(indZ),1);
  
  for(int i=0 ; i<NSbj ; i++){
    arma::mat C_star_sbj_temp = as<arma::mat>(C_star_sbj[i]);
    arma::mat C_star_sbj_temp_sel = C_star_sbj_temp.cols(arma::find(as<arma::fvec>(indZ)>0.5));
    arma::mat common_temp = trans(C_star_sbj_temp_sel) * as<arma::mat>(Inv_middle[i]);
    arma::mat temp1 = common_temp * C_star_sbj_temp_sel;
    arma::mat temp2 = common_temp * (as<arma::vec>(L_sbj[i])-as<arma::mat>(X_sbj[i])*beta);
    sumXUX_total += temp1;
    sumXUY_total += temp2;
    XXi_sbj[i] = temp1;
    Xi_sbj[i] = temp2;
  }
  
  arma::uvec pos_k = find(Clusters==small_k);
  int mat_size = as<arma::mat>(XXi_sbj[pos_k[0]]).n_cols;
  int pos_size = pos_k.size();
  
  arma::mat XXi_sum_temp = arma::zeros(mat_size, mat_size);
  arma::mat Xi_sum_temp = arma::zeros(mat_size, 1);
  
  for(int l=0 ; l<pos_size ; l++){
    XXi_sum_temp += as<arma::mat>(XXi_sbj[pos_k[l]]);
    Xi_sum_temp += as<arma::vec>(Xi_sbj[pos_k[l]]);
  }
  
  double val;
  double sign;
  arma::mat out_val;
  
  arma::mat temp_in = tau*inv_sympd((trans(sumXUX_total)+sumXUX_total)/2)*XXi_sum_temp + arma::eye(sum(indZ),sum(indZ));
  log_det(val, sign, temp_in);
  arma::mat temp_in2 = XXi_sum_temp + sumXUX_total/tau;
  arma::mat first_part = 0.5*trans(Xi_sum_temp)*inv_sympd((trans(temp_in2)+temp_in2)/2)*Xi_sum_temp;
  double second_part = (-0.5)*val;
  out_val = first_part + second_part;       

return( List::create(out_val) );

}

// [[Rcpp::export()]]
List MakeRk(int K, arma::vec tau, arma::vec active, List C_star_sbj, List inv_middle, List inv_middle0, List gamma, int NSbj, Function flatten_gamma){
  List Rk(K);
  List Inv_Phi_k(K);
  List Inv_XXi_k(K);
  
  for(int k=0 ; k<K ; k++){
    
    NumericVector indZ = flatten_gamma(gamma[k]);
    arma::mat sumXUX_total = arma::zeros(sum(indZ),sum(indZ));
    
    if(std::find(active.begin(),active.end(),(k+1)) != active.end()){
      
      for(int i=0 ; i<NSbj ; i++){
        arma::mat C_star_sbj_temp = as<arma::mat>(C_star_sbj[i]);
        arma::mat C_star_sbj_temp_sel = C_star_sbj_temp.cols(arma::find(as<arma::fvec>(indZ)>0.5));
        sumXUX_total += trans(C_star_sbj_temp_sel)*as<arma::mat>(inv_middle[i])*C_star_sbj_temp_sel;
        
      }
      
    }else{
      for(int i=0 ; i<NSbj ; i++){
        arma::mat C_star_sbj_temp = as<arma::mat>(C_star_sbj[i]);
        arma::mat C_star_sbj_temp_sel = C_star_sbj_temp.cols(arma::find(as<arma::fvec>(indZ)>0.5));
        sumXUX_total += trans(C_star_sbj_temp_sel)*as<arma::mat>(inv_middle0[i])*C_star_sbj_temp_sel;
      }
    }
    Rk[k] = sumXUX_total;
  }
  return(List::create(Rk));
}

// [[Rcpp::export()]]
List sample_randomeffects(int n, List Zi, List InvAdj, List Li, List C_star_sbj, List gamma, List X_sbj, arma::vec beta, arma::mat Psi, List theta_k, arma::vec Clusters, int r, Function flatten_gamma){
  int k = Li.size();
  List listb(k);
  arma::mat btb = arma::zeros(r,r);
  for(int i = 0; i < k; i++) {
    NumericVector indZ = flatten_gamma(gamma[Clusters[i]-1]);
    arma::mat C_star_sbj_temp = as<arma::mat>(C_star_sbj[i]);
    arma::mat C_star_sbj_temp_sel = C_star_sbj_temp.cols(arma::find(as<arma::fvec>(indZ)>0.5));
    arma::vec fixE = as<arma::vec>(theta_k[Clusters[i]-1]);
    
    arma::mat tempSigma = Psi-Psi*trans(as<arma::mat>(Zi[i]))*as<arma::mat>(InvAdj[i])*as<arma::mat>(Zi[i])*Psi;
    arma::mat Sigma = (trans(tempSigma)+tempSigma)/2;
    arma::vec mu = Psi*trans(as<arma::mat>(Zi[i]))*as<arma::mat>(InvAdj[i])*(as<arma::vec>(Li[i]) - C_star_sbj_temp_sel*fixE - as<arma::mat>(X_sbj[i])*beta );
    int p = Sigma.n_cols;
    
    arma::mat X = reshape(arma::vec(rnorm(p * n)), p, n);
    
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, Sigma);
    X = eigvec * diagmat(sqrt(eigval)) * X;
    X.each_col() += mu;
    listb[i] = trans(X);
    btb += X*trans(X);
    
  }
  return List::create(listb,btb);
}

// [[Rcpp::export()]]
List MakeRkPsi(List n_obs, int K, List Rk, arma::vec tau, arma::vec active, List C_star_sbj, List Z_sbj, arma::mat p_Psi, List gamma, List theta_k, int NSbj, Function flatten_gamma){
  List p_Rk(K);
  List log_psi_value(K);
  
  for(int k=0 ; k<K ; k++){
    
    NumericVector indZ = flatten_gamma(gamma[k]);
    arma::mat sumXUX_total = arma::zeros(sum(indZ),sum(indZ));
    arma::vec fixE = as<arma::vec>(theta_k[k]);
    if(std::find(active.begin(),active.end(),(k+1)) != active.end()){
      
      for(int i=0 ; i<NSbj ; i++){
        arma::mat temp_Mat = arma::eye(n_obs[i],n_obs[i])+as<arma::mat>(Z_sbj[i])*p_Psi*trans(as<arma::mat>(Z_sbj[i]));
        arma::mat IU = inv_sympd((trans(temp_Mat)+temp_Mat)/2);
        arma::mat C_star_sbj_temp = as<arma::mat>(C_star_sbj[i]);
        arma::mat C_star_sbj_temp_sel = C_star_sbj_temp.cols(arma::find(as<arma::fvec>(indZ)>0.5));
        sumXUX_total += trans(C_star_sbj_temp_sel)*IU*C_star_sbj_temp_sel;
      }
      
      arma::mat p_cov_mat = tau[k] * inv_sympd((trans(sumXUX_total)+sumXUX_total)/2);
      arma::mat c_cov_mat = tau[k] * inv_sympd((trans(as<arma::mat>(Rk[k]))+as<arma::mat>(Rk[k]))/2);
      
      double p_val;
      double p_sign;
      double c_val;
      double c_sign;
      
      log_det(p_val, p_sign, p_cov_mat);
      log_det(c_val, c_sign, c_cov_mat);
      
      log_psi_value[k] =  exp( (0.5*p_val-0.5*trans(fixE)*sumXUX_total*fixE) - (0.5*c_val-0.5*trans(fixE)*as<arma::mat>(Rk[k])*fixE) );
      
      
    }else{
      sumXUX_total += 0;
      log_psi_value[k] = 0;
    }
    
  }
  return(List::create(log_psi_value));
}

// [[Rcpp::export()]]
List ObtainOutput(int NSbj, int K, arma::vec Diri_p, List L_sbj, List C_star_sbj, List X_sbj, arma::vec beta, List Z_sbj, List bi, List gamma, List theta_k, Function flatten_gamma){
  
  List EachCategory(K);
  
  for(int k=0 ; k<K ; k++){
    
    NumericVector indZ = flatten_gamma(gamma[k]);
    arma::vec fixE = as<arma::vec>(theta_k[k]);
    double log_diri_p = log(Diri_p[k]);
    List Value_Sbj(NSbj); 
    for(int i=0 ; i<NSbj ; i++){
      arma::mat C_star_sbj_temp = as<arma::mat>(C_star_sbj[i]);
      arma::mat C_star_sbj_temp_sel = C_star_sbj_temp.cols(arma::find(as<arma::fvec>(indZ)>0.5));
      arma::vec temp_Li_mu = as<arma::vec>(L_sbj[i]) - C_star_sbj_temp_sel*fixE - as<arma::mat>(Z_sbj[i])*as<arma::vec>(bi[i]) - as<arma::mat>(X_sbj[i])*beta;
      Value_Sbj[i] = (-0.5)*trans(temp_Li_mu)*temp_Li_mu + log_diri_p;
    }
    EachCategory[k]=Value_Sbj;
  }
  
  return(List::create(EachCategory));
}


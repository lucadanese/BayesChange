#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_sf.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace RcppGSL;

//--------
// UTILS
//--------

int which_min_cpp(arma::vec x) {
  int n = x.size();
  int min_index = 0; // 0-based index
  double min_value = x[0];

  for (int i = 1; i < n; ++i) {
    if (x[i] < min_value) {
      min_value = x[i];
      min_index = i;
    }
  }

  return min_index; // Convert to 1-based index
}

double rshiftedgamma(double a, double b, double shift_param, gsl_rng *r){

  double res = gsl_ran_gamma(r,a,b) + shift_param;

  return res;
}

double compute_k_n(double k_0, double gamma, int n_gamma){
  double k_n = k_0 + (std::pow((1-gamma),2)/(1 - std::pow(gamma,2)))*(n_gamma-1) + 1 ;
  return k_n;
}

arma::mat compute_s_n(int d, arma::vec m_0, float k_0, int V_0, arma::mat phi_0, float gamma, int n_gamma, arma::mat gamma_k, int nbasis){

  float k_n;

  arma::mat vettore_m_n(nbasis, n_gamma - 1, arma::fill::zeros); // inizializziamo la matrice m_n
  arma::mat vettore_phi_n(nbasis, nbasis, arma::fill::zeros);     // inizializziamo la matrice Phi_n
  arma::vec m_n(d);
  arma::mat phi_n(d,d, arma::fill::zeros);
  k_n = k_0 + (std::pow((1-gamma),2)/(1 - std::pow(gamma,2)))*(n_gamma-1) + 1 ;

  // m_n
  if(n_gamma != 1){
    for(int i = 1; i < n_gamma; i++){
      vettore_m_n.col(i-1) = gamma_k.row(i).t() - (gamma * gamma_k.row(i-1).t());
    }
    m_n = (1/k_n)*( (m_0*k_0) + gamma_k.row(0).t() + ((1-gamma)/(1 - std::pow(gamma,2))) * sum(vettore_m_n,1));
  } else {
    m_n = (1/k_n)*((m_0*k_0) + gamma_k.row(0).t());
  }

  // Phi_n
  if(n_gamma != 1){
    for(int i = 1; i < n_gamma; i++){
      vettore_phi_n = vettore_phi_n + ((gamma_k.row(i).t() - gamma * gamma_k.row(i-1).t()) * (gamma_k.row(i).t() - gamma * gamma_k.row(i-1).t()).t())/(1-std::pow(gamma,2));
    }
    phi_n = phi_0 + (gamma_k.row(0).t() * gamma_k.row(0)) + vettore_phi_n + ((m_0*m_0.t())*k_0) - ((m_n * m_n.t())*k_n) ;
  }else{
    phi_n = phi_0 + (gamma_k.row(0).t() * gamma_k.row(0)) + ((m_0 * m_0.t())*k_0) - ((m_n * m_n.t())*k_n) ;
  }
  return phi_n;
}

float lgamma_multi(int d, int a){

  arma::vec vettore(d, arma::fill::zeros);

  for(int j = 0; j < d ; j++){
    float j_float = j;
    vettore(j) = std::lgamma(a + (-j_float/2));
  }
  return log(std::pow(M_PI,(d*(d-1)/4))) + arma::accu(vettore);
}

double indicator_1(int k, int n){
  double res = 0.0;
  if(1 < k && k < n){
    res = 1.0;
  }
  return res;
}

double indicator_2(int k, int n){
  double res = 0.0;
  if(k == n){
    res = 1.0;
  }
  return res;
}

double indicator_3(int k){
  double res = 0.0;
  if(k==1){
    res = 1.0;
  }
  return res;
}

int rint(arma::vec freqs){
  arma::vec probs = cumsum(freqs / sum(freqs));
  double u = arma::randu();
  for(arma::uword i = 0; i < freqs.n_elem; i++){
    if(u <= probs(i)) {
      return i;
    }
  }
  return -1;
}

arma::vec table_cpp(arma::vec vector){
  double k = max(vector) + 1.0;
  arma::vec table(k);

  table.fill(0);
  for(int i = 0; i < k; i++){
    for(arma::uword j = 0; j < vector.n_elem; j++){
      if(vector(j) == i){
        table(i) += 1;
      }
    }
  }
  return(table);
}

int rbinom(int n, double p){
  int out = 0;
  for(int i = 0; i < n; i++){
    if(arma::randu() < p){
      out += 1;
    }
  }
  return out;
}

arma::vec rmultin(int n, arma::vec probs){
  arma::vec out(probs);
  out.fill(0);
  int temp;
  for(int i = 0; i < n; i++){
    temp = rint(probs);
    out(temp) += 1;
  }
  return out;
}

double dmultinom_log_cpp(arma::vec x, arma::vec prob){

  arma::vec temp(x.n_elem);

  double s = sum(prob);

  double size = sum(x);

  for(arma::uword i = 0; i < x.n_elem; i++){

    prob(i) = prob(i)/s;

    temp(i) = x(i) * log(prob(i)) - lgamma(x(i) + 1);

  }

  double res = lgamma(size + 1) + sum(temp);

  return res;
}

double dbinom_log_cpp(double x,
                      double size,
                      double prob){

  double res = gsl_sf_lnchoose(size,x) + (x * log(prob)) + (size - x)*log(1-prob) ;

  return res;
}

double my_min(double a,
              double b){
  if(a <= b){
    return a;
  } else {
    return b;
  }
}

double my_max(double a,
              double b){
  if(a >= b){
    return a;
  } else {
    return b;
  }
}

double log_sum_exp(arma::vec log_vals){
  double M = max(log_vals);
  double lse = M + log(sum(exp(log_vals - M)));
  return lse;
}

arma::vec generate_random_partition(int n, gsl_rng *r ){
  arma::vec random_order(n);

  int num_groups = gsl_rng_uniform_int(r, n) + 1;
  int counter = 0;

  //for(int i = 0; i < (num_groups-1); i++){
  for(int i = 0; i < (num_groups+1); i++){
    random_order(i) = counter;
    counter = counter + 1;
  }

  if((num_groups - n) != 0){
    for(int i = (num_groups+1); i < n; i++){
      random_order(i) = randi(1,arma::distr_param(0,num_groups-1))(0);
    }
  }

  return random_order;

}

arma::vec clean_partition_cpp(arma::vec partition){
  arma::vec partition_temp = partition;
  double index, counter = 0;
  arma::uvec uni_clust = find_unique(partition, true);
  for(arma::uword i = 0; i < uni_clust.n_elem; i++){
    index = partition(uni_clust(i));
    arma::uvec sub = find(partition == index);

    for(arma::uword j = 0; j < sub.n_elem; j++){
      partition_temp(sub(j)) = counter;
    }
    counter += 1;
  }
  return partition_temp;
}

arma::vec generate_random_order(double t, double p, gsl_rng *r){
  arma::vec freqs, temp_probs, new_order, cfreq;
  int k;

  k = 1 + rbinom(t, p);
  temp_probs.resize(k);
  for(int l = 0; l < k; l++){
    temp_probs(l) = gsl_ran_gamma(r, 4, 2);
  }
  freqs = rmultin(t, temp_probs);
  new_order.resize(t);
  cfreq = cumsum(freqs);
  for(int i = 0; i < cfreq(0); i++){
    new_order(i) = 0;
  }
  for(arma::uword j = 1; j < cfreq.n_elem; j++){
    for(int i = cfreq(j-1); i < cfreq(j); i++){
      new_order(i) = j;
    }
  }
  while(new_order(0) > 0){
    new_order -= 1;
  }
  return new_order;
}

double AbsStirling1st(double r, double k){
  double n = my_max(r,k) + 1, res;
  if((k <= 0) | (k > r)){
    res = 0;
  } else {
    arma::mat matrix(n,n,arma::fill::zeros);
    matrix(0,0) = 1;
    matrix(1,1) = 1;
    for(int i = 2; i < n; i++){
      for(int j = 1; j <= i; j++){
        matrix(i,j) = - (i-1) * matrix(i-1,j) + matrix(i-1,j-1);
        matrix(j,i) = matrix(i,j);
      }
      matrix(i,i) = 1;
    }
    res = abs(matrix(r,k));
  }
  return(res);
}

arma::mat ExtractSubData(arma::mat data, arma::vec order, int index){

  arma::vec table_order = table_cpp(order);
  arma::mat res_mat(data.n_rows,table_order(index));
  arma::vec cumsum_vec = cumsum(table_order);

  if(index == 0){
    res_mat = data.cols(0, cumsum_vec(index) - 1);
  } else {
    res_mat = data.cols(cumsum_vec(index - 1),cumsum_vec(index) - 1);
  }

  return res_mat;

}

//---------------------------------
// LIKELIHOOD PRIOR AND POSTERIOR
//---------------------------------

double Likelihood_UniTS(arma::mat data, arma::vec order,
                        double phi, double a, double b, double c){

  double k = max(order) + 1;
  arma::vec res_vec, table_order = table_cpp(order), vec_likelihood(k);

  for(int i = 0; i < k; i++){

    double n_i = table_order(i);

    arma::mat gamma_k(1,n_i);

    if(i == 0){
      gamma_k = data.row(0).cols(0, n_i-1);
    } else if (i == k) {
      gamma_k = data.row(0).cols(data.n_cols - 1 - n_i, data.n_cols - 1);
    } else {
      arma::vec table_order_temp = table_order.subvec(0,i-1);
      gamma_k = data.row(0).cols(sum(table_order_temp), sum(table_order_temp) + n_i - 1);

    }

    if((n_i != 1) & (n_i != 2)){

      arma::mat S_i(n_i,n_i,arma::fill::zeros);
      S_i(0,0) = 1;
      for(int j = 1; j < n_i; j++){
        S_i(j,j) = 1 + std::pow(phi,2);
        S_i(j,j-1) = -phi;
        S_i(j-1,j) = -phi;
      }

      S_i(n_i-1,n_i-1) = 1;

      vec_likelihood(i) = a*std::log(2*b* ( 1 -std::pow(phi,2))) + gsl_sf_lngamma(n_i/2 + a) - n_i/2 * std::log(M_PI) - gsl_sf_lngamma(a) + 0.5*(std::log(c) + std::log(1+phi) + std::log(1+std::pow(phi,2)) - std::log(c) - std::log(n_i) + std::log(phi) + std::log(n_i-c-2)) - (n_i/2 + a) * log((gamma_k * S_i * gamma_k.t()).eval()(0,0) - (((1-phi)*std::pow(sum(gamma_k.row(0)) - phi * sum(gamma_k.row(0).cols(1,gamma_k.n_cols - 2)),2))/(c+n_i-phi*(n_i-c-2))) + 2*b*(1-std::pow(phi,2)));
    } else if (n_i == 1) {

      arma::mat S_i(1,1,arma::fill::zeros);
      S_i(0,0) = 1;

      vec_likelihood(i) = a*std::log(2*b* ( 1 -std::pow(phi,2))) + gsl_sf_lngamma(n_i/2 + a) - n_i/2 * std::log(M_PI) - gsl_sf_lngamma(a) + 0.5*(std::log(c) + std::log(1+phi) + std::log(1+std::pow(phi,2)) - std::log(c) - std::log(n_i) + std::log(phi) + std::log(n_i-c-2)) - (n_i/2 + a) * log((gamma_k * S_i * gamma_k.t()).eval()(0,0) - (((1-phi)*std::pow(sum(gamma_k.row(0)),2))/(c+n_i-phi*(n_i-c-2))) + 2*b*(1-std::pow(phi,2)));
    } else if (n_i == 2) {

      arma::mat S_i(n_i,n_i,arma::fill::zeros);
      S_i(0,0) = 1;
      for(int j = 1; j < n_i; j++){
        S_i(j,j) = 1 + std::pow(phi,2);
        S_i(j,j-1) = -phi;
        S_i(j-1,j) = -phi;
      }

      S_i(n_i-1,n_i-1) = 1;

      vec_likelihood(i) = a*std::log(2*b* ( 1 -std::pow(phi,2))) + gsl_sf_lngamma(n_i/2 + a) - n_i/2 * std::log(M_PI) - gsl_sf_lngamma(a) + 0.5*(std::log(c) + std::log(1+phi) + std::log(1+std::pow(phi,2)) - std::log(c) - std::log(n_i) + std::log(phi) + std::log(n_i-c-2)) - (n_i/2 + a) * log((gamma_k * S_i * gamma_k.t()).eval()(0,0) - (((1-phi)*std::pow(sum(gamma_k.row(0)) - phi * sum(gamma_k.row(0).col(1)),2))/(c+n_i-phi*(n_i-c-2))) + 2*b*(1-std::pow(phi,2)));

    }

  }

  return sum(vec_likelihood);

}


double LogLikelihood_TS(arma::mat data, arma::mat order,
                        double gamma_par, double a, double b, double c){

  arma::vec lkl(max(order.row(0))+1);
  double length_block;

  for(int j = 0; j < max(order.row(0)) + 1; j++){

    length_block = max(find(order.row(0) == j))+1 - min(find(order.row(0) == j));
    arma::mat data_block = data.cols(min(find(order.row(0) == j)),max(find(order.row(0) == j)));
    arma::mat Sj(length_block, length_block);

    if(length_block != 1.0){
      Sj.fill(0);
      Sj.diag().fill(1 + std::pow(gamma_par,2));
      Sj.row(0).col(0) = 1;
      Sj.row(Sj.n_rows - 1).col(Sj.n_cols - 1) = 1;
      for(arma::uword h = 0; h < Sj.n_rows - 2; h ++){
        Sj.row(h).col(h+1) = - gamma_par;
      }
      for(arma::uword h = 1; h < Sj.n_rows - 1; h ++){
        Sj.row(h).col(h-1) = - gamma_par;
      }
    } else {
      Sj.fill(1);
    }

    double r1 = a * log(2*b * (1-pow(gamma_par,2))) + lgamma(length_block/2 + a) - (log(M_PI) * (length_block/2) + lgamma(a));
    double r2 = 0.5 * (log(c*(1-gamma_par) * (1-pow(gamma_par,2))) - log(c + length_block - gamma_par * (length_block-c-2)) );

    arma::mat r3;

    if(length_block > 2.0){
      r3 = -(length_block/2 + a) * log( (data_block.row(0) * Sj * data_block.row(0).t()) - (( (1-gamma_par) * pow(sum(data_block.row(0)) - gamma_par * sum(data_block.row(0).cols(1,length_block-2)),2) / (c + length_block - gamma_par*(length_block-c-2)) )) + 2*b*(1-pow(gamma_par,2)));
    } else if (length_block == 2.0) {
      r3 = -(length_block/2 + a) * log( (data_block.row(0) * Sj * data_block.row(0).t()) - (( (1-gamma_par) * pow(sum(data_block.row(0)) - gamma_par * sum(data_block.row(0).col(1)),2) / (c + length_block - gamma_par*(length_block-c-2)) )) + 2*b*(1-pow(gamma_par,2)));
    } else {
      r3 = -(length_block/2 + a) * log( (data_block.row(0) * Sj * data_block.row(0).t()) - (( (1-0) * pow(sum(data_block.row(0)),2)) / (c+length_block - 0*(length_block-c-2))) + 2*b*(1-0));
    }
    lkl(j) = r1+r2+as_scalar(r3);
  }
  return(sum(lkl));
}


double Likelihood_MultiTS(arma::mat data, arma::vec order,
                          double gamma, double k_0, double nu_0,
                          arma::mat phi_0, arma::vec m_0){
  double d = data.n_rows, k = max(order) + 1;
  arma::vec res_vec, table_order = table_cpp(order), vec_likelihood(k);

  for(int j = 0; j < k; j++){

    double n_i = table_order(j);

    arma::mat vettore_m_n(d, n_i - 1, arma::fill::zeros);
    arma::mat vettore_phi_n(d, d, arma::fill::zeros);
    arma::vec m_n(d);
    arma::mat phi_n(d,d, arma::fill::zeros);
    arma::mat gamma_k(d,n_i);

    if(j == 0){
      gamma_k = data.rows(0,d-1).cols(0, n_i-1);
    } else if (j == k) {
      gamma_k = data.rows(0,d-1).cols(data.n_cols - 1 - n_i, data.n_cols - 1);
    } else {
      arma::vec table_order_temp = table_order.subvec(0,j-1);
      gamma_k = data.rows(0,d-1).cols(sum(table_order_temp), sum(table_order_temp) + n_i - 1);

    }

    double nu_n = nu_0 + (n_i);
    double k_n = k_0 + (std::pow((1-gamma),2)/(1 - std::pow(gamma,2)))*(n_i-1) + 1;

    // m_n

    if(n_i != 1){

      for(int j = 1; j < n_i; j++){
        vettore_m_n.col(j-1) = gamma_k.col(j) - (gamma * gamma_k.col(j-1));
      }

      m_n = (1/k_n)*( (m_0*k_0) + gamma_k.col(0) + ((1-gamma)/(1 - std::pow(gamma,2))) * sum(vettore_m_n,1));

    } else {

      m_n = (1/k_n)*((m_0*k_0) + gamma_k.col(0));

    }

    // Phi_n

    if(n_i != 1){

      for(int j = 1; j < n_i; j++){
        vettore_phi_n = vettore_phi_n + ((gamma_k.col(j) - gamma * gamma_k.col(j-1)) * (gamma_k.col(j) - gamma * gamma_k.col(j-1)).t())/(1-std::pow(gamma,2));
      }

      phi_n = phi_0 + (gamma_k.col(0) * gamma_k.col(0).t()) + vettore_phi_n + ((m_0*m_0.t())*k_0) - ((m_n * m_n.t())*k_n) ;

    } else {

      phi_n = phi_0 + (gamma_k.col(0) * gamma_k.col(0).t()) + ((m_0 * m_0.t())*k_0) - ((m_n * m_n.t())*k_n) ;

    }

    vec_likelihood(j) = log(pow(k_n,(-d/2))) - log(std::pow(k_0,(-d/2))) + (nu_0/2) * log(arma::det(phi_0)) + lgamma_multi(d, (nu_n/2)) - (nu_n/2)*log(arma::det(phi_n)) - lgamma_multi(d,(nu_0/2)) + log(std::pow(M_PI,(-n_i*d/2))) - log(pow((1-pow(gamma,2)),((n_i-1)/2)));

  }

  return sum(vec_likelihood);

}

double Prior_TS(arma::vec order,
                double theta, double sigma){

  arma::vec table_order = table_cpp(order), vec_a3;

  int k = max(order) + 1;
  int n = order.n_elem;

  double a1 = gsl_sf_lnfact(n) - gsl_sf_lnfact(k);

  double a2 = 0;

  for(int i = 1; i < k; i++){
    a2 += log(theta + i*sigma);
  }

  a2 = a2 - gsl_sf_lnpoch(theta + 1, n - 1);

  vec_a3.resize(k);

  for(int i = 0; i < k; i++){
    vec_a3(i) = gsl_sf_lnpoch(1 - sigma, table_order(i) - 1) - gsl_sf_lnfact(table_order(i));
  }

  double a3 = sum(vec_a3);

  double res = a1 + a2 + a3;

  return res;
}

double Posterior_MultiTS(arma::mat data, arma::vec order,
                         double gamma, double k_0, double nu_0,
                         double theta, double sigma,
                         arma::mat phi_0, arma::vec m_0){
  return Prior_TS(order, theta, sigma) + Likelihood_MultiTS(data, order, gamma, k_0, nu_0, phi_0, m_0);
}

double Posterior_UniTS(arma::mat data, arma::vec order,
                       double theta, double sigma, double phi,
                       double a, double b, double c){
  return Prior_TS(order, theta, sigma) + Likelihood_UniTS(data, order, phi, a, b, c);
}



arma::vec DSA_curve(double dt,
                    int T,
                    double rhoval,
                    double gamma,
                    arma::vec Betat,
                    double S0 = 1,
                    double R0 = 0){

  arma::vec out(T);
  double t = 0;
  int j = 0;

  double St = S0 + dt * (-Betat(j) * rhoval);
  double It = rhoval + dt * (Betat(j) * rhoval - gamma * rhoval);
  double Rt = R0 + dt * (gamma * rhoval);

  for(int i = 0; i < T/dt; i++){

    t = t + dt;
    St = St + dt * (-Betat(j) * St * It);
    It = It + dt * (Betat(j) * St * It - gamma*It);
    Rt = Rt + dt * (gamma * It);

    if(floor(t) == j + 1){
      out(j) = log(Betat(j)) + log(St) + log(It) ;
      j = j + 1;
    }
  }
  return(out - log(1 - St));
}

Rcpp::List SIR_curve(double dt,
                     int T,
                     arma::vec order,
                     double rhoval,
                     double gamma,
                     double a0,
                     double b0,
                     double S0 = 1,
                     double R0 = 0){

  arma::vec out(T);


  arma::vec t_beta(T);

  t_beta(0) = arma::randg(1, arma::distr_param(a0, 1.0 / b0) )[0];
  for(int i = 1; i < T; i++){
    if(order(i) == order(i - 1)){
      t_beta(i) = t_beta(i - 1);
    } else {
      t_beta(i) = arma::randg(1, arma::distr_param(a0, 1.0 / b0) )[0];
    }
  }

  double t = 0;
  int j = 0;

  arma::vec Betat = t_beta;

  double St = S0 + dt * (-Betat(j) * rhoval);
  double It = rhoval + dt * (Betat(j) * rhoval - gamma * rhoval);
  double Rt = R0 + dt * (gamma * rhoval);

  for(int i = 0; i < T/dt; i++){

    t = t + dt;
    St = St + dt * (-Betat(j) * St * It);
    It = It + dt * (Betat(j) * St * It - gamma*It);
    Rt = Rt + dt * (gamma * It);

    if(floor(t) == j + 1){
      out(j) = St;
      j = j + 1;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("st") = out,
    Rcpp::Named("beta") = Betat
  );

}

arma::mat integrated_curves_mat(double dt,
                                arma::vec order,
                                double a0,
                                double b0,
                                double gamma,
                                double rhoval,
                                int M,
                                double S0 = 1,
                                double R0 = 0){
  int T = order.n_elem;
  arma::mat out(M, T + 1);
  arma::vec t_beta(T), t_Imp(T);

  for(int m = 0; m < M; m++){

    t_beta(0) = arma::randg(1, arma::distr_param(a0, 1.0 / b0))[0];
    t_Imp(0) = R::dgamma(t_beta(0), a0, 1/b0, true);

    for(int i = 1; i < T; i++){
      if(order(i) == order(i - 1)){
        t_beta(i) = t_beta(i - 1);
      } else {
        t_beta(i) = arma::randg(1, arma::distr_param(a0, 1.0 / b0) )[0];
        t_Imp.resize(t_Imp.n_elem + 1);
        t_Imp(t_Imp.n_elem - 1) = R::dgamma(t_beta(i), a0, 1/ b0, true);
      }
    }

    out.row(m).cols(0,T-1) = DSA_curve(dt, T, rhoval, gamma, t_beta, S0, R0).t();
    out.row(m).col(T) = sum(t_Imp);

    t_Imp.resize(1);

  }
  return out;
}

// ------------------
// ACCEPTANCE RATIO
// ------------------

double AlphaSplitOrder_TS(arma::mat data, arma::vec new_order, arma::vec old_order,
                          int index, arma::vec freq_temp,
                          double a, double b, double c,
                          double gamma, double q){

  arma::vec lkl_vec(data.n_rows);

  for(arma::uword i = 0; i < data.n_rows; i++){
    lkl_vec(i) = LogLikelihood_TS(data.row(i), new_order.t(), gamma, a, b, c) - LogLikelihood_TS(data.row(i), old_order.t(), gamma, a, b, c);
  }

  double a11 = old_order.n_elem - std::count(freq_temp.begin(),freq_temp.end(),1);
  double a12 = std::count(old_order.begin(),old_order.end(),old_order(index));
  double a13 = freq_temp(index);

  double a1 = log((1-q)/q) + sum(lkl_vec) + log((a11 * (a12 - 1))/a13);
  double a2 = log(1.0);

  double res = my_min(a1,a2);

  return res;
}

double AlphaMergeOrder_TS(arma::mat data, arma::vec new_order, arma::vec old_order,
                          int index, arma::vec freq_temp,
                          double a, double b, double c,
                          double gamma, double q){

  arma::vec lkl_vec(data.n_rows);

  for(arma::uword i = 0; i < data.n_rows; i++){
    lkl_vec(i) = LogLikelihood_TS(data.row(i), new_order.t(), gamma, a, b, c) - LogLikelihood_TS(data.row(i), old_order.t(), gamma, a, b, c);
  }

  double a11 = max(old_order) + 1 - 1;
  double a12 = freq_temp.n_elem - std::count(freq_temp.begin(),freq_temp.end(),1.0);
  double a13 = freq_temp(index);
  double a14 = freq_temp(index+1);

  double a1 = log(q/(1-q)) + sum(lkl_vec) + log(a11) - log(a12 * (a13 + a14 - 1)) ;
  double a2 = log(1.0);

  double res = my_min(a1,a2);

  return res;
}

double AlphaShuffleOrder_TS(arma::mat data, arma::vec new_order, arma::vec old_order,
                            double a, double b, double c, double gamma, double q){

  arma::vec lkl_vec(data.n_rows);

  for(arma::uword i = 0; i < data.n_rows; i++){
    lkl_vec(i) = LogLikelihood_TS(data.row(i), new_order.t(), gamma, a, b, c) - LogLikelihood_TS(data.row(i), old_order.t(), gamma, a, b, c);
  }

  double a1 = sum(lkl_vec);
  double a2 = log(1.0);

  double res = my_min(a1,a2);

  return res;
}


double AlphaSplitOrder_MultiTS_Acc(arma::cube data, arma::vec new_order, arma::vec old_order,
                                   double q, double index, double gamma, double k_0, double nu_0,
                                   arma::mat phi_0, arma::vec m_0){

  double k = max(old_order) + 1, a11, a12 = 0, a13;
  arma::vec table_oldorder = table_cpp(old_order);

  if((k < data.n_cols) & (k > 1)){
    a11 = log((1-q)/q);

    for(arma::uword i = 0; i < data.n_slices; i++){

      a12 = a12 + Likelihood_MultiTS(data.slice(i), new_order, gamma, k_0, nu_0, phi_0, m_0) - Likelihood_MultiTS(data.slice(i), old_order, gamma, k_0, nu_0, phi_0, m_0);

    }
    //a12 = Likelihood_MultiTS(data, new_order, gamma, k_0, nu_0, phi_0, m_0) - Likelihood_MultiTS(data, old_order, gamma, k_0, nu_0, phi_0, m_0);
    a13 = log(( table_oldorder.n_elem - std::count(table_oldorder.begin(),table_oldorder.end(),1) *(table_oldorder(index) - 1))/k);
  } else {
    a11 = log(1-q);
    for(arma::uword i = 0; i < data.n_slices; i++){
      a12 = a12 + Likelihood_MultiTS(data.slice(i), new_order, gamma, k_0, nu_0, phi_0, m_0) - Likelihood_MultiTS(data.slice(i), old_order, gamma, k_0, nu_0, phi_0, m_0);
    }
    //a12 = Likelihood_MultiTS(data, new_order, gamma, k_0, nu_0, phi_0, m_0) - Likelihood_MultiTS(data, old_order, gamma, k_0, nu_0, phi_0, m_0);
    a13 = log(data.n_cols-1);
  }

  double a1 = a11 + a12 + a13;
  double a2 = log(1);
  double res = my_min(a1,a2);

  return res;
}



double AlphaMergeOrder_MultiTS_Acc(arma::cube data, arma::vec new_order, arma::vec old_order,
                                   double q, double index, double gamma, double k_0, double nu_0,
                                   arma::mat phi_0, arma::vec m_0){

  double k = max(old_order) + 1, a11, a12 = 0, a13;
  arma::vec table_oldorder = table_cpp(old_order);

  if((k < data.n_cols) & (k > 1)){
    a11 = log(q/(1-q));
    for(arma::uword i = 0; i < data.n_slices; i++){
      a12 = a12 + Likelihood_MultiTS(data.slice(i), new_order, gamma, k_0, nu_0, phi_0, m_0) - Likelihood_MultiTS(data.slice(i), old_order, gamma, k_0, nu_0, phi_0, m_0);
    }
    a13 = log((k-1) / ((table_oldorder.n_elem - std::count(table_oldorder.begin(),table_oldorder.end(),1) + 1) * (table_oldorder(index) + table_oldorder(index + 1) - 1)) );
  } else {
    a11 = log(q);
    a12 = Likelihood_MultiTS(data, new_order, gamma, k_0, nu_0, phi_0, m_0) - Likelihood_MultiTS(data, old_order, gamma, k_0, nu_0, phi_0, m_0);
    a13 = log(data.n_cols - 1);
  }

  double a1 = a11 + a12 + a13;
  double a2 = log(1);

  return  my_min(a1,a2);
}


double AlphaShuffleOrder_MultiTS(arma::cube data, arma::vec new_order, arma::vec old_order,
                                 double gamma, double k_0, double nu_0,
                                 arma::mat phi_0, arma::vec m_0){

  double a1 = 0;

  for(arma::uword i = 0; i < data.n_slices; i++){
    a1 = a1 + Likelihood_MultiTS(data.slice(i), new_order, gamma, k_0, nu_0, phi_0, m_0) - Likelihood_MultiTS(data.slice(i), old_order, gamma, k_0, nu_0, phi_0, m_0);
  }
  double a2 = log(1);

  return  my_min(a1,a2);
}



double AlphaMergePartition_cpp(arma::vec lkl_old_i,
                               arma::vec lkl_old_j,
                               arma::vec lkl_proposal,
                               double num_groups,
                               double num_times,
                               double num_obs,
                               double alpha_SM_,
                               arma::vec merge_i,
                               arma::vec merge_j,
                               arma::vec norm_const,
                               double coars = 1){

  double n_i_merge, n_j_merge, n_ij_merge, q_merge, p_merge,
  p_merge_1, p_merge_2, p_merge_2_1, p_merge_2_2, p_merge_2_3, l_merge,
  f_merge_n, f_merge_d,f_merge ;

  n_i_merge = std::count(merge_i.begin(),merge_i.end(),1.0);
  n_j_merge = std::count(merge_j.begin(),merge_j.end(),1.0);
  n_ij_merge = n_i_merge + n_j_merge;

  q_merge = (n_i_merge + n_j_merge - 2) * log(0.5);

  p_merge_1 = lgamma(alpha_SM_) + lgamma(alpha_SM_ + n_ij_merge) - lgamma(alpha_SM_ + n_i_merge) - lgamma(alpha_SM_ + n_j_merge);

  arma::vec vec_1_km1 = arma::linspace<arma::vec>(1.0, num_groups, num_groups - 1);
  arma::vec vec_0_km1 = arma::linspace<arma::vec>(0.0, num_groups - 1, num_groups - 1);

  arma::vec vec_temp(2);
  vec_temp(0) = log(1);
  vec_temp(1) = -sum(log(-vec_1_km1 + std::pow(2,num_times-1)) - log(-vec_0_km1 + std::pow(2,num_times-1)));

  p_merge_2_1 = log_sum_exp(vec_temp);

  arma::vec vec_1_kp1 = arma::linspace<arma::vec>(1.0, num_groups + 1, num_groups + 1);
  arma::vec vec_0_k = arma::linspace<arma::vec>(0.0, num_groups, num_groups + 1);

  vec_temp(0) = log(1);
  vec_temp(1) = -sum(log(-vec_1_kp1 + std::pow(2,num_times-1)) - log(-vec_0_k + std::pow(2,num_times-1)));

  p_merge_2_2 = log_sum_exp(vec_temp) + log_sum_exp(vec_temp);

  p_merge_2_3 = log_sum_exp(vec_temp) + log_sum_exp(vec_temp);

  p_merge_2 = p_merge_2_1 - p_merge_2_2 + p_merge_2_3;

  p_merge = p_merge_1 + p_merge_2;

  l_merge =  sum(lkl_proposal.elem(find(merge_i == 1.0))) + sum(lkl_proposal.elem(find(merge_j == 1.0))) - sum(lkl_old_i.elem(find(merge_i == 1.0))) - sum(lkl_old_j.elem(find(merge_j == 1.0)));

  f_merge_n = log_sum_exp(lkl_old_i - norm_const - log(num_obs)) +
    log_sum_exp(lkl_old_j - norm_const - log(num_obs));

  f_merge_d = log_sum_exp(lkl_proposal - norm_const - log(num_obs));

  f_merge = f_merge_n - f_merge_d;

  double res = q_merge + p_merge + (coars * l_merge) + (coars * f_merge);

  return(my_min(log(1),res));

}

double AlphaSplitPartition_cpp(arma::vec lkl_proposal_i,
                               arma::vec lkl_proposal_j,
                               arma::vec lkl_old,
                               double num_groups,
                               double num_times,
                               double num_obs,
                               double alpha_SM_,
                               arma::vec split_i,
                               arma::vec split_j,
                               arma::vec norm_const,
                               double coars){

  double n_i_split = std::count(split_i.begin(),split_i.end(),1.0);
  double n_j_split = std::count(split_j.begin(),split_j.end(),1.0);
  double n_ij_split = n_i_split + n_j_split;

  double q_split = -(n_i_split + n_j_split - 2) * log(0.5);

  double p_split_1 = -lgamma(alpha_SM_) + lgamma(alpha_SM_ + n_j_split) + lgamma(alpha_SM_ + n_i_split) - lgamma(alpha_SM_ + n_ij_split);

  arma::vec vec_temp(2);
  arma::vec vec_1_kp1 = arma::linspace<arma::vec>(1.0, num_groups + 1, num_groups + 1);
  arma::vec vec_0_k = arma::linspace<arma::vec>(0.0, num_groups, num_groups + 1);

  vec_temp(0) = log(1);
  vec_temp(1) = -sum(log(-vec_1_kp1 + std::pow(2,num_times-1)) - log(-vec_0_k + std::pow(2,num_times-1)));

  double p_split_2_1 = log_sum_exp(vec_temp);

  double p_split_2_2 = log_sum_exp(vec_temp);

  arma::vec vec_1_km1 = arma::linspace<arma::vec>(1.0, num_groups, num_groups - 1);
  arma::vec vec_0_km1 = arma::linspace<arma::vec>(0.0, num_groups - 1, num_groups - 1);

  vec_temp(0) = log(1);
  vec_temp(1) = -sum(log(-vec_1_km1 + std::pow(2,num_times-1)) - log(-vec_0_km1 + std::pow(2,num_times-1)));

  double p_merge_2_3 = log_sum_exp(vec_temp);

  double p_split_2 = p_split_2_1 + p_split_2_2 - p_merge_2_3;

  double p_split = p_split_1 + p_split_2;

  double l_split =  sum(lkl_proposal_i.elem(find(split_i == 1.0))) + sum(lkl_proposal_j.elem(find(split_j == 1.0))) - sum(lkl_old.elem(find(split_i == 1.0))) - sum(lkl_old.elem(find(split_j == 1.0)));

  double f_split_n_1 = log_sum_exp(lkl_old - norm_const - log(num_obs));

  double f_split_d_1 = log_sum_exp(lkl_proposal_i - norm_const - log(num_obs));

  double f_split_d_2 = log_sum_exp(lkl_proposal_j - norm_const - log(num_obs));

  double f_split = f_split_n_1 - f_split_d_1 - f_split_d_2;

  double res = q_split + p_split + (coars * l_split) + (coars * f_split);

  return(my_min(log(1),res));

}


double AlphaSplitOrder_UniTS(arma::mat data, arma::vec new_order, arma::vec old_order,
                             double q, double index, double theta, double sigma, double phi,
                             double a, double b, double c){

  double k = max(old_order) + 1, a11 = 0, a12 = 0, a13 = 0;
  arma::vec table_oldorder = table_cpp(old_order);

  if((k < data.n_cols) & (k > 1)){
    a11 = log((1-q)/q);
    a12 = Posterior_UniTS(data, new_order, theta, sigma, phi, a, b, c) - Posterior_UniTS(data, old_order, theta, sigma, phi, a, b, c);
    a13 = log(( table_oldorder.n_elem - std::count(table_oldorder.begin(),table_oldorder.end(),1) *(table_oldorder(index) - 1))/k);
  } else if (k == 1) {
    a11 = log(1-q);
    a12 = Posterior_UniTS(data, new_order, theta, sigma, phi, a, b, c) - Posterior_UniTS(data, old_order, theta, sigma, phi, a, b, c);
    a13 = log(data.n_cols-1);
  }

  double a1 = a11 + a12 + a13;
  double a2 = log(1);
  double res = my_min(a1,a2);

  return res;
}

double AlphaMergeOrder_UniTS(arma::mat data, arma::vec new_order, arma::vec old_order,
                             double q, double index, double theta, double sigma, double phi,
                             double a, double b, double c){

  double k = max(old_order) + 1, a11 = 0, a12 = 0, a13 = 0;
  arma::vec table_oldorder = table_cpp(old_order);

  if((k < data.n_cols) & (k > 1)){
    a11 = log(q/(1-q));
    a12 = Posterior_UniTS(data, new_order, theta, sigma, phi, a, b, c) - Posterior_UniTS(data, old_order, theta, sigma, phi, a, b, c);
    a13 = log((k-1) / ((table_oldorder.n_elem - std::count(table_oldorder.begin(),table_oldorder.end(),1) + 1) * (table_oldorder(index) + table_oldorder(index + 1) - 1)) );
  } else if (k == data.n_cols){
    a11 = log(q);
    a12 = Posterior_UniTS(data, new_order, theta, sigma, phi, a, b, c) - Posterior_UniTS(data, old_order, theta, sigma, phi, a, b, c);
    a13 = log(data.n_cols - 1);
  }

  double a1 = a11 + a12 + a13;
  double a2 = log(1);

  return  my_min(a1,a2);
}

double AlphaShuffleOrder_UniTS(arma::mat data, arma::vec new_order, arma::vec old_order,
                               double theta, double sigma, double phi,
                               double a, double b, double c){

  double a1 = Posterior_UniTS(data, new_order, theta, sigma, phi, a, b, c) - Posterior_UniTS(data, old_order, theta, sigma, phi, a, b, c);
  double a2 = log(1);

  return  my_min(a1,a2);
}

double AlphaSplitOrder_MultiTS(arma::mat data, arma::vec new_order, arma::vec old_order,
                               double q, double index, double gamma, double k_0, double nu_0,
                               double theta, double sigma,arma::mat phi_0, arma::vec m_0){

  double k = max(old_order) + 1, a11 = 0, a12 = 0, a13 = 0;
  arma::vec table_oldorder = table_cpp(old_order);

  if((k < data.n_cols) & (k > 1)){
    a11 = log((1-q)/q);
    a12 = Posterior_MultiTS(data, new_order, gamma, k_0, nu_0, theta, sigma, phi_0, m_0) - Posterior_MultiTS(data, old_order, gamma, k_0, nu_0, theta, sigma, phi_0, m_0);
    a13 = log(( table_oldorder.n_elem - std::count(table_oldorder.begin(),table_oldorder.end(),1) *(table_oldorder(index) - 1))/k);
  } else {
    a11 = log(1-q);
    a12 = Posterior_MultiTS(data, new_order, gamma, k_0, nu_0, theta, sigma, phi_0, m_0) - Posterior_MultiTS(data, old_order, gamma, k_0, nu_0, theta, sigma, phi_0, m_0);
    a13 = log(data.n_cols-1);
  }

  double a1 = a11 + a12 + a13;
  double a2 = log(1);
  double res = my_min(a1,a2);

  return res;
}

double AlphaMergeOrder_MultiTS(arma::mat data, arma::vec new_order, arma::vec old_order,
                               double q, double index, double gamma, double k_0, double nu_0,
                               double theta, double sigma,arma::mat phi_0, arma::vec m_0){

  double k = max(old_order) + 1, a11, a12, a13;
  arma::vec table_oldorder = table_cpp(old_order);

  if((k < data.n_cols) & (k > 1)){
    a11 = log(q/(1-q));
    a12 = Posterior_MultiTS(data, new_order, gamma, k_0, nu_0, theta, sigma, phi_0, m_0) - Posterior_MultiTS(data, old_order, gamma, k_0, nu_0, theta, sigma, phi_0, m_0);
    a13 = log((k-1) / ((table_oldorder.n_elem - std::count(table_oldorder.begin(),table_oldorder.end(),1) + 1) * (table_oldorder(index) + table_oldorder(index + 1) - 1)) );
  } else {
    a11 = log(q);
    a12 = Posterior_MultiTS(data, new_order, gamma, k_0, nu_0, theta, sigma, phi_0, m_0) - Posterior_MultiTS(data, old_order, gamma, k_0, nu_0, theta, sigma, phi_0, m_0);
    a13 = log(data.n_cols - 1);
  }

  double a1 = a11 + a12 + a13;
  double a2 = log(1);

  return  my_min(a1,a2);
}

double AlphaShuffleOrder_MultiTS(arma::mat data, arma::vec new_order, arma::vec old_order,
                                 double gamma, double k_0, double nu_0,
                                 double theta, double sigma,arma::mat phi_0, arma::vec m_0){

  double a1 = Posterior_MultiTS(data, new_order, gamma, k_0, nu_0, theta, sigma, phi_0, m_0) - Posterior_MultiTS(data, old_order, gamma, k_0, nu_0, theta, sigma, phi_0, m_0);
  double a2 = log(1);

  return  my_min(a1,a2);
}


// ------------------------------
// SPLIT-MERGE TO UPDATE ORDERS
// ------------------------------

Rcpp::List Split_cpp(arma::vec order){

  arma::vec freq_temp, res_order, temp_prob;
  double k, bound = 0, temp_obs;
  int temp_id;

  k = max(order) + 1;
  freq_temp.resize(k);

  for(arma::uword i = 0; i < freq_temp.n_elem; i++){
    if(std::count(order.begin(),order.end(),i) > 1){
      freq_temp(i) = 1;
    } else {
      freq_temp(i) = 0;
    }
  }

  temp_id = rint(freq_temp);

  temp_prob.resize(order.n_elem);
  temp_prob.fill(0.0);

  for(arma::uword i = 0; i < order.n_elem; i++){
    if(order(i) == temp_id){
      temp_prob(i) = 1.0;
      bound = i;
    }
  }

  temp_prob(bound) = 0.0;
  temp_obs = rint(temp_prob);

  res_order = order;
  for(arma::uword i = temp_obs + 1; i < order.n_elem; i++){
    res_order(i) += 1;
  }

  Rcpp::List out_list;
  out_list["split_index"] = temp_id;
  out_list["split_order"] = res_order;

  return out_list;
}

Rcpp::List Merge_cpp(arma::vec order){

  double k = max(order) + 1, temp_id;
  arma::vec freq_temp, merge_order;

  freq_temp.resize(k - 1);
  freq_temp.fill(1.0);
  temp_id = rint(freq_temp);
  merge_order = order;

  for(arma::uword i = 0; i < order.n_elem; i++){
    if(order(i) > temp_id){
      merge_order(i) -= 1;
    }
  }

  Rcpp::List out_list;
  out_list["merge_index"] = temp_id;
  out_list["merge_order"] = merge_order;

  return out_list;

}

Rcpp::List Shuffle_cpp(arma::vec order){

  arma::vec freq_temp, temp_prob;
  double k = max(order) + 1, bound = 0;

  freq_temp.resize(k);
  freq_temp.fill(1.0);
  freq_temp(k-1) = 0;

  double temp_id = rint(freq_temp);

  arma::vec shuffle_order = order;

  temp_prob.resize(order.n_elem);
  temp_prob.fill(0);

  double lower_bound = 0;

  for(arma::uword i = 0; i < order.n_elem; i++){
    if(shuffle_order(i) == temp_id || shuffle_order(i) == (temp_id + 1)){

      if(shuffle_order(i) == temp_id && lower_bound == 0){
        lower_bound = i;
      }

      temp_prob(i) = 1.0;
      bound = i;
      shuffle_order(i) = temp_id;
    }
  }

  temp_prob(bound) = 0;
  temp_prob(lower_bound) = 0;

  double temp_obs = rint(temp_prob);

  if(temp_obs <= 0){
    temp_obs = bound;
  }

  for(int i = temp_obs; i <= bound; i++){
    shuffle_order(i) = temp_id + 1;
  }

  Rcpp::List out_list;
  out_list["shuffle_index"] = temp_id;
  out_list["shuffle_order"] = shuffle_order;

  return out_list;

}

// -------------------
// ACCELERATION STEP
// -------------------

void SplitMergeAccUnivariateTS(arma::mat data,
                               arma::vec &order,
                               int iter, double q, double a, double b, double c, double gamma){

  arma::vec res_order, probs(2);

  for(int i = 0; i < iter; i++){

    double k = max(order) + 1;

    probs(0) = q * indicator_1(k,data.n_cols) + indicator_3(k);
    probs(1) = (1-q) * indicator_1(k,data.n_cols) + indicator_2(k, data.n_cols);

    probs(0) = probs(0)/(probs(0)+probs(1));
    probs(1) = probs(1)/(probs(0)+probs(1));

    double u = arma::randu();

    if(u <= probs(0)){

      // SPLIT

      /// propose a new order

      Rcpp::List split_list = Split_cpp(order);

      arma::vec split_order = split_list[1];
      int split_index = split_list[0];

      /// evaluate the proposed order

      double alpha_split = AlphaSplitOrder_TS(data, split_order, order, split_index, table_cpp(order), a, b, c, gamma, q);

      if(log(arma::randu()) <= alpha_split){
        res_order = split_order;
      } else {
        res_order = order;
      }

    } else {

      // MERGE

      /// propose a new order

      Rcpp::List merge_list = Merge_cpp(order);

      arma::vec merge_order = merge_list[1];
      int merge_index = merge_list[0];

      arma::vec freq_temp = table_cpp(order);

      /// evaluate the proposed order

      double alpha_merge = AlphaMergeOrder_TS(data, merge_order, order, merge_index, freq_temp, a, b, c, gamma, q);

      if(log(arma::randu()) <= alpha_merge){
        res_order = merge_order;
      } else {
        res_order = order;
      }

    }

    if(max(res_order) > 0){

      // SHUFFLE

      /// propose a new order

      Rcpp::List shuffle_list = Shuffle_cpp(res_order);

      arma::vec shuffle_order = shuffle_list[1];

      /// evaluate the proposed order

      double alpha_shuffle = AlphaShuffleOrder_TS(data, shuffle_order, order, a, b, c, gamma, q);

      if(log(arma::randu()) <= alpha_shuffle){
        res_order = shuffle_order;
      } else {
        res_order = order;
      }
    }

    order = res_order;
  }

}

void SplitMergeAccMultivariateTS(arma::cube data,
                                 arma::vec &order,
                                 int iter, double q, double k_0, double nu_0, arma::mat phi_0, arma::vec m_0, double gamma){

  arma::vec res_order, probs(2);

  for(int i = 0; i < iter; i++){

    double k = max(order) + 1;

    probs(0) = q * indicator_1(k,data.slice(0).n_cols) + indicator_3(k);
    probs(1) = (1-q) * indicator_1(k,data.slice(0).n_cols) + indicator_2(k, data.slice(0).n_cols);

    probs(0) = probs(0)/(probs(0)+probs(1));
    probs(1) = probs(1)/(probs(0)+probs(1));

    double u = arma::randu();

    if(u <= probs(0)){

      // SPLIT

      /// propose a new order

      Rcpp::List split_list = Split_cpp(order);

      arma::vec split_order = split_list[1];
      int split_index = split_list[0];

      /// evaluate the proposed order

      double alpha_split = AlphaSplitOrder_MultiTS_Acc(data, split_order, order, q, split_index, gamma, k_0, nu_0, phi_0, m_0);

      if(log(arma::randu()) <= alpha_split){
        res_order = split_order;
      } else {
        res_order = order;
      }

    } else {

      // MERGE

      /// propose a new order

      Rcpp::List merge_list = Merge_cpp(order);

      arma::vec merge_order = merge_list[1];
      int merge_index = merge_list[0];

      /// evaluate the proposed order

      double alpha_merge = AlphaMergeOrder_MultiTS_Acc(data, merge_order, order, merge_index, q, gamma, k_0, nu_0, phi_0, m_0);



      if(log(arma::randu()) <= alpha_merge){
        res_order = merge_order;
      } else {
        res_order = order;
      }

    }

    if(max(res_order) > 0){

      // SHUFFLE

      /// propose a new order

      Rcpp::List shuffle_list = Shuffle_cpp(res_order);

      arma::vec shuffle_order = shuffle_list[1];

      /// evaluate the proposed order

      double alpha_shuffle = AlphaShuffleOrder_MultiTS(data, shuffle_order, order, gamma, k_0, nu_0, phi_0, m_0);


      if(log(arma::randu()) <= alpha_shuffle){
        res_order = shuffle_order;
      } else {
        res_order = order;
      }
    }

    //Rcpp::Rcout << table_cpp(res_order).t() << "\t";

    order = res_order;
  }

}

// ------------------------
// NORMALISATION CONSTANT
// ------------------------


arma::vec norm_constant_uni(arma::mat data,
                            double gamma_par,
                            int R,
                            double a,
                            double b,
                            double c,
                            double p,
                            bool print_progress = true){
  arma::vec temp_llik_vec(data.n_rows), freqs, temp_probs, cfreq, new_order_vec;
  arma::mat temp_llik_mat(R,data.n_rows), curve_mat, new_order_mat(1,data.n_cols);
  int T = data.n_cols;
  double ord_lprob;

  int start_s = clock();
  int current_s;
  int nupd = round(R / 10);
  for(int r = 0; r < R; r++){

    double k = 1 + rbinom(T, p);
    temp_probs.resize(k);
    temp_probs.fill(1/k);

    freqs = rmultin(T-k, temp_probs);

    ord_lprob = dmultinom_log_cpp(freqs, temp_probs) + dbinom_log_cpp(freqs.n_elem,T-1,0.5);

    for(arma::uword i = 0; i < freqs.n_elem; i++){
      freqs(i) = freqs(i) + 1;
    }

    new_order_vec.resize(T);

    cfreq = cumsum(freqs);
    for(int i = 0; i < cfreq(0); i++){
      new_order_vec(i) = 0;
    }

    for(arma::uword j = 1; j < cfreq.n_elem; j++){
      for(int i = cfreq(j-1); i < cfreq(j); i++){

        new_order_vec(i) = j;
      }
    }

    while(new_order_vec(0) > 0){
      new_order_vec(0) -= 1;
    }

    new_order_mat.row(0) = new_order_vec.t();

    for(arma::uword i = 0; i < data.n_rows; i++){
      temp_llik_mat(r,i) = LogLikelihood_TS(data.row(i),new_order_mat.row(0),gamma_par,a,b,c) - ord_lprob;
    }

    // print time
    if(((r + 1) % nupd == 0) & (print_progress == true) ){
      current_s = clock();
      Rcpp::Rcout << "Normalization constant - completed:\t" << (r + 1) << "/" << R << " - in " <<
        double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
    }
    Rcpp::checkUserInterrupt();
  }

  for(arma::uword i = 0; i < data.n_rows; i++){
    temp_llik_vec(i) = log_sum_exp(temp_llik_mat.col(i)) + log(R) - (T-1)*log(2);
  }

  return temp_llik_vec;
}


arma::vec norm_constant_multi(arma::cube data,
                                double gamma_par,
                                int R,
                                double k_0,
                                double nu_0,
                                arma::mat phi_0,
                                arma::vec m_0,
                                double p,
                                bool print_progress = true){
  arma::vec temp_llik_vec(data.n_slices), freqs, temp_probs, cfreq, new_order_vec;
  arma::mat temp_llik_mat(R,data.n_slices), curve_mat, new_order_mat(1,data.slice(0).n_cols);
  int T = data.slice(0).n_cols;
  double ord_lprob;

  int start_s = clock();
  int current_s;
  int nupd = round(R / 10);
  for(int r = 0; r < R; r++){

    double k = 1 + rbinom(T, p);
    temp_probs.resize(k);
    temp_probs.fill(1/k);

    freqs = rmultin(T-k, temp_probs);

    ord_lprob = dmultinom_log_cpp(freqs, temp_probs) + dbinom_log_cpp(freqs.n_elem,T-1,0.5);

    for(arma::uword i = 0; i < freqs.n_elem; i++){
      freqs(i) = freqs(i) + 1;
    }

    new_order_vec.resize(T);

    cfreq = cumsum(freqs);
    for(int i = 0; i < cfreq(0); i++){
      new_order_vec(i) = 0;
    }

    for(arma::uword j = 1; j < cfreq.n_elem; j++){
      for(int i = cfreq(j-1); i < cfreq(j); i++){

        new_order_vec(i) = j;
      }
    }

    while(new_order_vec(0) > 0){
      new_order_vec(0) -= 1;
    }

    new_order_mat.row(0) = new_order_vec.t();

    for(arma::uword i = 0; i < data.n_slices; i++){
      temp_llik_mat(r,i) = Likelihood_MultiTS(data.slice(i),new_order_mat.row(0).t(),gamma_par,k_0,nu_0,phi_0,m_0) - ord_lprob;
    }

    // print time
    if(((r + 1) % nupd == 0) & (print_progress == true)){
      current_s = clock();
      Rcpp::Rcout << "Normalization constant - completed:\t" << (r + 1) << "/" << R << " - in " <<
        double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
    }
    Rcpp::checkUserInterrupt();
  }

  for(arma::uword i = 0; i < data.n_slices; i++){
    temp_llik_vec(i) = log_sum_exp(temp_llik_mat.col(i)) + log(R) - (T-1)*log(2);
  }

  return temp_llik_vec;
}



arma::vec norm_constant_epi(arma::mat data,
                             double gamma_par,
                             int num_orders,
                             double a0,
                             double b0,
                             arma::vec rho,
                             int M,
                             double dt,
                             gsl_rng *r,
                             bool print_progress = true,
                             double S0 = 1,
                             double R0 = 0,
                             double p = 0.03){
  arma::vec temp_llik_vec(data.n_rows), freqs, temp_probs, cfreq, new_order_vec;
  arma::mat temp_llik_mat(num_orders,data.n_rows), curve_mat, new_order_mat(1,data.n_cols);
  int T = data.n_cols;
  double ord_lprob;

  int start_s = clock();
  int current_s;
  int nupd = round(num_orders / 10);
  for(int r_iter = 0; r_iter < num_orders; r_iter++){

    double k = 1 + rbinom(data.n_cols, p);
    temp_probs.resize(k);
    temp_probs.fill(1/k);

    freqs = rmultin(data.n_cols-k, temp_probs);

    ord_lprob = dmultinom_log_cpp(freqs, temp_probs) + dbinom_log_cpp(freqs.n_elem,T-1,0.5);

    for(arma::uword i = 0; i < freqs.n_elem; i++){
      freqs(i) = freqs(i) + 1;
    }

    new_order_vec.resize(data.n_cols);

    cfreq = cumsum(freqs);
    for(int i = 0; i < cfreq(0); i++){
      new_order_vec(i) = 0;
    }

    for(arma::uword j = 1; j < cfreq.n_elem; j++){
      for(int i = cfreq(j-1); i < cfreq(j); i++){

        new_order_vec(i) = j;
      }
    }

    while(new_order_vec(0) > 0){
      new_order_vec(0) -= 1;
    }

    new_order_mat.row(0) = new_order_vec.t();

    arma::vec new_order = generate_random_order(data.n_cols, 2/data.n_cols, r);

    for(arma::uword i = 0; i < data.n_rows; i++){
      curve_mat = integrated_curves_mat(dt, new_order, a0, b0, gamma_par, rho(i), M, S0, R0);
      temp_llik_mat(r_iter,i) = log_sum_exp(curve_mat.cols(0,data.n_cols-1) * data.row(i).t() - curve_mat.col(data.n_cols)) - log(M)  - ord_lprob;
    }

    // print time
    if(((r_iter + 1) % nupd == 0) & (print_progress == true)){
      current_s = clock();
      Rcpp::Rcout << "Normalization constant - completed:\t" << (r_iter + 1) << "/" << num_orders << " - in " <<
        double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
    }
    Rcpp::checkUserInterrupt();
  }

  for(arma::uword i = 0; i < data.n_rows; i++){
    temp_llik_vec(i) = log_sum_exp(temp_llik_mat.col(i)) + log(num_orders) - (T-1)*log(2);
  }

  return temp_llik_vec;
}


// -------------------
// UPDATE PARAMETERS
// -------------------

double FullConditionalGamma(arma::mat data, arma::vec order, double gamma,
                            double theta, double sigma, double k_0, double nu_0,
                            arma::mat phi_0, arma::vec m_0){

  double k = max(order) + 1, n = data.n_cols, d = data.n_rows;

  arma::vec vec_num(k-1), vec_dx(k);
  arma::vec table_order = table_cpp(order);

  for(int i = 0; i < k-1; i++){
    vec_num(i) = theta + (i+1)*sigma;
  }

  double num = gsl_sf_lnfact(n) + (k*d/2)*log(k_0) + (k*nu_0/2) * log(arma::det(phi_0)) + sum(log(vec_num));

  double den = gsl_sf_lnfact(k) + gsl_sf_lnpoch((theta+1),(n-1)) + (n*d/2) * log(M_PI) + ((n-k)/2)*log(1-std::pow(gamma,2)) + k * lgamma_multi(d, nu_0/2);

  for(int i = 0; i < k; i++){

    double n_i = table_order(i);

    double k_ni = compute_k_n(k_0, gamma, n_i);

    arma::mat gamma_k = ExtractSubData(data, order, i);

    arma::mat S_ni = compute_s_n(d, m_0, k_0, nu_0, phi_0, gamma, n_i, gamma_k.t(), d);

    //
    double num_vec_dx = lgamma_multi(d,(nu_0 + n_i)/2) + gsl_sf_lnpoch((1-sigma),(n_i-1)) ;
    double den_vec_dx = gsl_sf_lnfact(n_i) + (d/2)*(k_ni) + ((nu_0 + n_i)/2) * log(det(S_ni));

    vec_dx(i) = num_vec_dx - den_vec_dx;

  }

  double res = num - den + sum(vec_dx);

  return res;

}

double FullConditionalSigma(arma::vec order, double theta, double sigma,
                            double a, double b, double c, double d){


  double k = max(order) + 1;
  arma::vec table_order = table_cpp(order);

  arma::vec vec_1(k-1);
  arma::vec vec_2(k);

  for(int i = 0; i < k - 1; i++){
    vec_1(i) = log(theta + (i+1)*sigma);
  }

  for(int i = 0; i < k; i++){
    vec_2(i) = gsl_sf_lnpoch(1 - sigma, table_order(i) - 1);
  }

  double res = (a-1) * log(sigma) + (b-1) * log(1-sigma) + (c-1) * log(theta + sigma) + log(exp(-d*sigma)) + sum(vec_1) + sum(vec_2);

  return res;

}

void UpdateGamma(double gamma_old, arma::mat data, arma::vec order,
                 double theta, double sigma, double k_0, double nu_0,
                 arma::mat phi_0, arma::vec m_0,
                 arma::vec &gamma_inf, arma::vec &gamma_inf_10, gsl_rng *r, double prior_var_gamma){

  gamma_inf.resize(gamma_inf.n_elem + 1);
  gamma_inf_10.resize(gamma_inf_10.n_elem + 1);


  double tau = log(gamma_old/(1-gamma_old));

  double tau_star = tau + gsl_ran_gaussian(r,prior_var_gamma); // How can I let the user set a different sd instead of 1?

  double gamma_new = exp(tau_star)/(1+exp(tau_star));

  double deriv_inv_star = abs(exp(tau_star) / std::pow(1+exp(tau_star),2));

  double deriv_inv = abs(exp(tau) / std::pow(1+exp(tau),2));

  double alpha_MH = FullConditionalGamma(data, order, gamma_new, theta, sigma, k_0,
                                         nu_0, phi_0, m_0) + log(deriv_inv_star) - FullConditionalGamma(data, order, gamma_old, theta, sigma, k_0,
                                         nu_0, phi_0, m_0) - log(deriv_inv);

  if(log(arma::randu()) <= my_min(alpha_MH,log(1))){
    gamma_inf(gamma_inf.n_elem - 1) = gamma_new;
    gamma_inf_10(gamma_inf_10.n_elem - 1) = 1;
  } else {
    gamma_inf(gamma_inf.n_elem - 1) = gamma_old;
    gamma_inf_10(gamma_inf_10.n_elem - 1) = 0;
  }

}

void UpdateSigma(arma::vec order, double theta, double sigma,
                 arma::vec &sigma_inf, arma::vec &sigma_inf_10, gsl_rng *r){

  sigma_inf.resize(sigma_inf.n_elem + 1);
  sigma_inf_10.resize(sigma_inf_10.n_elem + 1);

  double sigma_new = gsl_ran_beta(r, 1, 1);
  double alpha_MH = FullConditionalSigma(order, theta, sigma_new, 1, 1, 1, 1) - FullConditionalSigma(order, theta, sigma, 1, 1, 1, 1);

  if(log(arma::randu()) <= my_min(alpha_MH,log(1))){
    sigma_inf(sigma_inf.n_elem - 1) = sigma_new;
    sigma_inf_10(sigma_inf_10.n_elem - 1) = 1;
  } else {
    sigma_inf(sigma_inf.n_elem - 1) = sigma;
    sigma_inf_10(sigma_inf_10.n_elem - 1) = 0;
  }

}

void UpdateTheta(double theta,double sigma, arma::vec order,  arma::vec &theta_inf,  double prior_theta_c, double prior_theta_d, gsl_rng *r){

  theta_inf.resize(theta_inf.n_elem + 1);

  double k = max(order) + 1;
  double n = order.n_elem;
  arma::vec vec(k+1);

  double z = gsl_ran_beta(r, theta + 2, n);
  double f = gsl_ran_exponential(r, theta + 1);

  for(int i = 0; i < k + 1; i++){
    double omega_j_num = (n - sigma) * (n+1-sigma) * AbsStirling1st(k-1,i) + (2*n + 1 - 2*sigma) * sigma * AbsStirling1st(k-1,i-1) + std::pow(sigma,2) * AbsStirling1st(k-1,i-2) + gsl_sf_gamma(prior_theta_c + i);
    double omega_j_den = std::pow(sigma*(prior_theta_d + f - log(z)), i);
    vec(i) = omega_j_num / omega_j_den;

  }

  vec = vec/sum(vec);

  double u = gsl_rng_uniform(r);

  int component = 0;

  for(arma::uword i = 0; i < vec.n_elem; i++){

    if((vec(i) > u) & (component == 0)){
      component = i;
    }

  }

  theta_inf(theta_inf.n_elem - 1) = rshiftedgamma(prior_theta_c + (component - 1), prior_theta_d + f - std::log(z), sigma, r) ;

}



void update_rho(arma::mat data,
                arma::vec &rho,
                double a0,
                double b0,
                double c0,
                double d0,
                double MH_var,
                double gamma,
                double dt,
                int M,
                double S0,
                double R0,
                arma::vec &llik,
                arma::vec clust,
                arma::mat orders){

  for(arma::uword i = 0; i < data.n_rows; i++){

    double rho_temp = rho(i);
    double llik_temp = llik(i);
    double new_llik;
    double log_rho = log(rho_temp), log_rho_new, rho_new, acc_rate;
    int T = data.n_cols;
    arma::mat curve_mat;

    log_rho_new = log_rho + arma::randn() * sqrt(MH_var);
    rho_new = exp(log_rho_new);


    int clust_obs = clust(i);
    curve_mat = integrated_curves_mat(dt, orders.row(clust_obs).t(), a0, b0, gamma, rho_new, M, S0, R0);

    new_llik = log_sum_exp(curve_mat.cols(0,T-1) * data.row(i).t() - curve_mat.col(T)) - log(M);

    acc_rate = my_min(0, new_llik - llik_temp + (b0 - 1) * log(rho_new) - d0 * rho_new - (b0 - 1) * log(rho_temp) + d0 * rho_temp + rho_new - rho_temp);

    if(log(arma::randu()) < acc_rate){
      rho(i) = rho_new;
      llik(i) = new_llik;
    }
  }
}

//------------------------------------
// FUNCTIONS FOR CLUSTERING EPI DATA
//------------------------------------

void update_single_order(arma::mat data,
                         arma::vec clust,
                         int clust_id,
                         arma::mat &orders,
                         arma::vec &llik,
                         double q,
                         double dt,
                         double a0,
                         double b0,
                         double gamma,
                         arma::vec rho,
                         int M,
                         double S0 = 1,
                         double R0 = 0,
                         double coars = 1){

  int k = max(orders.row(clust_id)) + 1, temp_id, temp_obs,
    T = orders.n_cols, bound = 0, temp_count;
  arma::vec temp_llik = llik, temp_prob, new_order, freq_temp;
  arma::mat curve_mat;
  bool check;
  double u = arma::randu(), acc_rate;

  if(k == 1 || (u < q && k < T)){
    check = true;
  } else {
    check = false;
  }

  if(check == true){

    freq_temp.resize(k);
    for(arma::uword l = 0; l < freq_temp.n_elem; l++){
      if(std::count(orders.row(clust_id).begin(), orders.row(clust_id).end(), l) > 1){
        freq_temp(l) = 1;
      } else {
        freq_temp(l) = 0;
      }
    }

    temp_id = rint(freq_temp);
    temp_prob.resize(T);
    temp_prob.fill(0.0);
    for(int i = 0; i < T; i++){
      if(orders(clust_id,i) == temp_id){
        temp_prob(i) += 1.0;
        bound = i;
      }
    }
    temp_prob(bound) = 0.0;
    temp_obs = rint(temp_prob);
    new_order = orders.row(clust_id).t();
    for(int i = temp_obs + 1; i < T; i++){
      new_order(i) += 1;
    }

    for(arma::uword i = 0; i < clust.n_elem; i++){
      if(clust(i) == clust_id){
        curve_mat = integrated_curves_mat(dt, new_order, a0, b0, gamma, rho(i), M, S0 = 1, R0 = 0);
        temp_llik(i) = log_sum_exp(curve_mat.cols(0,T-1) * data.row(i).t() - curve_mat.col(T)) - log(M);
      }
    }

    if(k == 1){
      acc_rate = my_min(0, log(1 - q) - log(q) + (coars * sum(temp_llik)) - (coars * sum(llik)) +
        log(sum(freq_temp)) + log(std::count(orders.row(clust_id).begin(), orders.row(clust_id).end(), temp_id)) - log(k));
    } else {
      acc_rate = my_min(0, log(1 - q) + log(T - 1) + (coars * sum(temp_llik)) - (coars * sum(llik)));
    }

    u = arma::randu();
    if(u <= exp(acc_rate)){
      orders.row(clust_id) = new_order.t();
      llik = temp_llik;
    }

  } else {

    freq_temp.resize(k - 1);
    freq_temp.fill(1.0);
    temp_id = rint(freq_temp);
    new_order = orders.row(clust_id).t();
    for(int i = 0; i < T; i++){
      if(orders(clust_id, i) > temp_id){
        new_order(i) -= 1;
      }
    }

    for(arma::uword i = 0; i < clust.n_elem; i++){
      if(clust(i) == clust_id){
        curve_mat = integrated_curves_mat(dt, new_order, a0, b0, gamma, rho(i), M, S0 = 1, R0 = 0);
        temp_llik(i) = log_sum_exp(curve_mat.cols(0,T-1) * data.row(i).t() - curve_mat.col(T)) - log(M);

      }
    }

    temp_count = 0;
    for(int i = 0; i < k - 1; i++){
      if(std::count(orders.row(clust_id).begin(),orders.row(clust_id).end(), i) > 1){
        temp_count += 1;
      }
    }

    if(k == 1){
      acc_rate = my_min(0, log(q) - log(1 - q) + (coars * sum(temp_llik)) - (coars * sum(llik)) +
        log(k - 1) - log(temp_count) - log(std::count(new_order.begin(), new_order.end(), temp_id) - 1));
    } else {
      acc_rate = my_min(0, log(q) + log(T - 1) + (coars * sum(temp_llik)) - (coars * sum(llik)));
    }

    u = arma::randu();
    if(u <= exp(acc_rate)){
      orders.row(clust_id) = new_order.t();
      llik = temp_llik;
    }
  }

  k = max(orders.row(clust_id)) + 1;
  if(k > 1){
    freq_temp.resize(k - 1);
    freq_temp.fill(1.0);
    temp_id = rint(freq_temp);

    temp_prob.resize(T);
    temp_prob.fill(0.0);
    for(int i = 0; i < T; i++){
      if(orders(clust_id,i) == temp_id || orders(clust_id,i) == temp_id + 1){
        temp_prob(i) += 1.0;
        bound = i;
      }
    }
    temp_prob(bound) = 0.0;
    temp_obs = rint(temp_prob);
    new_order = orders.row(clust_id).t();
    for(int i = temp_obs + 1; i < T; i++){
      if(new_order(i) == temp_id || new_order(i) == temp_id + 1){
        new_order(i) = temp_id + 1;
      }
    }

    for(arma::uword i = 0; i < clust.n_elem; i++){
      if(clust(i) == clust_id){
        curve_mat = integrated_curves_mat(dt, new_order, a0, b0, gamma, rho(i), M, S0, R0);
        temp_llik(i) = log_sum_exp(curve_mat.cols(0,T-1) * data.row(i).t() - curve_mat.col(T)) - log(M);
      }
    }

    acc_rate = my_min(0, (coars * sum(temp_llik)) - (coars * sum(llik)));

    u = arma::randu();
    if(u <= exp(acc_rate)){
      orders.row(clust_id) = new_order.t();
      llik = temp_llik;
    }
  }
}

void update_partition(arma::mat data,
                      arma::vec &clust,
                      arma::mat &orders,
                      arma::vec &llik,
                      arma::vec norm_const,
                      double alpha,
                      double p,
                      double q,
                      double dt,
                      double a0,
                      double b0,
                      double gamma,
                      arma::vec rho,
                      int M,
                      int L,
                      double S0 = 1,
                      double R0 = 0,
                      double coars = 1){
  arma::vec temp_llik = llik, temp_clust = clust, freq_temp, prob_temp(clust.n_elem),
    temp_vec1(clust.n_elem), temp_vec2(clust.n_elem), temp_vec3(clust.n_elem);
  int id3, id4, k, u_bound;
  arma::uword id1, id2;
  double acc_rate;
  arma::mat temp_order(2, data.n_cols), curve_mat1, curve_mat2, curve_mat3;
  temp_order.fill(0);

  freq_temp.resize(clust.n_elem);
  freq_temp.fill(1.0);
  id1 = rint(freq_temp);
  freq_temp(id1) = 0.0;
  id2 = rint(freq_temp);
  k = max(clust) + 1;

  if(clust(id1) != clust(id2)){

    prob_temp.fill(0.0);
    for(arma::uword i = 0; i < clust.n_elem; i++){
      if((clust(i) == clust(id1)) || (clust(i) == clust(id2))){
        prob_temp(i) = 1.0;
      }
    }

    id3 = rint(prob_temp);

    temp_order.row(0) = orders.row(clust(id3));
    temp_clust.fill(1);

    temp_clust(id3) = 0;

    for(int l = 0; l < L; l++){
      update_single_order(data, temp_clust, 0, temp_order, temp_llik, q, dt,
                          a0, b0, gamma, rho, M, S0, R0);
    }

    for(arma::uword i = 0; i < clust.n_elem; i++){
      curve_mat1 = integrated_curves_mat(dt, temp_order.row(0).t(), a0, b0, gamma, rho(i), M, S0, R0);
      curve_mat2 = integrated_curves_mat(dt, orders.row(clust(id1)).t(), a0, b0, gamma, rho(i), M, S0, R0);
      curve_mat3 = integrated_curves_mat(dt, orders.row(clust(id2)).t(), a0, b0, gamma, rho(i), M, S0, R0);
      if((clust(i) == clust(id1)) || (clust(i) == clust(id2))){
        int T = data.n_cols;
        temp_llik(i) = log_sum_exp(curve_mat1.cols(0,T-1) * data.row(i).t() - curve_mat1.col(T)) - log(M);
        temp_vec1(i) = temp_llik(i) - norm_const(i);

        if(clust(i) == clust(id1)){
          llik(i) = log_sum_exp(curve_mat2.cols(0,T-1) * data.row(i).t() - curve_mat2.col(T)) - log(M);

          temp_vec2(i) = log_sum_exp(curve_mat2.cols(0,T-1) * data.row(i).t() - curve_mat2.col(T)) - log(M) - norm_const(i);
        } else if (clust(i) == clust(id2)){
          llik(i) = log_sum_exp(curve_mat3.cols(0,T-1) * data.row(i).t() - curve_mat3.col(T)) - log(M);
          temp_vec2(i) = log_sum_exp(curve_mat3.cols(0,T-1) * data.row(i).t() - curve_mat3.col(T)) - log(M) - norm_const(i);
        }

      } else {
        int T = data.n_cols;
        temp_vec1(i) = log_sum_exp(curve_mat1.cols(0,T-1) * data.row(i).t() - curve_mat1.col(T)) - log(M) - norm_const(i);
        temp_vec2(i) = log_sum_exp(curve_mat1.cols(0,T-1) * data.row(i).t() - curve_mat1.col(T)) - log(M) - norm_const(i);
      }
    }


    arma::vec vec1 = arma::regspace(1,k);
    arma::vec vec2 = arma::regspace(0,k-1);
    arma::vec vec3 = arma::regspace(1,k+1);
    arma::vec vec4 = arma::regspace(0,k);

    arma::vec vecNpar1(vec1.n_elem);
    vecNpar1.fill(std::pow(2,(data.n_cols-1)));

    arma::vec vecNpar2(vec2.n_elem);
    vecNpar2.fill(std::pow(2,(data.n_cols-1)));

    arma::vec vecNpar3(vec3.n_elem);
    vecNpar3.fill(std::pow(2,(data.n_cols-1)));

    arma::vec vecNpar4(vec4.n_elem);
    vecNpar4.fill(std::pow(2,(data.n_cols-1)));


    vecNpar1 = vecNpar1 - vec1;
    vecNpar2 = vecNpar2 - vec2;
    vecNpar3 = vecNpar3 - vec3;
    vecNpar4 = vecNpar4 - vec4;

    vecNpar1 = log(vecNpar1);
    vecNpar2 = log(vecNpar2);
    vecNpar3 = log(vecNpar3);
    vecNpar4 = log(vecNpar4);

    arma::vec diffvec1 = vecNpar1 - vecNpar2;
    arma::vec diffvec2 = vecNpar3 - vecNpar4;
    arma::vec diffvec3 = vecNpar3 - vecNpar4;

    int sum1 = sum(diffvec1);
    int sum2 = sum(diffvec2);
    int sum3 = sum(diffvec3);

    arma::vec sum1vec(2);
    sum1vec(0) = log(1);
    sum1vec(1) = -sum1;

    arma::vec sum2vec(2);
    sum2vec(0) = log(1);
    sum2vec(1) = -sum2;

    arma::vec sum3vec(2);
    sum3vec(0) = log(1);
    sum3vec(1) = -sum3;

    acc_rate = my_min(0, ((std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id1)) +  std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id2)) - 2) * log(0.5) +
      lgamma(alpha) +
      log_sum_exp(sum1vec) -
      log_sum_exp(sum2vec) +
      lgamma(alpha + std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id1)) +  std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id2))) -
      lgamma(alpha + std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id1))) - lgamma(alpha + std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id2))) +
      (coars * sum(temp_llik)) - (coars * sum(llik)) +
      log_sum_exp(temp_vec1) - log_sum_exp(temp_vec2)));


    if(log(arma::randu()) < acc_rate){
      clust.elem(find(clust == clust(id2))).fill(clust(id1));
      orders.row(clust(id1)) = temp_order.row(0);
      llik = temp_llik;
    }

  } else {


    prob_temp.fill(0.0);
    for(arma::uword i = 0; i < clust.n_elem; i++){
      if(clust(i) == clust(id1)){
        prob_temp(i) = 1.0;
      }
    }

    id3 = rint(prob_temp);
    prob_temp(id3) = 0;
    id4 = rint(prob_temp);

    prob_temp(id3) = 1;

    prob_temp(id1) = 0;
    prob_temp(id2) = 0;

    temp_clust(id2) = k;

    for(arma::uword i = 0; i < clust.n_elem; i++){
      if(prob_temp(i) == 1){

        if(arma::randu() < 0.5){
          temp_clust(i) = temp_clust(id1);
        } else {
          temp_clust(i) = temp_clust(id2);
        }
      }
    }


    temp_order.row(0)= orders.row(temp_clust(id1));
    temp_order.row(1) = orders.row(temp_clust(id1));

    arma::vec temp_clust_update(clust.n_elem);
    temp_clust_update.fill(2);

    temp_clust_update(id3) = 0;
    temp_clust_update(id4) = 1;


    for(int l = 0; l < L; l++){

      update_single_order(data, temp_clust_update, 0, temp_order, temp_llik, q, dt,
                          a0, b0, gamma, rho, M, S0, R0);



      update_single_order(data, temp_clust_update, 1, temp_order, temp_llik, q, dt,
                          a0, b0, gamma, rho, M, S0, R0);
    }

    for(arma::uword i = 0; i < temp_clust.n_elem; i++){

      if((temp_clust(i) == clust(id1)) && (i != id1) && (i != id2)){

        if(arma::randu() < 0.5){
          temp_clust(i) = clust(id1);
        } else {
          temp_clust(i) = k;
        }

      }


    }

    for(arma::uword i = 0; i < clust.n_elem; i++){
      curve_mat1 = integrated_curves_mat(dt, temp_order.row(0).t(), a0, b0, gamma, rho(i), M, S0, R0);
      curve_mat2 = integrated_curves_mat(dt, temp_order.row(1).t(), a0, b0, gamma, rho(i), M, S0, R0);
      curve_mat3 = integrated_curves_mat(dt, orders.row(clust(id1)).t(), a0, b0, gamma, rho(i), M, S0, R0);
      if(temp_clust(i) == clust(id1)){
        int T = data.n_cols;
        temp_llik(i) = log_sum_exp(curve_mat1.cols(0,T-1) * data.row(i).t() - curve_mat1.col(T)) - log(M);

        llik(i) = log_sum_exp(curve_mat3.cols(0,T-1) * data.row(i).t() - curve_mat3.col(T)) - log(M);
        temp_vec1(i) = temp_llik(i) - norm_const(i);
        temp_vec3(i) = log_sum_exp(curve_mat3.cols(0,T-1) * data.row(i).t() - curve_mat3.col(T)) - log(M) - norm_const(i);

      } else if(temp_clust(i) == k){

        int T = data.n_cols;
        temp_llik(i) = log_sum_exp(curve_mat2.cols(0,T-1) * data.row(i).t() - curve_mat2.col(T)) - log(M);

        llik(i) = log_sum_exp(curve_mat3.cols(0,T-1) * data.row(i).t() - curve_mat3.col(T)) - log(M);

        temp_vec1(i) = temp_llik(i) - norm_const(i);;
        temp_vec3(i) = log_sum_exp(curve_mat3.cols(0,T-1) * data.row(i).t() - curve_mat3.col(T)) - log(M) - norm_const(i);


      } else {
        temp_vec1(i) = temp_llik(i) - norm_const(i);
        temp_vec3(i) = temp_llik(i) - norm_const(i);
      }
    }

    arma::vec vec1 = arma::regspace(1,k+1);
    arma::vec vec2 = arma::regspace(0,k);
    arma::vec vec3 = arma::regspace(1,k);
    arma::vec vec4 = arma::regspace(0,k-1);

    arma::vec vecNpar1(vec1.n_elem);
    vecNpar1.fill(std::pow(2,(data.n_cols-1)));

    arma::vec vecNpar2(vec2.n_elem);
    vecNpar2.fill(std::pow(2,(data.n_cols-1)));

    arma::vec vecNpar3(vec3.n_elem);
    vecNpar3.fill(std::pow(2,(data.n_cols-1)));

    arma::vec vecNpar4(vec4.n_elem);
    vecNpar4.fill(std::pow(2,(data.n_cols-1)));

    vecNpar1 = vecNpar1 - vec1;
    vecNpar2 = vecNpar2 - vec2;
    vecNpar3 = vecNpar3 - vec3;
    vecNpar4 = vecNpar4 - vec4;

    vecNpar1 = log(vecNpar1);
    vecNpar2 = log(vecNpar2);
    vecNpar3 = log(vecNpar3);
    vecNpar4 = log(vecNpar4);

    arma::vec diffvec1 = vecNpar1 - vecNpar2;
    arma::vec diffvec2 = vecNpar3 - vecNpar4;
    arma::vec diffvec3 = vecNpar3 - vecNpar4;

    int sum1 = sum(diffvec1);
    int sum2 = sum(diffvec2);
    int sum3 = sum(diffvec3);

    arma::vec sum1vec(2);
    sum1vec(0) = log(1);
    sum1vec(1) = -sum1;

    arma::vec sum2vec(2);
    sum2vec(0) = log(1);
    sum2vec(1) = -sum2;

    arma::vec sum3vec(2);
    sum3vec(0) = log(1);
    sum3vec(1) = -sum3;

    acc_rate = my_min(0, (- (std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id1)) +  std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id2)) - 2) * log(0.5) +
      - lgamma(alpha + std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id1)) + std::count(temp_clust.begin(),temp_clust.end(),temp_clust(id2))) +
      log_sum_exp(sum1vec) -
      log_sum_exp(sum3vec) +
      (coars * sum(temp_llik)) - (coars * sum(llik)) + log_sum_exp(temp_vec1) -
      log_sum_exp(temp_vec3)));

    if(log(arma::randu()) < acc_rate){
      clust = temp_clust;
      orders.resize(k + 1, orders.n_cols);
      orders.row(clust(id1)) = temp_order.row(0);

      orders.row(k) = temp_order.row(1);
      llik = temp_llik;
    }


  }

  k = orders.n_rows;
  for(int i = 0; i < k; i++){

    if((int) std::count(clust.begin(), clust.end(), i) == 0){
      for(int j = k; j > i; j--){
        if((int) std::count(clust.begin(), clust.end(), j) != 0){
          clust(find(clust == j)).fill(i);
          orders.swap_rows(i,j);
          break;
        }
      }
    }
  }

  u_bound = 1;
  for(int i = 1; i < k + 1; i++){
    if(std::count(clust.begin(), clust.end(),i) > 0){
      u_bound += 1;
    }
  }
  orders.resize(u_bound, orders.n_cols);

}

Rcpp::List marginal_CP(arma::mat data,
                       int niter,
                       int nburn,
                       double alpha,
                       double q,
                       double dt,
                       double a0,
                       double b0,
                       double c0,
                       double d0,
                       double gamma,
                       double MH_var,
                       int M,
                       int R,
                       int L,
                       double S0 = 1,
                       double R0 = 0,
                       double p = 0.003,
                       int nupd = 0,
                       double coars = 1,
                       unsigned long user_seed = 1234){


  arma::vec rho(data.n_rows);
  rho.fill(0.001);

  arma::vec clust(data.n_rows), llik(data.n_rows);
  clust.fill(0);

  arma::mat orders(data.n_rows, data.n_cols);
  orders.fill(0);

  // set seed for gsl random distribution generator
  const gsl_rng_type * T;
  gsl_rng *r;
  gsl_rng_env_setup();
  T = gsl_rng_default; // Generator setup
  r = gsl_rng_alloc (T);
  gsl_rng_set(r, user_seed);
  //

  for(arma::uword i = 0; i < orders.n_rows; i++){

    orders.row(i) = generate_random_order(data.n_cols, p, r).t();

  }

  for(arma::uword i = 0; i < clust.n_elem; i++){

    arma::mat curve_mat = integrated_curves_mat(dt, orders.row(0).t(), a0, b0, gamma, rho(i), M, S0, R0);
    llik(i) = log_sum_exp(curve_mat.cols(0,data.n_cols-1) * data.row(i).t() - curve_mat.col(data.n_cols)) - log(M);

  }


  arma::mat res_clust(niter - nburn, data.n_rows);
  arma::cube res_orders(data.n_rows, data.n_cols, niter - nburn);
  arma::mat res_llik(niter - nburn, data.n_rows);
  arma::mat res_rho(niter - nburn, data.n_rows);

  //loop
  int start_s = clock();
  int current_s;
  if(nupd == 0){
    nupd = round(niter / 10);
  }

  Rcpp::Rcout << "\n------ MAIN LOOP ------\n\n";
  // start
  for(int iter = 0; iter < niter; iter++){

    update_rho(data, rho, a0, b0, c0, d0, MH_var, gamma, dt, M,
               S0, R0, llik, clust, orders);

    for(arma::uword j = 0; j < orders.n_rows; j++){
      update_single_order(data, clust, j, orders, llik,
                          q, dt, a0, b0, gamma, rho, M, S0, R0, coars);
    }

    if(iter >= nburn){
      res_clust.row(iter-nburn) = clust.t();
      res_orders.slice(iter-nburn) = orders;
      res_llik.row(iter-nburn).cols(0,data.n_rows-1) = llik.t();
      res_rho.row(iter-nburn).cols(0,data.n_rows-1) = rho.t();
    }
    // print time
    if((iter + 1) % nupd == 0){
      current_s = clock();
      Rcpp::Rcout << "Completed:\t" << (iter + 1) << "/" << niter << " - in " <<
        double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
    }
    Rcpp::checkUserInterrupt();
  }
  //double time = double(current_s-start_s)/CLOCKS_PER_SEC;

  Rcpp::List results;
  results["clust"] = res_clust;
  results["orders"] = res_orders;
  results["llik"] = res_llik;
  results["rho"] = res_rho;
  //results["time"] = time;
  return results;
}

//-----------------------------------------
// GENERATE SYNTHETIC TIMES OF INFECTIONS
//-----------------------------------------

Rcpp::List DoobGillespieAlg(double S0,
                            double I0,
                            double MaxTime,
                            arma::vec beta_vec,
                            double gamma_0,
                            gsl_rng *r,
                            double R0 = 0){

  arma::vec infection_time;
  arma::vec infection_01;

  double St = S0;
  double It = I0;
  double Rt = R0;
  double t = 1;
  double t_star;
  int flag;

  double beta_t, gamma_t = gamma_0;

  while((t < MaxTime) & (It > 0) & (St > 0)){

    beta_t = beta_vec(round(t)-1);

    double E1 = gsl_ran_exponential(r, (S0 /(beta_t * St * It)));
    double E2 = gsl_ran_exponential(r, 1/(gamma_t * It));

    if(E1 < E2){     // Infection event
      t_star = E1;
      St = St - 1;
      It = It + 1;
      flag = 1;
    } else {	     // Recovery event
      t_star = E2;
      It = It - 1;
      Rt = Rt + 1;
      flag = 0;
    }

    t += t_star;

    infection_time.resize(infection_time.n_elem + 1);
    infection_time(infection_time.n_elem - 1) = t;

    infection_01.resize(infection_01.n_elem + 1);
    infection_01(infection_01.n_elem - 1) = flag;

  }

  Rcpp::List results;
  results["TimeInfections"] = infection_time;
  results["FlagInfections"] = infection_01;

  return results;

}

//' Simulate epidemiological data
//'
//' @param S0 number of individuals in the population.
//' @param I0 number of infected individuals at time 0.
//' @param max_time maximum observed time.
//' @param beta_vec vector with the infection rate for each discrete time.
//' @param gamma_0 the recovery rate. for the population, must be in \eqn{(0,1)}.
//' @param user_seed seed for random distribution generation.
//' @return Function \code{sim_epi_data} returns a vector with the simulated infection times.
//'
//' @export
// [[Rcpp::export]]
arma::vec sim_epi_data(double S0,
                       double I0,
                       double max_time,
                       arma::vec beta_vec,
                       double gamma_0,
                       unsigned long user_seed = 1234){

// set seed for gsl random distribution generator
const gsl_rng_type * T;
gsl_rng *r;
gsl_rng_env_setup();
T = gsl_rng_default; // Generator setup
r = gsl_rng_alloc (T);
gsl_rng_set(r, user_seed);
//

// WARNINGS //
if(S0 < 1){
  Rcpp::stop("'S0' must be at least equal to 1.");
}

if(I0 > S0){
  Rcpp::stop("'I0' must be smaller than 'S0'.");
}

if((gamma_0 > 1) | (gamma_0 < 0)){
  Rcpp::stop("'gamma_0' must be in (0,1).");
}

if(beta_vec.n_elem != max_time){
  Rcpp::stop("number of elements in 'beta_vec' must be equal to 'max_time'.");
}

if(beta_vec.n_elem  != max_time){
  Rcpp::stop("number of elements in 'beta_vec' must be equal to 'max_time'.");
}


// ------- //

Rcpp::List list_simtimes = DoobGillespieAlg(S0,I0,max_time,beta_vec,gamma_0,r);

arma::vec list_times = list_simtimes[0];
arma::vec list_flags = list_simtimes[1];
arma::vec infection_times = list_times(arma::find(list_flags == 1));

return infection_times;

}

//---------------------
// PARTITION ESTIMATE
//---------------------

//' Compute the posterior similarity matrix
//'
//' @param M A matrix where each row corresponds to the output cluster of the corresponding iteration.
//' @return Function \code{psm} returns an \eqn{n}\eqn{\times}\eqn{n} posterior similarity matrix.
//'
//[[Rcpp::export]]
arma::mat psm(arma::mat M){
  // initialize results
  arma::mat result(M.n_cols, M.n_cols, arma::fill::zeros);

  for(arma::uword i = 0; i < M.n_cols; i++){
    for(arma::uword j = 0; j <= i; j++){
      result(i,j) = sum(M.col(i) == M.col(j));

      result(j,i) = result(i,j);
    }
    Rcpp::checkUserInterrupt();
  }
  return(result / M.n_rows);
}

//' Estimate order
//'
//' @param orders_mat A matrix where each row corresponds to the output cluster of the corresponding iteration.
//' @return Function \code{get_clust_VI} returns a point estimate for the clustering of the data.
//'
//' @export
//[[Rcpp::export]]
arma::vec get_clust_VI(arma::mat orders_mat){

  arma::vec out_res(orders_mat.n_cols);
  arma::vec result(orders_mat.n_rows);
  arma::mat psm_mat = psm(orders_mat);

  double f = 0.0;
  int n = psm_mat.n_cols;
  arma::vec tvec(n);

  for(arma::uword j = 0; j < orders_mat.n_rows; j++){
    f = 0.0;
    for(int i = 0; i < n; i++){
      tvec = psm_mat.col(i);
      f += (log2(sum(orders_mat.row(j) == orders_mat(j,i))) +
        log2(sum(tvec)) -
        2 * log2(sum(tvec.elem(find(orders_mat.row(j).t() == orders_mat(j,i))))))/n;
    }
    result(j) = f;
    Rcpp::checkUserInterrupt();
  }


  int index_min = which_min_cpp(result);

  out_res = orders_mat.row(index_min).t();

  return(out_res);
}

//-----------------
// MAIN FUNCTIONS
//-----------------

//' Detect Change Points on an univariate time series.
//'
//' @param data vector of observations.
//' @param n_iterations number of MCMC iteration.
//' @param q probability of performing a split at each iterations.
//' @param phi parameter \eqn{\phi} of the integrated likelihood function.
//' @param a,b,c parameters of the Normal-Gamma prior for \eqn{\mu} and \eqn{\lambda}.
//' @param par_theta_c,par_theta_d parameters of the shifted Gamma prior for \eqn{\theta}.
//' @param print_progress If TRUE (default) print the progress bar.
//' @param user_seed seed for random distribution generation.
//' @return Function \code{detect_cp_uni} returns a list containing the following components: \itemize{
//' \item{\code{$orders}} a matrix where each row corresponds to the output order of the corresponding iteration.
//' \item{\code{time}} computational time in seconds.
//' \item{\code{$sigma_MCMC}} traceplot for \eqn{\sigma}.
//' \item{\code{$sigma_MCMC_01}} a \eqn{0/1} vector, the \eqn{n}-th element is equal to \eqn{1} if the proposed \eqn{\sigma} was accepted, \eqn{0} otherwise.
//' \item{\code{$theta_MCMC}} traceplot for \eqn{\theta}.
//' }
//'
//' @examples
//'
//' data_vec <- as.numeric(c(rnorm(50,0,0.1), rnorm(50,1,0.25)))
//'
//' out <- detect_cp_uni(data = data_vec,
//'                             n_iterations = 2500,
//'                             q = 0.25,
//'                             phi = 0.1, a = 1, b = 1, c = 0.1)
//'
//' get_clust_VI(out$order)
//'
//' @export
//[[Rcpp::export]]
Rcpp::List detect_cp_uni(arma::vec data,
                                int n_iterations, double q, double phi, double a, double b, double c,
                                double par_theta_c = 1, double par_theta_d = 1, bool print_progress = true, unsigned long user_seed = 1234){


  // WARNINGS //
  if(n_iterations < 1){
    Rcpp::stop("number of iterations must be at least 1.");
  }

  if((q > 1) | (q < 0)){
    Rcpp::stop("'q' must be included in (0,1).");
  }

  if((phi > 1) | (phi < 0)){
    Rcpp::stop("'phi' must be included in (0,1).");
  }

  if(a < 0){
    Rcpp::stop("'a' must be positive.");
  }

  if(b < 0){
    Rcpp::stop("'b' must be positive.");
  }

  if(c < 0){
    Rcpp::stop("'c' must be positive.");
  }

  if(par_theta_c < 0){
    Rcpp::stop("'par_theta_c' must be positive.");
  }

  if(par_theta_d < 0){
    Rcpp::stop("'par_theta_d' must be positive.");
  }
  // ------- //

 int start_s = clock();
 int current_s = start_s;
 int nupd = round(n_iterations / 10);

 // set seed for gsl random distribution generator
 const gsl_rng_type * T;
 gsl_rng *r;
 gsl_rng_env_setup();
 T = gsl_rng_default; // Generator setup
 r = gsl_rng_alloc (T);
 gsl_rng_set(r, user_seed);
 //

 arma::vec res_order, probs(2), sigma_inf(1), sigma_inf_10(1), theta_inf(1);
 arma::mat data_mat(1,data.n_elem), res_mat(n_iterations, data_mat.n_cols);

 sigma_inf(0) = 0.1; // add a random initialisation for the gamma param
 sigma_inf_10(0) = 0;
 theta_inf(0) = 0.1;

 for(arma::uword i = 0; i < data.n_elem; i++){
   data_mat.row(0).col(i) = data(i);
 }

 //generate random starting order
 arma::vec order = generate_random_order(data.n_elem, 2/data.n_elem, r);

 for(int iter = 0; iter < n_iterations; iter++){

   double k = max(order) + 1;

   probs(0) = q * indicator_1(k,data_mat.n_cols) + indicator_3(k);
   probs(1) = (1-q) * indicator_1(k,data_mat.n_cols) + indicator_2(k, data_mat.n_cols);

   probs(0) = probs(0)/(probs(0)+probs(1));
   probs(1) = probs(1)/(probs(0)+probs(1));

   double u = arma::randu();

   if(u <= probs(0)){

     // SPLIT

     /// propose a new order

     Rcpp::List split_list = Split_cpp(order);

     arma::vec split_order = split_list[1];
     int split_index = split_list[0];

     /// evaluate the proposed order

     double alpha_split = AlphaSplitOrder_UniTS(data_mat, split_order, order, q, split_index, theta_inf(iter), sigma_inf(iter), phi, a, b, c);


     if(log(arma::randu()) <= alpha_split){
       res_order = split_order;
     } else {
       res_order = order;
     }

   } else {

     // MERGE

     /// propose a new order

     Rcpp::List merge_list = Merge_cpp(order);

     arma::vec merge_order = merge_list[1];
     int merge_index = merge_list[0];

     arma::vec freq_temp = table_cpp(order);

     /// evaluate the proposed order

     double alpha_merge = AlphaMergeOrder_UniTS(data_mat, merge_order, order, q, merge_index, theta_inf(iter), sigma_inf(iter), phi, a, b, c);

     if(log(arma::randu()) <= alpha_merge){
       res_order = merge_order;
     } else {
       res_order = order;
     }

   }

   if(max(res_order) > 0){

     // SHUFFLE

     /// propose a new order

     Rcpp::List shuffle_list = Shuffle_cpp(res_order);

     arma::vec shuffle_order = shuffle_list[1];

     /// evaluate the proposed order

     double alpha_shuffle = AlphaShuffleOrder_UniTS(data_mat, shuffle_order, order, theta_inf(iter), sigma_inf(iter), phi, a, b, c);

     if(log(arma::randu()) <= alpha_shuffle){
       res_order = shuffle_order;
     } else {
       res_order = order;
     }
   }

   order = res_order;

   // Posterior infererence on main parameters
   UpdateSigma(order, theta_inf(iter), sigma_inf(iter), sigma_inf, sigma_inf_10, r);

   UpdateTheta(theta_inf(iter), sigma_inf(iter+1), order, theta_inf, par_theta_c, par_theta_d, r);
   //

   res_mat.row(iter) = order.t();

   if(((iter + 1) % nupd == 0) & (print_progress == true)){
     current_s = clock();
     Rcpp::Rcout << "Completed:\t" << (iter + 1) << "/" << n_iterations << " - in " <<
       double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
   }
   Rcpp::checkUserInterrupt();

 }

 double time = double(current_s-start_s)/CLOCKS_PER_SEC;

 Rcpp::List out_list;
 out_list["orders"] = res_mat;
 out_list["time"] = time;
 out_list["sigma_MCMC"] = sigma_inf;
 out_list["sigma_MCMC_01"] = sigma_inf_10;
 out_list["theta_MCMC"] = theta_inf;

 gsl_rng_free (r);

 return(out_list);

}


//' Detect Change Points on multivariate time series
//'
//' @param data a matrix where each row is a component of the time series and the columns correpospond to the times.
//' @param n_iterations number of MCMC iterations.
//' @param q probability of performing a split at each iteration.
//' @param k_0,nu_0,phi_0,m_0 parameters for the Normal-Inverse-Wishart prior for \eqn{(\mu,\lambda)}.
//' @param par_theta_c,par_theta_d parameters for the shifted Gamma prior for \eqn{\theta}.
//' @param prior_var_gamma parameters for the Gamma prior for \eqn{\gamma}.
//' @param print_progress If TRUE (default) print the progress bar.
//' @param user_seed seed for random distribution generation.
//' @return Function \code{detect_cp_multi} returns a list containing the following components: \itemize{
//' \item{\code{$orders}} a matrix where each row corresponds to the output order of the corresponding iteration.
//' \item{\code{time}} computational time in seconds.
//' \item{\code{$gamma_MCMC}} traceplot for \eqn{\gamma}.
//' \item{\code{$gamma_MCMC_01}} a \eqn{0/1} vector, the \eqn{n}-th element is equal to \eqn{1} if the proposed \eqn{\gamma} was accepted, \eqn{0} otherwise.
//' \item{\code{$sigma_MCMC}} traceplot for \eqn{\sigma}.
//' \item{\code{$sigma_MCMC_01}} a \eqn{0/1} vector, the \eqn{n}-th element is equal to \eqn{1} if the proposed \eqn{\sigma} was accepted, \eqn{0} otherwise.
//' \item{\code{$theta_MCMC}} traceplot for \eqn{\theta}.
//' }
//'
//'@examples
//'
//' data_mat <- matrix(NA, nrow = 3, ncol = 100)
//'
//' data_mat[1,] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
//' data_mat[2,] <- as.numeric(c(rnorm(50,0,0.125), rnorm(50,1,0.225)))
//' data_mat[3,] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))
//'
//' out <- detect_cp_multi(data = data_mat,
//'                               n_iterations = 2500,
//'                               q = 0.25,k_0 = 0.25, nu_0 = 4, phi_0 = diag(1,3,3), m_0 = rep(0,3),
//'                               par_theta_c = 2, par_theta_d = 0.2, prior_var_gamma = 0.1)
//'
//' get_clust_VI(out$order)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List detect_cp_multi(arma::mat data,
                                  int n_iterations, double q, double k_0, double nu_0,
                                  arma::mat phi_0, arma::vec m_0,
                                  double par_theta_c = 1, double par_theta_d = 1, double prior_var_gamma = 0.1,
                                  bool print_progress = true, unsigned long user_seed = 1234){

    // WARNINGS //
    if(n_iterations < 1){
      Rcpp::stop("number of iterations must be at least 1.");
    }

    if((q > 1) | (q < 0)){
      Rcpp::stop("'q' must be included in (0,1).");
    }

    if(k_0 < 0){
      Rcpp::stop("'k_0' must be positive.");
    }

    if(nu_0 < 0){
      Rcpp::stop("'nu_0' must be positive.");
    }

    if(phi_0.n_rows != data.n_rows){
      Rcpp::stop("number of rows in 'phi_0' must equal number of observations.");
    }

    if(phi_0.n_cols != data.n_rows){
      Rcpp::stop("number of columns in 'phi_0' must equal number of observations.");
    }

    if(m_0.n_elem != data.n_rows){
      Rcpp::stop("number of elements in 'm_0' must equal number of observations.");
    }

    if(par_theta_c < 0){
      Rcpp::stop("'par_theta_c' must be positive.");
    }

    if(par_theta_d < 0){
      Rcpp::stop("'par_theta_d' must be positive.");
    }

    if(prior_var_gamma < 0){
      Rcpp::stop("'prior_var_gamma' must be positive.");
    }
    // ------- //

   // set seed for gsl random distribution generator
   const gsl_rng_type * T;
   gsl_rng *r;
   gsl_rng_env_setup();
   T = gsl_rng_default; // Generator setup
   r = gsl_rng_alloc (T);
   gsl_rng_set(r, user_seed);
   //

   arma::vec res_order, probs(2), gamma_inf(1), gamma_inf_10(1), sigma_inf(1), sigma_inf_10(1), theta_inf(1);
   arma::mat res_mat(n_iterations, data.n_cols);

   gamma_inf(0) = 0.5; // add a random initialisation for the gamma param
   gamma_inf_10(0) = 0;
   sigma_inf(0) = 0.1;
   sigma_inf_10(0) = 0;
   theta_inf(0) = 0.1;

   //generate random starting order
   arma::vec order = generate_random_order(data.n_cols, 2/data.n_cols, r);

   int start_s = clock();
   int current_s = start_s;
   int nupd = round(n_iterations / 10);

   for(int iter = 0; iter < n_iterations; iter++){

     double k = max(order) + 1;

     probs(0) = q * indicator_1(k,data.n_cols) + indicator_3(k);
     probs(1) = (1-q) * indicator_1(k,data.n_cols) + indicator_2(k, data.n_cols);

     probs(0) = probs(0)/(probs(0)+probs(1));
     probs(1) = probs(1)/(probs(0)+probs(1));

     double u = arma::randu();

     if(u <= probs(0)){

       // SPLIT

       /// propose a new order

       Rcpp::List split_list = Split_cpp(order);

       arma::vec split_order = split_list[1];
       int split_index = split_list[0];

       /// evaluate the proposed order

       double alpha_split = AlphaSplitOrder_MultiTS(data, split_order, order, q, split_index, gamma_inf(iter), k_0, nu_0, theta_inf(iter), sigma_inf(iter), phi_0, m_0);

       if(log(arma::randu()) <= alpha_split){
         res_order = split_order;
       } else {
         res_order = order;
       }

     } else {

       // MERGE

       /// propose a new order

       Rcpp::List merge_list = Merge_cpp(order);

       arma::vec merge_order = merge_list[1];
       int merge_index = merge_list[0];

       arma::vec freq_temp = table_cpp(order);

       /// evaluate the proposed order

       double alpha_merge = AlphaMergeOrder_MultiTS(data, merge_order, order, q, merge_index, gamma_inf(iter), k_0, nu_0, theta_inf(iter), sigma_inf(iter), phi_0, m_0);

       if(log(arma::randu()) <= alpha_merge){
         res_order = merge_order;
       } else {
         res_order = order;
       }

     }

     if(max(res_order) > 0){

       // SHUFFLE

       /// propose a new order

       Rcpp::List shuffle_list = Shuffle_cpp(res_order);

       arma::vec shuffle_order = shuffle_list[1];

       /// evaluate the proposed order

       double alpha_shuffle = AlphaShuffleOrder_MultiTS(data, shuffle_order, order, gamma_inf(iter), k_0, nu_0, theta_inf(iter), sigma_inf(iter), phi_0, m_0);

       if(log(arma::randu()) <= alpha_shuffle){
         res_order = shuffle_order;
       } else {
         res_order = order;
       }
     }

     order = res_order;

     // Posterior infererence on main parameters
     UpdateGamma(gamma_inf(iter), data, order, theta_inf(iter), sigma_inf(iter), k_0, nu_0,
                 phi_0, m_0, gamma_inf, gamma_inf_10, r, prior_var_gamma);

     UpdateSigma(order, theta_inf(iter), sigma_inf(iter), sigma_inf, sigma_inf_10, r);

     UpdateTheta(theta_inf(iter), sigma_inf(iter+1), order, theta_inf, par_theta_c, par_theta_d, r);
     //

     res_mat.row(iter) = order.t();

     if(((iter + 1) % nupd == 0) & (print_progress == true)){
       current_s = clock();
       Rcpp::Rcout << "Completed:\t" << (iter + 1) << "/" << n_iterations << " - in " <<
         double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
     }
     Rcpp::checkUserInterrupt();

   }

   double time = double(current_s-start_s)/CLOCKS_PER_SEC;

   Rcpp::List out_list;
   out_list["orders"] = res_mat;
   out_list["time"] = time;
   out_list["gamma_MCMC"] = gamma_inf;
   out_list["gamma_MCMC_01"] = gamma_inf_10;
   out_list["sigma_MCMC"] = sigma_inf;
   out_list["sigma_MCMC_01"] = sigma_inf_10;
   out_list["theta_MCMC"] = theta_inf;

   gsl_rng_free (r);

   return(out_list);

 }



//' Clustering Epidemiological survival functions with common changes in time
//'
//' @param data a matrix where each entry is the number of infected for a population (row) at a specific discrete time (column).
//' @param n_iterations Second value
//' @param M number of Monte Carlo iterations when computing the likelihood of the survival function.
//' @param B number of orders for the normalisation constant.
//' @param L number of split-merge steps for the proposal step.
//' @param gamma recovery rate fixed constant for each population at each time.
//' @param alpha \eqn{\alpha} for the acceptance ration in the split-merge procedure.
//' @param q probability of performing a split when updating the single order for the proposal procedure.
//' @param dt,a0,b0,c0,d0 parameters for the computation of the integrated likelihood of the survival functions.
//' @param MH_var variance for the Metropolis-Hastings estimation of the proportion of infected at time 0.
//' @param S0,R0 parameters for the SDE solver.
//' @param p prior average number of change points for each order.
//' @param coars coarsening parameter.
//' @param print_progress If TRUE (default) print the progress bar.
//' @param user_seed seed for random distribution generation.
//' @return Function \code{clust_cp_epi} returns a list containing the following components: \itemize{
//' \item{\code{$clust}} a matrix where each row corresponds to the output cluster of the corresponding iteration.
//' \item{\code{$orders}} a multidimensional matrix where each slice is a matrix with the orders associated to the output cluster of that iteration.
//' \item{\code{time}} computational time in seconds.
//' \item{\code{$llik}} a matrix containing the log-likelihood of each population at each iteration.
//' \item{\code{$rho}} traceplot for the proportion of infected individuals at time 0.
//' }
//'
//'@examples
//'\donttest{
//' data_mat <- matrix(NA, nrow = 5, ncol = 50)
//'
//' betas <- list(c(rep(0.45, 25),rep(0.14,25)),
//'               c(rep(0.55, 25),rep(0.11,25)),
//'               c(rep(0.50, 25),rep(0.12,25)),
//'               c(rep(0.52, 10),rep(0.15,40)),
//'               c(rep(0.53, 10),rep(0.13,40)))
//'
//'  inf_times <- list()
//'
//'  for(i in 1:5){
//'
//'    inf_times[[i]] <- sim_epi_data(10000, 10, 50, betas[[i]], 1/8)
//'
//'    vec <- rep(0,50)
//'    names(vec) <- as.character(1:50)
//'
//'    for(j in 1:50){
//'      if(as.character(j) %in% names(table(floor(inf_times[[i]])))){
//'        vec[j] = table(floor(inf_times[[i]]))[which(names(table(floor(inf_times[[i]]))) == j)]
//'      }
//'    }
//'    data_mat[i,] <- vec
//'  }
//'
//'  out <- clust_cp_epi(data = data_mat, n_iterations = 3000, M = 250, B = 1000, L = 1)
//'
//'  get_clust_VI(out$clust[1000:3000,])
//'}
//' @export
// [[Rcpp::export]]
Rcpp::List clust_cp_epi(arma::mat data,
                          int n_iterations,
                          int M,
                          int B,
                          int L,
                          double gamma = 1/8,
                          double alpha = 1,
                          double q = 0.1,
                          double dt = 0.1,
                          double a0 = 4,
                          double b0 = 10,
                          double c0 = 1,
                          double d0 = 1,
                          double MH_var = 0.01,
                          double S0 = 1,
                          double R0 = 0,
                          double p = 0.003,
                          double coars = 1,
                          bool print_progress = true,
                          unsigned long user_seed = 1234){

  // WARNINGS //
  if(n_iterations < 1){
    Rcpp::stop("number of iterations must be at least 1.");
  }

  if(M < 1){
    Rcpp::stop("'M' must be at least equal to 1.");
  }

  if(B < 1){
    Rcpp::stop("'B' must be at least equal to 1.");
  }

  if(L < 1){
    Rcpp::stop("'L' must be at least equal to 1.");
  }

  if((gamma < 0)| (gamma > 1)){
    Rcpp::stop("'gamma' must be in the interval (0,1).");
  }

  if(alpha < 0){
    Rcpp::stop("'alpha' must be positive.");
  }

  if((q < 0) | (q > 1)){
    Rcpp::stop("'q' must be in the interval (0,1).");
  }

  if((dt < 0) | (dt > 1)){
    Rcpp::stop("'dt' must be in the interval (0,1).");
  }

  if(a0 < 0){
    Rcpp::stop("'a0' must be positive.");
  }

  if(b0 < 0){
    Rcpp::stop("'b0' must be positive.");
  }

  if(c0 < 0){
    Rcpp::stop("'c0' must be positive.");
  }

  if(d0 < 0){
    Rcpp::stop("'d0' must be positive.");
  }

  if(MH_var < 0){
    Rcpp::stop("'MH_var' must be positive.");
  }

  if(S0 < 0){
    Rcpp::stop("'S0' must be positive.");
  }

  if(R0 < 0){
    Rcpp::stop("'R0' must be positive.");
  }

  if((p < 0) | (p > 1)){
    Rcpp::stop("'p' must be in the interval (0,1).");
  }

  if((coars < 0) | (coars > 1)){
    Rcpp::stop("'coars' must be in the interval (0,1).");
  }


  // ------- //


 // set seed for gsl random distribution generator
 const gsl_rng_type * T;
 gsl_rng *r;
 gsl_rng_env_setup();
 T = gsl_rng_default; // Generator setup
 r = gsl_rng_alloc (T);
 gsl_rng_set(r, user_seed);
 //

 arma::vec rho(data.n_rows);
 rho.fill(0.0045);

 arma::vec clust(data.n_rows), llik(data.n_rows);
 clust = arma::regspace(0, data.n_rows-1);
 arma::mat orders(data.n_rows, data.n_cols);
 orders.fill(0);

 arma::vec norm_vec = norm_constant_epi(data, gamma, B, a0, b0, rho, M, dt,
                                        r, print_progress);

 for(arma::uword i = 0; i < clust.n_elem; i++){
   arma::mat curve_mat = integrated_curves_mat(dt, orders.row(0).t(),
                                               a0, b0, gamma, rho(i), M, S0, R0);
   llik(i) = log_sum_exp(curve_mat.cols(0,data.n_cols-1) * data.row(i).t() - curve_mat.col(data.n_cols)) - log(M);
 }

 for(int l = 0; l < 1; l++){
   for(arma::uword j = 0; j < orders.n_rows; j++){
     update_single_order(data, clust, j, orders, llik,
                         q, dt, a0, b0, gamma, rho, M, S0, R0, coars);
   }
 }
 arma::mat res_clust(n_iterations, data.n_rows);
 arma::cube res_orders(data.n_rows, data.n_cols, n_iterations);
 arma::mat res_llik(n_iterations, data.n_rows);
 arma::mat res_rho(n_iterations, data.n_rows);

 //loop
 int start_s = clock();
 int current_s = start_s;
 int nupd = round(n_iterations / 10);

 if(print_progress == true){
   Rcpp::Rcout << "\n------ MAIN LOOP ------\n\n";
 }

 // start
 for(int iter = 0; iter < n_iterations; iter++){

   update_rho(data, rho, a0, b0, c0, d0, MH_var, gamma, dt, M,
              S0, R0, llik, clust, orders);


   update_partition(data, clust, orders, llik, norm_vec, alpha, p, q,
                    dt, a0, b0, gamma, rho, M, L, S0, R0,coars);


   for(arma::uword j = 0; j < orders.n_rows; j++){
     update_single_order(data, clust, j, orders, llik,
                         q, dt, a0, b0, gamma, rho, M, S0, R0, coars);
   }


   res_clust.row(iter) = clust.t();
   res_llik.row(iter).cols(0,data.n_rows-1) = llik.t();
   res_rho.row(iter).cols(0,data.n_rows-1) = rho.t();
   res_orders.slice(iter).rows(0, orders.n_rows-1) = orders;

   // print time
   if(((iter + 1) % nupd == 0) & (print_progress == true)){
     current_s = clock();
     Rcpp::Rcout << "Completed:\t" << (iter + 1) << "/" << n_iterations << " - in " <<
       double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
   }
   Rcpp::checkUserInterrupt();



 }

 double time = double(current_s-start_s)/CLOCKS_PER_SEC;

 Rcpp::List results;
 results["clust"] = res_clust;
 results["orders"] = res_orders;
 results["time"] = time;
 results["llik"] = res_llik;
 results["rho"] = res_rho;
 return results;
}



//' Clustering univariate times series with common changes in time
//'
//' @param data a matrix where each row is an observation and each column corresponds to a discrete time.
//' @param n_iterations number of MCMC iterations.
//' @param B number of orders for the normalisation constant.
//' @param L number of split-merge steps for the proposal step.
//' @param gamma,a,b,c parameters \eqn{\gamma} of the integrated likelihood.
//' @param q probability of a split in the split-merge proposal and acceleration step.
//' @param alpha_SM \eqn{\alpha} for the split-merge proposal and acceleration step.
//' @param coars coarsening coefficient, must be in (0,1].
//' @param print_progress If TRUE (default) print the progress bar.
//' @param user_seed seed for random distribution generation.
//' @return Function \code{clust_cp_uni} returns a list containing the following components: \itemize{
//' \item{\code{$clust}} a matrix where each row corresponds to the output cluster of the corresponding iteration.
//' \item{\code{$orders}} a multidimensional array where each slice is a matrix and represent an iteration. The row of each matrix correspond the order associated to the corresponding cluster.
//' \item{\code{time}} computational time in seconds.
//' \item{\code{$norm_vec}} a vector containing the normalisation constant computed at the beginning of the algorithm.
//' }
//'
//' @examples
//'
//' data_mat <- matrix(NA, nrow = 5, ncol = 100)
//'
//' data_mat[1,] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
//' data_mat[2,] <- as.numeric(c(rnorm(50,0,0.125), rnorm(50,1,0.225)))
//' data_mat[3,] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))
//' data_mat[4,] <- as.numeric(c(rnorm(25,0,0.135), rnorm(75,1,0.225)))
//' data_mat[5,] <- as.numeric(c(rnorm(25,0,0.155), rnorm(75,1,0.280)))
//'
//' out <- clust_cp_uni(data = data_mat, n_iterations = 5000, B = 1000, L = 1, gamma = 0.5)
//'
//' get_clust_VI(out$clust[2500:5000,])
//'
//' @export
// [[Rcpp::export]]
Rcpp::List clust_cp_uni(arma::mat data,
                          int n_iterations,
                          int B,
                          int L,
                          double gamma,
                          double a = 1,
                          double b = 1,
                          double c = 1,
                          double q = 0.5,
                          double alpha_SM = 0.1, double coars = 1,
                          bool print_progress = true,
                          unsigned long user_seed = 1234){

// WARNINGS //
if(n_iterations < 1){
  Rcpp::stop("number of iterations must be at least 1.");
}

if(B < 1){
  Rcpp::stop("'B' must be at least equal to 1.");
}

if(L < 1){
  Rcpp::stop("'L' must be at least equal to 1.");
}

if((gamma > 1) | (gamma < 0)){
  Rcpp::stop("'gamma' must be included in (0,1).");
}

if(a < 0){
  Rcpp::stop("'a' must be positive.");
}

if(b < 0){
  Rcpp::stop("'b' must be positive.");
}

if(c < 0){
  Rcpp::stop("'c' must be positive.");
}


if((q > 1) | (q < 0)){
  Rcpp::stop("'q' must be included in (0,1).");
}

if(alpha_SM <= 0){
  Rcpp::stop("'alpha_SM' must be positive.");
}

if((coars <= 0) | (coars > 1)){
  Rcpp::stop("'coars' must be in (0,1].");
}

// ------- //

arma::mat res_clust(n_iterations, data.n_rows), res_lkl(n_iterations, data.n_rows), orders_temp_clean;
arma::cube res_orders(data.n_rows, data.n_cols, n_iterations);
arma::vec freq_temp(data.n_rows), prob_temp(data.n_rows), proposed_order_i(data.n_cols), proposed_order_j(data.n_cols), order_i(data.n_cols), order_j(data.n_cols),
proposed_order(data.n_cols), proposed_partition(data.n_rows), prob_temp_j(data.n_rows), prob_temp_i(data.n_rows), proposed_partition_clean(data.n_cols),
merge_i(data.n_rows), merge_j(data.n_rows), old_order, order_0, lkl_proposal_m, lkl_old_i_m, lkl_old_j_m, lkl_old_s, lkl_proposal_i_s,lkl_proposal_j_s;
int id1, id2, id3, id4;
double alpha;

// set seed for gsl random distribution generator
const gsl_rng_type * T;
gsl_rng *r;
gsl_rng_env_setup();
T = gsl_rng_default; // Generator setup
r = gsl_rng_alloc (T);
gsl_rng_set(r, user_seed);
//

// generate starting partition and starting orders
arma::vec partition_temp = generate_random_partition(data.n_rows, r);
arma::mat orders_temp(data.n_rows, data.n_cols);
for(int i = 0; i < (max(partition_temp) + 1); i++){
  orders_temp.row(i) = generate_random_order(data.n_cols, 2.0/data.n_cols, r).t();
}
// computing starting likelihood
arma::vec lkl_temp(data.n_rows);

for(arma::uword i = 0; i < data.n_rows; i++){

  lkl_temp(i) = LogLikelihood_TS(data.row(i), orders_temp.row(partition_temp(i)), gamma,a,b,c);

}

// COMPUTE NORMALISATION CONSTANT

arma::vec norm_const = norm_constant_uni(data, gamma, B, a, b, c, 2.0/data.n_cols,print_progress);

if(print_progress == true){
  Rcpp::Rcout << "\n------ MAIN LOOP ------\n\n";
}

// MAIN LOOP

int start_s = clock();
int current_s = start_s;
int nupd = round(n_iterations / 10);

for(int iter = 0; iter < n_iterations; iter++){

  //define k and n
  int k = max(partition_temp) + 1;
  int n = data.n_rows;
  //

  // select two random obs
  freq_temp.fill(1.0);
  id1 = rint(freq_temp);
  freq_temp(id1) = 0;
  id2 = rint(freq_temp);
  //


  if(partition_temp(id1) != partition_temp(id2)){
    // MERGE

    proposed_partition = partition_temp;

    merge_i.fill(0.0);
    merge_j.fill(0.0);

    merge_i.elem(find(partition_temp == partition_temp(id1))).fill(1.0);
    merge_j.elem(find(partition_temp == partition_temp(id2))).fill(1.0);

    proposed_partition.elem(find(proposed_partition == proposed_partition(id2))).fill(proposed_partition(id1));

    prob_temp.fill(0.0);
    for(int i = 0; i < n; i++){
      if(partition_temp(i) == partition_temp(id1) || partition_temp(i) == partition_temp(id2)){
        prob_temp(i) = 1.0;
      }
    }

    order_i = orders_temp.row(partition_temp(id1)).t();
    order_j = orders_temp.row(partition_temp(id2)).t();

    id3 = rint(prob_temp);

    proposed_order = orders_temp.row(partition_temp(id3)).t();

    // update the proposed order with a split-and-merge procedure
    SplitMergeAccUnivariateTS(data.row(id3), proposed_order, L, q, a, b, c, gamma);
    //

    lkl_proposal_m = lkl_temp;
    lkl_old_i_m = lkl_temp;
    lkl_old_j_m = lkl_temp;

    for(int i = 0; i < n; i++){
      lkl_old_i_m(i) = LogLikelihood_TS(data.row(i), order_i.t(), gamma,a,b,c);
      lkl_old_j_m(i) = LogLikelihood_TS(data.row(i), order_j.t(), gamma,a,b,c);
      lkl_proposal_m(i) = LogLikelihood_TS(data.row(i), proposed_order.t(), gamma,a,b,c);
    }


    // evaluate the proposed new partition

    alpha = AlphaMergePartition_cpp(lkl_old_i_m,
                                    lkl_old_j_m,
                                    lkl_proposal_m,
                                    k, data.n_cols, n, alpha_SM,
                                    merge_i,
                                    merge_j,
                                    norm_const, coars);

    if(log(arma::randu()) <= alpha){

      proposed_partition_clean = clean_partition_cpp(proposed_partition);

      // fixing orders' matrix

      orders_temp_clean = orders_temp;

      for(arma::uword i = 0; i < orders_temp.n_rows; i++){

        if(prob_temp(i) == 1.0){
          orders_temp_clean.row(proposed_partition_clean(i)) = proposed_order.t();
        } else {
          orders_temp_clean.row(proposed_partition_clean(i)) = orders_temp.row(proposed_partition(i));
        }

      }

      //

      orders_temp = orders_temp_clean;
      partition_temp = proposed_partition_clean;

      for(int i = 0; i < n; i++){
        lkl_temp(i) = LogLikelihood_TS(data.row(i), orders_temp.row(partition_temp(i)).t(),gamma,a,b,c);
      }

    }

    //


  } else {

    // SPLIT

    proposed_partition = partition_temp;

    prob_temp.fill(0.0);
    for(arma::uword i = 0; i < proposed_partition.n_elem; i++){
      if(proposed_partition(i) == proposed_partition(id1)){
        prob_temp(i) = 1.0;
      }
    }

    proposed_partition(id2) = k;

    prob_temp_i.fill(0.0);
    prob_temp_i(id1) = 1.0;
    prob_temp_j.fill(0.0);
    prob_temp_j(id2) = 1.0;

    for(arma::uword i = 0; i < proposed_partition.n_elem; i++){
      if(prob_temp(i) == 1){
        if(arma::randu() < 0.5){
          proposed_partition(i) = proposed_partition(id1);
          prob_temp_i(i) = 1.0;
        } else {
          proposed_partition(i) = proposed_partition(id2);
          prob_temp_j(i) = 1.0;
        }
      }
    }

    old_order = orders_temp.row(proposed_partition(id1)).t();
    proposed_order_i = orders_temp.row(proposed_partition(id1)).t();
    proposed_order_j = orders_temp.row(proposed_partition(id2)).t();

    // Update the proposed orders with a split-and-merge procedure

    id3 = rint(prob_temp_i);
    id4 = rint(prob_temp_j);

    SplitMergeAccUnivariateTS(data.row(id3), proposed_order_i, L, q, a, b, c, gamma);
    SplitMergeAccUnivariateTS(data.row(id4), proposed_order_j, L, q, a, b, c, gamma);

    //

    lkl_old_s = lkl_temp;
    lkl_proposal_i_s = lkl_temp;
    lkl_proposal_j_s = lkl_temp;

    for(int i = 0; i < n; i++){
      lkl_proposal_i_s(i) = LogLikelihood_TS(data.row(i), proposed_order_i.t(), gamma,a,b,c);
      lkl_proposal_j_s(i) = LogLikelihood_TS(data.row(i), proposed_order_j.t(), gamma,a,b,c);
      lkl_old_s(i) = LogLikelihood_TS(data.row(i), old_order.t(), gamma,a,b,c);
    }

    // evaluate the proposed new partition

    alpha = AlphaSplitPartition_cpp(lkl_proposal_i_s,
                                    lkl_proposal_j_s,
                                    lkl_old_s,
                                    k, data.n_cols, n, alpha_SM,
                                    prob_temp_i,
                                    prob_temp_j,
                                    norm_const, coars);

    if(log(arma::randu()) <= alpha){

      proposed_partition_clean = clean_partition_cpp(proposed_partition);

      // qua fixare la matrice degli ordini

      orders_temp_clean = orders_temp;

      for(arma::uword i = 0; i < orders_temp.n_rows; i++){

        if(prob_temp_i(i) == 1.0){
          orders_temp_clean.row(proposed_partition_clean(i)) = proposed_order_i.t();
        } else if (prob_temp_j(i) == 1.0){
          orders_temp_clean.row(proposed_partition_clean(i)) = proposed_order_j.t();
        } else {
          orders_temp_clean.row(proposed_partition_clean(i)) = orders_temp.row(proposed_partition(i));
        }

      }

      //

      partition_temp = proposed_partition_clean;
      orders_temp = orders_temp_clean;

      for(int i = 0; i < n; i++){
        lkl_temp(i) = LogLikelihood_TS(data.row(i), orders_temp.row(partition_temp(i)).t(),gamma,a,b,c);
      }

    }


  }

  // ACCELERATION STEP

  for(int i = 0; i <= max(partition_temp); i++){

    proposed_order = orders_temp.row(min(find(partition_temp == i))).t();

    arma::uvec obs = find(partition_temp == i);

    SplitMergeAccUnivariateTS(data.rows(obs), proposed_order, 1, q, a, b, c, gamma);

    orders_temp.row(min(find(partition_temp == i))) = proposed_order.t();

  }

  for(int i = 0; i < n; i++){
        lkl_temp(i) = LogLikelihood_TS(data.row(i), orders_temp.row(partition_temp(i)).t(),gamma,a,b,c);
      }

  res_clust.row(iter) = partition_temp.t();
  res_orders.slice(iter) = orders_temp;
  res_lkl.row(iter) = lkl_temp.t();

  // print time
  if(((iter + 1) % nupd == 0) & (print_progress == true)){
    current_s = clock();
    Rcpp::Rcout << "Completed:\t" << (iter + 1) << "/" << n_iterations << " - in " <<
      double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
  }
  Rcpp::checkUserInterrupt();

}

double time = double(current_s-start_s)/CLOCKS_PER_SEC;

Rcpp::List out_list;
out_list["clust"] = res_clust;
out_list["orders"] = res_orders;
out_list["time"] = time;
out_list["lkl"] = res_lkl;
out_list["norm_vec"] = norm_const;

return out_list;

}

//' Clustering multivariate times series with common changes in time
//'
//' @param data a multidimensional matrix where each element is a matrix whose rows are the observations and columns the dimensions.
//' @param n_iterations number of MCMC iterations.
//' @param B number of orders for the normalisation constant.
//' @param L number of split-merge steps for the proposal step.
//' @param gamma,k_0,nu_0,phi_0,m_0 parameters of the integrated likelihood.
//' @param q probability of a split in the split-merge proposal and acceleration step.
//' @param alpha_SM \eqn{\alpha} for the split-merge proposal and acceleration step.
//' @param coars coarsening coefficient, must be in (0,1].
//' @param print_progress If TRUE (default) print the progress bar.
//' @param user_seed seed for random distribution generation.
//' @return Function \code{clust_cp_multi} returns a list containing the following components: \itemize{
//' \item{\code{$clust}} a matrix where each row corresponds to the output cluster of the corresponding iteration.
//' \item{\code{$orders}} a multidimensional array where each slice is a matrix and represent an iteration. The row of each matrix correspond the order associated to the corresponding cluster.
//' \item{\code{time}} computational time in seconds.
//' \item{\code{$norm_vec}} a vector containing the normalisation constant computed at the beginning of the algorithm.
//' }
//'
//' @examples
//'
//' data_array <- array(data = NA, dim = c(3,100,5))
//'
//' data_array[1,,1] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
//' data_array[2,,1] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
//' data_array[3,,1] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
//'
//' data_array[1,,2] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
//' data_array[2,,2] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
//' data_array[3,,2] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
//'
//' data_array[1,,3] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))
//' data_array[2,,3] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))
//' data_array[3,,3] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))
//'
//' data_array[1,,4] <- as.numeric(c(rnorm(25,0,0.135), rnorm(75,1,0.225)))
//' data_array[2,,4] <- as.numeric(c(rnorm(25,0,0.135), rnorm(75,1,0.225)))
//' data_array[3,,4] <- as.numeric(c(rnorm(25,0,0.135), rnorm(75,1,0.225)))
//'
//' data_array[1,,5] <- as.numeric(c(rnorm(25,0,0.155), rnorm(75,1,0.280)))
//' data_array[2,,5] <- as.numeric(c(rnorm(25,0,0.155), rnorm(75,1,0.280)))
//' data_array[3,,5] <- as.numeric(c(rnorm(25,0,0.155), rnorm(75,1,0.280)))
//'
//' out <- clust_cp_multi(data = data_array, n_iterations = 3000, B = 1000, L = 1,
//'                         gamma = 0.1, k_0 = 0.25, nu_0 = 5, phi_0 = diag(0.1,3,3), m_0 = rep(0,3))
//'
//' get_clust_VI(out$clust[1000:3000,])
//'
//' @export
// [[Rcpp::export]]
Rcpp::List clust_cp_multi(arma::cube data,
                            int n_iterations,
                            int B,
                            int L,
                            double gamma,
                            double k_0,
                            double nu_0,
                            arma::mat phi_0,
                            arma::vec m_0,
                            double q = 0.5,
                            double alpha_SM = 0.1, double coars = 1, bool print_progress = true, unsigned long user_seed = 1234){

  arma::mat res_clust(n_iterations, data.n_slices), res_lkl(n_iterations,data.n_slices), orders_temp_clean;
  arma::cube res_orders(data.n_slices, data.slice(0).n_cols, n_iterations);
  arma::vec freq_temp(data.n_slices), prob_temp(data.n_slices), proposed_order_i(data.slice(0).n_cols), proposed_order_j(data.slice(0).n_cols), order_i(data.slice(0).n_cols), order_j(data.slice(0).n_cols),
  proposed_order(data.slice(0).n_cols), proposed_partition(data.n_slices), prob_temp_j(data.n_slices), prob_temp_i(data.n_slices), proposed_partition_clean(data.slice(0).n_cols),
  merge_i(data.n_slices), merge_j(data.n_slices), old_order, order_0, lkl_proposal_m, lkl_old_i_m, lkl_old_j_m, lkl_old_s, lkl_proposal_i_s,lkl_proposal_j_s;
  int id1, id2, id3, id4;
  double alpha;

  // WARNINGS //
  if(n_iterations < 1){
    Rcpp::stop("number of iterations must be at least 1.");
  }

  if(B < 1){
    Rcpp::stop("'B' must be at least equal to 1.");
  }

  if(L < 1){
    Rcpp::stop("'L' must be at least equal to 1.");
  }

  if((gamma > 1) | (gamma < 0)){
    Rcpp::stop("'gamma' must be included in (0,1).");
  }

  if(k_0 < 0){
    Rcpp::stop("'k_0' must be positive.");
  }

  if(nu_0 < 0){
    Rcpp::stop("'nu_0' must be positive.");
  }

  if(phi_0.n_rows != data.slice(0).n_rows){
    Rcpp::stop("number of rows and columns in 'phi_0' must correspond to the number of dimensions of the time series.");
  }

  if(m_0.n_elem != data.slice(0).n_rows){
    Rcpp::stop("number of elements in 'm_0' must correspond to the number of dimensions of the time series.");
  }

  if((q > 1) | (q < 0)){
    Rcpp::stop("'q' must be included in (0,1).");
  }

  if(alpha_SM <= 0){
    Rcpp::stop("'alpha_SM' must be positive.");
  }

  if((coars <= 0) | (coars > 1)){
    Rcpp::stop("'coars' must be in (0,1].");
  }

  // ------- //


  // set seed for gsl random distribution generator
  const gsl_rng_type * T;
  gsl_rng *r;
  gsl_rng_env_setup();
  T = gsl_rng_default; // Generator setup
  r = gsl_rng_alloc (T);
  gsl_rng_set(r, user_seed);
  //

  // generate starting partition and starting orders
  arma::vec partition_temp = generate_random_partition(data.n_slices, r);
  arma::mat orders_temp(data.n_slices, data.slice(0).n_cols);

  double num_groups_temp = max(partition_temp) + 1;

  for(int i = 0; i < num_groups_temp; i++){
    orders_temp.row(i) = generate_random_order(data.slice(0).n_cols, 2.0/data.slice(0).n_cols, r).t();
  }

  // computing starting likelihood
  arma::vec lkl_temp(data.n_slices);

  for(arma::uword i = 0; i < data.n_slices; i++){
    lkl_temp(i) = Likelihood_MultiTS(data.slice(i), orders_temp.row(partition_temp(i)).t(), gamma,k_0,nu_0,phi_0,m_0);
  }
  // COMPUTE NORMALISATION CONSTANT

  arma::vec norm_const = norm_constant_multi(data, gamma, B, k_0, nu_0, phi_0, m_0, 2.0/data.slice(0).n_cols,print_progress);

  if(print_progress == true){
    Rcpp::Rcout << "\n------ MAIN LOOP ------\n\n";
  }

  // MAIN LOOP

  int start_s = clock();
  int current_s = start_s;
  int nupd = round(n_iterations / 10);

  for(int iter = 0; iter < n_iterations; iter++){

    //define k
    //int k = max(partition_temp) + 1;
    double k = max(partition_temp) + 1;
    //

    // select two random obs
    freq_temp.fill(1.0);
    id1 = rint(freq_temp);
    freq_temp(id1) = 0;
    id2 = rint(freq_temp);
    //

    if(partition_temp(id1) != partition_temp(id2)){
      // MERGE

      proposed_partition = partition_temp;

      merge_i.fill(0.0);
      merge_j.fill(0.0);

      merge_i.elem(find(partition_temp == partition_temp(id1))).fill(1.0);
      merge_j.elem(find(partition_temp == partition_temp(id2))).fill(1.0);

      proposed_partition.elem(find(proposed_partition == proposed_partition(id2))).fill(proposed_partition(id1));

      prob_temp.fill(0.0);
      for(arma::uword i = 0; i < data.n_slices; i++){
        if(partition_temp(i) == partition_temp(id1) || partition_temp(i) == partition_temp(id2)){
          prob_temp(i) = 1.0;
        }
      }

      order_i = orders_temp.row(partition_temp(id1)).t();
      order_j = orders_temp.row(partition_temp(id2)).t();

      id3 = rint(prob_temp);

      proposed_order = orders_temp.row(partition_temp(id3)).t();

      // update the proposed order with a split-and-merge procedure
      SplitMergeAccMultivariateTS(data.slices(id3, id3), proposed_order, L, q, k_0, nu_0, phi_0, m_0, gamma);
      //

      lkl_proposal_m = lkl_temp;
      lkl_old_i_m = lkl_temp;
      lkl_old_j_m = lkl_temp;

      for(arma::uword i = 0; i < data.n_slices; i++){
        lkl_old_i_m(i) = Likelihood_MultiTS(data.slice(i), order_i, gamma,k_0,nu_0,phi_0,m_0);
        lkl_old_j_m(i) = Likelihood_MultiTS(data.slice(i), order_j, gamma,k_0,nu_0,phi_0,m_0);
        lkl_proposal_m(i) = Likelihood_MultiTS(data.slice(i), proposed_order, gamma,k_0,nu_0,phi_0,m_0);
      }

      // evaluate the proposed new partition

      alpha = AlphaMergePartition_cpp(lkl_old_i_m,
                                      lkl_old_j_m,
                                      lkl_proposal_m,
                                      k, data.n_cols, data.n_slices, alpha_SM,
                                      merge_i,
                                      merge_j,
                                      norm_const, coars);

      if(log(arma::randu()) <= alpha){

        proposed_partition_clean = clean_partition_cpp(proposed_partition);

        // fixing orders' matrix

        orders_temp_clean = orders_temp;

        for(arma::uword i = 0; i < orders_temp.n_rows; i++){

          if(prob_temp(i) == 1.0){
            orders_temp_clean.row(proposed_partition_clean(i)) = proposed_order.t();
          } else {
            orders_temp_clean.row(proposed_partition_clean(i)) = orders_temp.row(proposed_partition(i));
          }

        }

        //

        orders_temp = orders_temp_clean;
        partition_temp = proposed_partition_clean;

        for(arma::uword i = 0; i < data.n_slices; i++){
          lkl_temp(i) = Likelihood_MultiTS(data.slice(i), orders_temp.row(partition_temp(i)).t(), gamma,k_0,nu_0,phi_0,m_0);
        }

      }

      //


    } else {

      // SPLIT

      proposed_partition = partition_temp;

      prob_temp.fill(0.0);
      for(arma::uword i = 0; i < proposed_partition.n_elem; i++){
        if(proposed_partition(i) == proposed_partition(id1)){
          prob_temp(i) = 1.0;
        }
      }

      proposed_partition(id2) = k;

      prob_temp_i.fill(0.0);
      prob_temp_i(id1) = 1.0;
      prob_temp_j.fill(0.0);
      prob_temp_j(id2) = 1.0;

      for(arma::uword i = 0; i < proposed_partition.n_elem; i++){
        if(prob_temp(i) == 1){
          if(arma::randu() < 0.5){
            proposed_partition(i) = proposed_partition(id1);
            prob_temp_i(i) = 1.0;
          } else {
            proposed_partition(i) = proposed_partition(id2);
            prob_temp_j(i) = 1.0;
          }
        }
      }

      old_order = orders_temp.row(proposed_partition(id1)).t();
      proposed_order_i = orders_temp.row(proposed_partition(id1)).t();
      proposed_order_j = orders_temp.row(proposed_partition(id2)).t();

      // Update the proposed orders with a split-and-merge procedure

      id3 = rint(prob_temp_i);
      id4 = rint(prob_temp_j);

      SplitMergeAccMultivariateTS(data.slices(id3,id3), proposed_order_i, L, q, k_0, nu_0, phi_0, m_0, gamma);
      SplitMergeAccMultivariateTS(data.slices(id4,id4), proposed_order_j, L, q, k_0, nu_0, phi_0, m_0, gamma);

      //

      lkl_old_s = lkl_temp;
      lkl_proposal_i_s = lkl_temp;
      lkl_proposal_j_s = lkl_temp;

      for(arma::uword i = 0; i < data.n_slices; i++){
        lkl_proposal_i_s(i) = Likelihood_MultiTS(data.slice(i), proposed_order_i, gamma, k_0, nu_0, phi_0, m_0);
        lkl_proposal_j_s(i) = Likelihood_MultiTS(data.slice(i), proposed_order_j, gamma, k_0, nu_0, phi_0, m_0);
        lkl_old_s(i) = Likelihood_MultiTS(data.slice(i), old_order, gamma, k_0, nu_0, phi_0, m_0);
      }

      // evaluate the proposed new partition

      alpha = AlphaSplitPartition_cpp(lkl_proposal_i_s,
                                      lkl_proposal_j_s,
                                      lkl_old_s,
                                      k, data.n_cols, data.n_slices, alpha_SM,
                                      prob_temp_i,
                                      prob_temp_j,
                                      norm_const, coars);

      if(log(arma::randu()) <= alpha){

        proposed_partition_clean = clean_partition_cpp(proposed_partition);

        // qua fixare la matrice degli ordini

        orders_temp_clean = orders_temp;

        for(arma::uword i = 0; i < orders_temp.n_rows; i++){

          if(prob_temp_i(i) == 1.0){
            orders_temp_clean.row(proposed_partition_clean(i)) = proposed_order_i.t();
          } else if (prob_temp_j(i) == 1.0){
            orders_temp_clean.row(proposed_partition_clean(i)) = proposed_order_j.t();
          } else {
            orders_temp_clean.row(proposed_partition_clean(i)) = orders_temp.row(proposed_partition(i));
          }

        }

        //

        partition_temp = proposed_partition_clean;
        orders_temp = orders_temp_clean;

        for(arma::uword i = 0; i < data.n_slices; i++){
          lkl_temp(i) = Likelihood_MultiTS(data.slice(i), orders_temp.row(partition_temp(i)).t(), gamma, k_0, nu_0, phi_0, m_0);
        }

      }

    }

    // ACCELERATION STEP

    for(int i = 0; i <= max(partition_temp); i++){
      proposed_order = orders_temp.row(min(find(partition_temp == i))).t();
      arma::uvec obs = find(partition_temp == i);
      SplitMergeAccMultivariateTS(data.slices(obs), proposed_order, 1, q, k_0, nu_0, phi_0, m_0, gamma);
      orders_temp.row(min(find(partition_temp == i))) = proposed_order.t();
    }

    for(arma::uword i = 0; i < data.n_slices; i++){
          lkl_temp(i) = Likelihood_MultiTS(data.slice(i), orders_temp.row(partition_temp(i)).t(), gamma, k_0, nu_0, phi_0, m_0);
        }

    res_clust.row(iter) = partition_temp.t();
    res_orders.slice(iter) = orders_temp;
    res_lkl.row(iter) = lkl_temp.t();

    // print time
    if(((iter + 1) % nupd == 0) & (print_progress == true)){
      current_s = clock();
      Rcpp::Rcout << "Completed:\t" << (iter + 1) << "/" << n_iterations << " - in " <<
        double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
    }
    Rcpp::checkUserInterrupt();

  }

  double time = double(current_s-start_s)/CLOCKS_PER_SEC;



  Rcpp::List out_list;
  out_list["clust"] = res_clust;
  out_list["orders"] = res_orders;
  out_list["time"] = time;
  out_list["lkl"] = res_lkl;
  out_list["norm_vec"] = norm_const;

  return out_list;

}




#include <Rcpp.h>
#include <math.h>
#include <R.h>
#include <Rmath.h> 
#include <stdlib.h> 
#include <stdio.h>

using namespace Rcpp;

double awert(double x){

   /* This function provides an absolute value of a double x */
 
   if(x < 0)
     return(-x);
   else
    return(x);
}


double epanc(double x){

   /* This function calculates the value of the Epanechnikov kernel at point x */

   if(x > 1 || x < -1)
     return(0.0);
   else
     return (3.0 / 4.0) * (1 - x * x);
}


double s_t_0(double u, double h, int t_len){
 
   /* This function computes the kernel constant S_T0 needed for the local linear weights 
      of the kernel average hat(psi)_T(u,h) */ 

   int t, t_min, t_max;
   double result;

   result = 0.0;

   t_min = (int)(ceil((u - h) * (float)t_len));      
   if(t_min < 1) t_min = 1;
   t_max = (int)(floor((u + h) * (float)t_len));      
   if(t_max > t_len) t_max = t_len;

   for(t = t_min; t < (t_max+1); t++){
      result += epanc((t / (float)t_len - u) / h);
   }
   return(result / (t_len * h));
}


double s_t_1(double u, double h, int t_len){
 
   /* t_lenhis function computes the kernel constant S_T1 needed for the local linear weights 
      of the kernel average hat(psi)_T(u,h) */ 

   int t, t_min, t_max;
   double result;

   result = 0.0;

   t_min = (int)(ceil((u - h) * (float)t_len));      
   if(t_min < 1) t_min = 1;
   t_max = (int)(floor((u + h) * (float)t_len));      
   if(t_max > t_len) t_max = t_len;

   for(t = t_min; t < (t_max+1); t++){
      result += epanc((t / (float)t_len - u) / h) * ((t / (float)t_len - u) / h);
   }
   return(result / (t_len * h));
}

double s_t_2(double u, double h, int t_len){
   
   /* This function computes the kernel constant S_T2 needed for the local linear weights 
    of the kernel average hat(psi)_T(u,h) */ 
   
   int t, t_min, t_max;
   double result;
   
   result = 0.0;
   
   t_min = (int)(ceil((u - h) * (float)t_len));      
   if(t_min < 1) t_min = 1;
   t_max = (int)(floor((u + h) * (float)t_len));      
   if(t_max > t_len) t_max = t_len;
   
   for(t = t_min; t < (t_max+1); t++){
      result += epanc((t / (float)t_len - u) / h) * ((t / (float)t_len - u) / h) * ((t / (float)t_len - u) / h);
   }
   return(result / (t_len * h));
}

NumericVector kernel_averages(int t_len, NumericVector gset, NumericVector correct,
                              NumericVector data, double sigmahat, int n,
                              int deriv_order = 0){
   /* Inputs: 
      t_len       Integer, length of time series
      gset        vector of location-bandwidth points (u_1,...,u_n,h_1,...h_n)
      correct     vector of corrections of lenghth n (lambda(h_1), ..., lambda(h_n))
      data        time series of length t_len (y_1, ..., y_t_len)
      sigmahat    value by which we normalise our test statstics
      n           number of location-bandwidth points in gset
      deriv_order Order of the derivative of the trend that is investigated
    
           
      Output:
      result      vector of normalized kernel averages without corrections 
                  (\phi(u_1,h_1)/sigmahat,...,\phi(u_n,h_n)/sigmahat) and vector of absolute values with
                  corrections  (abs(\phi(u_1,h_1)/sigmahat) - \lambda(h_1),...,abs(\phi(u_n,h_n)/sigmahat) - \lambda(h_n))
                  stacked together
   */

   int i, t, t_min, t_max;
   double u, h, x, temp, result_temp, weight_norm;

   NumericVector result(2 * n);

   for(i = 0; i < n; i++){
      u = gset[i];
      h = gset[i + n];
      t_min = (int)(ceil((u - h) * (float)t_len));      
      if(t_min < 1) t_min = 1;
      t_max = (int)(floor((u + h) * (float)t_len));      
      if(t_max > t_len) t_max = t_len;

      result_temp = 0;
      weight_norm = 0;
      
      if(u > h && u < 1.0-h){
        for(t = t_min; t < (t_max+1); t++){
           x = ((t / (float)t_len - u) / h);           
           temp = epanc(x) * x;
           result_temp = result_temp + temp * data[t-1];
           weight_norm = weight_norm + temp * temp;
	      }
      } else {
        for(t = t_min; t < (t_max+1); t++){
           x = ((t / (float)t_len - u) / h);
           if(deriv_order == 1){
              temp = epanc(x) * (s_t_0(u, h, t_len) * x - s_t_1(u, h, t_len));
           } else if(deriv_order == 0){
              temp = epanc(x) * (s_t_2(u, h, t_len) - s_t_1(u, h, t_len) * x);
           } else {temp = 0;}
           result_temp = result_temp + temp * data[t-1];
           weight_norm = weight_norm + temp * temp;
         }
      }
      result[i] = result_temp / (sqrt(weight_norm) * sigmahat);
      result[i + n] = awert(result[i]) - correct[i];
   }
  return (result);
}

// [[Rcpp::export]]
NumericVector simulate_gaussian(int t_len, int n_ts, int sim_runs,
                                Rcpp::NumericVector gset, Rcpp::IntegerVector ijset,
                                double sigma, int deriv_order = 0,
                                bool correction = true, bool epidem = false){
   /* Inputs: 
    t_len       Integer, length of time series
    n_ts        Integer. Number of time series.
    sim_runs    Integer. Number of simulations needed to produce the quantiles.
    gset        A double vector of location-bandwidth points (u_1,...,u_n,h_1,...h_n).
    ijset       An integer vector of indices of the countries for the comparison
                (i_1, i_2, ..., j_{n_comparisons}).
    sigma       A double that is equal to sqrt{long-run varaince}
                in case of n_ts = 1, or the overdispersion in case of
                n_ts > 1. If not given, then the default is 1.
    deriv_order Integer. Order of the derivative of the trend that is investigated.
                Default is 0 => we analyse whether the trend itself >< than 0.
    correction  Boolean variable, true if we adjust the critical values
                in order to balance the significance of different hypotheses
    epidem      Boolean variable, true if we are working with epidemic time trends. Default is TRUE. 

        Output:
    Phi_vec     A double vector of length sim_runs that consists of the calculated
                Gaussian statistic for each simulation run.
    */
   
   int n = gset.length() / 2;
   Rcpp::NumericVector Phi_vec(sim_runs);
   Rcpp::NumericVector correct_a(n);
   Rcpp::NumericVector correct_b(n);
   Rcpp::NumericMatrix Phi(n_ts, n_ts);
   NumericVector vals(2 * n);
   NumericVector vals_cor(n);
   int pos, i, j, k, t, l, n_comparisons;
   double u, h, len, result_temp;
   
   for (k = 0; k < n; k++){
      len = 2 * gset[k + n];
      /* Here we need the whole length of the interval! */
      correct_b[k] = sqrt(2 * log(1 / len));
      correct_a[k] = sqrt(log(exp(1) / len)) / log(log(exp(exp(1)) / len));
   }
   
   if (n_ts == 1){
      for (pos = 0; pos < sim_runs; pos++) {
         NumericVector Z(t_len);
         Z = rnorm(t_len);
         NumericVector data_(t_len);
         data_ = sigma * Z;
         vals = kernel_averages(t_len, gset, correct_b, data_, sigma, n, deriv_order);
         vals_cor = vals[Rcpp::Range(n, 2 * n - 1)];
         Phi_vec(pos) = max(vals_cor);
      }
   } else {
      
      n_comparisons = ijset.length() / 2;

      for (pos = 0; pos < sim_runs; pos++) {
         NumericMatrix Z(t_len, n_ts);
         for (i = 0; i < n_ts; i++) {
            Z(_, i) = rnorm(t_len);
         }
         for (l = 0; l < n_comparisons; l++) {
            i = ijset[l] - 1;
            j = ijset[l + n_comparisons] - 1;
            NumericVector result(n);
            NumericVector result_cor(n);
            if (epidem) {
               for(k = 0; k < n; k++){
                  u = gset[k];
                  h = gset[k + n];
                  result_temp = 0;
                  for(t = 1; t < (t_len + 1); t++){
                     if (t / (float)t_len >= (u - h) && t / (float)t_len <= (u + h)){
                        result_temp += Z(t - 1, i) - Z(t - 1, j);
                     }
                  }
                  result[k] = awert(result_temp) / (sqrt(2 * t_len * 2 * h));
                  result_cor[k] = correct_a[k] * (awert(result_temp) / (sqrt(2 * t_len * 2 * h)) - correct_b[k]);
               }
            }
            else {
               NumericVector data_temp(t_len);
               double z_i_mean, z_j_mean;
               z_i_mean = mean(Z(_, i));
               z_j_mean = mean(Z(_, j));
               data_temp = (Z(_, i) - z_i_mean) - (Z(_, j) - z_j_mean);
               Rcpp::NumericVector result_temp(2 * n);
               result_temp = kernel_averages(t_len, gset, correct_b, data_temp, sqrt(2), n, 0);
               result = result_temp[Rcpp::Range(0, n - 1)];
               result_cor = result_temp[Rcpp::Range(n, 2 * n - 1)];
            }
            if (correction)
               Phi(i, j) = max(result_cor);
            else
               Phi(i, j) = max(result);
         }
         Phi_vec[pos] = max(Phi);
      }
   }
   return(Phi_vec);
}

// [[Rcpp::export]]
NumericVector simulate_gaussian_2(int t_len, int n_ts, Rcpp::NumericVector gset,
                                  Rcpp::IntegerVector ijset, double sigma = 1,
                                  int deriv_order = 0, bool correction = true,
                                  bool epidem = false){
   /* Inputs: 
    t_len       Integer, length of time series
    n_ts        Integer. Number of time series.
    gset        A double vector of location-bandwidth points
                (u_1,...,u_n,h_1,...h_n).
    ijset       An integer vector of indices of the countries for the comparison
                (i_1, i_2, ..., j_{n_comparisons}).
    sigma       A double that is equal to sqrt{long-run varaince}
                in case of n_ts = 1, or the overdispersion in case of
                n_ts > 1. If not given, then the default is 1.
    deriv_order Integer. Order of the derivative of the trend that is investigated.
    Default is 0 => we analyse whether the trend itself >< than 0.
    correction  Boolean variable, true if we adjust the critical values
    in order to balance the significance of different hypotheses
    epidem      Boolean variable, true if we are working with epidemic time trends. Default is TRUE. 
    
    Output:
    Phi_vec     A double vector of length sim_runs that consists of the calculated
    Gaussian statistic for each simulation run.
    */
   
   int n = gset.length() / 2;
   double Phi;
   Rcpp::NumericVector correct_a(n);
   Rcpp::NumericVector correct_b(n);
   NumericVector vals(2 * n);
   NumericVector vals_cor(n);
   int i, j, k, t, l, n_comparisons;
   double u, h, len, result_temp;
   
   for (k = 0; k < n; k++){
      len = 2 * gset[k + n];
      /* Here we need the whole length of the interval! */
      correct_b[k] = sqrt(2 * log(1 / len));
      correct_a[k] = sqrt(log(exp(1) / len)) / log(log(exp(exp(1)) / len));
   }
   
   if (n_ts == 1){
      NumericVector Z(t_len);
      Z = rnorm(t_len);
      NumericVector data_(t_len);
      data_ = sigma * Z;
      vals = kernel_averages(t_len, gset, correct_b, data_, sigma, n, deriv_order);
      vals_cor = vals[Rcpp::Range(n, 2 * n - 1)];
      if (correction)
         Phi = max(vals_cor);
      else
         Phi = max(vals[Rcpp::Range(0, n - 1)]);
      }
   } else {
      NumericMatrix Phi_mat(n_ts, n_ts);      
      n_comparisons = ijset.length() / 2;
      
      NumericMatrix Z(t_len, n_ts);
      for (i = 0; i < n_ts; i++) {
         Z(_, i) = rnorm(t_len);
      }
      for (l = 0; l < n_comparisons; l++) {
         i = ijset[l] - 1;
         j = ijset[l + n_comparisons] - 1;
         NumericVector result(n);
         NumericVector result_cor(n);
         if (epidem) {
            for(k = 0; k < n; k++){
               u = gset[k];
               h = gset[k + n];
               result_temp = 0;
               for(t = 1; t < (t_len + 1); t++){
                  if (t / (float)t_len >= (u - h) && t / (float)t_len <= (u + h)){
                     result_temp += Z(t - 1, i) - Z(t - 1, j);
                  }
               }
               result[k] = awert(result_temp) / (sqrt(2 * t_len * 2 * h));
               result_cor[k] = correct_a[k] * (awert(result_temp) / (sqrt(2 * t_len * 2 * h)) - correct_b[k]);
            }
         }
         else {
            NumericVector data_temp(t_len);
            double z_i_mean, z_j_mean;
            z_i_mean = mean(Z(_, i));
            z_j_mean = mean(Z(_, j));
            data_temp = (Z(_, i) - z_i_mean) - (Z(_, j) - z_j_mean);
            Rcpp::NumericVector result_temp(2 * n);
            result_temp = kernel_averages(t_len, gset, correct_b, data_temp, sqrt(2), n, 0);
            result = result_temp[Rcpp::Range(0, n - 1)];
            result_cor = result_temp[Rcpp::Range(n, 2 * n - 1)];
         }
         if (correction)
            Phi_mat(i, j) = max(result_cor);
         else
            Phi_mat(i, j) = max(result);
      }
      Phi = max(Phi_mat);
   }
   return(Phi);
}


// [[Rcpp::export]]
Rcpp::List compute_multiple_statistics(int t_len, int n_ts, Rcpp::NumericMatrix data,
                                       Rcpp::NumericVector gset, Rcpp::IntegerVector ijset,
                                       double sigma){
   /* Calculates the statistics (pairwise and overall) in case of multiple time series

        Inputs: 
    t_len       Integer, length of time series
    n_ts        Integer. Number of time series.
    data        Matrix of doubles, each column contains one time series.
    gset        A double vector of location-bandwidth points (u_1,...,u_n,h_1,...h_n).
    ijset       An integer vector of indices of the countries for the comparison
                (i_1, i_2, ..., j_{n_comparisons}).
    sigma       A double, estimator of the overdispersion parameter.

         Output:
    stat            A double matrix of pairwise multiscale statistics Psi_ij.
    vals_matrix     Matrix with n rows and n_ts * (n_ts - 1) / 2 columns where each column
                    contains values of the normalised Kernel averages in order to be able to perform 
                    the test on every separate interval.
    vals_cor_matrix Matrix with n rows and n_ts * (n_ts - 1) / 2 columns where each column
                    contains values of the normalised and corrected Kernel averages\
                    (that are adjusted such that we balance the significance of different null hypotheses)
                    in order to be able to perform the test on every separate interval.
    */
   
   int n = gset.length() / 2;
   Rcpp::NumericVector correct_b(n);
   Rcpp::NumericVector correct_a(n);
   Rcpp::NumericMatrix stat(n_ts, n_ts);
   NumericVector vals(n);
   NumericVector vals_cor(n);
   int i, j, k, l, t;
   int n_comparisons = ijset.length() / 2;
   Rcpp::NumericMatrix vals_cor_mat(n, n_comparisons);
   Rcpp::NumericMatrix vals_mat(n, n_comparisons);
   double len;
   double u, h, result_temp, weight_norm;

   for (k = 0; k < n; k++){
      len = 2 * gset[k + n];
      /* Here we need the whole length of the interval! */
      correct_b[k] = sqrt(2 * log(1 / len));
      correct_a[k] = sqrt(log(exp(1) / len)) / log(log(exp(exp(1)) / len));
   }
   
   for (l = 0; l < n_comparisons; l++) {
      i = ijset[l] - 1;
      j = ijset[l + n_comparisons] - 1;
      for(k = 0; k < n; k++){
         u = gset[k];
         h = gset[k + n];
               
         result_temp = 0;
         weight_norm = 0;
               
         for(t = 1; t < (t_len + 1); t++){
            if (t / (float)t_len >= (u - h) && t / (float)t_len <= (u + h)){
               result_temp += data(t-1, i) - data(t-1, j);
               weight_norm += data(t-1, i) + data(t-1, j);
            }
            vals[k] = awert(result_temp) / (sqrt(weight_norm) * sigma);
            vals_cor[k] = correct_a[k] * (awert(result_temp) / (sqrt(weight_norm) * sigma) - correct_b[k]);
         }
      }
      stat(i, j) = max(vals_cor);
      vals_mat(_, l) = vals;
      vals_cor_mat(_, l) = vals_cor;
   }
   return(Rcpp::List::create(Rcpp::Named("vals_matrix") = vals_mat,
                             Rcpp::Named("vals_cor_matrix") = vals_cor_mat,
                             Rcpp::Named("stat") = stat));
}

// [[Rcpp::export]]
Rcpp::List compute_multiple_statistics_2(int t_len, int n_ts, Rcpp::NumericMatrix data,
                                         Rcpp::NumericVector gset, Rcpp::IntegerVector ijset,
                                         Rcpp::NumericVector sigma_vec){
   /* Calculates the statistics (pairwise and overall) in case of multiple time series
    
    Inputs: 
    t_len       Integer, length of time series
    n_ts        Integer. Number of time series.
    data        Matrix of doubles, each column contains one time series.
    gset        A double vector of location-bandwidth points (u_1,...,u_n,h_1,...h_n).
    ijset       An integer vector of indices of the countries for the comparison
    (i_1, i_2, ..., j_{n_comparisons}).
    sigma_vec   A vector of the estimates of the long-run variances.
    
    Output:
    stat            A double matrix of pairwise multiscale statistics Psi_ij.
    vals_matrix     Matrix with n rows and n_ts * (n_ts - 1) / 2 columns where each column
                    contains values of the normalised Kernel averages in order to be able to perform 
                    the test on every separate interval.
    vals_cor_matrix Matrix with n rows and n_ts * (n_ts - 1) / 2 columns where each column
                    contains values of the normalised and corrected Kernel averages\
                    (that are adjusted such that we balance the significance of different null hypotheses)
                    in order to be able to perform the test on every separate interval.
    */
   int n = gset.length() / 2;
   Rcpp::NumericVector correct_b(n);
   Rcpp::NumericMatrix stat(n_ts, n_ts);
   NumericVector data_temp(t_len);
   int i, j, k, l;
   int n_comparisons = ijset.length() / 2;
   Rcpp::NumericMatrix vals_cor_mat(n, n_comparisons);
   Rcpp::NumericMatrix vals_mat(n, n_comparisons);
   double len, sigma;

   for (k = 0; k < n; k++){
      len = 2 * gset[k + n];
      /* Here we need the whole length of the interval! */
      correct_b[k] = sqrt(2 * log(1 / len));
   }
   
   for (l = 0; l < n_comparisons; l++) {
      i = ijset[l] - 1;
      j = ijset[l + n_comparisons] - 1;
      sigma = sqrt(sigma_vec[i] * sigma_vec[i] + sigma_vec[j] * sigma_vec[j]);
      data_temp = data(_, i) - data(_, j);
      Rcpp::NumericVector result(2 * n);
      result = kernel_averages(t_len, gset, correct_b, data_temp, sigma, n, 0);
      stat(i, j) = max(result[Rcpp::Range(n, 2 * n - 1)]);
      vals_mat(_, l) = result[Rcpp::Range(0, n - 1)];
      vals_cor_mat(_, l) = result[Rcpp::Range(n, 2 * n - 1)];
   }
   return(Rcpp::List::create(Rcpp::Named("vals_matrix") = vals_mat,
                             Rcpp::Named("vals_cor_matrix") = vals_cor_mat,
                             Rcpp::Named("stat") = stat));
}

// [[Rcpp::export]]
Rcpp::List compute_single_statistics(int t_len, Rcpp::NumericVector data, Rcpp::NumericVector gset,
                              double sigma, int deriv_order) {
   /* Calculates the statistics in case of a single time series.
    
    Inputs: 
    t_len       Integer, length of time series
    data        Vector of doubles of length t_len that contains the time series for analysis.
    gset        A double vector of location-bandwidth points (u_1,...,u_n,h_1,...h_n).
    sigma       A double, estimator of the overdispersion parameter.
    deriv_order Integer. Order of the derivative of the trend that is investigated.
                Default is 0 => we analyse whether the trend itself >< than 0.
    
    Output:
    stat        A double, value of the test statistics calculated as max(vals_cor).
    vals        A double vector of length n of kernel averages (sign included).
    vals_cor    A double vector of absolute value of normalized kernel averages with
                correction \eqn{(abs(\phi(u_1,h_1)/sigmahat) - \lambda(h_1), ...,
                abs(\phi(u_n,h_n)/sigmahat) - \lambda(h_n))}
    */
   
   int n = gset.length() / 2;
   Rcpp::NumericVector correct(n);
   Rcpp::NumericVector vals(2 * n);
   int k;
   double stat;

   for (k = 0; k < n; k++){
      correct[k] = sqrt(2 * log(1 / (2 * gset[k + n])));
   }
   
   vals = kernel_averages(t_len, gset, correct, data, sigma, n, deriv_order);
   stat  = max(vals[Rcpp::Range(n, 2 * n - 1)]);
   return Rcpp::List::create(Rcpp::Named("vals") = vals[Rcpp::Range(0, n - 1)],
                             Rcpp::Named("vals_cor") = vals[Rcpp::Range(n, 2 * n - 1)],
                             Rcpp::Named("stat") = stat);
}
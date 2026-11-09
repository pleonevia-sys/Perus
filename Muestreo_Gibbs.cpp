#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector actualizar_z_cpp(NumericVector datos, 
                               NumericVector mu, 
                               NumericVector sigma2, 
                               NumericVector pi_pesos, 
                               int K_max) {
  int n = datos.size();
  IntegerVector z(n);
  
  // Bucle optimizado en C++ para las 53,000 observaciones
  for (int i = 0; i < n; i++) {
    NumericVector probs(K_max);
    double suma_probs = 0.0;
    
    for (int k = 0; k < K_max; k++) {
      double diff = datos[i] - mu[k];
      double var_k = sigma2[k];
      // Densidad Normal evaluada de forma directa
      double dens = (1.0 / sqrt(2.0 * M_PI * var_k)) * exp(-0.5 * (diff * diff) / var_k);
      probs[k] = pi_pesos[k] * dens;
      suma_probs += probs[k];
    }
    
    // Muestreo categorico ponderado seguro
    double u = R::runif(0.0, 1.0) * (suma_probs + 1e-300);
    double acum = 0.0;
    int cluster_elegido = K_max - 1;
    
    for (int k = 0; k < K_max; k++) {
      acum += probs[k];
      if (u <= acum) {
        cluster_elegido = k;
        break;
      }
    }
    
    // Retornar indice en base 1 para R (1, 2, ..., K_max)
    z[i] = cluster_elegido + 1;
  }
  
  return z;
}
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double HCalC(NumericVector dt_vector) {
  int n = dt_vector.size();
  double sum_abs_diff = 0.0;
  
  for (int i = 0; i < n; ++i) {
    double abs_diff = 0.0;
    for (int j = 0; j < n; ++j) {
      abs_diff += std::abs(dt_vector[i] - dt_vector[j]);
    }
    sum_abs_diff += abs_diff;
  }
  
  double mean_dt = mean(dt_vector);
  return sum_abs_diff / (2.0 * n * n * mean_dt);
}

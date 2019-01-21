#include <math.h>
// #include "helper_functions.h"

double weighted_entropy(double *v, int dim) {
  double ans = 0, v_sum = sum(v, dim);
  if (v_sum != 0) {
    int i;
    for (i = 0; i < dim; i++) {
      if (v[i] != 0) {
        ans -= v[i] * log2(v[i] / v_sum);
      }
    }
    // ans *= v_sum;
  }
  return ans;
}

double* clusters_entropy(double *clusters, int k, int dim) {
  int i;
  // printf("initialize results array of size %d\n", k);
  double* results = (double *)malloc(sizeof(double) * k);
  // printf("initialized array\n");
  // printf("calculate entropies\n");
  for (i = 0; i < k; i++) {
    results[i] = weighted_entropy(&clusters[i*dim], dim);
  }
  // printf("return results\n");
  return results;
}

double partition_entropy(double *clusters, int k, int dim)
  {
    double* results = clusters_entropy(clusters, k, dim);
    double ans = sum(results, k);
    free(results);
    return ans;
  }

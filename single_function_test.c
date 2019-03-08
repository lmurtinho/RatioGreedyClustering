#include "helper/helper_functions.h"
#include "helper/random_functions.h"
#include "dhillon/dhillon.h"
#include "dominance/dominance.h"
#include "coresets/coresets.h"
#include <stdio.h>

int main () {
  int i, k, /* n = 170946, dim=103, */ n = 51840, dim = 20,
    max_iter=100, iter[] = {0};
  int size_k = 40;
  int k_list[] = {1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
                  1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
                  2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000,
                  2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000};

  clock_t start, end;

  double *data = data_from_csv("ng20.csv", n, dim);
  printf("got data\n");
  k = k_list[size_k-1];
  double *clusters = (double *)malloc(sizeof(double)*k*dim);
  int *assigned = (int *)malloc(sizeof(int)*n);

  FILE *results, *f;

  results = fopen("co_ng20_results.csv", "w");
  fprintf(results, "n_clusters,time,entropy,iterations");
  //fclose(times);

  f = fopen("rg.csv", "w");

  for (i = 0; i < size_k; i++) {
    k = k_list[i];
    printf("\n\nWith %d clusters\n", k);
    start = clock();
    // printf("started\n");

    assigned =
      co_clustering(data, n, k, dim, 2500, 100, iter);
      //rg_clustering(data, n, k, dim, 0, iter, cluster_assign_ext);
      // di_clustering(data, n, k, dim, 100, iter,
      // f, f, f, start);
    end = clock();
    // printf("ended\n");
    check_assignment(assigned, 0, k, n);
    //times = fopen("new_times.csv", "w");
    maximization(data, assigned, clusters, n, k, dim);
    fprintf(results, "\n%d,%f,%f,%d", k,
      ( (double)(end - start) ) / CLOCKS_PER_SEC,
      partition_entropy(clusters, k, dim), iter[0]);
    //fclose(times);
    printf("time: %f seconds\t", ( (double)(end - start) ) / CLOCKS_PER_SEC);
    printf("entropy: %f\t", partition_entropy(clusters, k, dim));
    printf("iters: %d\n", iter[0]);

    // double *ce = clusters_entropy(clusters, k, dim);
    // printf("clusters entropy\n");
    // print_array(ce, 1, k);
    // free(ce);
  }

  fclose(results);
  fclose(f);
  free(data);
  free(clusters);
  free(assigned);
}

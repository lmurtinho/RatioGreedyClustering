#include "helper_functions.h"
#include "random_functions.h"
#include "dominance.h"
#include "dhillon.h"
#include <stdio.h>

int main () {
  int i, k, n = 170946, dim=103, /* n = 51840, dim = 20, */
    max_iter=100, iter[] = {0};
  int size_k = 29;
  int k_list[] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                  17, 18, 19, 20, 30, 40, 50, 100, 200, 300, 400, 500,
                  1000, 2000};

  clock_t start, end;

  double *data = data_from_csv("rcv1.csv", n, dim);
  printf("got data\n");
  k = k_list[size_k-1];
  double *clusters = (double *)malloc(sizeof(double)*k*dim);
  int *assigned = (int *)malloc(sizeof(int)*n);

  FILE *results, *f;

  results = fopen("dom_heap_rcv1_results_sun2.csv", "w");
  fprintf(results, "n_clusters,time,entropy,iterations");
  //fclose(times);

  f = fopen("dom_heap.csv", "w");

  for (i = 0; i < size_k; i++) {
    k = k_list[i];
    printf("\n\nWith %d clusters\n", k);
    start = clock();
    // printf("started\n");

    assigned =
      rg_clustering(data, n, k, dim, 0, iter, cluster_assign_ext);
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

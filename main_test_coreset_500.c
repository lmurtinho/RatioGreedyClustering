#include "helper/helper_functions.h"
#include "helper/random_functions.h"
#include "dhillon/dhillon.h"
#include "dominance/dominance.h"
#include "coresets/coresets.h"
#include <stdio.h>

int main (int argc, char *argv[]) {
  int i, k, n = atoi(argv[2]), dim = atoi(argv[3]),
    max_iter=100, iter[] = {0}, size_k = 29,
    k_list[] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                  17, 18, 19, 20, 30, 40, 50, 100, 200, 300, 400, 500,
                  1000, 2000};
  char *name = argv[1];
  char *data_file_name  = concat(name, ".csv");
  char *times_file_name = concat(name, "_CO_500_times.csv");
  char *entrs_file_name = concat(name, "_CO_500_entrs.csv");
  char *iters_file_name = concat(name, "_CO_500_iters.csv");

  clock_t start, end;

  double *data = data_from_csv(data_file_name, n, dim);
  printf("got data\n");
  k = k_list[size_k-1];
  double *clusters = (double *)malloc(sizeof(double)*k*dim);
  int *assigned = (int *)malloc(sizeof(int)*n);

  FILE *times, *entrs, *iters;

  times = fopen(times_file_name, "w");
  fprintf(times, "n_clusters,CO_500_init,CO_500_iter1,CO_500_iter5,CO_500_iter10,CO_500_final");
  //fclose(times);

  entrs = fopen(entrs_file_name, "w");
  fprintf(entrs, "n_clusters,CO_500_init,CO_500_iter1,CO_500_iter5,CO_500_iter10,CO_500_final");
  //fclose(entrs);

  iters = fopen(iters_file_name, "w");
  fprintf(iters, "n_clusters,CO_500_init,CO_500_iter1,CO_500_iter5,CO_500_iter10,CO_500_final");
  //fclose(iters);

  for (i = 0; i < size_k; i++) {
    k = k_list[i];
    printf("\n\nWith %d clusters\n", k);

    printf("CO clustering with 500 samples (init)\n");
    start = clock();
    assigned = co_clustering(data, n, k, dim, 500, 0, iter);
    end = clock();
    check_assignment(assigned, 0, k, n);
    fprintf(times, "\n%d,%f", k, ( (double)(end - start) ) / CLOCKS_PER_SEC);
    maximization(data, assigned, clusters, n, k, dim);
    fprintf(entrs, "\n%d,%f", k, partition_entropy(clusters, k, dim));
    fprintf(iters, "\n%d,%d", k, iter[0]);
    printf("time: %f seconds\t", ( (double)(end - start) ) / CLOCKS_PER_SEC);
    printf("entropy: %f\t", partition_entropy(clusters, k, dim));
    printf("iters: %d\n", iter[0]);

    printf("CO clustering with 500 samples (iter1)\n");
    start = clock();
    assigned = co_clustering(data, n, k, dim, 500, 1, iter);
    end = clock();
    check_assignment(assigned, 0, k, n);
    fprintf(times, ",%f", ( (double)(end - start) ) / CLOCKS_PER_SEC);
    maximization(data, assigned, clusters, n, k, dim);
    fprintf(entrs, ",%f", partition_entropy(clusters, k, dim));
    fprintf(iters, ",%d", iter[0]);
    printf("time: %f seconds\t", ( (double)(end - start) ) / CLOCKS_PER_SEC);
    printf("entropy: %f\t", partition_entropy(clusters, k, dim));
    printf("iters: %d\n", iter[0]);

    printf("CO clustering with 500 samples (iter5)\n");
    start = clock();
    assigned = co_clustering(data, n, k, dim, 500, 5, iter);
    end = clock();
    check_assignment(assigned, 0, k, n);
    fprintf(times, ",%f", ( (double)(end - start) ) / CLOCKS_PER_SEC);
    maximization(data, assigned, clusters, n, k, dim);
    fprintf(entrs, ",%f", partition_entropy(clusters, k, dim));
    fprintf(iters, ",%d", iter[0]);
    printf("time: %f seconds\t", ( (double)(end - start) ) / CLOCKS_PER_SEC);
    printf("entropy: %f\t", partition_entropy(clusters, k, dim));
    printf("iters: %d\n", iter[0]);

    printf("CO clustering with 500 samples (iter10)\n");
    start = clock();
    assigned = co_clustering(data, n, k, dim, 500, 10, iter);
    end = clock();
    check_assignment(assigned, 0, k, n);
    fprintf(times, ",%f", ( (double)(end - start) ) / CLOCKS_PER_SEC);
    maximization(data, assigned, clusters, n, k, dim);
    fprintf(entrs, ",%f", partition_entropy(clusters, k, dim));
    fprintf(iters, ",%d", iter[0]);
    printf("time: %f seconds\t", ( (double)(end - start) ) / CLOCKS_PER_SEC);
    printf("entropy: %f\t", partition_entropy(clusters, k, dim));
    printf("iters: %d\n", iter[0]);

    printf("CO clustering with 500 samples (maxiter100)\n");
    start = clock();
    assigned = co_clustering(data, n, k, dim, 500, 100, iter);
    end = clock();
    check_assignment(assigned, 0, k, n);
    fprintf(times, ",%f", ( (double)(end - start) ) / CLOCKS_PER_SEC);
    maximization(data, assigned, clusters, n, k, dim);
    fprintf(entrs, ",%f", partition_entropy(clusters, k, dim));
    fprintf(iters, ",%d", iter[0]);
    printf("time: %f seconds\t", ( (double)(end - start) ) / CLOCKS_PER_SEC);
    printf("entropy: %f\t", partition_entropy(clusters, k, dim));
    printf("iters: %d\n", iter[0]);
  }

  fclose(times);
  fclose(iters);
  fclose(entrs);
  free(data);
  free(clusters);
  free(assigned);
  free(iters_file_name);
  free(times_file_name);
  free(entrs_file_name);
  free(data_file_name);
}

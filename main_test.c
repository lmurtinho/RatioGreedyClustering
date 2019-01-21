#include "helper/helper_functions.h"
#include "helper/random_functions.h"
#include "dhillon/dhillon.h"
#include "dominance/dominance.h"
#include <stdio.h>

int main (int argc, char *argv[]) {
  int i, k, n = atoi(argv[2]), dim = atoi(argv[3]),
    max_iter=100, iter[] = {0}, size_k = 29,
    k_list[] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                  17, 18, 19, 20, 30, 40, 50, 100, 200, 300, 400, 500,
                  1000, 2000};
  char *name = argv[1];
  char *data_file_name  = concat(name, ".csv");
  char *times_file_name = concat(name, "_times.csv");
  char *entrs_file_name = concat(name, "_entrs.csv");
  char *iters_file_name = concat(name, "_iters.csv");

  clock_t start, end;

  double *data = data_from_csv(data_file_name, n, dim);
  printf("got data\n");
  k = k_list[size_k-1];
  double *clusters = (double *)malloc(sizeof(double)*k*dim);
  int *assigned = (int *)malloc(sizeof(int)*n);

  FILE *times, *entrs, *iters;

  times = fopen(times_file_name, "w");
  fprintf(times, "n_clusters,Random,RG,DI_init,DI_iter1,DI_iter5,DI_iter10,DI_final");
  //fclose(times);

  entrs = fopen(entrs_file_name, "w");
  fprintf(entrs, "n_clusters,Random,RG,DI_init,DI_iter1,DI_iter5,DI_iter10,DI_final");
  //fclose(entrs);

  iters = fopen(iters_file_name, "w");
  fprintf(iters, "n_clusters,Random,RG,DI_init,DI_iter1,DI_iter5,DI_iter10,DI_final");
  //fclose(iters);

  for (i = 0; i < size_k; i++) {
    k = k_list[i];
    printf("\n\nWith %d clusters\n", k);

    printf("Random clustering\n");
    start = clock();
    assigned = rd_clustering(data, n, k, dim, 0, iter);
    end = clock();
    check_assignment(assigned, 0, k, n);
    //times = fopen("new_times.csv", "w");
    fprintf(times, "\n%d,%f", k, ( (double)(end - start) ) / CLOCKS_PER_SEC);
    //fclose(times);
    maximization(data, assigned, clusters, n, k, dim);
    fprintf(entrs, "\n%d,%f", k, partition_entropy(clusters, k, dim));
    fprintf(iters, "\n%d,0", k);

    printf("time: %f seconds\t", ( (double)(end - start) ) / CLOCKS_PER_SEC);
    printf("entropy: %f\t", partition_entropy(clusters, k, dim));
    printf("iters: 0\n");

    printf("RG clustering\n");
    start = clock();
    assigned = rg_clustering(data, n, k, dim, 0, iter, cluster_assign_ext);
    end = clock();
    check_assignment(assigned, 0, k, n);
    fprintf(times, ",%f", ( (double)(end - start) ) / CLOCKS_PER_SEC);
    maximization(data, assigned, clusters, n, k, dim);
    fprintf(entrs, ",%f", partition_entropy(clusters, k, dim));
    fprintf(iters, ",0");
    printf("time: %f seconds\t", ( (double)(end - start) ) / CLOCKS_PER_SEC);
    printf("entropy: %f\t", partition_entropy(clusters, k, dim));
    printf("iters: 0\n");

    printf("DI clustering\n");
    start = clock();
    assigned = di_clustering(data, n, k, dim, 100, iter, times, entrs, iters, start);
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

#include "dhillon_functions.h"

typedef bool (*InitClusterFunc) (double *data, double *clusters,
  double *norm_data, double *data_logs, int n, int k, int dim);

bool expectation(double *norm_data, double *data_logs,
                  int *assigned, double *clusters,
                  int n, int k, int dim) {
  int i, j, curr, new;
  double dist, d, c_sum;

  // for (i = 0; i < k; i++) {
  //   c_sum = sum(&clusters[i*dim], dim);
  //   if ( c_sum <= 0) {
  //     printf("sum of cluster %d = %f\n", i, c_sum);
  //   }
  // }

  double* cluster_logs = get_logs(clusters, k, dim);
  int *vecs_per_cluster = (int *)malloc(sizeof(int)*k);
  fill_array_int(vecs_per_cluster, 0, k);

  bool changed = false;                      // no change made so far

  for (i = 0; i < n; i++) {              // for each vector
    dist = BIG_DOUBLE;                 // current distance is largest possible
    curr = assigned[i];                 // get vector's current assignment
    new = -1;                            // new assignment is non-existent
    // dist = kl_div(norm_data, data_logs, cluster_logs, i, new, dim);
    for (j = 0; j < k; j++) {            // for each cluster
      // CHANGE LATER
      d = kl_div(norm_data, data_logs, cluster_logs, i, j, dim);  // get distance between cluster and vector
      if (d < dist)  {                    // if smaller than smallest distance
        new = j;                         // new assignment is current cluster
        vecs_per_cluster[j]++;
        dist = d;                        // new distance is current distance
      }
      // if (j == 288) {
      //   printf("dist to 288: %f\tdist to %d: %f\n", d, new, dist);
      //   printf("sum of cluster log: %f\n", sum(&cluster_logs[j*dim], dim));
      // }
      if (d >= BIG_DOUBLE) {
        printf("infinite distance!\n");
      }
    }
    if (new != curr) {                   // if assignment changed
      changed = true;                   // set change to true
      assigned[i] = new;                // change assignment to new cluster
    }
    if (new == -1) {
      printf("not assigned!\n");
    }
  }

  for (i = 0; i < k; i++) {
    if (vecs_per_cluster[i] == 0) {
      j = (int) (rand() % n);
      assigned[j] = i;
      changed = true;
    }
  }

  free(vecs_per_cluster);
  free(cluster_logs);
  return changed;
}

// K-MEANS WITH DHILLON INITIALIZATION
int* di_clustering(double *data, int n, int k, int dim, int max_iter, int *iter,
            FILE *times, FILE *entrs, FILE *iters, clock_t start) {

  // create arrays for cluster assignment and cluster vectors
  clock_t end;
  int *assigned = (int *)malloc(sizeof(int) * n);
  int i, j, c;
  double *clusters = (double *)malloc(sizeof(double) * k * dim), pe;
  double timer;
  iter[0] = 0;
  // CHANGED HERE
  // normalize_array(data, n, dim);

  // assign initial centers (++ implementation)
  // printf("initializing clusters\n");
  initial_di(data, clusters, n, k, dim);
  // timer = ( (double)(end - start) ) / CLOCKS_PER_SEC;
  // printf("initialized in %.4f seconds\n", timer);
  // expectation(norm_data, data_logs, assigned, clusters, n, k, dim);
  end = clock();

  pe = partition_entropy(clusters, k, dim);
  fprintf(times, ",%f", ( (double)(end - start) ) / CLOCKS_PER_SEC);
  fprintf(entrs, ",%f", pe);
  fprintf(iters, ",0");

  double *norm_data = (double *)malloc(sizeof(double) * n * dim);
  copy_to(data, norm_data, n * dim);
  normalize_array(norm_data, n, dim);

  double *data_logs = get_logs_from_normal(norm_data, n * dim);

  // maximization(data, assigned, clusters, n, k, dim);

  // printf("first assignment\n");
  bool changed = expectation(norm_data, data_logs, assigned, clusters, n, k, dim);
  // check_assignment(assigned, 0, k, n);
  // printf("\n\n\n");

  // start = clock();
  // double old_pe;
  // double *ce = (double *)malloc(sizeof(double)*k);
  while(changed && (iter[0] < max_iter) ) {
    // old_pe = partition_entropy(clusters, k, dim);
    end = clock();
    iter[0]++;
    // printf("iteration %d\t", iter[0]);
    maximization(data, assigned, clusters, n, k, dim);
    // CHANGED HERE
    // normalize_array(clusters, k, dim);
    pe = partition_entropy(clusters, k, dim);
    if ( (iter[0] == 1) || (iter[0] == 5) || (iter[0] == 10) ) {
      fprintf(times, ",%f", ( (double)(end - start) ) / CLOCKS_PER_SEC);
      fprintf(entrs, ",%f", pe);
      fprintf(iters, ",%d", iter[0]);
    }// printf("entropy difference: %f", old_pe/pe);
    // if (pe > old_pe) {
    //   printf("worse partition at iteration %d!\n", iter[0]);
    //   printf("%.4f vs. %.4f\n", old_pe, pe);
    // }
    // ce = clusters_entropy(clusters, k, dim);
    // printf("clusters entropies:\n");
    // print_array(ce, 1, k);
    // start = clock();
    changed = expectation(norm_data, data_logs, assigned, clusters, n, k, dim);
    // printf("iter %d\n", iter[0]);
    // check_assignment(assigned, 0, k, n);
    // printf("\tpartition entropy: %f\n", pe);
  }
  if (iter[0] < 10) {
    fprintf(times, ",NA");
    fprintf(entrs, ",NA");
    fprintf(iters, ",NA");
  }
  if (iter[0] < 5) {
    fprintf(times, ",NA");
    fprintf(entrs, ",NA");
    fprintf(iters, ",NA");
  }
  if (iter[0] < 1) {
    fprintf(times, ",NA");
    fprintf(entrs, ",NA");
    fprintf(iters, ",NA");
  }
  // end = clock();
  // timer = ( (double)(end - start) ) / CLOCKS_PER_SEC;
  // printf("Average per iteration: %.4f seconds\n", timer / iter[0]);
  // printf("done in %d iterations\n", iter[0]);
  free(clusters);
  free(norm_data);
  free(data_logs);
  return assigned;
}

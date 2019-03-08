#include <math.h>

// initial centers found according to ++ implementation
void initial_pp(double *data, double *clusters, double *norm_data,
    double *data_logs, double *cluster_logs, double *distances,
    int n, int k, int dim) {
  double d, sum_distances;
  int curr_cluster, i, j;

  double *probs         = (double *)malloc(sizeof(double) * n);
  double *norm_clusters = (double *)malloc(sizeof(double) * k * n);

  fill_array(distances, (double)INFINITY, n);

  curr_cluster = rand() % n;
  copy_to(&data[curr_cluster * dim], clusters, dim);
  copy_to(clusters, norm_clusters, dim);
  normalize_array(norm_clusters, 1, dim);

  for (i = 0; i < dim; i++) {
    cluster_logs[i] = log2(norm_clusters[i]);
  }

  for (i = 1; i < k; i++) {
    for (j = 0; j < n; j++) {
      d = kl_div(norm_data, data_logs, cluster_logs, j, i-1, dim);
      d *= d;
      if (d < distances[j]) {
        distances[j] = d;
      }
    }

    copy_to(distances, probs, n);
    sum_distances = sum(distances, n);
    normalize_vector(probs, sum_distances, n);

    curr_cluster = find_next_number(probs, n);

    copy_to(&data[curr_cluster * dim], &clusters[i * dim], dim);
    copy_to(&clusters[i * dim], &norm_clusters[i * dim], dim);
    normalize_array(&norm_clusters[i * dim], 1, dim);

    for (j = i*dim; j < (i*dim + dim); j++) {
      cluster_logs[j] = log2(norm_clusters[j]);
    }

  }
  free(probs);
  free(norm_clusters);
  return;
}

int* get_vecs_per_cluster(int *assigned, int n, int k) {
  int i, *ans = (int *)malloc(sizeof(int) * k);
  fill_array_int(ans, 0, n);
  for (i = 0; i < n; i++) {
    ans[assigned[i]]++;
  }
  return ans;
}

double *get_norm_data(double *data, int n, int dim) {
  double *ans = (double *)malloc(sizeof(double) * n * dim);
  copy_to(data, ans, n * dim);
  normalize_array(ans, n, dim);
  return ans;
}

double* get_distances(double *norm_data, double *data_logs,
  double *cluster_logs, int *assigned, int n, int dim) {
    int i;
    double *ans = (double *)malloc(sizeof(double) * n);
    for (i = 0; i < n; i++) {
      ans[i] = kl_div(norm_data, data_logs, cluster_logs, i, assigned[i], dim);
    }
    return ans;
  }

double* get_distances_per_cluster(double *distances, int *assigned,
  int n, int k) {
  int i;
  double *ans = (double *)malloc(sizeof(double) * k);
  fill_array(ans, 0, k);
  for (i = 0; i < n; i++) {
    ans[assigned[i]] += distances[i];
  }
  return ans;
}

double* get_sensitivities(double *data, double *norm_data, double *data_logs,
    double *clusters, int n, int k, int dim) {

  clock_t start = clock(), end;
  int i, c;
  double v, alpha = 16 * (log2(k) + 2), // Step 1
    *ans = (double *)malloc(sizeof(double) * n),
    *cluster_logs = (double *)malloc(sizeof(double) * k * dim),
    *distances = (double *)malloc(sizeof(double) * n);

    end = clock();
    printf("Set up in %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);
    start = clock();

  // Steps 2-3
  initial_pp(data, clusters, norm_data, data_logs, cluster_logs, distances,
    n, k, dim);

  end = clock();
  printf("++ initialization in %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);
  start = clock();

  int *assigned = (int *)malloc(sizeof(int) * n);
  expectation(norm_data, data_logs, assigned, clusters, n, k, dim);

  end = clock();
  printf("First assignment in %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);
  start = clock();

  double *distances_per_cluster = get_distances_per_cluster(distances, assigned,
      n, k);
    // Step 4
  end = clock();
  printf("Distances in %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);
  start = clock();

  double distance_avg = sum(distances, n) / n;

  end = clock();
  printf("Average distance in %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);
  start = clock();

  // Steps 5-6
  for (i = 0; i < n; i ++) {
    c = assigned[i];
    v = alpha * distances[i] / distance_avg;
    v += 2 * alpha * distances_per_cluster[c] / (k * distance_avg);
    v += (4 * n) / k;
    ans[i] = v;
  }

  end = clock();
  printf("Sensitivities in %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);
  start = clock();

  free(assigned);
  free(cluster_logs);
  free(distances);
  free(distances_per_cluster);

  end = clock();
  printf("Freed arrays in %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);

  return ans;
}

int* coreset_construction(double *data, double *norm_data, double *data_logs,
    double *clusters, int n, int k, int dim, int m) {
  clock_t start = clock(), end;
  double *sensitivities = get_sensitivities(data, norm_data, data_logs,
    clusters, n, k, dim);
  end = clock();
  printf("Calculated sensitivities in %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);

  start = clock();
  int *ans = random_sampling(sensitivities, n, m, false);
  end = clock();
  printf("Random sampling in %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);

  return ans;
}

void kmeans(double *data, double *norm_data, double *data_logs,
    double *clusters, int n, int k, int dim, int max_iter, int *iter) {

  int i, j, c, *assigned = (int *)malloc(sizeof(int) * n);
  fill_array_int(assigned, -1, n);
  iter[0] = 0;
  // printf("\niter %d\n", iter[0]);
  bool changed = expectation(norm_data, data_logs, assigned, clusters,
    n, k, dim);
  // check_assignment(assigned, 0, k, n);
  while(changed && (iter[0] < max_iter) ) {
    iter[0]++;
    // printf("\niter %d\n", iter[0]);
    maximization(data, assigned, clusters, n, k, dim);
    changed = expectation(norm_data, data_logs, assigned, clusters, n, k, dim);
    // check_assignment(assigned, 0, k, n);
  }
  //printf("done\n\n");
}

int* co_clustering(double *data, int n, int k, int dim, int m,
  int max_iter, int *iter) {
  clock_t start = clock(), end;
  double *norm_data = get_norm_data(data, n, dim),
         *data_logs = get_logs_from_normal(norm_data, n * dim),
         *clusters  = (double *)malloc(sizeof(double) * k * dim);
  end = clock();
  printf("Set up in %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);
  start = clock();

  int *coreset_indices = coreset_construction(data, norm_data, data_logs,
    clusters, n, k, dim, m);

  end = clock();
  printf("Found coreset indices in %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);
  start = clock();

  double *coreset_data = get_indices(data, coreset_indices, n, dim, m);
  double *coreset_norm = get_indices(norm_data, coreset_indices, n, dim, m);
  double *coreset_logs = get_indices(data_logs, coreset_indices, n, dim, m);

  end = clock();
  printf("Built coreset arrays in %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);
  start = clock();

  kmeans(coreset_data, coreset_norm, coreset_logs, clusters, m, k, dim,
    max_iter, iter);

  end = clock();
  printf("Ran k-means in coreset in %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);
  start = clock();

  // printf("kmeans done\n");
  int *assigned = (int *)malloc(sizeof(int) * n);
  expectation(norm_data, data_logs, assigned, clusters, n, k, dim);

  end = clock();
  printf("Assigned whole set in %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);
  start = clock();

  free(norm_data);
  free(data_logs);
  free(clusters);
  free(coreset_indices);
  free(coreset_data);
  free(coreset_norm);
  free(coreset_logs);

  end = clock();
  printf("Freed data in %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);


  return assigned;
}

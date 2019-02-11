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

  int i, c;
  double v, alpha = 16 * (log2(k) + 2), // Step 1
    *ans = (double *)malloc(sizeof(double) * n),
    *cluster_logs = (double *)malloc(sizeof(double) * k * dim),
    *distances = (double *)malloc(sizeof(double) * n);

  // Steps 2-3
  initial_pp(data, clusters, norm_data, data_logs, cluster_logs, distances,
    n, k, dim);
  int *assigned = (int *)malloc(sizeof(int) * n);
  expectation(norm_data, data_logs, assigned, clusters, n, k, dim);

  double *distances_per_cluster = get_distances_per_cluster(distances, assigned,
      n, k),
    // Step 4
    distance_avg = sum(distances) / n;

  // Steps 5-6
  for (i = 0; i < n; i ++) {
    c = assigned[i];
    v = alpha * distances[i] / distance_avg;
    v += 2 * alpha * distances_per_cluster[c] / (k * distance_avg);
    v += (4 * n) / k;
    ans[i] = v;
  }
  free(norm_data);
  free(data_logs);
  free(clusters);
  free(assigned);
  free(cluster_logs);
  free(distances);
  free(distances_per_cluster);
  return ans;
}

int* coreset_construction(double *data, double *norm_data, double *data_logs,
    double *clusters, int n, int k, int dim, int m) {
  double *sensitivities = get_sensitivities(data, norm_data, data_logs,
    clusters, n, k, dim);
  return random_sampling(sensitivities, n, m, false);
}

int* kmeans(double *data, double *norm_data, double *data_logs,
    double *clusters, int n, int k, int dim, int max_iter, int *iter) {

  int i, j, c, *assigned = (int *)malloc(sizeof(int) * n);
  fill_array_int(assigned, -1, n);
  iter[0] = 0;
  bool changed = expectation(norm_data, data_logs, assigned, clusters,
    n, k, dim);
  while(changed && (iter[0] < max_iter) ) {
    iter[0]++;
    maximization(data, assigned, clusters, n, k, dim);
    changed = expectation(norm_data, data_logs, assigned, clusters, n, k, dim);
  }
  return assigned;
}

int* coreset_clustering(double *data, int n, int k, int dim, int m,
  int max_iter, int *iter) {
  clock_t start = clock();
  double *norm_data = get_norm_data(data, n, dim),
         *data_logs = get_logs_from_normal(norm_data, n * dim),
         *clusters  = (double *)malloc(sizeof(double) * k * dim);
  int *coreset_indices = coreset_construction(data, norm_data, data_logs,
    clusters, n, k, dim, m);
  double *coreset_data = get_indices(data, coreset_indices, n, dim, m);
  double *coreset_norm = get_indices(norm_data, coreset_indices, n, dim, m);
  double *coreset_logs = get_indices(data_logs, coreset_indices, n, dim, m);

  int *coreset_assigned = kmeans(coreset_data, coreset_norm, coreset_logs,
    clusters, m, k, dim, max_iter, iter);
  int *assigned = (int *)malloc(sizeof(int) * n);
  maximization(data, assigned, clusters, n, k, dim);

  free(norm_data);
  free(data_logs);
  free(clusters);
  free(coreset_indices);
  free(coreset_data);
  free(coreset_norm);
  free(coreset_logs);
  free(coreset_assigned);

  return assigned;
}

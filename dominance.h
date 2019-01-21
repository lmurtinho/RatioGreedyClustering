#include "dominance_functions.h"
#include "dominance_heap.h"

typedef void (*ClusterAssignFunc) (double *cluster_data, int *cluster_assignment,
  int first_cluster, int avail_clusters, int n, int dim);

double *transform_data(double *data, int n, int k, int dim) {
  double *ans = (double *)malloc(sizeof(double)*n*k),
    *data_sum = array_sum(data, n, dim), rest, val;
    // printf("data sum:\n");
    // print_array(data_sum, 1, dim);
  int* le = largest_elements(data_sum, dim, k-1), i, j, p;

  for (i = 0; i < n; i++) {
    rest = sum(&data[i*dim], dim);

    for (j = 0; j < k-1; j++) {
      p = le[j];
      val = data[i*dim + p];
      rest -= val;
      ans[i*k+j] = val;
    }
    ans[i*k+(k-1)] = rest;
  }

  return ans;
}

int* rg_main(double *data, int n, int k, int dim, int max_iter,
  int *iter, ClusterAssignFunc func) {
    int i, j, d;

    int *assigned = (int *)malloc(sizeof(int) * n);



    // FIND CLUSTER ASSIGNMENTS FOR K = DIM AND VECTOR RATIOS
    double *clusters = (double *)malloc(sizeof(double) * k * dim);
    double *ratio = (double *)malloc(sizeof(double) * n);
    dom_first_pass(data, clusters, assigned, ratio, n, k, dim);

    // for (i = 0; i < n; i++) {
    //   if ( (ratio[i] < 0) || (ratio[i] > 1) ) {
    //     printf("wrong ratio %d: %f\n", i, ratio[i]);
    //   }
    // }

    // SORT INDICES OF CLUSTERS BY ENTROPY
    double *original_results = clusters_entropy(clusters, dim, dim);
    double o_sum = sum(original_results, dim);
    normalize_vector(original_results, o_sum, dim);
    // printf("normalized entropies\n");
    // print_array(original_results, 1, dim);
    int *c_indices = largest_elements(original_results, dim, dim);
    // printf("components ordered by entropy\n");
    // print_int_array(c_indices, 1, dim);
    // SORT INDICES OF VECTORS BY ENTROPY
    int *v_indices = largest_elements(ratio, n, n);

    // PARTITION VECTOR INDICES BY FIRST ASSIGNMENT CLUSTER
    int *vinds_per_cluster = (int *)malloc(sizeof(int) * dim * n);
    int *vecs_per_cluster = (int *)malloc(sizeof(int) * dim);
    first_partition(vinds_per_cluster, vecs_per_cluster, v_indices, assigned,
      n, dim);
    // printf("vecs per component:\n");
    // print_int_array(vecs_per_cluster, 1, dim);

    // FIND ADDITIONAL CLUSTERS PER COMPONENT
    int *additional_clusters = assign_clusters(c_indices,
      original_results, vecs_per_cluster, k, dim);
    // printf("clusters per component:\n");
    // print_int_array(additional_clusters, 1, dim);

    //
    double *cluster_data = (double *)malloc(sizeof(double)*n*dim);
    int *cluster_assignment = (int *)malloc(sizeof(int)*n);
    int first_cluster = 0;

    for (i = 0; i < dim; i++) {
      d = c_indices[i];
      if (vecs_per_cluster[d] > 0) {
        fill_array_int(cluster_assignment, -1, n);
        get_vecs_for_cluster(data, cluster_data, &vinds_per_cluster[d * n], dim,
          vecs_per_cluster[d]);
        // printf("i = %d\tcomponent = %d\tn_vecs=%d\tn_clusters = %d\t",
        //   i, d, vecs_per_cluster[d], additional_clusters[d]);
        // printf("entropy: %f\n", original_results[d]);
        func(cluster_data, cluster_assignment, first_cluster,
          additional_clusters[d], vecs_per_cluster[d], dim);

        for (j = 0; j < vecs_per_cluster[d]; j++) {
          assigned[vinds_per_cluster[d*n + j]] = cluster_assignment[j];
        }
        first_cluster += additional_clusters[d];
      }
    }

    free(clusters);
    free(ratio);
    free(original_results);
    free(c_indices);
    free(v_indices);
    free(vinds_per_cluster);
    free(vecs_per_cluster);
    free(additional_clusters);
    free(cluster_data);
    free(cluster_assignment);
    return assigned;
}

  int* rg_clustering(double *data, int n, int k, int dim, int max_iter,
                       int *iter, ClusterAssignFunc func) {
    iter[0] = 0;

    int *assigned = (int *)malloc(sizeof(int)*n);

    if (k < dim) {
      double *new_data = transform_data(data, n, k, dim);
      // printf("data sum: %f\t", sum(data, n*dim));
      // printf("new data sum: %f\n", sum(new_data, n*k));
      assigned = rg_main(new_data, n, k, k, max_iter, iter, func);
      free(new_data);
    }
    else {
      assigned = rg_main(data, n, k, dim, max_iter, iter, func);
    }

    return assigned;
  }

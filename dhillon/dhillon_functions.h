#ifndef DHILLON_FUNCTIONS
#define DHILLON_FUNCTIONS

#include "../helper/helper_functions.h"

// find cluster of vector v in dominance algorithm
int get_cluster(double *v, int *le, int dim, int k) {
  // calculate sum of v
  double v_sum = sum(v, dim);

  double max = -INFINITY;
  int idx, i, max_idx = -1;

  for (i = 0; i < k-1; i++) {
    idx = le[i];
    v_sum -= v[idx];
    if (v[idx] > max) {
      max = v[idx];
      max_idx = i;
    }
  }

  if (v_sum > max) {
    max_idx = k-1;
  }

  return max_idx;
}

int* dominance_small_k(double *data, int n, int k, int dim) {
  // printf("\nStarting dominance\n");
  int i;
  int l = (k > dim)? dim: k;

  int* assigned = (int *)malloc(sizeof(int) * n);
  double* data_sum = array_sum(data, n, dim);
  int* le = largest_elements(data_sum, dim, l);

  /*
  printf("largest elements:\t");
  for (i = 0; i < l; i++) {
    printf("%d\t", le[i]);
  }
  printf("\n");
  */
  // printf("\ngetting clusters\n");
  for (i = 0; i < n; i++) {
    assigned[i] = get_cluster(&data[i * dim], le, dim, l);
    // printf("%d\t", i);
  }
  // printf("\ngot clusters\n");

  free(data_sum);
  free(le);

  // printf("\nEnding dominance\n");
  return assigned;
}

double kl_div(double *norm_data, double *data_logs, double *cluster_logs,
              int n, int k, int dim) {
  double ans = 0, log;
  int i;
  int v = n * dim;
  int c = k * dim;

  for (i = 0; i < dim; i++) {
    log = data_logs[v + i] - cluster_logs[c + i];
    ans += norm_data[v + i] * log;
  }

  return ans;
}

void initial_di(double *data, double *clusters,
                     int n, int k, int dim) {
  int iter[] = {0};
  int *assigned = dominance_small_k(data, n, k, dim);

  // printf("checking initial assigned:\n");
  // check_assignment(assigned, 0, dim, n);
  // printf("\n\n\n");


  if (k > dim) {

    int *clusters_per_component = (int *)malloc(sizeof(int)*dim);
    fill_array_int(clusters_per_component, k / dim, dim);
    int addtl_clusters = k - int_sum(clusters_per_component, dim);
    int i = 0;

	// TODO: TURN INTO FUNCTION
    while (addtl_clusters > 0) {
      clusters_per_component[i]++;
      i++;
      addtl_clusters--;
    }

    // printf("clusters per component:\n");
    // print_int_array(clusters_per_component, 1, dim);
    // printf("total clusters assigned: %d\n", int_sum(clusters_per_component, dim));

    int *initial_cluster = (int *)malloc(sizeof(int) * dim);
    int *current_cluster = (int *)malloc(sizeof(int) * dim);
    initial_cluster[0] = 0;
    current_cluster[0] = 0;

	// TODO: TURN INTO FUNCTION
    for (i = 1; i < dim; i++) {
      initial_cluster[i] = initial_cluster[i-1] + clusters_per_component[i-1];
      current_cluster[i] = initial_cluster[i];
    }

    // printf("initial cluster for each component:\n");
    // print_int_array(initial_cluster, 1, dim);
    // printf("current cluster for each component:\n");
    // print_int_array(initial_cluster, 1, dim);

    // int *v_indices = (int *)malloc(sizeof(int)*n);
    //
    // for (i = 0; i < n; i++) {
    //   v_indices[i] = i;
    // }

    // int *vecs_per_cluster = get_vecs_per_cluster(v_indices, assigned, n, dim);
    //
    // printf("vecs per cluster\n");
    // print_int_array(vecs_per_cluster, 1, dim);

    int c;

    // for (i = 0; i < n; i++) {
    //   c = assigned[i];
    //   n_clusters = clusters_per_component[c];
    //   i_c = initial_cluster[c];
    //   assigned[i] = i_c + (rand() % n_clusters);
    // }

    for (i = 0; i < n; i++) {
      c = assigned[i];
      assigned[i] = current_cluster[c];
      // if ( (k == 300) && (current_cluster[c] == 288) ) {
      //   printf("element %d asigned to cluster %d\n", i, current_cluster[c]);
      // }
      current_cluster[c]++;
      if (current_cluster[c] >=
          (initial_cluster[c] + clusters_per_component[c]) ) {
        // printf("reached maximum for component %d\n", c);
        current_cluster[c] = initial_cluster[c];
      }
    }

    // printf("checking assigned after first pass:\n");
    // check_assignment(assigned, 0, k, n);
    // printf("\n\n\n");

    // free(v_indices);
    // free(vecs_per_cluster);
    free(clusters_per_component);
    free(initial_cluster);
    free(current_cluster);
  }

  maximization(data, assigned, clusters, n, k, dim);
  // printf("checking assigned after maximization:\n");
  // check_assignment(assigned, 0, k, n);
  // printf("cluster 288\n");
  // print_array(&clusters[288*dim], 1, dim);
  // printf("assigned to cluster 288:\n");
  // int i;
  // for (i = 0; i < n; i++) {
  //   if (assigned[i] == 288) {
  //     printf("%d\t", i);
  //   }
  // }
  // printf("\n\n\n");
  free(assigned);
}

#endif

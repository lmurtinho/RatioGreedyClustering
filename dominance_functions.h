#ifndef DOMINANCE_FUNCTIONS
#define DOMINANCE_FUNCTIONS

#include "helper_functions.h"

// Makes first pass at data for dominance clustering when k >= d
void dom_first_pass(double *data, double *clusters, int *assigned,
                    double *ratio, int n, int k, int dim) {
  int i, j;
  double v_sum, v_max, val;

  // INITIALIZE CLUSTER ARRAY
  fill_array(clusters, 0, k*dim);

  // FOR EACH VECTOR
  for (i = 0; i < n; i++) {
    v_max = -1;
    v_sum = 0;

    // FOR EACH COMPONENT
    for (j = 0; j < dim; j++) {
      val = data[i*dim + j];
      v_sum += val;
      if (val > v_max) {
        assigned[i] = j;
        v_max = val;
      }
    }

    if (v_sum > 0) {
      ratio[i] = v_max / v_sum;
      if (ratio[i] > 1) {
        ratio[i] = 1;
      }
    }
    else {
      ratio[i] = 0;
    }

    add_array(&clusters[assigned[i]*dim], &data[i*dim], dim);

  }

}

// Partition vector indices (sorted by ratios) into the vectors' clusters
void first_partition (int *vinds_per_cluster, int *vecs_per_cluster,
                      int *v_indices, int *assigned, int n, int dim) {

  int i, v, c, p;

  // INITIALIZE ARRAY OF NUMBER OF VECTORS PER CLUSTER
  for (i = 0; i < dim; i++) {
    vecs_per_cluster[i] = 0;
  }

  // INITIALIZE ARRAY OF SORTED VECTOR INDICES PER CLUSTER
  for (i = 0; i < dim * n; i++) {
    vinds_per_cluster[i] = -1;
  }

  for (i=0; i < n; i++) {
    v = v_indices[i]; // v = i-th vector sorted by ratio
    c = assigned[v]; // c = cluster to which v is assigned
    p = vecs_per_cluster[c]; // p = current number of vectors already in c
    vinds_per_cluster[c*n + p] = v; // add v to c
    vecs_per_cluster[c]++; // increment count of vectors in c
  }
}

int* assign_clusters(int *c_indices, double *original_results,
  int *vecs_per_cluster, int k, int dim) {

  int *additional_clusters = (int *)malloc(sizeof(int) * dim);
  int idx, nxt, i, initial_rem_clusters = k, rem_clusters;

  fill_array_int(additional_clusters, 0, dim);
  for (i = 0; i < dim; i++) {
    idx = c_indices[i];
    if (vecs_per_cluster[idx] > 0) {
      additional_clusters[idx]++;
      initial_rem_clusters--;
    }
  }

  rem_clusters = initial_rem_clusters;

  for (i = 0; i < dim; i++) {

    idx = c_indices[i];
    nxt = (int) ( ( initial_rem_clusters * original_results[idx] ) + 1);

    if (nxt > (vecs_per_cluster[idx]) ) {
      nxt = (vecs_per_cluster[idx]);
    }

    if (nxt > rem_clusters) {
      nxt = rem_clusters;
    }

    additional_clusters[idx] += nxt;
    rem_clusters -= nxt;

    if (rem_clusters <= 0) {
      break;
    }
  }

  i = 0;
  while(rem_clusters > 0) {
    idx = c_indices[i];
    additional_clusters[idx]++;
    rem_clusters--;
    i++;
  }

  return additional_clusters;
}

void get_vecs_for_cluster(double *data, double *cluster_data,
	int *vinds_per_cluster, int dim, int vecs) {
	int i, idx;
	for (i = 0; i < vecs; i++) {
		idx = vinds_per_cluster[i];
		copy_to(&data[idx*dim], &cluster_data[i*dim], dim);
	}
}

#endif

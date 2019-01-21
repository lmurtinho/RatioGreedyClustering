#ifndef DOM_4_FUNCTIONS
#define DOM_4_FUNCTIONS

#include "heap_functions.h"

HEAP create_heap_element (double *cluster_vecs, double *heap_vecs,
    double *cluster_ents, double *heap_ents, int i, int n, int dim) {
  HEAP ans;

  ans.id = i;
  ans.ini = i;
  ans.end = i+1;

  ans.ini_vec = &cluster_vecs[i*dim];
  ans.end_vec = &cluster_vecs[(i+1)*dim];
  ans.heap_vec = &heap_vecs[i*dim];
  ans.val = heap_ents[i] - cluster_ents[i] - cluster_ents[i+1];

  ans.prev = (i == 0) ? -1 : i-1;
  ans.next = (i == (n-1) ) ? -1 : i+1;
  return ans;
}

void adjust_heap_element(HEAP *heap, double *heap_ents, double *cluster_ents,
    int el, int n, int dim) {
  copy_to(heap[el].ini_vec, heap[el].heap_vec, dim);
  add_array(heap[el].heap_vec, heap[el].end_vec, dim);
  heap_ents[heap[el].ini] = weighted_entropy(heap[el].heap_vec, dim);
  heap[el].val = heap_ents[heap[el].ini] - cluster_ents[heap[el].ini] -
    cluster_ents[heap[el].end];
  heap_restore(heap, el, n);
}

void cluster_assign_ext (double *cluster_data, int *cluster_assignment,
  int first_cluster, int avail_clusters, int n, int dim) {

  if (avail_clusters == 1) {
    fill_array_int(cluster_assignment, first_cluster, n);
    return;
  }

  // FILE *f = fopen("dominance_greedy.txt", "a");

  int i, j, ini, end, c = n-1;

  HEAP *heap = (HEAP *)malloc(sizeof(HEAP)*(n-1));
  HEAP root;

  double *cluster_vecs = (double *)malloc(sizeof(double)*n*dim);
  double *heap_vecs = (double *)malloc(sizeof(double)*(n-1)*dim);

  double *cluster_ents = (double *)malloc(sizeof(double)*n);
  double *heap_ents = (double *)malloc(sizeof(double)*(n-1));

  // VECTOR ENTROPIES
  for (i = 0; i < n; i++) {
    copy_to(&cluster_data[i*dim], &cluster_vecs[i*dim], dim);
    cluster_ents[i] = weighted_entropy(&cluster_vecs[i*dim], dim);
    cluster_assignment[i] = i;
  }

  for (i = 0; i < c; i++) {

    // CLUSTER ENTROPIES
    copy_to(&cluster_vecs[i*dim], &heap_vecs[i*dim], dim);
    add_array(&heap_vecs[i*dim], &cluster_vecs[(i+1)*dim], dim);
    heap_ents[i] = weighted_entropy(&heap_vecs[i*dim], dim);

    // HEAP ELEMENTS
    heap[i] = create_heap_element(cluster_vecs, heap_vecs, cluster_ents,
      heap_ents, i, c, dim);
  }

  heap_initialize(heap, c);

  // for (i = 0; i < c; i++) {
  //   heap_check_position(heap, i, c);
  // }

  while (c > avail_clusters-1) {

    root = heap_remove_root(heap, c);
    c--;

    // for (i = 0; i < c; i++) {
    //   heap_check_position(heap, i, c);
    // }

    cluster_assignment[root.end] = -1;
    add_array(root.ini_vec, root.end_vec, dim);
    cluster_ents[root.ini] = weighted_entropy(root.ini_vec, dim);

    if (root.next >= 0) {
      adjust_heap_element(heap, heap_ents, cluster_ents, root.next, c, dim);
    }

    if (root.prev >= 0) {
      adjust_heap_element(heap, heap_ents, cluster_ents, root.prev, c, dim);
    }
  }

  cluster_assignment[0] = first_cluster;
  for (i = 1; i < n; i++) {
    if (cluster_assignment[i] == -1) {
      cluster_assignment[i] = cluster_assignment[i-1];
    }
    else {
      cluster_assignment[i] = cluster_assignment[i-1] + 1;
    }
  }

  free(heap);
  free(cluster_vecs);
  free(heap_vecs);
  free(cluster_ents);
  free(heap_ents);
}

#endif

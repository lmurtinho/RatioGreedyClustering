#include <stdlib.h>
#include "helper_functions.h"
#include "printing_functions.h"

void quicksort_indices(int *indices, double *array, int dim) {

  // printf("\n\n\nSorting indices\n");

  if (dim < 2) return;

  int i, s = 0, l = 0, p = 0, comp, pos = rand() % dim;

  // printf("pivot at position %d\n", pos);

  int *smaller = (int*)malloc(sizeof(int)*dim);
  int *larger  = (int*)malloc(sizeof(int)*dim);
  int *pivots  = (int*)malloc(sizeof(int)*dim);

  for (i = 0; i < dim; i++) {

    comp = compare(array, indices[i], indices[pos]);

    if (comp > 0) {
      larger[l] = indices[i];
      l++;
    }
    else if (comp < 0) {
      smaller[s] = indices[i];
      s++;
    }
    else {
      pivots[p] = indices[i];
      p++;
    }

  }

  // printf("smaller:\n");
  // print_int_array(smaller, 1, s);
  // printf("pivots:\n");
  // print_int_array(pivots, 1, p);
  // printf("larger:\n");
  // print_int_array(larger, 1, l);



  quicksort_indices(smaller, array, s);
  quicksort_indices(larger, array, l);
  copy_to_int(smaller, indices, s);
  copy_to_int(pivots, &indices[s], p);
  copy_to_int(larger, &indices[s + p], l);

  free(smaller);
  free(pivots);
  free(larger);

  return;
}

// Returns positions of k largest elements in vector v
int* largest_elements(double *v, int dim, int k) {

  int i, pos;

  int *indices = (int *)malloc(sizeof(int) * dim);
  for (i = 0; i < dim; i++) {
    indices[i] = i;
  }
  quicksort_indices(indices, v, dim);

  reverse_int(indices, dim);

  int *ans = (int *)malloc(sizeof(int) * k);

  for (i = 0; i < k; i++) {
    ans[i] = indices[i];
  }

  free(indices);
  return ans;
}

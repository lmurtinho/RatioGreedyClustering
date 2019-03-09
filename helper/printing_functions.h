#ifndef PRINTING_FUNCTIONS
#define PRINTING_FUNCTIONS

#include <stdio.h>

// Prints double array with n rows and dim columns
void print_array(double *a, int n, int dim) {
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < dim; j++) {
      printf("%.12f\t", a[i * dim + j]);
    }
    printf("\t\t%.12f\n", sum(&a[i * dim], dim));
  }
  return;
}

// Prints int array with n rows and dim columns
void print_int_array(int *a, int n, int dim) {
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < dim; j++) {
      printf("%d\t", a[i * dim + j]);
    }
  }
  printf("\n");
  return;
}

#endif

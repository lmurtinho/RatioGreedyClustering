#include "helper/helper_functions.h"
#include "helper/random_functions.h"
#include "dhillon/dhillon.h"
#include "coresets/coresets.h"

int main (int argc, char *argv[]) {
  char *file_name = concat(argv[1], ".csv");
  int n = atoi(argv[2]), dim = atoi(argv[3]), k = atoi(argv[4]),
    m = atoi(argv[5]), i = atoi(argv[6]), seed = atoi(argv[7]), j = 0,
    iter[] = {0};
  double *data = data_from_csv(file_name, n, dim);
  srand(seed);
  clock_t start = clock(), end;
  for (j = 0; j < i; j++) {
    printf("round number %d\n", j);
    co_clustering(data, n, k, dim, m, 100, iter, argv[1]);
  }
  end = clock();
  printf("Ran in %f seconds\n", (double) (end - start) / CLOCKS_PER_SEC);
  free(file_name);
  free(data);
}

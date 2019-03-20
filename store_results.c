#include "helper/helper_functions.h"
#include "helper/random_functions.h"
#include "dhillon/dhillon.h"
#include "coresets/coresets.h"

int main(int argc, char *argv[]) {

  FILE *f;
  char *file_name;
  int i, n, k, dim, m, n_iter, exp, iter[] = {0};

  if (access("results.csv", F_OK) == -1 ) {
    f = fopen("results.csv", "w");
    fprintf(f, "experiment,dataset,method,n_clusters,n_set,iter,entropy,time\n");
  }
  else {
    f = fopen("results.csv", "a");
  }
  exp = last_experiment() + 1;
  printf("%d\n", exp);
  fclose(f);

  if (strncmp(argv[1], "ng20", 4) == 0) {
    n = 51840;
    dim = 20;
    file_name = "ng20.csv";
  }
  else if (strncmp(argv[1], "rcv1", 4) == 0) {
    n = 170946;
    dim = 103;
    file_name = "rcv1.csv";
  }
  else {
    printf("incorrect data set!\n");
    return -1;
  }

  double *data = data_from_csv(file_name, n, dim);

  k = atoi(argv[2]);
  n_iter = atoi(argv[3]);
  srand(atoi(argv[4]));
  m = atoi(argv[5]);

  printf("Running %d iterations of the coresets method ", n_iter);
  printf("for dataset %s with %d clusters and %d subsamples ", argv[1], k, m);
  printf("using seed %s\n\n", argv[4]);

  for (i = 0; i < n_iter; i++) {
    printf("round number %d, experiment %d\n", i, exp + i);
    iter[0] = 0;
    co_clustering(data, n, k, dim, m, 100, iter, argv[1], exp + i);
  }

  free(data);
}

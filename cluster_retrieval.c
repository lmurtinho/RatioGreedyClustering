#include "helper/helper_functions.h"

int main(int argc, char *argv[]) {
  int k = 2, dim = 20;
  char iter = *argv[2];
  char *name = argv[1];
  double *cluster = cluster_from_csv(name, k, dim, iter);
  print_array(cluster, k, dim);
}

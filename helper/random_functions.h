int random_number(int max_number) {
    return rand() % max_number;
  }

int *rd_clustering(double *data, int n, int k, int dim, int max_iter, int *iter) {
    iter[0] = 0;
    int i, *assigned = (int *)malloc(sizeof(int) * n);
    for (i = 0; i < n; i++) {
      assigned[i] = random_number(k);
    }
    return assigned;
  }

double* create_random_array(int n, int dim, int max) {
    int i;
    double *ans = (double *)malloc(sizeof(double) * n * dim);
    for (i = 0; i < n * dim; i++) {
      ans[i] = rand() % max + 1;
    }
    return ans;
  }

int* create_random_int_array(int n, int dim, int max) {
    int i;
    int *ans = (int *)malloc(sizeof(int) * n * dim);
    for (i = 0; i < n * dim; i++) {
      ans[i] = rand() % max + 1;
    }
    return ans;
  }

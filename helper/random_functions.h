int find_next_number(double *probs, int n) {
  double r = (double)rand() / (double)(RAND_MAX);
  double sum_probs = 0;
  int i;

  for (i = 0; i < n; i++) {
    sum_probs += probs[i];
    if (sum_probs > r) {
      break;
    }
  }
  return i;
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

int* random_sampling(double *values, int n, int k, bool replace) {
  int i, j, *ans = (int *)malloc(sizeof(int) * k);
  double *probs = normalize_vector(values, sum(values, n), n);
  for (i = 0; i < k; i++) {
    ans[i] = find_next_number(probs, n);
    if (!replace) {
      probs[ans[i]] = 0;
      probs = normalize_vector(probs, sum(probs, n), n);
    }
  }
  free(probs);
  return ans;
}

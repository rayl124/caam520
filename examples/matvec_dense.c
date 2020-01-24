#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "time_diff.h"

int main(int argc, char **argv)
{
  const int n = atoi(argv[1]);
  const int reps = atoi(argv[2]);

  double *A = malloc(n*n*sizeof(double));
  double *x = malloc(n*sizeof(double));
  double *y = malloc(n*sizeof(double));

  struct timespec start, end;

  clock_gettime(CLOCK_MONOTONIC, &start);
  for (int r = 0; r < reps; r++) {
    // Compute matrix-vector product.
    for (int i = 0; i < n; i++) {
      y[i] = 0.0;
      for (int j = 0; j < n; j++) {
        y[i] += A[i*n + j]*x[j];
      }
    }
  }
  clock_gettime(CLOCK_MONOTONIC, &end);
  const double delta = time_diff(&end, &start);

  printf("n = %d, time = %fs (%d repetitions)\n", n, delta, reps);

  free(A);
  free(x);
  free(y);
  return 0;
}

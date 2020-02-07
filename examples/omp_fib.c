// omp_fib.c - Compute the Fibonacci sequence using OpenMP tasks.

#include <stdio.h>
#include <omp.h>

int fib(int n)
{
  if (n < 2) return n;

  int l, r;

  printf("fib(%d) on thread #%d\n", n, omp_get_thread_num());

  // Offload the computation of the previous two Fibonacci numbers
  // into tasks.
  // NOTE: This will *not* accelerate the computation, and it will
  //       in fact slow it down.
  //       The purpose of this example is only to demonstrate how
  //       OpenMP tasks can be used to implement recursive and
  //       non-loop-based algorithms with multi-threading.
  #pragma omp task shared(l)
  l = fib(n - 1);

  #pragma omp task shared(r)
  r = fib(n - 2);

  #pragma omp taskwait

  return l + r;
}

int main()
{
  const int n = 8;
  int result;

  #pragma omp parallel
  {
    #pragma omp single
    result = fib(n);
  }

  printf("Final result: fib(%d) = %d\n", n, result);
  return 0;
}

#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <omp.h>


// Scan A array and write result into prefix_sum array;
// use long data type to avoid overflow
void scan_seq(long* prefix_sum, const long* A, long n) {
  if (n == 0) return;
  prefix_sum[0] = 0;
  for (long i = 1; i < n; i++) {
    prefix_sum[i] = prefix_sum[i-1] + A[i-1];
  }
}

void scan_omp(long* prefix_sum, const long* A, long n) {
  // TODO: implement multi-threaded OpenMP scan
  const int N_THREAD = 8;
  prefix_sum[0] = 0;
  long offsets[N_THREAD + 1];
  long len = n / N_THREAD;
  printf("start\n");
  #pragma omp parallel num_threads(N_THREAD)
  {
    int id = omp_get_thread_num();
    long startIdx = id * len + 1;
    long endIdx = id == N_THREAD - 1 ? n - 1 : startIdx + len - 1;
    prefix_sum[startIdx] = A[startIdx - 1];
    for(long i = startIdx + 1; i <= endIdx; i++){
      prefix_sum[i] = prefix_sum[i - 1] + A[i - 1];
    }
    offsets[id + 1] = prefix_sum[endIdx];
    printf("id :  %d\n", id);
  }
  #pragma omp barrier
  offsets[0] = 0;
  for(int i = 1; i <= N_THREAD; i++){
    offsets[i] += offsets[i - 1];
    printf("offset %d is %ld", i, offsets[i]);
  }

  int curOff = 0;
  int cnt = 0;
  for(int i = 1; i < n ; i++){
    cnt++;
    if(cnt > len){
      curOff++;
      cnt = 1;
      printf("curOff is %d", curOff);
    }
    prefix_sum[i] += offsets[curOff];
  }

}

int main() {
  long N = 100000000;
  long* A = (long*) malloc(N * sizeof(long));
  long* B0 = (long*) malloc(N * sizeof(long));
  long* B1 = (long*) malloc(N * sizeof(long));
  for (long i = 0; i < N; i++) A[i] = rand();

  double tt = omp_get_wtime();
  scan_seq(B0, A, N);
  printf("sequential-scan = %fs\n", omp_get_wtime() - tt);

  tt = omp_get_wtime();
  scan_omp(B1, A, N);
  printf("parallel-scan   = %fs\n", omp_get_wtime() - tt);

  long err = 0;
  for (long i = 0; i < N; i++){
    err = std::max(err, std::abs(B0[i] - B1[i]));
    if(std::abs(B0[i] - B1[i]) > 1){
      printf("error idx : %ld\n", i);
      printf("correct number is : %ld", B0[i]);
      printf("my number is : %ld", B1[i]);
    }

  }

  printf("error = %ld\n", err);

  for(long i = N - 1; i >= N - 10; i--){
    printf("B0 %ld\n", B0[i]);
  }

  for(long i = N - 1; i >= N - 10; i--){
    printf("B1 %ld\n", B1[i]);
  }

  free(A);
  free(B0);
  free(B1);
  return 0;
}

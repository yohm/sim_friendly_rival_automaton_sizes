#include <iostream>
#include <cstdio>
#include <mpi.h>

int GetAutomatonSize(const char strategy[64]) {
  // [TODO] implement me
  return 0;
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  std::cout << "Hello, MPI World!" << std::endl;

  int my_rank, total_size;
  MPI_Comm_size(MPI_COMM_WORLD, &total_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  std::cout << "my_rank, total_size: " << my_rank << ", " << total_size << std::endl;

  char filename[30];
  sprintf(filename, "input_file_%04d", my_rank);
  FILE* fp = fopen(filename, "r");
  if (fp == NULL) { MPI_Abort(MPI_COMM_WORLD, 1); }

  long automaton_sizes[65] = {0};
  char* line = NULL;
  size_t len = 0;
  while ((getline(&line, &len, fp)) != -1) {
    printf("%s", line);
    int atm_size = GetAutomatonSize(line);
    automaton_sizes[atm_size]++;
  }

  long results[65] = {0};
  MPI_Reduce(automaton_sizes, results, 65, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

  if (my_rank == 0) {
    for (int i = 0; i < 65; i++) {
      std::cout << i << ' ' << results[i] << std::endl;
    }
  }

  MPI_Finalize();
  return 0;
}

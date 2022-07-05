#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <chrono>
#include <mpi.h>
#include "icecream.hpp"
#include "StrategyM3.hpp"

int GetAutomatonSize(const std::string& strategy_str) {
  StrategyM3 strategy(strategy_str.c_str());
  return strategy.MinimizeDFA(true).to_map().size();
}

using automaton_sizes_t = std::array<size_t, 65>;

void ExpandWildcard(const std::string& line, automaton_sizes_t& automaton_sizes) {
  for (int i = 0; i < line.size(); i++) {
    if (line[i] == '*') {
      std::string temp = line;
      temp[i] = 'c';
      ExpandWildcard(temp, automaton_sizes);
      temp[i] = 'd';
      ExpandWildcard(temp, automaton_sizes);
      return;
    }
  }
  int size = GetAutomatonSize(line);
  automaton_sizes[size]++;
}

automaton_sizes_t AutomatonSizes(const std::string& line) {
  automaton_sizes_t sizes;
  sizes.fill(0);

  ExpandWildcard(line, sizes);
  return sizes;
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  std::cout << "Hello, MPI World!" << std::endl;

  int my_rank, total_size;
  MPI_Comm_size(MPI_COMM_WORLD, &total_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  std::cout << "my_rank, total_size: " << my_rank << ", " << total_size << std::endl;

  std::ifstream fin("../friendly_rival_10k");
  std::vector<std::string> lines;
  while (fin) {
    std::string s;
    fin >> s;
    if (fin) {
      lines.push_back(s);
    }
  }
  // IC(lines);

  std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

  automaton_sizes_t sizes;
  sizes.fill(0);
  for (const auto& line : lines) {
    automaton_sizes_t r = AutomatonSizes(line);
    // IC(r);
    for (size_t i = 0; i < r.size(); ++i) {
      sizes[i] += r[i];
    }
  }

  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;

  size_t sum = 0;
  double total = 0.0;
  if (my_rank == 0) {
    for (int i = 0; i < 65; i++) {
      std::cout << i << ' ' << sizes[i] << std::endl;
      sum += sizes[i];
      total += i * sizes[i];
    }
    IC(sum, total/sum, elapsed_seconds.count(), elapsed_seconds.count()/sum);
  }

  MPI_Finalize();
  return 0;
}

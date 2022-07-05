#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <chrono>
#include <mpi.h>
#include "icecream.hpp"
#include "StrategyM3.hpp"


template <int N, int Z>  // N: number of states, Z: number of alphabets
class DFA {
  public:
  DFA(std::array<size_t,N> f, const std::array<std::array<size_t,Z>, N>& _delta) :
    F(std::move(f)), delta(_delta) {
    // validity check
    for (auto i: F) { assert(i == 0 || i == 1); }
    for (const auto& a: delta) {
      for (auto i: a) { assert(i < N); }
    }
  };
  std::array<size_t,N> F;  // set of final states
  const std::array< std::array<size_t,Z>,N >& delta; // state-transition function.
  // represented by NxZ two-dim vector. delta[q][s] = transition from state-q by input s

  using partition_t = std::array<int, N>;

  bool Equivalent(size_t i, size_t j, const partition_t& partition_0) {
    for (int z = 0; z < Z; z++) {  // for input z
      int dest_i = delta[i][z];
      int dest_j = delta[j][z];
      if (partition_0[dest_i] != partition_0[dest_j]) return false;
    }
    return true;
  }

  partition_t Minimize() {
    partition_t partition_0;
    partition_0.fill(-1);

    // initialize grouping by the final state
    int accept_idx = -1, reject_idx = -1;
    for (int i = 0; i < N; i++) {
      if (F[i]) {
        if (accept_idx == -1) {
          accept_idx = i;
        }
        partition_0[i] = accept_idx;
      } else {
        if (reject_idx == -1) {
          reject_idx = i;
        }
        partition_0[i] = reject_idx;
      }
    }

    while (true) {
      partition_t partition;
      partition.fill(-1);
      int i = 0;
      while (i < N) {
        partition[i] = i;
        int i_next = N;
        for (int j = i+1; j < N; j++) {
          if (partition[j] >= 0) continue;  // skip if j is already merged
          if (partition_0[i] == partition_0[j] && Equivalent(i, j, partition_0)) {
            partition[j] = i;   // merge i & j
          }
          else if (i_next == N) { i_next = j; }  // keep the first unmerged node
        }
        i = i_next;
      }
      if (partition_0 == partition) break;
      partition_0 = partition;
    }

    return partition_0;
  }

  size_t MinimizedSize() {
    partition_t partition_0 = Minimize();
    size_t size = 0;
    for (int i = 0; i < N; i++) {
      if (partition_0[i] == i) size++;
    }
    return size;
  }

  std::map<size_t,std::vector<size_t>> PartitionToMap(const partition_t& partition_0) {
    std::map<size_t,std::vector<size_t>> m;
    for (int i = 0; i < N; i++) {
      m[partition_0[i]].emplace_back(i);
    }
    return std::move(m);
  }
};

std::array<std::array<size_t,4>, 64> ConstructDelta() {
  std::array<std::array<size_t,4>, 64> delta;
  using bt = std::bitset<6>;
  for (int i = 0; i < 64; i++) {
    const bt b = bt(i) << 1;
    const bt n0 = b & bt(0b110110) | bt(0b000000);
    const bt n1 = b & bt(0b110110) | bt(0b000001);
    const bt n2 = b & bt(0b110110) | bt(0b001000);
    const bt n3 = b & bt(0b110110) | bt(0b001001);
    delta[i] = {
      n0.to_ulong(),
      n1.to_ulong(),
      n2.to_ulong(),
      n3.to_ulong()
    };
  }
  return delta;
}
const std::array<std::array<size_t,4>, 64> g_delta = ConstructDelta();

int MinimizedDFASize(const char str[64]) {
  std::array<size_t,64> F;
  for (size_t i = 0; i < 64; i++) {
    F[i] = (str[i] == 'c') ? 1 : 0;
  }
  DFA<64,4> dfa(F, g_delta);
  return dfa.MinimizedSize();
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
  int size = MinimizedDFASize(line.c_str());
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

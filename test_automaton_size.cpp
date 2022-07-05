#include <iostream>
#include <vector>
#include <array>
#include <bitset>
#include <set>
#include <chrono>
#include <cassert>
#include "icecream.hpp"
#include "StrategyM3.hpp"


template <int N, int Z>  // N: number of states, Z: number of alphabets
class DFA {
  public:
  DFA(std::set<size_t> f, std::array<std::array<size_t,Z>, N> _delta) :
  F(std::move(f)), delta(std::move(_delta)) {
    // validity check
    for (auto i: F) { assert(i < N); }
    for (const auto& a: delta) {
      for (auto i: a) { assert(i < N); }
    }
  };
  std::set<size_t> F;  // set of final states
  std::array< std::array<size_t,Z>,N > delta; // state-transition function.
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
    int accept_idx = *F.begin(), reject_idx = -1;
    for (int i: F) {
      partition_0[i] = accept_idx;
    }
    for (int i = 0; i < N; i++) {
      if (partition_0[i] < 0) {
        if (reject_idx < 0) reject_idx = i;
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

  std::map<size_t,std::vector<size_t>> PartitionToMap(const partition_t& partition_0) {
    std::map<size_t,std::vector<size_t>> m;
    for (int i = 0; i < N; i++) {
      m[partition_0[i]].emplace_back(i);
    }
    return std::move(m);
  }
};


void ConvertDFA(const char str[64]) {
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

  std::set<size_t> F;
  for (size_t i = 0; i < 64; i++) {
    if (str[i] == 'c') {
      F.insert(i);
    }
  }
  DFA<64,4> dfa(F, delta);
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  for (int i = 0; i < 10000; i++) {
    dfa.Minimize();
  }
  end = std::chrono::system_clock::now();
  std::cout << "dfa.Minimize(): " << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << "ms" << std::endl;

  //auto partition = dfa.Minimize();
  // auto m = dfa.PartitionToMap(partition);
  // IC(m, m.size());
}


int main(int argc, char* argv[]) {

  const std::string line = "cddcccddcdcdddddcccddddcdddcddcdcccccdddddcddcddddccdddddcdccccd";

  StrategyM3 str(line.c_str());

  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  for (int i = 0; i < 10000; i++) {
    auto dfa = str.MinimizeDFA(true);
  }
  end = std::chrono::system_clock::now();
  std::cout << "MinimizeDFA: " << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << "ms" << std::endl;
  // auto dfa = str.MinimizeDFA(true);
  // IC(dfa.to_map(), dfa.to_map().size());

  ConvertDFA(line.c_str());

  return 0;
}
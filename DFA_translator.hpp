#include <iostream>
#include <array>
#include <map>

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

  std::map<size_t,std::vector<size_t>> MinimizedPartitionMap() {
    partition_t partition = Minimize();
    std::map<size_t,std::vector<size_t>> m;
    for (int i = 0; i < N; i++) {
      m[partition[i]].emplace_back(i);
    }
    return std::move(m);
  }
};


namespace DFA_translator {
  std::array<std::array<size_t,4>, 64> ConstructDeltaForFullDFA() {
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

  constexpr std::array<std::array<size_t,4>, 64> delta_for_full_DFA = {{
                                                                         {0, 1, 8, 9},
                                                                         {2, 3, 10, 11},
                                                                         {4, 5, 12, 13},
                                                                         {6, 7, 14, 15},
                                                                         {0, 1, 8, 9},
                                                                         {2, 3, 10, 11},
                                                                         {4, 5, 12, 13},
                                                                         {6, 7, 14, 15},
                                                                         {16, 17, 24, 25},
                                                                         {18, 19, 26, 27},
                                                                         {20, 21, 28, 29},
                                                                         {22, 23, 30, 31},
                                                                         {16, 17, 24, 25},
                                                                         {18, 19, 26, 27},
                                                                         {20, 21, 28, 29},
                                                                         {22, 23, 30, 31},
                                                                         {32, 33, 40, 41},
                                                                         {34, 35, 42, 43},
                                                                         {36, 37, 44, 45},
                                                                         {38, 39, 46, 47},
                                                                         {32, 33, 40, 41},
                                                                         {34, 35, 42, 43},
                                                                         {36, 37, 44, 45},
                                                                         {38, 39, 46, 47},
                                                                         {48, 49, 56, 57},
                                                                         {50, 51, 58, 59},
                                                                         {52, 53, 60, 61},
                                                                         {54, 55, 62, 63},
                                                                         {48, 49, 56, 57},
                                                                         {50, 51, 58, 59},
                                                                         {52, 53, 60, 61},
                                                                         {54, 55, 62, 63},
                                                                         {0, 1, 8, 9},
                                                                         {2, 3, 10, 11},
                                                                         {4, 5, 12, 13},
                                                                         {6, 7, 14, 15},
                                                                         {0, 1, 8, 9},
                                                                         {2, 3, 10, 11},
                                                                         {4, 5, 12, 13},
                                                                         {6, 7, 14, 15},
                                                                         {16, 17, 24, 25},
                                                                         {18, 19, 26, 27},
                                                                         {20, 21, 28, 29},
                                                                         {22, 23, 30, 31},
                                                                         {16, 17, 24, 25},
                                                                         {18, 19, 26, 27},
                                                                         {20, 21, 28, 29},
                                                                         {22, 23, 30, 31},
                                                                         {32, 33, 40, 41},
                                                                         {34, 35, 42, 43},
                                                                         {36, 37, 44, 45},
                                                                         {38, 39, 46, 47},
                                                                         {32, 33, 40, 41},
                                                                         {34, 35, 42, 43},
                                                                         {36, 37, 44, 45},
                                                                         {38, 39, 46, 47},
                                                                         {48, 49, 56, 57},
                                                                         {50, 51, 58, 59},
                                                                         {52, 53, 60, 61},
                                                                         {54, 55, 62, 63},
                                                                         {48, 49, 56, 57},
                                                                         {50, 51, 58, 59},
                                                                         {52, 53, 60, 61},
                                                                         {54, 55, 62, 63}
                                                                       }};

  std::array<std::array<size_t,2>, 64> ConstructDeltaForSimplifiedDFA(const char str[64]) {
    std::array<std::array<size_t,2>, 64> delta;
    using bt = std::bitset<6>;
    for (int i = 0; i < 64; i++) {
      const bt b = bt(i) << 1;
      const bt n0 = b & bt(0b110110) | bt(0b000000);
      const bt n1 = b & bt(0b110110) | bt(0b000001);
      const bt n2 = b & bt(0b110110) | bt(0b001000);
      const bt n3 = b & bt(0b110110) | bt(0b001001);
      if (str[i] == 'c') {
        delta[i] = {
          n0.to_ulong(),
          n1.to_ulong()
        };
      }
      else {
        delta[i] = {
          n2.to_ulong(),
          n3.to_ulong()
        };
      }
    }
    return delta;
  }

  int MinimizedDFASize(const char str[64]) {
    std::array<size_t,64> F;
    for (size_t i = 0; i < 64; i++) {
      F[i] = (str[i] == 'c') ? 1 : 0;
    }
    DFA<64,4> dfa(F, delta_for_full_DFA);
    return dfa.MinimizedSize();
  }

  int MinimizedDFASizeSimple(const char str[64]) {
    std::array<size_t,64> F;
    for (size_t i = 0; i < 64; i++) {
      F[i] = (str[i] == 'c') ? 1 : 0;
    }
    DFA<64,2> dfa(F, ConstructDeltaForSimplifiedDFA(str));
    return dfa.MinimizedSize();
  }

  std::map<size_t,std::vector<size_t>> MinimizedPartition(const char str[64]) {
    std::array<size_t,64> F;
    for (size_t i = 0; i < 64; i++) {
      F[i] = (str[i] == 'c') ? 1 : 0;
    }
    DFA<64,4> dfa(F, delta_for_full_DFA);
    return dfa.MinimizedPartitionMap();
  }

  std::map<size_t,std::vector<size_t>> MinimizedPartitionSimple(const char str[64]) {
    std::array<size_t,64> F;
    for (size_t i = 0; i < 64; i++) {
      F[i] = (str[i] == 'c') ? 1 : 0;
    }
    DFA<64,2> dfa(F, ConstructDeltaForSimplifiedDFA(str));
    return dfa.MinimizedPartitionMap();
  }
}

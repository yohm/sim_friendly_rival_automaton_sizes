#include <iostream>
#include <array>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <numeric>
#include <algorithm>


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
      const bt b = bt(i & 0b011011) << 1;
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
    auto delta = ConstructDeltaForSimplifiedDFA(str);
    DFA<64,2> dfa(F, delta);
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

  class SimpAutomGraph {
    public:
    SimpAutomGraph() {};
    SimpAutomGraph(const char str[64]) {
      const std::map<size_t,std::vector<size_t>> autom = MinimizedPartitionSimple(str);

      auto find_root = [&autom](size_t i) -> size_t {
        for (auto& [k,v] : autom) {
          if (std::find(v.begin(), v.end(), i) != v.end()) {
            return k;
          }
        }
        return 100ul; // must not happen
      };
      auto next_states = [&find_root, &str](size_t n) -> std::vector<size_t> {
        size_t nb = (n & 0b011011) << 1;
        if (str[n] == 'c') {
          return {find_root(nb), find_root(nb | 0b000001)};
        } else {
          return {find_root(nb | 0b001000), find_root(nb | 0b001001)};
        }
      };

      auto add_node = [this,&find_root,&next_states,&str](size_t n) ->size_t {
        size_t node_idx = find_root(n);
        nodes.push_back( std::make_pair(node_idx, str[node_idx]));
        edges[node_idx] = next_states(node_idx);
        return node_idx;
      };

      // add node 0 & 63
      size_t root_0 = add_node(0);
      edges[root_0].push_back(find_root(8));  // additionally, append 0->8

      size_t root_63 = add_node(63);
      edges[root_63].push_back(find_root(55));  // additionally, append 63->55

      for (auto& [k,v] : autom) {
        if (k == root_0 || k == root_63) {
          continue;
        }
        add_node(k);
      }

      // calculate path from 8 & 55
      std::vector<size_t> path_8 = _TracePath(str, 8);
      for (auto& n : path_8) {
        size_t n_co = ((n & 0b000111) << 3) | ((n & 0b111000) >> 3);
        auto p = std::make_pair(find_root(n), find_root(n_co));
        bool found = (std::find(path_from_8.begin(), path_from_8.end(), p) != path_from_8.end());
        path_from_8.push_back(p);
        if (found) break;
      }
      std::vector<size_t> path_55 = _TracePath(str, 55);
      for (auto& n : path_55) {
        size_t n_co = ((n & 0b000111) << 3) | ((n & 0b111000) >> 3);
        auto p = std::make_pair(find_root(n), find_root(n_co));
        bool found = (std::find(path_from_55.begin(), path_from_55.end(), p) != path_from_55.end());
        path_from_55.push_back(p);
        if (found) break;
      }
    }

    std::vector<size_t> _TracePath(const char str[64], size_t init) {
      std::vector<size_t> path;
      size_t node_idx = init;
      auto NextITGState = [&str](size_t n) -> size_t {
        size_t n_co = ((n & 0b000111) << 3) | ((n & 0b111000) >> 3); // state for the co-player
        size_t nb = (n & 0b011011) << 1;
        if (str[n] == 'c' && str[n_co] == 'c') {
          return nb;
        } else if (str[n] == 'c' && str[n_co] == 'd') {
          return nb | 0b000001;
        } else if (str[n] == 'd' && str[n_co] == 'c') {
          return nb | 0b001000;
        } else {
          return nb | 0b001001;
        }
      };
      // if node_idx is not included in path, add it
      path.push_back(node_idx);
      while (true) {
        node_idx = NextITGState(node_idx);
        auto it = std::find(path.begin(), path.end(), node_idx);
        if (it == path.end()) {
          path.push_back(node_idx);
        } else {
          break;
        }
      }
      return path;
    }

    std::vector<std::pair<size_t,char>> nodes; // nodes[i] = (node_idx, node_label)
    std::map<size_t,std::vector<size_t>> edges;
    std::vector<std::pair<size_t,size_t> > path_from_8, path_from_55;

    std::string ToString() const {
      std::stringstream ss;
      for (size_t i = 0; i < nodes.size(); i++) {
        ss << nodes[i].first << "," << nodes[i].second;
        for (size_t j : edges.at(nodes[i].first)) {
          ss << "," << j;
        }
        ss <<";";
      }
      ss << ' ';
      for (auto& [k,v] : path_from_8) {
        ss << k << "," << v << ";";
      }
      ss << ' ';
      for (auto& [k,v] : path_from_55) {
        ss << k << "," << v << ";";
      }
      return ss.str();
    }

    std::string ToRegularizedString() const {
      // make vector from 0 to nodes.size()
      std::vector<size_t> v(nodes.size());
      std::iota(v.begin(), v.end(), 0);

      std::string s = "999";
      do {
        SimpAutomGraph g = *this;
        g.NodeRemap(v);
        std::string s2 = g.ToString();
        if (s2 < s) {
          s = s2;
        }
      } while (std::next_permutation(v.begin()+2, v.end()));
      return s;
    }

    void NodeRemap(const std::vector<size_t>& remap) {
      if (remap.size() != nodes.size()) {
        throw std::runtime_error("NodeRemap: remap size does not match");
      }
      std::map<size_t,size_t> m;
      decltype(nodes) new_nodes;
      for (size_t i = 0; i < nodes.size(); i++) {
        m[nodes[i].first] = remap[i];
        new_nodes.push_back(std::make_pair(remap[i], nodes[i].second));
      }

      nodes = new_nodes;
      decltype(edges) new_edges;
      for (auto& [k,v] : edges) {
        new_edges[m[k]] = std::vector<size_t>();
        for (size_t& i : v) {
          new_edges[m[k]].push_back(m[i]);
        }
      }
      edges = new_edges;

      decltype(path_from_8) new_path_from_8;
      for (auto& [k,v] : path_from_8) {
        new_path_from_8.push_back(std::make_pair(m[k], m[v]));
      }
      path_from_8 = new_path_from_8;

      decltype(path_from_55) new_path_from_55;
      for (auto& [k,v] : path_from_55) {
        new_path_from_55.push_back(std::make_pair(m[k], m[v]));
      }
      path_from_55 = new_path_from_55;

      // sort nodes by index
      std::sort(nodes.begin(), nodes.end(),
        [](const std::pair<size_t,char>& a, const std::pair<size_t,char>& b) {
          return a.first < b.first;
        }
      );
    }
  };

  std::string SerializeSimpleAutom(const char str[64]) {
    SimpAutomGraph g(str);
    return g.ToRegularizedString();
  }
}

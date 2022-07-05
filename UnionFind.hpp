//
// Created by Yohsuke Murase on 2020/05/22.
//

#ifndef CPP__UNIONFIND_HPP_
#define CPP__UNIONFIND_HPP_

#include <cstdlib>
#include <vector>
#include <set>
#include <map>

class UnionFind {
 public:
  explicit UnionFind(size_t n) : N(n), parent(n) {
    for (size_t i = 0; i < n; i++) { parent[i] = i; }
  }
  size_t org_size() const { return N; }
  size_t root(size_t i) {
    if (parent[i] != i) {
      size_t r = root(parent[i]);
      parent[i] = r;
    }
    return parent[i];
  }
  bool is_root(size_t i) const {
    return (parent[i] == i);
  }
  bool merge(size_t i, size_t j) {
    size_t ri = root(i);
    size_t rj = root(j);
    if (ri == rj) return false;
    else if (ri > rj) { parent[ri] = rj; }
    else if (ri < rj) { parent[rj] = ri; }
    return true;
  }
  std::map<size_t, std::vector<size_t> > to_map() {
    std::map<size_t, std::vector<size_t> > m;
    for (size_t i = 0; i < parent.size(); i++) {
      size_t r = root(i);
      m[r].emplace_back(i);
    }
    return std::move(m);
  }
  bool operator==(const UnionFind& rhs) const {
    return (N == rhs.N && parent == rhs.parent);
  }
 private:
  size_t N;
  std::vector<size_t> parent;
};

#endif //CPP__UNIONFIND_HPP_
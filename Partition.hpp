//
// Created by Yohsuke Murase on 2020/06/17.
//

#ifndef CPP__PARTITION_HPP_
#define CPP__PARTITION_HPP_

#include <vector>
#include <array>
#include <set>
#include <map>

class Partition {
 public:
  explicit Partition(size_t n) : N(n) {
    groups[0] = std::set<size_t>();
    for (size_t i = 0; i < N; i++) {
      groups[0].insert(i);
    }
  }
  size_t size() const { return groups.size(); }
  size_t org_size() const { return N; }
  // split set i into sub & si-sub, returns ids of divided sets sorted by their sizes
  std::array<size_t,2> split(size_t i, const std::set<size_t> &sub) {
    std::set<size_t> sub2;
    std::set_difference(groups[i].cbegin(), groups[i].cend(), sub.cbegin(), sub.cend(), std::inserter(sub2, sub2.end()));
    if(sub.size() + sub2.size() != groups[i].size()) { throw std::runtime_error("unknown member is given"); }
    if (sub.empty() || sub2.empty()) { throw std::runtime_error("cannot divide"); }
    groups.erase(i);
    size_t i1 = *sub.begin();
    groups[i1] = sub;
    size_t i2 = *sub2.begin();
    groups[i2] = sub2;
    if (sub.size() <= sub2.size()) { return {i1, i2}; }
    else { return {i2, i1}; }
  }
  const std::set<size_t>& group(size_t i) const { return groups.at(i); }
  const size_t group_id(size_t si) const {
    for (const auto& kv: groups) {
      if (kv.second.find(si) != kv.second.end()) {
        return kv.first;
      }
    }
    std::cerr << "unknown member" << std::endl;
    throw std::runtime_error("unknown member");
  }
  std::vector<size_t> group_ids() const {
    std::vector<size_t> ans;
    for (const auto& it : groups) {
      ans.emplace_back(it.first);
    }
    return ans;
  }
  const std::map<size_t, std::set<size_t>>& to_map() const { return groups; }
  std::map<size_t, std::vector<size_t>> to_vec_map() const {
    std::map<size_t, std::vector<size_t>> m;
    for (const auto& kv: groups) {
      m[kv.first] = std::vector<size_t>(kv.second.begin(), kv.second.end());
    }
    return m;
  }
  friend std::ostream &operator<<(std::ostream &os, const Partition &p) {
    for (const auto &kv: p.groups) {
      os << kv.first << " => " << kv.second.size() << " [\n  ";
      for (const auto &x: kv.second) {
        os << x << ", ";
      }
      os << "]," << std::endl;
    }
    return os;
  };
 private:
  size_t N;
  std::map<size_t, std::set<size_t> > groups;
  std::vector<size_t> parent;
};


#endif //CPP__PARTITION_HPP_
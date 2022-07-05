#include <iostream>
#include <set>
#include <map>
#include "StrategyM3.hpp"

StrategyM3::StrategyM3(const std::array<Action, 64> &acts) : actions(acts) {}

StrategyM3::StrategyM3(const char acts[64]) {
  for (size_t i = 0; i < 64; i++) {
    actions[i] = C2A(acts[i]);
  }
}

StrategyM3::StrategyM3(uint64_t acts) {
  for(size_t i = 0; i <64; i++) {
    if( (acts & (1ul << i)) == 0 ) { actions[i] = C; }
    else { actions[i] = D; }
  }
}

std::string StrategyM3::Name() const {
  const static std::map<uint64_t,std::string> m = {
    {ALLC().ID(), "ALLC"},
    {ALLD().ID(), "ALLD"},
    {TFT().ID(), "TFT"},
    {WSLS().ID(), "WSLS"},
    {GRIM().ID(), "GRIM"},
    {TF2T().ID(), "TF2T"},
    {TFT_ATFT().ID(), "TFT-ATFT"},
    {CAPRI().ID(), "CAPRI"},
    {AON(2).ID(), "AON2"},
    {AON(3).ID(), "AON3"}
  };

  if (m.find(ID()) != m.end()) {
    return m.at(ID());
  }
  else {
    return "";
  }
}

StrategyM3 StrategyM3::ConstructFromName(const std::string& name) {
  const static std::map<std::string,StrategyM3> m = {
    {"ALLC", ALLC()},
    {"ALLD", ALLD()},
    {"TFT", TFT()},
    {"WSLS", WSLS()},
    {"GRIM", GRIM()},
    {"TF2T", TF2T()},
    {"TFT-ATFT", TFT_ATFT()},
    {"CAPRI", CAPRI()},
    {"AON2", AON(2)},
    {"AON3", AON(3)}
  };

  if (m.find(name) != m.end()) {
    return m.at(name);
  }
  else {
    std::cerr << "Error: unknown strategy " << name << std::endl;
    std::cerr << "  supported strategies are [";
    for (const auto& kv: m) { std::cerr << kv.first << ", "; }
    std::cerr << "]" << std::endl;
    throw std::runtime_error("unknown strategy name");
  }
}

void StrategyM3::Inspect(std::ostream& out) const {
  out << "ID: " << ID() << "\n"
      << "String: " << ToString() << "\n"
      << "Name: " << Name() << "\n";

  for (size_t i = 0; i < 64; i++) {
    out << StateM3(i) << "| " << actions[i] << "    ";
    if (i % 8 == 7) out << "\n";
  }
  bool is_efficient = IsEfficientTopo();
  bool is_defensible = IsDefensible();
  out << "IsEfficientTopo: " << is_efficient << "\n"
      << "IsEfficient: " << IsEfficient() << "\n"
      << "IsDefensible: " << is_defensible << "\n"
      << "IsDefenisbleDFA: " << IsDefensibleDFA() << "\n"
      << "IsDistinguishable: " << IsDistinguishableTopo() << "\n"
      << "CooperationLevel: " << CooperationLevel() << "\n";
  if (is_efficient) {
    std::vector<int> path1 = {1};
    while (true) {
      int next = NextITGState( *path1.rbegin() );
      bool found = (std::find(path1.begin(), path1.end(), next) != path1.end());
      path1.push_back(next);
      if (found) break;
    }
    out << "Path in g(S,S) from (ccc,ccd): ";
    for (int n: path1) { out << StateM3(n) << " -> "; }
    out << "\n";

    std::vector<int> path62 = {62};
    while (true) {
      int next = NextITGState( *path62.rbegin() );
      bool found = (std::find(path62.begin(), path62.end(), next) != path62.end());
      path62.push_back(next);
      if (found) break;
    }
    out << "Path in g(S,S) from (ddd,ddc): ";
    for (int n: path62) { out << StateM3(n) << " -> "; }
    out << "\n";
  }
  if (!is_defensible) {
    auto n_cycles = ShortestNegativeCycles();
    out << "Shortest negative cycles in g(S):\n";
    for (auto cycle: n_cycles) {
      out << "  ";
      for (size_t n: cycle) {
        out << StateM3(n) << ' ';
      }
      out << "\n";
    }
    out << "\n";
  }
  auto p0 = MinimizeDFAHopcroft(false), p1 = MinimizeDFAHopcroft(true);
  out << "DFT minimized states: (" << p0.size() << " / " << p1.size() << ")\n";
  out << "--- automaton without noise:\n" << p0
      << "--- automaton with noise:\n" << p1;
}

std::vector<StateM3> StrategyM3::NextPossibleStates(StateM3 current) const {
  std::vector<StateM3> next_states;
  Action act_a = ActionAt(current);
  next_states.push_back(current.NextState(act_a, C));
  next_states.push_back(current.NextState(act_a, D));
  return std::move(next_states);
}

std::ostream &operator<<(std::ostream &os, const StrategyM3 &strategy) {
  for (size_t i = 0; i < 64; i++) {
    os << strategy.actions[i] << '|' << StateM3(i) << "  ";
    if (i % 8 == 7) { os << std::endl; }
  }
  return os;
}

std::string StrategyM3::ToString() const {
  char c[65];
  for (size_t i = 0; i < 64; i++) {
    c[i] = A2C(actions[i]);
  }
  c[64] = '\0';
  return std::string(c);
}

inline int8_t MIN(int8_t a, int8_t b) { return (a < b) ? a : b; }

bool StrategyM3::IsDefensible() const {
  const size_t N = 64;

  d_matrix_t d;

  // construct adjacency matrix
  const int INF = 32; // 32 is large enough since the path length is between -16 to 16.
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      d[i][j] = INF;
    }
  }

  for (size_t i = 0; i < N; i++) {
    StateM3 si(i);
    std::vector<StateM3> sjs = NextPossibleStates(si);
    for (auto sj: sjs) {
      size_t j = sj.ID();
      d[i][j] = si.RelativePayoff();
    }
    if (d[i][i] < 0) { return false; }
  }

  for (size_t k = 0; k < N; k++) {
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        d[i][j] = MIN(d[i][j], d[i][k] + d[k][j]);
      }
      if (d[i][i] < 0) { return false; }
    }
  }
  return true;
}

std::vector<std::vector<size_t>> StrategyM3::ShortestNegativeCycles() const {
  const size_t N = 64;
  std::vector<std::map<size_t,int>> links(N);
  std::vector<size_t> starts;
  // construct adjacency matrix
  for (size_t i = 0; i < N; i++) {
    StateM3 si(i);
    std::vector<StateM3> sjs = NextPossibleStates(si);
    for (auto sj: sjs) {
      size_t j = sj.ID();
      links[i][j] = si.RelativePayoff();
    }
    if (si.RelativePayoff() < 0) starts.push_back(i);
  }

  struct path_t {
    std::vector<size_t> p;
    int weight;
    path_t(size_t start) : weight(0) { p.push_back(start); };
  };
  std::vector<path_t> paths;
  for (size_t s: starts) { paths.push_back(s); }
  std::vector<path_t> negative_cycles;

  for (int d = 0; d < 10; d++) {
    std::vector<path_t> new_paths;
    for (path_t& path: paths) {
      size_t last = path.p[path.p.size()-1];
      for (const auto pair: links[last]) {
        path_t new_path = path;
        new_path.p.push_back(pair.first);
        new_path.weight += pair.second;
        if (new_path.p[0] == pair.first && new_path.weight < 0) {  //negative cycle found
          negative_cycles.push_back(new_path);
        }
        else {
          new_paths.push_back(new_path);
        }
      }
    }
    if (!negative_cycles.empty()) break;
    paths = new_paths;
  }

  std::vector<std::vector<size_t>> ans;
  for (path_t& c: negative_cycles) {
    ans.emplace_back(c.p);
  }
  return ans;
}

bool StrategyM3::IsDefensibleDFA() const {
  const auto autom = MinimizeDFA(false).to_map();
  std::vector<std::vector<size_t> > groups;
  groups.reserve(autom.size());
  for (const auto &kv : autom) { groups.emplace_back(kv.second); }
  const size_t AN = groups.size(); // automaton size

  // initialize d_matrix
  const int INF = 32; // 32 is large enough since the path length is between -16 to 16.
  std::vector<std::vector<int> > d(AN);
  for (size_t i = 0; i < AN; i++) { d[i].resize(AN, INF); }

  auto group_index_of_state = [&groups](const StateM3 &ns) -> size_t {
    auto found = std::find_if(groups.cbegin(), groups.cend(), [&ns](std::vector<size_t> g) {
      return std::find(g.begin(), g.end(), ns.ID()) != g.end();
    });
    assert(found != groups.cend());
    return std::distance(groups.cbegin(), found);
  };

  // set distance matrix
  for (size_t i = 0; i < AN; i++) {
    StateM3 sa = StateM3(*groups[i].cbegin());  // get a first state
    Action act_a = ActionAt(sa);  // A's action is same for all states in this group
    for (const auto act_b: std::array<Action, 2>({C, D})) {
      StateM3 ns = sa.NextState(act_a, act_b);
      size_t j = group_index_of_state(ns);
      d[i][j] = MIN(d[i][j], ns.RelativePayoff());
    }
    if (d[i][i] < 0) { return false; }
  }

  for (size_t k = 0; k < AN; k++) {
    for (size_t i = 0; i < AN; i++) {
      for (size_t j = 0; j < AN; j++) {
        d[i][j] = MIN(d[i][j], d[i][k] + d[k][j]);
      }
      if (d[i][i] < 0) { return false; }
    }
  }
  return true;
}

std::array<int, 64> StrategyM3::DestsOfITG() const {
  std::array<int, 64> dests = {};
  std::array<bool, 64> fixed = {false};

  for (int i = 0; i < 64; i++) {
    std::array<bool, 64> visited = {false}; // initialize by false
    visited[i] = true;
    StateM3 init(i);
    int next = NextITGState(init);
    while (next >= 0) {
      if (visited[next] || fixed[next]) { break; }
      visited[next] = true;
      next = NextITGState(StateM3(next));
    }
    int d = next;
    if (next >= 0) {
      d = fixed[next] ? dests[next] : next;
    }
    for (uint64_t j = 0; j < 64; j++) {
      if (visited[j]) {
        dests[j] = d;
        fixed[j] = true;
      }
    }
  }
  return dests;
}

int StrategyM3::NextITGState(const StateM3 &s) const {
  Action move_a = ActionAt(s);
  Action move_b = ActionAt(s.SwapAB());
  if ((move_a == C || move_a == D) && (move_b == C || move_b == D)) {
    return s.NextState(move_a, move_b).ID();
  }
  return -1;
}

std::array<double, 64> StrategyM3::StationaryState(double e, const StrategyM3 *coplayer) const {
  return StationaryStateLapack(e, coplayer);
}

std::array<double, 64> StrategyM3::StationaryStateLapack(double e, const StrategyM3 *coplayer) const {
  if (coplayer == nullptr) { coplayer = this; }

  constexpr lapack_int N = 64, NRHS = 1;
  std::array<double,4096> A = {0.0};
  std::array<double,64> b = {0.0};
  // double A[4096], b[64];  // column-major matrix A[i][j] = A_{ji}

  // Eigen::Matrix<double, 64, 64> A;

  for (int i = 0; i < 64; i++) {
    const StateM3 si(i);
    for (int j = 0; j < 64; j++) {
      // calculate transition probability from j to i
      const StateM3 sj(j);
      // StateM3 next = NextITGState(sj);
      Action act_a = ActionAt(sj);
      Action act_b = coplayer->ActionAt(sj.SwapAB());
      StateM3 next = sj.NextState(act_a, act_b);
      int d = next.NumDiffInT1(si);
      // row: i, col: j
      if (d < 0) {
        A[i+j*N] = 0.0;
      } else if (d == 0) {
        A[i+j*N] = (1.0 - e) * (1.0 - e);
      } else if (d == 1) {
        A[i+j*N] = (1.0 - e) * e;
      } else if (d == 2) {
        A[i+j*N] = e * e;
      } else {
        assert(false);
      }
    }
    A[i+i*N] = A[i+i*N] - 1.0;  // subtract unit matrix
  }
  for (int j = 0; j < 64; j++) { A[63+j*N] += 1.0; }  // normalization condition

  for (int i = 0; i < 63; i++) { b[i] = 0.0; }
  b[63] = 1.0;

  const lapack_int LDA = 64, LDB = 64;
  lapack_int IPIV[64];
  lapack_int info;
  info = LAPACKE_dgesv(LAPACK_COL_MAJOR, N, NRHS, A.data(), LDA, IPIV, b.data(), LDB);

  return b;
}

std::array<double, 64> StrategyM3::StationaryStateEigenDense(double e, const StrategyM3 *coplayer) const {
  if (coplayer == nullptr) { coplayer = this; }
  Eigen::Matrix<double, 64, 64> A;

  for (int i = 0; i < 64; i++) {
    const StateM3 si(i);
    for (int j = 0; j < 64; j++) {
      // calculate transition probability from j to i
      const StateM3 sj(j);
      // StateM3 next = NextITGState(sj);
      Action act_a = ActionAt(sj);
      Action act_b = coplayer->ActionAt(sj.SwapAB());
      StateM3 next = sj.NextState(act_a, act_b);
      int d = next.NumDiffInT1(si);
      if (d < 0) {
        A(i, j) = 0.0;
      } else if (d == 0) {
        A(i, j) = (1.0 - e) * (1.0 - e);
      } else if (d == 1) {
        A(i, j) = (1.0 - e) * e;
      } else if (d == 2) {
        A(i, j) = e * e;
      } else {
        assert(false);
      }
    }
    A(i, i) = A(i, i) - 1.0;  // subtract unit matrix
  }
  for (int i = 0; i < 64; i++) { A(63, i) += 1.0; }  // normalization condition

  Eigen::VectorXd b(64);
  for (int i = 0; i < 63; i++) { b(i) = 0.0; }
  b(63) = 1.0;

  Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

  std::array<double, 64> ans = {0};
  for (int i = 0; i < 64; i++) {
    ans[i] = x(i);
  }
  return ans;
}

std::array<double, 64> StrategyM3::StationaryStateEigenSparse(double e, const StrategyM3 *coplayer) const {
  if (coplayer == nullptr) { coplayer = this; }

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletVec;

  for (int i = 0; i < 64; i++) {
    const StateM3 si(i);
    for (int j = 0; j < 64; j++) {
      // calculate transition probability from j to i
      const StateM3 sj(j);
      Action act_a = ActionAt(sj);
      Action act_b = coplayer->ActionAt(sj.SwapAB());
      StateM3 next = sj.NextState(act_a, act_b);
      int d = next.NumDiffInT1(si);
      double aij = 0.0;
      if (d < 0) {
        // A(i, j) = 0.0;
      } else if (d == 0) {
        // A(i, j) = (1.0 - e) * (1.0 - e);
        aij = (1.0 - e) * (1.0 - e);
      } else if (d == 1) {
        // A(i, j) = (1.0 - e) * e;
        aij = (1.0 - e) * e;
      } else if (d == 2) {
        // A(i, j) = e * e;
        aij = e * e;
      } else {
        assert(false);
      }
      if (i == j) { aij -= 1.0; }  // subtract unix matrix
      if (i == 63) { aij += 1.0; } // normalization condition
      if (aij != 0.0) {
        tripletVec.emplace_back(i, j, aij);
      }
    }
  }
  Eigen::SparseMatrix<double> A(64, 64);
  A.setFromTriplets(tripletVec.cbegin(), tripletVec.cend());

  Eigen::VectorXd b(64);
  for (int i = 0; i < 63; i++) { b(i) = 0.0; }
  b(63) = 1.0;

  // Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
  solver.compute(A);
  Eigen::VectorXd x = solver.solve(b);

  // std::cerr << "#iterations:     " << solver.iterations() << std::endl;
  // std::cerr << "estimated error: " << solver.error() << std::endl;

  std::array<double, 64> ans = {0};
  for (int i = 0; i < 64; i++) {
    ans[i] = x(i);
  }
  return ans;
}

double StrategyM3::CooperationLevel(double e) const {
  auto p = StationaryState(e);
  double c_level = 0.0;
  for (size_t n = 0; n < 64; n++) {
    StateM3 state(n);
    if (state.a_1 == C) { c_level += p[n]; }
  }
  return c_level;
}

std::array<double,2> StrategyM3::Payoffs(const StrategyM3 &coplayer, double benefit, double e) const {
  auto p = StationaryState(e, &coplayer);
  double c_ij = 0.0, c_ji = 0.0;  // cooperation level from i to j and vice versa
  for (size_t n = 0; n < 64; n++) {
    StateM3 s(n);
    if (s.a_1 == C) {
      c_ij += p[n];
    }
    if (s.b_1 == C) {
      c_ji += p[n];
    }
  }
  constexpr double cost = 1.0;
  return { c_ji * benefit - c_ij * cost, c_ij * benefit - c_ji * cost };
}

DirectedGraph StrategyM3::ITG() const {
  DirectedGraph g(64);
  for (int i = 0; i < 64; i++) {
    StateM3 sa(i);
    Action act_a = ActionAt(sa);
    Action act_b = ActionAt(sa.SwapAB());
    int n = sa.NextState(act_a, act_b).ID();
    g.AddLink(i, n);
  }
  return std::move(g);
}

bool StrategyM3::IsEfficientTopo() const {
  if (actions[0] != C) { return false; }

  auto UpdateGn = [](DirectedGraph &gn) {
    components_t sinks = gn.SinkSCCs();
    for (const comp_t &sink: sinks) {
      for (long from: sink) {
        StateM3 s_from(from);
        for (const auto &s_to : s_from.NoisedStates()) {
          uint64_t to = s_to.ID();
          if (!gn.HasLink(from, to)) {
            gn.AddLink(from, to);
          }
        }
      }
    }
  };

  std::vector<int> checked(64, 0);
  checked[0] = 1;
  auto complete = [&checked]() {
    for (int i: checked) {
      if (i == 0) { return false; }
    }
    return true;
  };

  DirectedGraph gn = ITG();
  for (int n = 0; !complete(); n++) {
    if (n > 0) {
      UpdateGn(gn);
    }
    for (int i = 1; i < 64; i++) {
      if (checked[i] == 1) { continue; }
      if (gn.Reachable(i, 0)) {
        if (gn.Reachable(0, i)) {
          return false;   // inefficient
        } else {
          checked[i] = 1;
        }
      }
    }
  }

  return true;
}
bool StrategyM3::IsDistinguishableTopo() const {
  const StrategyM3 allc = StrategyM3::ALLC();
  if (actions[0] != C) { return true; }

  auto UpdateGn = [](DirectedGraph &gn) {
    components_t sinks = gn.SinkSCCs();
    for (const comp_t &sink: sinks) {
      for (long from: sink) {
        StateM3 s_from(from);
        for (const auto &s_to : s_from.NoisedStates()) {
          uint64_t to = s_to.ID();
          if (!gn.HasLink(from, to)) {
            gn.AddLink(from, to);
          }
        }
      }
    }
  };

  std::vector<int> checked(64, 0);
  checked[0] = 1;
  auto complete = [&checked]() {
    for (int i: checked) {
      if (i == 0) { return false; }
    }
    return true;
  };

  DirectedGraph gn(64);
  for (int i = 0; i < 64; i++) {
    StateM3 sa(i);
    StateM3 sb = sa.SwapAB();
    Action act_a = ActionAt(sa);
    Action act_b = allc.ActionAt(sb);
    assert(act_b == C);  // assert AllC
    int j = sa.NextState(act_a, act_b).ID();
    gn.AddLink(i, j);
  }

  for (int n = 0; !complete(); n++) {
    if (n > 0) {
      UpdateGn(gn);
    }
    for (int i = 1; i < 64; i++) {
      if (checked[i] == 1) { continue; }
      if (gn.Reachable(i, 0)) {
        if (gn.Reachable(0, i)) {
          return true;   // inefficient
        } else {
          checked[i] = 1;
        }
      }
    }
  }

  return false;
}

UnionFind StrategyM3::MinimizeDFA(bool noisy) const {
  const size_t N = 64;

  UnionFind uf_0(N);
  // initialize grouping by the action c/d
  size_t c_rep = N, d_rep = N;
  for (size_t i = 0; i < N; i++) {
    if (actions[i] == C) {
      c_rep = i;
      break;
    }
  }
  for (size_t i = 0; i < N; i++) {
    if (actions[i] == D) {
      d_rep = i;
      break;
    }
  }
  for (size_t i = 0; i < N; i++) {
    size_t target = (actions[i] == C) ? c_rep : d_rep;
    uf_0.merge(i, target);
  }

  while (true) {
    UnionFind uf(N);
    size_t i = 0;
    while (i < N) {
      size_t i_next = N;
      for (size_t j = i+1; j < N; j++) {
        if (!uf.is_root(j)) continue;  // skip if it is already merged
        if (uf_0.root(i) == uf_0.root(j) && _Equivalent(i, j, uf_0, noisy)) {
          uf.merge(i, j);
        }
        else if (i_next == N) { i_next = j; }  // keep the first unmerged node
      }
      i = i_next;
    }
    if (uf_0 == uf) break;
    uf_0 = uf;
  }
  return uf_0;
}

bool StrategyM3::_Equivalent(size_t i, size_t j, UnionFind &uf_0, bool noisy) const {
  assert(actions[i] == actions[j]);
  Action act_a = actions[i];
  Action err_a = (act_a == C) ? D : C;
  std::array<Action, 2> acts_b = {C, D};
  for (const Action &act_b : acts_b) {
    size_t ni = StateM3(i).NextState(act_a, act_b).ID();
    size_t nj = StateM3(j).NextState(act_a, act_b).ID();
    if (uf_0.root(ni) != uf_0.root(nj)) { return false; }
    if (noisy) {
      size_t ni2 = StateM3(i).NextState(err_a, act_b).ID();
      size_t nj2 = StateM3(j).NextState(err_a, act_b).ID();
      if (uf_0.root(ni2) != uf_0.root(nj2)) { return false; }
    }
  }
  return true;
}

Partition StrategyM3::MinimizeDFAHopcroft(bool noisy) const {
  Partition partition(64);

  // initialize partition by the action c/d
  std::set<size_t> c_set, d_set;
  for (size_t i = 0; i < 64; i++) {
    if (actions[i] == C) { c_set.insert(i); }
    else { d_set.insert(i); }
  }
  if (c_set.empty() || d_set.empty()) { return std::move(partition); }
  partition.split(0, c_set);
  size_t smaller = (c_set.size() < d_set.size()) ? (*c_set.begin()) : (*d_set.begin());

  std::set<splitter_t> waiting;  // initialize waiting set
  std::vector<int> inputs = {0, 1};  // 0: c, 1: d, (0: cc, 1: cd, 2: dc, 3: dd)
  if (noisy) { inputs.push_back(2); inputs.push_back(3); }
  for(int b: inputs) {
    waiting.insert({smaller, b});
  }

  while( !waiting.empty() ) {
    splitter_t splitter = *waiting.cbegin();  // take a splitter from a waiting list
    const std::set<size_t> Q = partition.group(splitter.first);  // copy splitter because this must remain same during this iteration
    const int b = splitter.second;
    waiting.erase(waiting.begin());
    auto groups = partition.group_ids();
    for(size_t p: groups) {  // for each P in partition
      // check if group i is splittable or not
      auto p1p2 = _SplitBySplitter(partition, p, Q, b, noisy);
      const std::set<size_t> &p1 = p1p2.at(0), &p2 = p1p2.at(1);
      if (p1.empty() || p2.empty() ) { continue; }  // P is not split by this splitter
      partition.split(p, p1);
      for (int b: inputs) {
        auto found = waiting.find(splitter_t(p, b) );
        if (found != waiting.end()) {
          // replace (P, b) by (P1, b) and (P2, b) in W
          waiting.erase(found);
          waiting.insert({*p1.cbegin(), b});
          waiting.insert({*p2.cbegin(), b});
        }
        else {
          if (p1.size() < p2.size()) {
            waiting.insert({*p1.cbegin(), b});
          }
          else {
            waiting.insert({*p2.cbegin(), b});
          }
        }
      }
    }
  }
  return std::move(partition);
}

std::array<std::set<size_t>,2> StrategyM3::_SplitBySplitter(const Partition &partition, size_t p, const std::set<size_t> &Q, int b, bool noisy) const {
  const std::set<size_t> &P = partition.group(p);
  Action act_b = (b & 1) ? D : C;
  // get members of P which go to a member of Q by b
  std::set<size_t> P1;
  for (size_t si: P) {
    StateM3 s(si);
    Action act_a = actions[si];
    if (noisy) { act_a = (b & 2) ? D : C; }
    size_t next = s.NextState(act_a, act_b).ID();
    if (Q.find(next) != Q.end()) {
      P1.insert(si);
    }
  }
  std::set<size_t> P2;
  std::set_difference(P.begin(), P.end(), P1.begin(), P1.end(), std::inserter(P2, P2.end()));
  return {P1, P2};
}

Action CAPRIn_action_at(size_t i, bool s_capri) {
  typedef std::bitset<6> B;
  const B bit(i);
  std::string s = bit.to_string('c', 'D');
  std::reverse(s.begin(), s.end());
  const std::string Istr = s.substr(0,3) + '-' + s.substr(3,3);

  const B oldest = 0b100'100ul, latest = 0b001'001ul;
  // In this implementation, upper/lower bits correspond to A/B, respectively
  const B a_mask = 0b111'000ul;
  const B b_mask = 0b000'111ul;
  const size_t na = (bit & a_mask).count(), nb = (bit & b_mask).count();

  // E: Exploit others if payoff difference is n or greater
  if (!s_capri) {
    int d = std::max(na, nb) - std::min(na, nb);
    if (d >= 2) { return D; }
  }

  size_t last_ccc = 3;
  for (size_t t = 0; t < 3; t++) {
    if ((bit & (latest<<t)) == B(0ul)) {  // CCC is found at t step before
      last_ccc = t;
      break;
    }
  }

  // C: cooperate if the mutual cooperation is formed at last two rounds
  if ((bit & latest) == 0ul) {
    return C;
  }
  else if (last_ccc > 0 && last_ccc < 3) {
    B mask = latest;
    for (size_t t = 0; t < last_ccc; t++) { mask = ((mask << 1ul) | latest); }
    // A: Accept punishment by prescribing *C* if all your relative payoffs are at least zero.
    size_t pa = (bit & mask & a_mask).count();
    size_t pb = (bit & mask & b_mask).count();
    if (pa >= pb) {
      return C;
    }
      // P: Punish by *D* if any of your relative payoffs is negative.
    else {
      return D;
    }
  }
  // R: grab the chance to recover
  if (bit == 0b111'110 || bit == 0b110'111) {
    // R: If payoff profile is (+1,+1,-1), prescribe *C*.
    return C;
  }
  // In all other cases, *D*
  return D;

}

StrategyM3 StrategyM3::CAPRI2() {
  const size_t N = 64;
  std::array<Action,64> acts{};
  for (size_t i = 0; i < N; i++) {
    acts[i] = CAPRIn_action_at(i, false);
  }
  return std::move(StrategyM3(acts));
}

StrategyM3 StrategyM3::sCAPRI2() {
  const size_t N = 64;
  std::array<Action,64> acts{};
  for (size_t i = 0; i < N; i++) {
    acts[i] = CAPRIn_action_at(i, true);
  }
  return std::move(StrategyM3(acts));
}

StrategyM3 StrategyM3::AON(size_t n) {
  if (n > 3 || n <= 0) { throw std::runtime_error("n must be less than 4"); }
  const size_t N = 64;
  std::array<Action,64> acts{};
  size_t mask = 0b001;
  if (n == 2) mask = 0b011;
  else if (n == 3) mask = 0b111;
  for (size_t i = 0; i < N; i++) {
    size_t a_histo = i & 0b111 & mask;
    size_t b_histo = (i>>3) & 0b111 & mask;
    if (a_histo == b_histo) acts[i] = C;
    else acts[i] = D;
  }
  return std::move(StrategyM3(acts));
}

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <map>
#include <chrono>
#include <bitset>
#include <mpi.h>
#include <omp.h>
#include "icecream.hpp"
#include "caravan.hpp"


using json = nlohmann::json;

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
size_t output_automaton_size_max = 0;

int MinimizedDFASize(const char str[64]) {
  std::array<size_t,64> F;
  for (size_t i = 0; i < 64; i++) {
    F[i] = (str[i] == 'c') ? 1 : 0;
  }
  DFA<64,4> dfa(F, g_delta);
  return dfa.MinimizedSize();
}

using automaton_sizes_t = std::array<size_t, 65>;

void ExpandWildcard(const std::string& line, automaton_sizes_t& automaton_sizes, std::vector<std::string>& found) {
  for (int i = 0; i < line.size(); i++) {
    if (line[i] == '*') {
      std::string temp = line;
      temp[i] = 'c';
      ExpandWildcard(temp, automaton_sizes, found);
      temp[i] = 'd';
      ExpandWildcard(temp, automaton_sizes, found);
      return;
    }
  }
  int size = MinimizedDFASize(line.c_str());
  automaton_sizes[size]++;
  if (size <= output_automaton_size_max) {
    found.emplace_back(line);
  }
}

automaton_sizes_t AutomatonSizes(const std::string& line, std::vector<std::string>& found) {
  automaton_sizes_t sizes;
  sizes.fill(0);

  ExpandWildcard(line, sizes, found);
  return sizes;
}

void AddExpandedLines(const std::string &line, std::vector<std::string>& lines, int level) {
  if (level == 0) {
    lines.push_back(line);
    return;
  }
  for (int i = 0; i < line.size(); i++) {
    if (line[i] == '*') {
      std::string temp = line;
      temp[i] = 'c';
      AddExpandedLines(temp, lines, level-1);
      temp[i] = 'd';
      AddExpandedLines(temp, lines, level-1);
      return;
    }
  }
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  int my_rank, total_size;
  MPI_Comm_size(MPI_COMM_WORLD, &total_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0) {
    int num_threads = omp_get_max_threads();
    std::cerr << "total_size = " << total_size << std::endl;
    std::cerr << "num_threads = " << num_threads << std::endl;
  }

  if (argc != 3) {
    if (my_rank == 0) {
      std::cerr << "Usage: " << argv[0] << " <input_file> <output_automaton_size_max>" << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }

  std::ifstream fin;
  std::ofstream fout;
  if (my_rank == 0) {
    fin.open(argv[1]);
    if (!fin) {
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    fout.open("automatons");
  }
  output_automaton_size_max = std::stoul(argv[2]);

  std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

  automaton_sizes_t total_sizes;
  total_sizes.fill(0);

  auto read_lines_and_push_task = [&fin,&start](caravan::Queue& q) -> uint64_t {
    json lines;
    size_t num_strategies = 0ul;
    const size_t max_line_size = 8192;
    const size_t max_num_strategies = 1ul << 22;

    std::string line;
    while (std::getline(fin, line)) {
      if (line.empty()) break;
      lines.emplace_back(line);
      size_t n = 1ul;
      for (size_t i = 0; i < line.size(); i++) {
        if (line[i] == '*') { n = n << 1; }
      }
      num_strategies += n;
      if (lines.size() >= max_line_size || num_strategies >= max_num_strategies) {
        break;
      }
    }
    // std::cerr << num_strategies << ' ' << lines.size() << std::endl;
    uint64_t task_id = q.Push(lines);
    if (task_id % 100000 == 0) {
      std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;
      std::cerr << "task_id:" << task_id << ' ' << elapsed_seconds.count() << std::endl;
    }
    return task_id;
  };

  // define a pre-process: create json object that contains parameters of tasks
  // This function is called only at the master process.
  std::function<void(caravan::Queue&)> on_init = [&fin,&read_lines_and_push_task,total_size](caravan::Queue& q) {
    while(fin) {
      uint64_t task_id = read_lines_and_push_task(q);
      // IC(task_id);
      if (task_id >= total_size) {
        break;
      }
    }
  };

  // After the task was executed at a worker process, its result is returned to the master process.
  // When the master process receives the result, this callback function is called at the master process.
  std::function<void(int64_t, const json&, const json&, caravan::Queue&)> on_result_receive = [&fin,&read_lines_and_push_task,&total_sizes,&fout](int64_t task_id, const json& input, const json& output, caravan::Queue& q) {
    auto v = output["sizes"].get<std::vector<size_t>>();
    for (size_t i = 0; i < v.size(); i++) {
      total_sizes[i] += v[i];
    }
    if (!output["found"].empty()) {
      for (auto& line : output["found"]) {
        fout << line.get<std::string>() << std::endl;
      }
    }
    if (fin) {
      uint64_t task_id = read_lines_and_push_task(q);
    }
  };

  // Define the function which is executed at a worker process.
  // The input parameter for the task is given as the argument.
  std::function<json(const json& input)> do_task = [](const json& input) {
    std::vector<std::string> lines;
    for (const json& j : input) {
      std::string line = j.get<std::string>();
      int n_ast = 0;
      for (size_t i = 0; i < line.size(); i++) {
        if (line[i] == '*') { n_ast++; }
      }

      int max = 20;
      if (n_ast > max) {
        AddExpandedLines(line, lines, n_ast-max);
        // IC(line, n_ast);
      }
      else {
        lines.push_back(line);
      }
    }

    automaton_sizes_t sizes;
    sizes.fill(0);
    std::vector<std::string> found;
    #pragma omp parallel for shared(lines,sizes) schedule(dynamic,1)
    for (size_t i = 0; i < lines.size(); i++) {
      automaton_sizes_t r = AutomatonSizes(lines[i], found);
      for (size_t i = 0; i < r.size(); ++i) {
        #pragma omp atomic
        sizes[i] += r[i];
      }
    }

    json output;
    output["sizes"] = sizes;
    output["found"] = found;
    return output;
  };


  caravan::Start(on_init, on_result_receive, do_task, MPI_COMM_WORLD);

  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;

  if (my_rank == 0) {
    size_t sum = 0;
    double total = 0.0;
    for (int i = 0; i < 65; i++) {
      std::cout << i << ' ' << total_sizes[i] << std::endl;
      sum += total_sizes[i];
      total += i * total_sizes[i];
    }
    IC(sum, total/sum, elapsed_seconds.count(), elapsed_seconds.count()/sum);
  }

  MPI_Finalize();
  return 0;
}

# sim_friendly_rival_automaton_sizes

Comprehensive calculation of DFA minimization on memory-3 strategies.
It is parallelized using MPI and OpenMP. Clone this repository with submodules.

```
git submodule --recursive git@github.com:yohm/sim_friendly_rival_automaton_sizes.git
```


## Executables

Use cmake to build executables. On Fugaku, run `./build_fugaku.sh` instead of cmake.
Some of the executalbes require lapack and eigen.

### main

Execution:

```shell
mpiexec -np <num process> main <input file> <output size max> <output size max for simplified automata>
```

For instance, the following code reads input from the file `friendly_rival_10k` and print out the histogram of the automaton sizes.
When the full automaton size is less than 13 or the simplified automaton size is less than 8, the strategies are printed to the file `automatons`.

```shell
mpiexec -np 4 main friendly_rival_10k 13 8
```

To run on Fugaku,
1. split the huge input list (that is 540GB) into smaller pieces using `split_lines.rb`.
  `ruby split_lines.rb strategy_list 1000000000`
2. You'll find the directories `split_***`. In each directory, submit the job.
  `cd split_000; pjsub fugaku_job.sh`
3. Repeat the above for all directories.

### sort_by_size

Sort the list of strategies according to the size of the automata. Serialized strings of simplified automata are printed to standard error.

```shell
g++ -O3 sort_by_size.cpp -o sort_by_size
./sort_by_size automatons > automatons_sorted 2> serialized
```

### test_automaton_size

Unit test program. It is also used for printing the automaton for a strategy like the following:

```shell
./test_automaton_size cdcdcdcdcdcdcdcddccccdcccdddcdddcdcdcdcdcdcdcdcdcdcccdcccdddcddd 1
```

### minimum

A minimum code for tuning & debugging of the code.

```shell
g++ -O3 minimum.cpp -o minimum
./minimum
```



#%%
import graphviz

# %%
line = "5 cdcdcdcdcccdcccddccdcdcdcccdcccdcdcdcdcdcccdcccddccdcdcdcccdcccd 0,2,4,6,10,14,18,20,22,26,30,32,34,36,38,42,46,50,52,54,58,62;1,3,5,7,11,15,19,21,23,27,31,33,35,37,39,43,47,51,53,55,59,63;8,12,24,28,40,44,56,60;9,13,17,25,29,41,45,49,57,61;16,48 0,2,4,6,9,10,13,14,17,18,20,22,25,26,29,30,32,34,36,38,41,42,45,46,49,50,52,54,57,58,61,62;1,3,5,7,11,15,19,21,23,27,31,33,35,37,39,43,47,51,53,55,59,63;8,12,24,28,40,44,56,60;16,48"
# %%
def line_to_partition(line):
  a = line.split()
  size = int(a[0])
  strategy = a[1]
  def string_to_partition(s):
    partition = list(range(64))
    for n_list in s.split(';'):
      nodes = [int(n) for n in n_list.split(',')]
      idx = nodes[0]
      for n in nodes:
        partition[n] = idx
    return partition

  partition_full = string_to_partition(a[2])
  partition_simple = string_to_partition(a[3])

  return size,strategy,partition_full,partition_simple

size,strategy,partition_f,partition_s = line_to_partition(line)
size,strategy,partition_f,partition_s
" ".join(map(str, partition_s)), " ".join(map(str, partition_f))
# %%
dests_graph = []
for i in range(64):
  b = (i << 1) & 0b111111
  dests = [
    (b & 0b110110) | 0b000000,
    (b & 0b110110) | 0b000001,
    (b & 0b110110) | 0b001000,
    (b & 0b110110) | 0b001001
  ]
  dests_graph.append(dests)

dests_graph


# %%
def make_dot(strategy, partition, simplified=False):
  color_list = ["skyblue", "salmon"]
  dot = graphviz.Digraph(strategy)
  # set nodes
  for idx,n in enumerate(partition):
    if idx == n:  # root node
      color = color_list[0] if strategy[n] == 'c' else color_list[1]
      dot.node(str(n), str(n), color=color, style='filled')
  # set links
  for idx,n in enumerate(partition):
    if idx == n:      # if n is a root node
      for dest_idx,d in enumerate(dests_graph[n]):
        label = ['cc', 'cd', 'dc', 'dd'][dest_idx]
        edge_colors = {'cc': 'gray', 'cd': 'salmon', 'dc': 'skyblue', 'dd': 'gray'}
        color = edge_colors[label]
        style = 'solid' if strategy[n] == label[0] else 'dashed'
        draw_edge = True
        if simplified:
          draw_edge = False if style == 'dashed' else True
          if n == 0 and label == 'dc':
            draw_edge = True
          if n == partition[63] and label == 'cd':
            draw_edge = True
        if draw_edge:
          dot.edge(str(n), str(partition[d]), label=label, color=color, style=style)
  return dot

make_dot(strategy, partition_f, simplified=False)

# %%
make_dot(strategy, partition_s, simplified=True)

# %%
# open file and read the first line
with open('../job/automatons_sorted', 'r') as f:
  for line_idx,line in enumerate(f):
    size,strategy,partition_f,partition_s = line_to_partition(line)
    print(strategy)
    if size > 5:
      break
    dot_f = make_dot(strategy, partition_f,simplified=False)
    dot_f.render(f"automaton_{line_idx}_full")
    dot_s = make_dot(strategy, partition_s,simplified=True)
    dot_s.render(f"automaton_{line_idx}_simple")
# %%

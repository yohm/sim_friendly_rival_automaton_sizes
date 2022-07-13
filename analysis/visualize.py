#%%
import graphviz

# %%
line = "5 cdcdcdcdcccdcccddccdcdcdcccdcccdcdcdcdcdcccdcccddccdcdcdcccdcccd 0,2,4,6,10,14,18,20,22,26,30,32,34,36,38,42,46,50,52,54,58,62, 1,3,5,7,11,15,19,21,23,27,31,33,35,37,39,43,47,51,53,55,59,63, 8,12,24,28,40,44,56,60, 9,13,17,25,29,41,45,49,57,61, 16,48, "
# %%
def line_to_partition(line):
  a = line.split()
  print(line)
  size = int(a[0])
  strategy = a[1]
  partition = list(range(64))
  for x in a[2:]:
    nodes = [int(n) for n in x.rstrip(',').split(',')]
    idx = nodes[0]
    for n in nodes:
      partition[n] = idx
  return size,strategy,partition

size,strategy,partition = line_to_partition(line)
size,strategy,partition
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
def make_dot(strategy, partition):
  color_list = ["skyblue", "salmon"]
  dot = graphviz.Digraph(strategy)
  # set nodes
  for idx,n in enumerate(partition):
    if idx == n:  # root node
      color = color_list[0] if strategy[n] == 'c' else color_list[1]
      dot.node(str(n), str(n), color=color, style='filled')
  # set links
  for idx,n in enumerate(partition):
    if idx == n:
      for dest_idx,d in enumerate(dests_graph[n]):
        label = ['cc', 'cd', 'dc', 'dd'][dest_idx]
        color = color_list[0] if label[1] == 'c' else color_list[1]
        style = 'solid' if strategy[n] == label[0] else 'dashed'
        dot.edge(str(n), str(partition[d]), label=label, color=color, style=style)
  return dot

make_dot(strategy, partition)

# %%
# open file and read the first line
with open('../job/automatons_sorted', 'r') as f:
  for line_idx,line in enumerate(f):
    size,strategy,partition = line_to_partition(line)
    if size > 6:
      break
    dot = make_dot(strategy, partition)
    dot.render(f"automaton_{line_idx}")
# %%

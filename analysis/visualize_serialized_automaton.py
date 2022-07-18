#%%
from collections import defaultdict
import graphviz
from collections import defaultdict

# %%
autom = "0,c,0,1,2;1,d,0,1,0;2,d,2,3;3,c,0,0;"
recov_path = "2,1;3,1;0,0;0,0;"
escape_path = "0,0;0,0;"
# %%
def visualize_automaton(autom, attr=None, graph_label=None):
  if not attr:
    attr = autom
  dot = graphviz.Digraph(attr)
  dot.graph_attr = {'label': graph_label, 'labelloc': 't', 'fontsize': '16'}
  color_list = ["skyblue", "salmon"]
  nodes = [node_s.split(',') for node_s in autom.split(';')]
  # remove last node
  nodes = nodes[:-1]
  for n in nodes:
    color = color_list[0] if n[1] == 'c' else color_list[1]
    dot.node(n[0], color=color, style='filled')
  for n in nodes:
    label_color = {'cc': 'gray', 'cd': 'salmon', 'dc': 'skyblue', 'dd': 'gray'}
    label = 'cc' if n[1] == 'c' else 'dc'
    dot.edge(n[0], n[2], label=label, color=label_color[label], style='solid')
    label = 'cd' if n[1] == 'c' else 'dd'
    dot.edge(n[0], n[3], label=label, color=label_color[label], style='solid')
    if len(n) == 5 and n[1] == 'c':
      label = 'dc'
      dot.edge(n[0], n[4], label=label, color=label_color[label], style='dashed')
    if len(n) == 5 and n[1] == 'd':
      label = 'cd'
      dot.edge(n[0], n[4], label=label, color=label_color[label], style='dashed')
  return dot

label = f"recovery path: {recov_path}\nescape path: {escape_path}"
visualize_automaton(autom, graph_label=label)
# %%
recovery_path_count = defaultdict(int)
escape_path_count = defaultdict(int)
with open("../job/serialized") as f:
  lines = f.readlines()
  for idx,line in enumerate(lines):
    autom,recov_path,escape_path,strategy,freq = line.split(' ')
    recovery_path_count[recov_path] += 1
    escape_path_count[(recov_path,escape_path)] += 1
    label = f"recovery path: {recov_path}\nescape path: {escape_path}"
    dot = visualize_automaton(autom, attr=f"{autom}_{strategy}", graph_label=label)
    dot.render(f"simplified_automaton_{idx}")
# %%
recovery_path_count, escape_path_count
# %%

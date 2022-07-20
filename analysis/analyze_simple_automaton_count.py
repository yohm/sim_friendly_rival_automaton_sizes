#%%
import graphviz
from collections import defaultdict
# %%
# load file
file = open("../job/simple_automaton_count_merged", "r")

histo = defaultdict(int)
automs = defaultdict(set)

# for each line in file, split by spaces
for line in file:
  autom,recov,escape,count = line.split()
  key = recov + ' ' + escape
  histo[key] += int(count)
  automs[key].add(autom)

histo, automs
# %%
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

visualize_automaton("0,c,0,1,2;1,d,3,1,2;2,c,0,3;3,d,3,0;", graph_label="2,1;3,3;0,0;0,0; 2,3;3,3;0,0;0,0;")
# %%
from IPython import display
sorted_histo = sorted(histo.items(), key=lambda x: x[0])

for key, value in sorted_histo:
  autom = list(automs[key])[0]
  dot = visualize_automaton(autom, graph_label=f"{key} {value}")
  display.display_svg(dot)
  #display.SVG(dot.pipe('svg'))
# %%
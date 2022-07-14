#%%
import graphviz

# %%
autom = "0,c,0,1,2;1,d,3,1,1;2,c,0,4;3,d,3,0;4,d,2,0"
# %%
def visualize_automaton(autom, attr=None):
  if not attr:
    attr = autom
  dot = graphviz.Digraph(attr)
  color_list = ["skyblue", "salmon"]
  nodes = [node_s.split(',') for node_s in autom.split(';')]
  for n in nodes:
    color = color_list[0] if n[1] == 'c' else color_list[1]
    dot.node(n[0], color=color, style='filled')
  for n in nodes:
    color = color_list[0] if n[1] == 'c' else color_list[1]
    label = 'cc' if n[1] == 'c' else 'dc'
    dot.edge(n[0], n[2], label=label, color=color, style='solid')
    label = 'cd' if n[1] == 'c' else 'dd'
    dot.edge(n[0], n[3], label=label, color=color, style='solid')
    if len(n) == 5 and n[1] == 'c':
      dot.edge(n[0], n[4], label='dc', color=color_list[1], style='dashed')
    if len(n) == 5 and n[1] == 'd':
      dot.edge(n[0], n[4], label='cd', color=color_list[0], style='dashed')
  return dot

visualize_automaton(autom)
# %%
with open("../job/serialized") as f:
  lines = f.readlines()
  for idx,line in enumerate(lines):
    autom,strategy,freq = line.split(' ')
    dot = visualize_automaton(autom, f"{autom}_{strategy}")
    dot.render(f"simplified_automaton_{idx}")
    if idx >= 30:
      break
# %%

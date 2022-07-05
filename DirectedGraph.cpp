#include <iostream>
#include "DirectedGraph.hpp"

DirectedGraph::DirectedGraph(size_t num_nodes) : m_num_nodes(num_nodes) {
  m_links.resize(m_num_nodes);
}

void DirectedGraph::AddLink(long from, long to) {
  m_links[from].push_back(to);
}

std::ostream &operator<<(std::ostream &os, const DirectedGraph &graph) {
  os << "m_num_nodes: " << graph.m_num_nodes << "\n";
  os << "m_links: \n";
  for( long from=0; from < graph.m_num_nodes; from++) {
    for( auto to : graph.m_links[from] ) {
      os << "  " << from << " " << to << "\n";
    }
  }
  return os;
}

std::set<long> DirectedGraph::TransitionNodes() const {
  std::set<long> transition_nodes;

  components_t components;
  SCCs(components);
  for( const std::vector<long>& component: components) {
    if( component.size() == 1 ) {
      long n = component[0];
      if( !HasSelfLoop(n) ) { transition_nodes.insert(n); }
    }
  }

  return std::move(transition_nodes);
}

components_t DirectedGraph::NonTransitionComponents() const {
  components_t sccs;
  SCCs(sccs);
  components_t ans;
  for(const auto& scc: sccs) {
    if(scc.size() > 1 ) {
      ans.push_back(scc);
    }
    else {
      long n = scc[0];
      if( HasSelfLoop(n) ) { ans.push_back(scc); }
    }
  }
  return std::move(ans);
}

DirectedGraph::ComponentFinder::ComponentFinder(const DirectedGraph &m_g) : m_g(m_g) {
  size_t n = m_g.m_num_nodes;
  m_t = 0;
  desc.resize(n,-1);
  low.resize(n,-1);
  on_stack.resize(n,false);
}

void DirectedGraph::ComponentFinder::SCCs( components_t& components ) {
  components.clear();

  for( long v=0; v < m_g.m_num_nodes; v++ ) {
    if( desc[v] < 0 ) {
      StrongConnect(v, components);
    }
  }
}

void DirectedGraph::ComponentFinder::StrongConnect(long v, components_t &components) {
  desc[v] = m_t;
  low[v] = m_t;
  m_t++;

  stack.push(v);
  on_stack[v] = true;

  for( long w : m_g.m_links.at(v) ) {
    if( desc[w] < 0 ) {
      StrongConnect(w, components);
      if( low[w] < low[v] ) { low[v] = low[w]; }
    }
    else if( on_stack[w] ) {
      if( desc[w] < low[v] ) { low[v] = desc[w]; }
    }
  }

  std::vector<long> comp;
  if( low[v] == desc[v] ) {
    while(true) {
      long w = stack.top();
      stack.pop();
      on_stack[w] = false;
      comp.push_back(w);
      if( v==w ) { break; }
    }
    components.push_back(comp);
  }
}

bool DirectedGraph::HasSelfLoop(long n) const {
  bool b = false;
  for( long j: m_links[n]) {
    if( j == n ) {
      b = true;
    }
  }
  return b;
}

components_t DirectedGraph::SinkSCCs() const {
  components_t sccs;
  SCCs(sccs);
  components_t ans;
  for(const auto& scc: sccs) {
    bool no_out = true;
    for(const long i: scc) {
      for(const long j: m_links[i]) {
        if(std::find(scc.begin(),scc.end(),j) == scc.end() ) {
          no_out = false;
          break;
        }
      }
      if(!no_out) { break; }
    }
    if(no_out) {
      ans.push_back(scc);
    }
  }
  return std::move(ans);
}

bool DirectedGraph::HasLink(long from, long to) const {
  return ( std::find(m_links[from].begin(), m_links[from].end(), to) != m_links[from].end() );
}

bool DirectedGraph::Reachable(long from, long to) const {
  std::vector<int> visited(m_num_nodes, 0);
  std::function<bool(long)> dfs = [to,&visited,this,&dfs] (long i) {
    visited[i] = 1;
    for(long j: this->m_links[i]) {
      if(j == to) { return true; }
      if(visited[j]) { continue; }
      bool b = dfs(j);
      if(b) { return true; }
    }
    return false;
  };
  return dfs(from);
}
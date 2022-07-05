#ifndef DIRECTED_GRAPH_HPP
#define DIRECTED_GRAPH_HPP

#include <iostream>
#include <vector>
#include <map>
#include <stack>
#include <set>
#include <functional>
#include <queue>
#include <algorithm>

typedef std::vector<long> comp_t;
typedef std::vector<comp_t> components_t;

class DirectedGraph {
public:
  DirectedGraph(size_t num_nodes);

  void AddLink(long from, long to);
  bool HasLink(long from, long to) const;
  bool Reachable(long from, long to) const;
  friend std::ostream &operator<<(std::ostream &os, const DirectedGraph &graph);
  void SCCs(components_t& components) const {
    ComponentFinder cf(*this);
    cf.SCCs(components);
  }
  std::set<long> TransitionNodes() const;
  components_t NonTransitionComponents() const;
  components_t SinkSCCs() const; // SCCs that has no outgoing link
  template <class T>
  void ForEachLink( const T& f) const {
    for( long i=0; i<m_num_nodes; i++) {
      for( long j: m_links[i]) {
        f(i,j);
      }
    }
  }
  template <class T> void BFS(long init, const T& f) const {
    std::vector<int> visited(m_num_nodes, 0);
    std::queue<long> q;
    q.push(init);
    visited[init] = 1;
    while( q.size() > 0 ) {
      long i = q.front();
      q.pop();
      f(i);
      for(long j: m_links[i]) {
        if(visited[j] == 0) {
          visited[j] = 1;
          q.push(j);
        }
      }
    }
  }
  template <class T> void DFS(long init, const T& f) const {
    std::vector<int> visited(m_num_nodes, 0);
    std::stack<long> s;
    s.push(init);
    while( s.size() > 0 ) {
      long i = s.top();
      s.pop();
      if( visited[i] > 0 ) { continue; }
      visited[i] = 1;
      f(i);
      for(long j: m_links[i]) {
        if(visited[j] == 0) { s.push(j); }
      }
    }
  }

  const size_t m_num_nodes;
  std::vector<std::vector<long> > m_links;
private:
  bool HasSelfLoop(long n) const; // return true if node n has a self-loop

  class ComponentFinder {
  public:
    ComponentFinder(const DirectedGraph &m_g);
    void SCCs( components_t& components );
  private:
    const DirectedGraph& m_g;
    long m_t;
    std::vector<long> desc;
    std::vector<long> low;
    std::stack<long> stack;
    std::vector<bool> on_stack;

    void StrongConnect( long v, components_t& components);

  };
};

#endif //DIRECTED_GRAPH_HPP

#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>

using namespace boost;

namespace socialgraph{
	
	typedef std::pair<int, int> Edge;

	typedef std::vector<Edge> Edges;

	typedef adjacency_list_traits<vecS, vecS, bidirectionalS> Traits;

	typedef adjacency_list<vecS, vecS, bidirectionalS, property<vertex_name_t, std::string> > Graph;
	
	//define residual graph
	typedef adjacency_list<vecS, vecS, bidirectionalS,
		no_property,
		property<edge_capacity_t, long,
		property<edge_residual_capacity_t, long,
		property<edge_reverse_t,Traits::edge_descriptor> > > > ResidualGraph;

	typedef graph_traits<Graph>::vertex_iterator Graph_vertex_iter;

	typedef graph_traits<Graph>::edge_iterator Graph_edge_iter;
	
	typedef graph_traits<Graph>::vertex_descriptor Vertex;
	
	typedef graph_traits<Graph>::vertices_size_type Size;
	
	void printAllVertices(Graph g);

	void printAllEdges(Graph g);

	ResidualGraph createResidualGraph(Graph g);
	
	ResidualGraph createResidualGraph2(Graph g);
	
	typedef property_map<ResidualGraph, edge_capacity_t>::type Capacity;
	typedef property_map<ResidualGraph, edge_reverse_t>::type Reverse;
	typedef property_map<ResidualGraph, edge_residual_capacity_t>::type Residual_Capacity;
	typedef std::vector<default_color_type> Color;
	typedef std::vector<Traits::edge_descriptor> Pred;
	
	int maxFlow(ResidualGraph rg, int s, int t, Capacity &capacity, Reverse &reverse, Residual_Capacity &residual_capacity, Color &color, Pred &pred);
	
	template < typename TimeMap > class bfs_time_visitor:public default_bfs_visitor {
		typedef typename property_traits < TimeMap >::value_type T;
	public:
		bfs_time_visitor(TimeMap tmap, T & t):m_timemap(tmap), m_time(t) { }
		template < typename Vertex, typename Graph >
		  void discover_vertex(Vertex u, const Graph & g) const
		{
		  put(m_timemap, u, m_time++);
		}
		TimeMap m_timemap;
		T & m_time;
	};
	
	extern Graph global_graph;
	
	//for centrality maps
	typedef property_map<Graph, vertex_index_t>::type VertexIndexMap;
	
	typedef std::map<Edge, int> StdEdgeIndexMap;
	typedef associative_property_map<StdEdgeIndexMap> EdgeIndexMap;
	
		
	//weighted graphs
	typedef adjacency_list<vecS, vecS, bidirectionalS, no_property, property<edge_weight_t, int> > WeightedGraph;
	
	typedef graph_traits<WeightedGraph>::vertex_descriptor WeightedVertex;
	typedef graph_traits<WeightedGraph>::edge_descriptor WeightedEdge;
	
	WeightedGraph createWeightedGraph(Graph g);
	void ShortestPaths(WeightedGraph wg, int s, std::vector<double> &centrality);

}

#endif

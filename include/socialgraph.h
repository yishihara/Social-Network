#ifndef SOCIALGRAPH_H
#define SOCIALGRAPH_H

#include <vector>
#include <string>
#include <iostream>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <filehandler.h>

#define PI (3.14159)

using namespace boost;

//definitions for edges
typedef std::pair<int, int> Edge;
typedef std::vector<Edge> Edges;

typedef adjacency_list<vecS, vecS, bidirectionalS, property<vertex_name_t, std::string> > Graph;
typedef graph_traits<Graph>::edge_iterator EdgeIterator;

//vertex edge descriptors
typedef Graph::vertex_descriptor VertexDescriptor;
typedef Graph::edge_descriptor EdgeDescriptor;

//index maps
typedef std::map<EdgeDescriptor, int> StdEdgeIndexMap;
typedef boost::associative_property_map<StdEdgeIndexMap> EdgeIndexMap;
typedef property_map<Graph, vertex_index_t>::type VertexIndexMap;

//weighted graphs
typedef adjacency_list<vecS, vecS, bidirectionalS, no_property, property<edge_weight_t, int> > WeightedGraph;
typedef graph_traits<WeightedGraph>::vertex_descriptor WeightedVertexDescriptor;
typedef graph_traits<WeightedGraph>::edge_descriptor WeightedEdgeDescriptor;

//undirected weighted graphs
typedef adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, int> > UndirectedGraph;
typedef graph_traits<UndirectedGraph>::vertex_descriptor UndirectedVertexDescriptor;
typedef graph_traits<UndirectedGraph>::edge_descriptor UndirectedEdgeDescriptor;


//definitions for visual struct
struct Node{
	Node() : x(0), y(0), dx(0), dy(0), r(0), g(0), b(0), delta(0.0), mtheta(0){}
	int cluster;
	float x, y, dx, dy;
	float r, g, b;
	float delta;
	float mtheta;
};
typedef std::vector<Node> Nodes;

struct Force{
	Force() : x(0), y(0), delta(0.0){}
	float x;
	float y;
	float delta;
};

typedef std::vector<Force> Forces;

typedef std::vector<std::vector<double> > MatrixScore;
typedef std::vector<std::vector<int> > Matrix;
typedef std::vector<Matrix> MatrixList;

typedef std::pair<double, Edge> Score;
typedef std::vector<Score> Scores;

Graph createGraph(std::ifstream &fin);
WeightedGraph createWeightedGraph(Graph g);
Matrix createAllPairsShortestPaths(WeightedGraph wg);
Edges newEdgeList(Graph g, Graph g2);
Matrix newList2Matrix(Graph g, Edges edges);
Scores SortScores(MatrixScore scores);
Scores SortIndividualScores(MatrixScore scores, int node);
MatrixScore initializeScore(Graph g);
MatrixScore ExtAdamicAdar(Graph g);
MatrixScore AdamicAdar(Graph g);
void normalize(MatrixScore scores, Graph g);
MatrixScore MultiTimeAdamicAdar(Graph g, Matrix distance, Matrix new_matrix);
bool comparisonSort(Score i, Score j);
float AvgDegreeAA(Graph g, Edges new_edges, Matrix distance, int d);
MatrixScore BCBasedLP(Graph g, Matrix distance, std::vector<double> v_centrality_vec, float beta);
Matrix multiplyMatrix(Matrix A, Matrix B);
std::vector<double> VertexBetweennessCentrality(Graph g, std::vector<double> &v_centrality_vec, std::vector<double> &e_centrality_map);
float BCBasedLPst(Graph g, Matrix distance, std::vector<double> v_centrality_vec, float beta, int s, int t);
Graph ClusterEdgeCentrality(Graph cg, double max_centrality);
Graph EdgesMST(Graph g, std::vector<double> ebc);
MatrixScore EdgeCentrality(Graph g, std::vector<double> ebc);
void NormalizePr(MatrixScore &scores, Matrix distance);
bool comparisonDouble(double i, double j);
std::vector<double> sortDouble(std::vector<double> scores);

class SocialGraph{

public:
	//set functions
	void initialize(int argc, char *argv[]);
	void setGraph(std::ifstream &fin);
	void setGraph2(std::ifstream &fin);
	void setWeightedGraph();
	void setNodes();
	void setNode(int i, int x, int y);

	//get functions
	Graph getGraph();
	Graph getGraph2();
	Graph getClusterGraph();
	Graph getClusterGraph2();
	Nodes getNodes();
	facebook::Name getName();
	WeightedGraph getWeightedGraph();
	Graph getSG();
	Graph getSG2();
	Edges getSGEdges();
	std::vector<double> getVBC();
	std::vector<double> getVBC2();
	MatrixScore getEBC();
	MatrixScore getEBC2();
	MatrixScore getProbability();
	MatrixScore getProbabilityAA();

	//spring force functions
	void repF(int i);
	void repulsiveForce();
	void attractiveForce(Graph);
	float limitForce();
	float springModel(Graph);

	//betweenness centrality based clustering
	void clusterColor(Graph tempg);

	void DepthFirst(std::vector<int> &visited, int max, std::vector<int> &hits, int &max_short, int t);
	void ShortestPaths(int s, std::vector<double> &centrality);

	void setMatrix();
	void addMatrix();
	void printMatrix();

	std::vector<double> Dependency(int s);
	MatrixScore AllDependency();

	int convertNode(int s);
	void AllScores();

	void allPairsShortestPaths();
	Matrix getDistance();

	void setAvgDegree();
	float getAvgDegree();

	Edges newEdgesList();
	void setNewEdges();
	std::vector<Edges> getNewEdges();
	int getNumNewEdges();

	void test();
	void showDetected();

	void selectSG(int s);
	void selectFutureSG(int s);
	float SpringSG(int s, char flag);


private:
	facebook::Name name;
	facebook::Map map, map2;
	Graph g, g2, g3, cg, cg2;
	WeightedGraph wg;
	Nodes nodes;
	std::vector<Edges> newEdges;
	int numNewEdges;
	Edges new_edges, new_edges2;
	Matrix new_matrix, new_matrix2;
	Matrix distance, distance2;
	MatrixList matrixlist;
	Forces forces;
	float avgDegree;
	MatrixScore probability, probabilityAA;
	std::vector<int> vizN;
	int countd2;
	Edges detected_edges;

	//for visualization
	Graph sg;
	Graph sg2;
	Edges sgEdges;
	std::vector<double> vertex_bc, vertex_bc2;
	MatrixScore edge_bc, edge_bc2;
};

#endif

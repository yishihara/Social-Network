#include <math.h>
#include <algorithm>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/bc_clustering.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/bind.hpp>
#include <iomanip>
#include <time.h>
#include <socialgraph.h>


void SocialGraph::setGraph(std::ifstream &fin){
	//read file
	facebook::Data data = facebook::readFile(fin, name);

	//create map for name and id
	map = facebook::createMap(data);

	//create graph
	g = facebook::createGraph(data, map);

	for(int i=0; i<num_vertices(g); i++){
		put(vertex_name, g, i, name[i]);
	}

	forces.resize(num_vertices(g));

	//create centrality graph
	cg = g;
}

void SocialGraph::setGraph2(std::ifstream &fin){
	facebook::Name t_name;
	facebook::Data data = facebook::readFile(fin, t_name);

	map2 = facebook::createMap(data);

	g2 = facebook::createGraph(data, map2);
}

void SocialGraph::setWeightedGraph(){
	WeightedGraph twg(num_vertices(g));
	
	property_map<WeightedGraph, edge_weight_t>::type weightmap = get(edge_weight, twg);
	EdgeIterator e, e_end;
	for(tie(e, e_end) = edges(cg); e != e_end; e++){
		WeightedEdgeDescriptor e1; bool inserted;
		tie(e1, inserted) = add_edge(source(*e, cg), target(*e, cg), twg);
		weightmap[e1] = 1;
	}
	
	wg = twg;
}

void SocialGraph::setNodes(){
	nodes.resize(num_vertices(g));
	for(int i=0; i<num_vertices(g); i++){
		nodes[i].x = (rand() % 200) -100;
		nodes[i].y = (rand() % 200) -100;
		nodes[i].dx = 0;
		nodes[i].dy = 0;
	}
}

void SocialGraph::setNode(int i, int x, int y){
	nodes[i].x = x;
	nodes[i].y = y;
}

Graph SocialGraph::getGraph(){
	return g;
}

Graph SocialGraph::getGraph2(){
	return g2;
}

WeightedGraph SocialGraph::getWeightedGraph(){
	return wg;
}

facebook::Name SocialGraph::getName(){
	return name;
}

Graph SocialGraph::getClusterGraph(){
	return cg;
}

Graph SocialGraph::getClusterGraph2(){
	return cg2;
}

Nodes SocialGraph::getNodes(){
	return nodes;
}

void SocialGraph::repF(int i){
	for(int j=i+1; j<nodes.size(); j++){
		float nx = nodes[i].x - nodes[j].x;
		float ny = nodes[i].y - nodes[j].y;
		float r = sqrt(pow(nx, 2) + pow(ny, 2));
		if(r > 0){
			float f = 900 / (r * r);
			forces[i].x += f * nx / r;
			forces[i].y += f * ny / r;
			forces[j].x -= f * nx / r;
			forces[j].y -= f * ny / r;
		}
	}
}

void SocialGraph::repulsiveForce(){
	for(int i=0; i<nodes.size(); i++){
		repF(i);
	}
}

void SocialGraph::attractiveForce(Graph temp_g){
	EdgeIterator e, e_end;
	for(tie(e, e_end) = edges(temp_g); e != e_end; e++){
		int s = source(*e, temp_g);
		int t = target(*e, temp_g);
		float nx = nodes[s].x - nodes[t].x;
		float ny = nodes[s].y - nodes[t].y;
		float r = sqrt(pow(nx, 2) + pow(ny, 2));
		if(r > 0){
			float f = -(r * r) / 30;
			forces[s].x += f * nx / r;
			forces[s].y += f * ny / r;
			forces[t].x -= f * nx / r;
			forces[t].y -= f * ny / r;

			/*
			forces[s].x += f * nx / (r * out_degree(s,temp_g));
			forces[s].y += f * ny / (r * out_degree(s,temp_g));
			forces[t].x -= f * nx / (r * out_degree(t,temp_g));
			forces[t].y -= f * ny / (r * out_degree(t,temp_g));
			*/
		}
	}
}

float SocialGraph::limitForce(){
	float dt = .2;
	float damping = 0.01;
	float kinetic = 0.0;

	for(int i=0; i<nodes.size(); i++){
		nodes[i].dx = (nodes[i].dx + dt * forces[i].x) * damping;
		nodes[i].dy = (nodes[i].dy + dt * forces[i].y) * damping;
		
		forces[i].x = 0;
		forces[i].y = 0;

		if(nodes[i].dx > 1) nodes[i].dx = 1;
		if(nodes[i].dx < -1) nodes[i].dx = -1;
		if(nodes[i].dy > 1) nodes[i].dy = 1;
		if(nodes[i].dy < -1) nodes[i].dy = -1;

/*
		if(nodes[i].dx > 10) nodes[i].dx = 10;
		if(nodes[i].dx < -10) nodes[i].dx = -10;
		if(nodes[i].dy > 10) nodes[i].dy = 10;
		if(nodes[i].dy < -10) nodes[i].dy = -10;
*/

		nodes[i].x = nodes[i].x + dt * nodes[i].dx;
		nodes[i].y = nodes[i].y + dt * nodes[i].dy;

		kinetic += sqrt(pow(nodes[i].dx, 2) + pow(nodes[i].dy, 2));
	}

	return kinetic;
}

float SocialGraph::springModel(Graph temp_g){
	float kinetic = 0.0;

	repulsiveForce();
	attractiveForce(temp_g);
	
	kinetic = limitForce();

	return kinetic;
}

void changeColor(float &r, float &g, float &b){
	if(r < 1.0){
		r += .5;
	}
	else{
		r = 0.0;
		if(g < 1.0){
			g += .5;
		}
		else{
			g = 0.0;
			if(b < .5){
				b += .5;
			}
			else
				b = 0.0;
		}
	}
}

void SocialGraph::clusterColor(Graph tempg){
	float r = 0.0;
	float g = 0.0;
	float b = 0.0;

	wg = ::createWeightedGraph(tempg);

	for(int i=0; i<num_vertices(wg); i++){
		for(int j=i+1; j<num_vertices(wg); j++){
			if(edge(i,j,wg).second == 0 && edge(j,i,wg).second != 0)
				remove_edge(j,i,wg);
			else if(edge(i,j,wg).second != 0 && edge(j,i,wg).second == 0)
				remove_edge(i,j,wg);
		}
	}

	std::vector<int> visited(num_vertices(wg), 0);
	int cluster = 0;
	for(int i=0; i<num_vertices(wg); i++){
		if(visited[i] == 0){
			cluster++;

			std::vector<WeightedVertexDescriptor> p(num_vertices(wg));
			std::vector<int> d(num_vertices(wg));
			dijkstra_shortest_paths(wg, i, predecessor_map(&p[0]).distance_map(&d[0]));

			changeColor(r,g,b);
			//std::cout << r << " " << g << " " << b << std::endl;
			for(int j=0; j<num_vertices(wg); j++){
				if(d[j] != INT_MAX){
					visited[j] = 1;
					nodes[j].r = r;
					nodes[j].g = g;
					nodes[j].b = b;
					nodes[j].cluster = cluster;
				}
			}
		}
	}
	std::cout << "number of clusters: " << cluster << std::endl;
}

bool ContainsNode(std::vector<int> visited, int node){
	for(int i=0; i<visited.size(); i++){
		if(visited[i] == node)
			return true;
	}
	return false;
}

void SocialGraph::DepthFirst(std::vector<int> &visited, int max, std::vector<int> &hits, int &max_short, int t){
	int back = visited.back();
	graph_traits<WeightedGraph>::out_edge_iterator e, e_end;

	for(tie(e, e_end) = out_edges(back, wg); e != e_end; e++){
		int node = target(*e, wg);
		
		if(ContainsNode(visited, node)) continue;

		if(node == t){
			visited.push_back(node);
			for(int i=1; i<visited.size()-1; i++){
				hits[visited[i]]++;
			}
			max_short++;
			
			int n = (int) visited.size() -1;
			visited.erase(visited.begin() + n);

			break;
		}

		visited.push_back(node);

		if(visited.size() <= max)
			DepthFirst(visited, max, hits, max_short, t);

		int n = (int)visited.size() - 1;
		visited.erase(visited.begin() + n);
	}
}

void SocialGraph::ShortestPaths(int s, std::vector<double> &centrality){
	std::vector<WeightedVertexDescriptor> p(num_vertices(wg));
	std::vector<int> d(num_vertices(wg));
	dijkstra_shortest_paths(wg, s, predecessor_map(&p[0]).distance_map(&d[0]));

	for(int i=0; i<num_vertices(wg); i++){
		//std::cout << i << " " << d[i] << std::endl;
		if(d[i] > 1 && d[i] != INT_MAX){
			std::vector<int> hits(num_vertices(wg), 0);
			int max_short = 0;

			std::vector<int> visited;
			visited.push_back(s);
			DepthFirst(visited, d[i], hits, max_short, i);

			for(int j=0; j<num_vertices(wg); j++){
				if(max_short > 0)
					centrality[j] += (double)hits[j]/(double)max_short;
			}
		}
	}
}

void SocialGraph::setMatrix(){
	Matrix matrix(num_vertices(g));

	for(int i=0; i<num_vertices(g); i++){
		std::vector<int> list(num_vertices(g));
		for(int j=0; j<num_vertices(g); j++){
			if(edge(i, j, g).second == 1)
				list[j] = 1;
			else
				list[j] = 0;
		}
		matrix[i] = list;
	}

	matrixlist.push_back(matrix);
}

void SocialGraph::addMatrix(){
	Matrix matrix1 = matrixlist[0];
	Matrix matrix2 = matrixlist[matrixlist.size() -1];
	Matrix temp(num_vertices(g));
	
	for(int i=0; i<num_vertices(g); i++){
		std::vector<int> list(num_vertices(g), 0);
		for(int j=0; j<num_vertices(g); j++){
			for(int k=0; k<num_vertices(g); k++){
				list[j] += matrix1[i][k] * matrix2[k][j];
			}
		}
		temp[i] = list;
	}

	matrixlist.push_back(temp);
}

void SocialGraph::printMatrix(){
	Matrix matrix = matrixlist[matrixlist.size() -1];

	for(int i=0; i<num_vertices(g); i++){
		for(int j=0; j<num_vertices(g); j++)
			std::cout << matrix[i][j] << " ";
		std::cout << std::endl;
	}
}

std::vector<double> SocialGraph::Dependency(int s){
	std::vector<WeightedVertexDescriptor> p(num_vertices(wg));
	std::vector<int> d(num_vertices(wg));
	dijkstra_shortest_paths(wg, s, predecessor_map(&p[0]).distance_map(&d[0]));

	std::vector<double> centrality(num_vertices(g), 0.0);

	for(int i=0; i<num_vertices(g); i++){
		//if(d[i] > 1 && d[i] != INT_MAX){
		if(d[i] > 1 && d[i] < 4){
			while(d[i] > matrixlist.size())
				addMatrix();

			//number of shortest paths
			int num = matrixlist[d[i]-1][s][i];

			std::vector<double> temp(num_vertices(g), 0.0);
			for(int j=1; j<d[i]; j++){
				for(int k=0; k<num_vertices(g); k++)
					temp[k] += matrixlist[j-1][s][k] * matrixlist[d[i]-1-j][k][i];
			}
			for(int k=0; k<num_vertices(g); k++)
				centrality[k] += (double)temp[k] / (double)num;
		}
	}

	return centrality;
}

MatrixScore SocialGraph::AllDependency(){
	MatrixScore all_centrality;
	
	for(int i=0; i<num_vertices(g); i++){
		std::vector<double> centrality = Dependency(i);
		all_centrality.push_back(centrality);
	}

	for(int i=0; i<num_vertices(g); i++){
		for(int j=0; j<num_vertices(g); j++){
			if(edge(i, j, g).second == 1)
				all_centrality[i][j] = 0;
		}
	}

	return all_centrality;
}

int SocialGraph::convertNode(int s){
	std::string source = map.left.find(s)->second;
	int t = map2.right.find(source)->second;
	return t;
}

Matrix SocialGraph::getDistance(){
	return distance;
}

void SocialGraph::setAvgDegree(){
	float avg=0.0;

	for(int i=0; i<num_vertices(g); i++){
		avg += out_degree(i, g);
	}
	avg = avg / num_vertices(g);

	avgDegree = avg;
}

float SocialGraph::getAvgDegree(){
	return avgDegree;
}

Edges SocialGraph::newEdgesList(){
	Edges edges;
	for(int i=0; i<num_vertices(g); i++){
		for(int j=i+1; j<num_vertices(g); j++){
			if(edge(i, j, g).second == 0){
				std::string source = map.left.find(i)->second;
				std::string target = map.left.find(j)->second;
				int s = map2.right.find(source)->second;
				int t = map2.right.find(target)->second;
				if(s > 0 && t > 0 && s < num_vertices(g2) + 1 && t < num_vertices(g2)){
					if(edge(s, t, g2).second == 1){
						Edge e(i, j);
						edges.push_back(e);
					}
				}
			}
		}
	}
	return edges;
}

void SocialGraph::setNewEdges(){
	Edges edges = newEdgesList();
	
	for(int i=0; i<edges.size(); i++){
		int s = edges[i].first;
		int t = edges[i].second;

		Edge e; e.first = s; e.second = t;
		if(distance[s][t] < num_vertices(g) && distance[s][t] > 0){
			while(newEdges.size() <= distance[s][t]){
				Edges es;
				newEdges.push_back(es);
			}
			newEdges[distance[s][t]].push_back(e);
		}
		else
			newEdges[0].push_back(e);
	}

	int count = 0;
	float allavg=0.0;
	for(int i=0; i<newEdges.size(); i++){
		count += newEdges[i].size();
		for(int j=0; j<newEdges[i].size(); j++){
			int s = newEdges[i][j].first;
			int t = newEdges[i][j].second;
			allavg += out_degree(s, g) + out_degree(t, g);
		}
	}
	std::cout << "new edges(all_new) = " << count << std::endl;
	numNewEdges = count;
	allavg = allavg/(count*2);
	std::cout << "avg degree(all_new) = " << allavg << std::endl;
	std::cout << std::endl;

	for(int i=0; i<newEdges.size(); i++){
		std::cout << "new edges (" << i << ") = " << newEdges[i].size() << std::endl;
		float avgnew = 0.0;
		for(int j=0; j<newEdges[i].size(); j++){
			int s = newEdges[i][j].first;
			int t = newEdges[i][j].second;
			avgnew += out_degree(s, g) + out_degree(t, g);
		}
		avgnew = avgnew/(newEdges[i].size()*2);
		std::cout << "avg degree (" << i << ") = " << avgnew << std::endl;
	}
	std::cout << std::endl;


	count = 0;
	int count2 = 0;
	float avg=0.0;
	float avg2=0.0;
	for(int i=0; i<newEdges[3].size(); i++){
		int flag = 0;
		int s = newEdges[3][i].first;
		int t = newEdges[3][i].second;
		for(int j=0; j<newEdges[2].size(); j++){
			int s2 = newEdges[2][j].first;
			int t2 = newEdges[2][j].second;
			if(s2 == s){
				if(distance[t2][t] == 1){
					flag = 1;
				}
			}
			else if(t2 == s){
				if(distance[s2][t] == 1){
					flag = 1;
				}
			}
			else if(s2 == t){
				if(distance[t2][s] == 1){
					flag = 1;
				}
			}
			else if(t2 == t){
				if(distance[s2][s] == 1){
					flag = 1;
				}
			}
			if(flag == 1){
				count++;
				avg += out_degree(s,g) + out_degree(t, g);
				break;
			}
		}
		if(flag == 0){
			count2++;
			avg2 += out_degree(s, g) + out_degree(t, g);
		}
	}
	std::cout << "new edges(3)<-(2) = " << count << std::endl;
	avg = avg / (count * 2);
	std::cout << "avg degree of (3)<-(2) = " << avg << std::endl;

	std::cout << "new edges(3)!(2) = " << count2 << std::endl;
	avg2 = avg2 / (count2 * 2);
	std::cout << "avg degree of (3)!(2) = " << avg2 << std::endl;
}

std::vector<Edges> SocialGraph::getNewEdges(){
	return newEdges;
}


//new functions
void SocialGraph::initialize(int argc, char *argv[]){
	std::ifstream fin;
	fin.open(argv[1]);
	facebook::Data data = facebook::readFile(fin, name);
	map = facebook::createMap(data);
	g = facebook::createGraph(data, map);
	for(int i=0; i<num_vertices(g); i++){
		put(vertex_name, g, i, name[i]);
	}
	forces.resize(num_vertices(g));
	fin.close();

	distance = ::createAllPairsShortestPaths(::createWeightedGraph(g));
	setAvgDegree();


	//initialize graph for t+1
	fin.open(argv[2]);
	Graph temp_g = ::createGraph(fin);
	g2 = g;
	fin.close();

	new_edges = ::newEdgeList(g, temp_g);
	new_matrix = ::newList2Matrix(g, new_edges);
	countd2=0;
	for(int i=0; i<new_edges.size(); i++){
		if(distance[new_edges[i].first][new_edges[i].second] < num_vertices(g)){
			if(distance[new_edges[i].first][new_edges[i].second] == 2){
				countd2++;
			}
			add_edge(new_edges[i].first, new_edges[i].second, g2);
			add_edge(new_edges[i].second, new_edges[i].first, g2);
		}
	}

	distance2 = ::createAllPairsShortestPaths(::createWeightedGraph(g2));

	setNodes();
	setAvgDegree();

	std::cout << "number of vertices: " << num_vertices(g) << std::endl;
	std::cout << "number of edges: " << num_edges(g)/2 << std::endl;
	std::cout << "number of new edges: " << new_edges.size() << std::endl;
	std::cout << "number of new edges at disance 2: " << countd2 << std::endl;

	#pragma omp parallel for
	for(int i=0; i<2; i++){
		if(i == 0)
			cg = ::ClusterEdgeCentrality(g, num_vertices(g));
		if(i == 1)
			cg2 = ::ClusterEdgeCentrality(g2, num_vertices(g2));
	}

	sg = g;
}

void SocialGraph::test(){
	Scores q;

	probabilityAA = ::AdamicAdar(g);
	q = ::SortScores(probabilityAA);

	//use this to count how many were the same as
	//betweenness centrality based method
	Edges daa;
	Edges dbc;
	
	int a=0;
	for(int i=0; i<new_edges.size(); i++){
	//for(int i=0; i<countd2; i++){
		int s = q[i].second.first;
		int t = q[i].second.second;
		if(edge(s, t, g2).second == 1){
			a++;
			daa.push_back(Edge(s,t));
		}
	}
	std::cout << "adamic-adar:" << a << std::endl;

	std::cout << "new edges of distance > 2: " << a << std::endl;
	std::cout << "average node degree: " << avgDegree << std::endl;

	std::vector<double> vbc(num_vertices(g), 0.0);
	std::vector<double> ebc(num_edges(g), 0.0);
	::VertexBetweennessCentrality(g, vbc, ebc);


	float beta = 0;
	probability = ::BCBasedLP(g, distance, vbc, beta);
	Scores sbc = ::SortScores(probability);
	a=0;
	int b=0;
	for(int i=0; i<new_edges.size(); i++){
	//for(int i=0; i<countd2; i++){
		int s = sbc[i].second.first;
		int t = sbc[i].second.second;
		if(edge(s, t, g2).second == 1 && (probability[s][t] > 0 || probability[t][s] > 0)){
			detected_edges.push_back(Edge(s,t));
			a++;
			dbc.push_back(Edge(s,t));
			if(distance[s][t] > 2)
				b++;
		}
	}
	std::cout << a << ":" << b << std::endl;

	//check how many were the same
	int overlap=0;
	for(int i=0; i<dbc.size(); i++){
		int s = dbc[i].first;
		int t = dbc[i].second;
		for(int j=0; j<daa.size(); j++){
			int s1 = daa[j].first;
			int t1 = daa[j].second;
			if((s1 == s && t1 == t) || (s1 == t && t1 == s))
				overlap++;
		}
	}
	std::cout << "overlap: " << overlap << "/" << dbc.size() << " = " << (float)overlap/dbc.size() << std::endl;

	vertex_bc = vbc;

	//calculate number of paths at what distance
	std::vector<int> count;
	for(int x=0; x<num_vertices(g); x++){
		for(int y=x+1; y<num_vertices(g); y++){
			if(distance[x][y] < num_vertices(g)){
				while(count.size() <= distance[x][y]){
					count.push_back(0);
				}
				count[distance[x][y]]++;
			}
		}
	}
	for(int i=0; i<count.size(); i++){
		std::cout << i << ":" << count[i] << std::endl;
	}

	::NormalizePr(probability, distance);
	::NormalizePr(probabilityAA, distance);

	edge_bc = ::EdgeCentrality(g, ebc);

	//finds probability of correct new edges detected
	//based on how many new edges detected of one individual
	std::vector<int> nEdges(num_vertices(g),0);
	for(int i=0; i<new_edges.size(); i++){
		int s = new_edges[i].first;
		int t = new_edges[i].second;
		if(distance[s][t] == 2){
			nEdges[s]++;
			nEdges[t]++;
		}
	}
	int necount=0; int necorrect=0;
	for(int i=0; i<num_vertices(g); i++){
		necount += nEdges[i];
		if(nEdges[i] > 0){
			Scores nscores = SortIndividualScores(probability, i);

			for(int j=0; j<nEdges[i]; j++){
				int s = nscores[j].second.first;
				int t = nscores[j].second.second;
				if(edge(s,t,g2).second == 1)
					necorrect++;
			}
		}
	}
	std::cout << necorrect << "/" << necount << std::endl;


	//change size of vertex and edge betweenness centrality
	//for visualization
	vbc = sortDouble(vertex_bc);
	double minbc = vbc[0.05*vbc.size()];
	double maxbc = vbc[0.95*vbc.size()];
	double diff = maxbc - minbc;
	for(int i=0; i<num_vertices(g); i++){
		if(vertex_bc[i] < minbc)
			vertex_bc[i] = 10;
		else if(vertex_bc[i] > maxbc)
			vertex_bc[i] = 25;
		else
			vertex_bc[i] = 10 + 15 * (vertex_bc[i] - minbc)/diff;
	}

	ebc = sortDouble(ebc);
	minbc = ebc[0];
	maxbc = ebc[ebc.size()*.5];
	diff = maxbc - minbc;

	for(int i=0; i<num_vertices(g); i++){
		for(int j=i+1; j<num_vertices(g); j++){
			if(edge(i,j,g).second == 1){
			if(edge_bc[i][j] < minbc){
				edge_bc[i][j] = 5;
				edge_bc[j][i] = 5;
			}
			else if(edge_bc[i][j] > maxbc){
				edge_bc[i][j] = 1;
				edge_bc[j][i] = 1;
			}
			else{
				edge_bc[i][j] = 5.0 - 4.0 * (edge_bc[i][j] - minbc)/diff;
				edge_bc[j][i] = edge_bc[i][j];
			}
			}
		}
	}

	//copy of above for future graph
	std::vector<double> vbc2(num_vertices(g2), 0.0);
	std::vector<double> ebc2(num_edges(g2), 0.0);
	::VertexBetweennessCentrality(g2, vbc2, ebc2);
	vertex_bc2 = vbc2;
	edge_bc2 = ::EdgeCentrality(g2, ebc2);

	vbc2 = sortDouble(vertex_bc2);
	minbc = vbc2[0.05*vbc2.size()];
	maxbc = vbc2[0.95*vbc2.size()];
	diff = maxbc - minbc;
	for(int i=0; i<num_vertices(g); i++){
		if(vertex_bc2[i] < minbc)
			vertex_bc2[i] = 10;
		else if(vertex_bc2[i] > maxbc)
			vertex_bc2[i] = 25;
		else
			vertex_bc2[i] = 10 + 15 * (vertex_bc2[i] - minbc)/diff;
	}

	ebc2 = sortDouble(ebc2);
	minbc = ebc2[0];
	maxbc = ebc2[ebc2.size()*.5];
	diff = maxbc - minbc;

	for(int i=0; i<num_vertices(g); i++){
		for(int j=i+1; j<num_vertices(g); j++){
			if(edge(i,j,g2).second == 1){
			if(edge_bc2[i][j] < minbc){
				edge_bc2[i][j] = 5;
				edge_bc2[j][i] = 5;
			}
			else if(edge_bc2[i][j] > maxbc){
				edge_bc2[i][j] = 1;
				edge_bc2[j][i] = 1;
			}
			else{
				edge_bc2[i][j] = 5.0 - 4.0 * (edge_bc2[i][j] - minbc)/diff;
				edge_bc2[j][i] = edge_bc2[i][j];
			}
			}
		}
	}

}

void SocialGraph::showDetected(){
	static int k=0;

	for(int i=0; i<num_vertices(g); i++){
		nodes[i].r = 0; nodes[i].g = 0; nodes[i].b = 0;
	}
	vizN.clear();

	k++;
	if(k == detected_edges.size())
		k = 0;
	
	int s = detected_edges[k].first;
	int t = detected_edges[k].second;
	nodes[s].b = 1.0; nodes[t].b = 1.0;

	vizN.push_back(s); vizN.push_back(t);

	for(int i=0; i<num_vertices(g); i++){
		if(distance[s][i] + distance[i][t] == distance[s][t]
		&& s != i && t != i){
			nodes[i].g = 1.0;
			vizN.push_back(i);
		}
	}

	std::cout << "changed " << k << ", normalized bc score = " << probability[s][t] << std::endl;
}

Graph SocialGraph::getSG(){
	return sg;
}

Graph SocialGraph::getSG2(){
	return sg2;
}

void SocialGraph::selectSG(int s){
	if(s != -1){

	Graph tempg(num_vertices(g));

	std::vector<int> nhbr;
	nhbr.push_back(s);
	for(int t=0; t<num_vertices(g); t++){
		if(distance[s][t] == 2 || distance[s][t] == 1){
			nhbr.push_back(t);

			//make nodes closer to selected node
			//for visualization
			/*
			nodes[t].x = nodes[s].x + (20.0*rand()/RAND_MAX - 10.0);
			nodes[t].y = nodes[s].y + (20.0*rand()/RAND_MAX - 10.0);
			*/
		}
	}
	
	for(int x=0; x<nhbr.size()-1; x++){
		for(int y=x+1; y<nhbr.size(); y++){
			int i = nhbr[x];
			int j = nhbr[y];
			if(edge(i,j,g).second == 1 && edge(i,j,tempg).second != 1){
				add_edge(i,j,tempg);
				add_edge(j,i,tempg);
			}
		}
	}

	sg = tempg;

	/*
	for(int i=0; i<num_vertices(g); i++){
		if(i == s){
			nodes[i].r = 1; nodes[i].g = 0; nodes[i].b = 0;
		}
		else{
			nodes[i].r = 0; nodes[i].g = 0; nodes[i].b = 0;
		}
	}
	*/

	}
}

void SocialGraph::selectFutureSG(int s){
	sg2 = sg;
	sgEdges.clear();

	int count=0;
	for(int i=0; i<num_vertices(g); i++){
		if(distance[s][i] == 2 && edge(s,i,g2).second == 1){
			add_edge(s,i,sg2);
			add_edge(i,s,sg2);
			count++;
			sgEdges.push_back(Edge(s,i));
		}
	}
	std::cout << count << " new edges" << std::endl;
}

float SocialGraph::SpringSG(int s, char flag){
	float kinetic = 0.0;

	if(s !=  -1){

	repulsiveForce();
	if(flag == 0)
		attractiveForce(sg);
	else if(flag == 1 || flag == 2){
		attractiveForce(sg);
		for(int i=0; i<num_vertices(g); i++){
			if(distance[s][i] == 2){
				float nx = nodes[s].x - nodes[i].x;
				float ny = nodes[s].y - nodes[i].y;
				float r = sqrt(pow(nx,2) + pow(ny,2));
				nx = nx / r;
				ny = ny / r;
				if(r>0 &&
				((flag==1)?(probability[s][i]):(probabilityAA[s][i])) >= 1){
				//if(r > 0 && (probability[s][i] >= 1)){
					float f = -2 * (r*r)/30;
					forces[s].x += f * nx;
					forces[s].y += f * ny;
					forces[i].x -= f * nx;
					forces[i].y -= f * ny;
				}
			}
		}
	}
	else if(flag == 3)
		attractiveForce(sg2);

	
	float dt = 1.0;
	float damping = 0.01;

	for(int i=0; i<nodes.size(); i++){
		nodes[i].dx = (nodes[i].dx + dt * forces[i].x) * damping;
		nodes[i].dy = (nodes[i].dy + dt * forces[i].y) * damping;
		
		forces[i].x = 0;
		forces[i].y = 0;

		if(nodes[i].dx > 1) nodes[i].dx = 1;
		if(nodes[i].dx < -1) nodes[i].dx = -1;
		if(nodes[i].dy > 1) nodes[i].dy = 1;
		if(nodes[i].dy < -1) nodes[i].dy = -1;


		if(i != s){
			nodes[i].x = nodes[i].x + dt * nodes[i].dx;
			nodes[i].y = nodes[i].y + dt * nodes[i].dy;
		}

		kinetic += sqrt(pow(nodes[i].dx, 2) + pow(nodes[i].dy, 2));
	}

	}


	return kinetic;

}

Edges SocialGraph::getSGEdges(){
	return sgEdges;
}

std::vector<double> SocialGraph::getVBC(){
	return vertex_bc;
}

std::vector<double> SocialGraph::getVBC2(){
	return vertex_bc2;
}

MatrixScore SocialGraph::getEBC(){
	return edge_bc;
}

MatrixScore SocialGraph::getEBC2(){
	return edge_bc2;
}

MatrixScore SocialGraph::getProbability(){
	return probability;
}

MatrixScore SocialGraph::getProbabilityAA(){
	return probabilityAA;
}

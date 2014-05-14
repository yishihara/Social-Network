#include <socialgraph.h>
#include <math.h>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/bc_clustering.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

Graph createGraph(std::ifstream &fin){
	facebook::Name name;
	facebook::Data data = facebook::readFile(fin, name);

	facebook::Map map = facebook::createMap(data);

	Graph g = facebook::createGraph(data, map);
	
	return g;
}

WeightedGraph createWeightedGraph(Graph g){
	WeightedGraph wg(num_vertices(g));
	
	property_map<WeightedGraph, edge_weight_t>::type weightmap = get(edge_weight, wg);
	EdgeIterator e, e_end;
	for(tie(e, e_end) = edges(g); e != e_end; e++){
		WeightedEdgeDescriptor e1; bool inserted;
		tie(e1, inserted) = add_edge(source(*e, g), target(*e, g), wg);
		weightmap[e1] = 1;
	}
	
	return wg;
}

Matrix createAllPairsShortestPaths(WeightedGraph wg){
	int V = num_vertices(wg);
	Matrix D(V, std::vector<int>(V));
	johnson_all_pairs_shortest_paths(wg, D);
	
	return D;
}

//have to corresponding vertices
Edges newEdgeList(Graph g, Graph g2){
	Edges edges;
	for(int i=0; i<num_vertices(g); i++){
		for(int j=i+1; j<num_vertices(g); j++){
			if(edge(i, j, g).second == 0){
				if(edge(i, j, g2).second == 1){
					Edge e(i, j);
					edges.push_back(e);
				}
			}
		}
	}
	return edges;
}

Matrix newList2Matrix(Graph g, Edges edges){
	Matrix matrix(num_vertices(g));
	for(int i=0; i<num_vertices(g); i++){
		std::vector<int> temp(num_vertices(g), 0);
		matrix[i] = temp;
	}
	for(int i=0; i<edges.size(); i++){
		matrix[edges[i].first][edges[i].second] = 1;
		matrix[edges[i].second][edges[i].first] = 1;
	}
	return matrix;
}

MatrixScore initializeScore(Graph g){
	MatrixScore scores(num_vertices(g));

	for(int i=0; i<num_vertices(g); i++){
		std::vector<double> score(num_vertices(g), 0.0);
		scores[i] = score;
	}
	return scores;
}

std::vector<double> sortDouble(std::vector<double> scores){
	std::vector<double> temps;
	for(int i=0; i<scores.size(); i++){
		if(scores[i] > 0)
			temps.push_back(scores[i]);
	}
	std::sort(temps.begin(), temps.end());
	return temps;
}

bool comparisonSort(Score i, Score j){
	if(i.first < j.first) return false;
	if(j.first < i.first) return true;
	return j.second.first < i.second.first;
}

Scores SortScores(MatrixScore mscores){
	Scores scores;

	for(int i=0; i<mscores.size(); i++){
		for(int j=i+1; j<mscores.size(); j++){
			scores.push_back(Score(mscores[i][j], Edge(i,j)));
		}
	}
	std::sort(scores.begin(), scores.end(), comparisonSort);

	return scores;
}

Scores SortIndividualScores(MatrixScore mscores, int node){
	Scores scores;

	for(int i=0; i<mscores.size(); i++){
		scores.push_back(Score(mscores[node][i], Edge(node, i)));
	}

	std::sort(scores.begin(), scores.end(), comparisonSort);

	return scores;
}

float AvgDegree(Graph g){
	float avg=0.0;

	for(int i=0; i<num_vertices(g); i++)
		avg += out_degree(i,g);

	avg = avg / num_vertices(g);
	
	return avg;
}

MatrixScore ExtAdamicAdar(Graph g){
	Graph temp_g = g;
	WeightedGraph wg = createWeightedGraph(temp_g);
	Matrix distance = createAllPairsShortestPaths(wg);
	MatrixScore scores = initializeScore(g);

	for(int i=0; i<num_vertices(g); i++){
		for(int j=0; j<num_vertices(g); j++){
			if(edge(i,j,g).second == 1)
				scores[i][j] = 1;
		}
	}

	float avg = AvgDegree(g);
	float thresh = (float)avg/log(2);
	
	for(int k=0; k<2; k++){

	for(int x=0; x<num_vertices(temp_g); x++){
		for(int y=0; y<num_vertices(temp_g); y++){
			if(distance[x][y] == 2){
				for(int i=0; i<num_vertices(g); i++){
					if(distance[x][i] == 1 && distance[i][y] == 1){
						scores[x][y] += (double)(scores[x][i] * scores[y][i])/log(out_degree(i,g));
					}
				}
				if(scores[x][y] > thresh)
					scores[x][y] = 1.0;
				else
					scores[x][y] = scores[x][y] / thresh;
				add_edge(x,y,temp_g);
			}
		}
	}
	wg = createWeightedGraph(temp_g);
	distance = createAllPairsShortestPaths(wg);

	}

	return scores;
}

MatrixScore AdamicAdar(Graph g){
	Matrix distance = createAllPairsShortestPaths(createWeightedGraph(g));
	MatrixScore scores = initializeScore(g);

	for(int x=0; x<num_vertices(g); x++){
		#pragma omp parallel for
		for(int y=x+1; y<num_vertices(g); y++){
			if(distance[x][y] == 2){
				for(int i=0; i<num_vertices(g); i++){
					if(x!=i && y!=i && out_degree(i,g) > 1){
						int d = distance[x][i] + distance[y][i];
						if(d == distance[x][y])
							scores[x][y] += (double)1/log(out_degree(i,g));
					}
				}
			}
		}
	}

	return scores;
}

void normalize(MatrixScore scores, Graph g){
	//mean
	int count=0;
	double avg=0.0;
	for(int x=0; x<num_vertices(g); x++){
		for(int y=x+1; y<num_vertices(g); y++){
			if(edge(x,y,g).second == 0){
				avg += scores[x][y];
				count++;
			}
		}
	}
	avg = avg/count;

	double dev=0.0;
	for(int x=0; x<num_vertices(g); x++){
		for(int y=x+1; y<num_vertices(g); y++){
			if(edge(x,y,g).second == 0){
				dev += pow(scores[x][y] - avg, 2);
			}
		}
	}
	dev = sqrt(dev/count);
	std::cout << avg << " " << dev << std::endl;
}

MatrixScore MultiTimeAdamicAdar(Graph g, Matrix distance, Matrix new_matrix){
	MatrixScore scores(num_vertices(g));
	for(int i=0; i<num_vertices(g); i++){
		std::vector<double> score(num_vertices(g), 0.0);
		scores[i] = score;
	}

	std::vector<int> is_new(num_vertices(g), 0);
	for(int i=0; i<num_vertices(g); i++){
		for(int j=0; j<num_vertices(g); j++){
			if(new_matrix[i][j] == 1){
				is_new[i] = 1;
				break;
			}
		}
	}

	for(int s=0; s<num_vertices(g); s++){
		for(int t=s+1; t<num_vertices(g); t++){
			if(distance[s][t] == 2){
				for(int i=0; i<num_vertices(g); i++){
					if(distance[s][i] == 1 && distance[t][i] == 1){
						float si = (new_matrix[s][i] == 1)?1.5:1;
						float ti = (new_matrix[t][i] == 1)?1.5:1;
						scores[s][t] += (si*ti)/log10(out_degree(i,g));
					}
				}
			}
		}
	}

	return scores;
}

float AvgDegreeAA(Graph g, Edges new_edges, Matrix distance){
	int count=0; float temp_avg=0.0;
	for(int i=0; i<new_edges.size(); i++){
		int s = new_edges[i].first;
		int t = new_edges[i].second;

		if(distance[s][t] == 2){
			for(int j=0; j<num_vertices(g); j++){
				if(distance[s][j] + distance[j][t] == distance[s][t]){
					temp_avg += out_degree(j,g);
					count++;
				}
			}
		}
	}

	return temp_avg/count;
}

Matrix createMatrix(Graph g){
	Matrix matrix(num_vertices(g));

	#pragma omp parallel for
	for(int i=0; i<num_vertices(g); i++){
		std::vector<int> list(num_vertices(g));
		for(int j=0; j<num_vertices(g); j++){
			if(edge(i,j,g).second == 1)
				list[j] = 1;
			else
				list[j] = 0;
		}
		matrix[i] = list;
	}

	return matrix;
}

Matrix multiplyMatrix(Matrix A, Matrix B){
	Matrix matrix(A.size());

	#pragma omp parallel for
	for(int i=0; i<A.size(); i++){
		std::vector<int> list(A.size(), 0);
		for(int j=0; j<A.size(); j++){
			for(int k=0; k<A.size(); k++){
				list[j] += A[i][k] * B[k][j];
			}
		}
		matrix[i] = list;
	}

	return matrix;
}

std::vector<double> VertexBetweennessCentrality(Graph g, std::vector<double> &v_centrality_vec, std::vector<double> &e_centrality_vec){
	StdEdgeIndexMap my_e_index;
	EdgeIndexMap e_index(my_e_index);
	int i=0;
	BGL_FORALL_EDGES(edge, g, Graph){
		my_e_index.insert(std::pair<EdgeDescriptor, int>(edge,i));
		++i;
	}

	//std::vector<double>e_centrality_vec(num_edges(g), 0.0);
	iterator_property_map<std::vector<double>::iterator, EdgeIndexMap> e_centrality_map(e_centrality_vec.begin(), e_index);

	VertexIndexMap v_index = get(vertex_index, g);
	//std::vector<double> v_centrality_vec(num_vertices(g), 0.0);
	iterator_property_map<std::vector<double>::iterator, VertexIndexMap> v_centrality_map(v_centrality_vec.begin(), v_index);
	brandes_betweenness_centrality(g, v_centrality_map, e_centrality_map);
	//relative_betweenness_centrality(g, v_centrality_map);

	return v_centrality_vec;
}

MatrixScore BCBasedLP(Graph g, Matrix distance, std::vector<double> v_centrality_vec, float beta){
	MatrixScore score = initializeScore(g);
	//Matrix distance = createAllPairsShortestPaths(createWeightedGraph(g));

	#pragma omp parallel for
	for(int x=0; x<num_vertices(g); x++){
		for(int y=x+1; y<num_vertices(g); y++){
		/*
			if(distance[x][y] == 3){
		
			std::vector<float> z;
			double avg=0.0;
			for(int i=0; i<num_vertices(g); i++){
				if((distance[x][i] + distance[i][y] == distance[x][y]) && x!=i && y!=i){
					z.push_back(v_centrality_vec[i]);
					avg += v_centrality_vec[i];
				}
			}

			avg = avg/z.size();

			for(int i=0; i<z.size(); i++){
				score[x][y] += pow(avg - z[i],2);
			}

			if(score[x][y] != 0)
				score[x][y] = sqrt(score[x][y]);
			else{
				double tavg = (v_centrality_vec[x] + v_centrality_vec[y])/2;
				tavg = pow(v_centrality_vec[x]-tavg,2) + pow(v_centrality_vec[y]-tavg,2);
				score[x][y] = sqrt(tavg);
			}

			}
		*/

			for(int i=0; i<num_vertices(g); i++){
				if((distance[x][i] + distance[i][y] == distance[x][y]) 
				&& x!=i && y!=i){
					//score[x][y] += 1/v_centrality_vec[i];
					//score[x][y] += 1/pow(v_centrality_vec[i],2);
					//score[x][y] += 1/log(v_centrality_vec[i]);
					//score[x][y] += 1/log(1+v_centrality_vec[i]);
					score[x][y] += 1/(pow(log(1+v_centrality_vec[i]),2));
				}
			}

			if(distance[x][y] >= 2)
				score[x][y] *= pow(beta, distance[x][y]-2);
		}
	}

	return score;
}

float BCBasedLPst(Graph g, Matrix distance, std::vector<double> v_centrality_vec, float beta, int s, int t){

	float score=0.0;
	int count=0;
	std::vector<float> z;
	double avg=0.0;
	for(int i=0; i<num_vertices(g); i++){
		if((distance[s][i] + distance[i][t] == distance[s][t]) && s!=i && t!=i){
			//score += v_centrality_vec[i];
			z.push_back(v_centrality_vec[i]);
			avg += v_centrality_vec[i];
			//count++;
		}
	}

	avg = avg/z.size();

	return avg;

	/*
	for(int i=0; i<z.size(); i++){
		score += pow(avg - z[i],2);
	}

	score = sqrt(score/z.size());

	return score;
	*/
}

/*
	Matrix distance = createAllPairsShortestPaths(createWeightedGraph(g));

	MatrixScore score = initializeScore(g);

	MatrixList ml;
	ml.push_back(createMatrix(g));

	for(int i=0; i<num_vertices(g); i++){
		for(int j=i+1; j<num_vertices(g); j++){
			if(distance[i][j] > 1 && distance[i][j] < num_vertices(g)){
				while(distance[i][j] > ml.size())
					ml.push_back(multiplyMatrix(ml[0], ml[ml.size()-1]));

				int num = ml[distance[i][j]-1][i][j];

				for(int k=0; k<num_vertices(g); k++){
					if(distance[i][k] + distance[k][j] == distance[i][j] && k!=i && k!=j){
						score[i][j] += (float)(ml[distance[i][k]-1][i][k] * ml[distance[k][j]-1][k][j])/num;
					}
				}
			}
		}
	}

	return score;
*/

Graph ClusterEdgeCentrality(Graph g, double max_centrality){
	Graph cg = g;
	StdEdgeIndexMap my_e_index;
	EdgeIndexMap e_index(my_e_index);
	int i=0;
	BGL_FORALL_EDGES(edge, cg, Graph){
		my_e_index.insert(std::pair<EdgeDescriptor, int>(edge, i));
		++i;
	}

	std::vector<double>e_centrality_vec(num_edges(cg), 0.0);
	iterator_property_map<std::vector<double>::iterator, EdgeIndexMap> e_centrality_map(e_centrality_vec.begin(), e_index);

	bc_clustering_threshold<double>terminate(max_centrality, cg, false);
	betweenness_centrality_clustering(cg, terminate, e_centrality_map);

	for(int i=0; i<num_vertices(cg); i++){
		for(int j=0; j<num_vertices(cg); j++){
			if(edge(i, j, cg).second == 1){
			}
		}
	}
	
	return cg;
}

Graph EdgesMST(Graph g, std::vector<double> ebc){
	MatrixScore weights = initializeScore(g);

	int k=0;
	BGL_FORALL_EDGES(e, g, Graph){
		int s = source(e, g);
		int t = target(e, g);
		if(s < t)
			weights[s][t] += ebc[k];
		else if(t < s)
			weights[t][s] += ebc[k];
		k++;
	}
	
	for(int i=0; i<num_vertices(g); i++){
		for(int j=i+1; j<num_vertices(g); j++){
			weights[i][j] = weights[i][j]/2;
			weights[j][i] = weights[i][j];
		}
	}

	UndirectedGraph ug(num_vertices(g));

	property_map<UndirectedGraph, edge_weight_t>::type weightmap = get(edge_weight,ug);

	BGL_FORALL_EDGES(e, g, Graph){
		int s = source(e, g);
		int t = target(e, g);
		if(edge(s,t,ug).second != 1){
			UndirectedEdgeDescriptor ue;
			bool inserted;
			tie(ue, inserted) = add_edge(s,t,ug);
			weightmap[ue] = (size_t)weights[s][t];
		}
	}

	std::vector<UndirectedEdgeDescriptor> spanning_tree;
	kruskal_minimum_spanning_tree(ug, std::back_inserter(spanning_tree));

	Graph tempg(num_vertices(g));
	for(int i=0; i<spanning_tree.size(); i++){
		int s = source(spanning_tree[i], ug);
		int t = target(spanning_tree[i], ug);

		add_edge(s,t,tempg);
		add_edge(t,s,tempg);
	}

	return tempg;
}

MatrixScore EdgeCentrality(Graph g, std::vector<double> ebc){
	MatrixScore score = initializeScore(g);

	int k=0;
	BGL_FORALL_EDGES(e, g, Graph){
		int s = source(e, g);
		int t = target(e, g);
		score[s][t] = ebc[k];
		score[t][s] = ebc[k];
		k++;
	}

	return score;
}

void NormalizePr(MatrixScore &scores, Matrix distance){
	//mean
	int count=0;
	double avg=0.0;
	for(int x=0; x<scores.size(); x++){
		for(int y=x+1; y<scores[x].size(); y++){
			if(distance[x][y] == 2){
				avg += scores[x][y];
				count++;
			}
		}
	}

	avg = avg/count;

	double dev=0.0;
	for(int x=0; x<scores.size(); x++){
		for(int y=x+1; y<scores[x].size(); y++){
			if(distance[x][y] == 2){
				dev += pow(scores[x][y] - avg, 2);
			}
		}
	}
	dev = sqrt(dev/count);

	count =0;

	for(int x=0; x<scores.size(); x++){
		for(int y=x+1; y<scores[x].size(); y++){
			if(distance[x][y] == 2){
				scores[x][y] = (scores[x][y] - avg)/dev;
				scores[y][x] = scores[x][y];
				if(scores[x][y] >= 1.0)
					count++;
			}
		}
	}
	std::cout << "number of scores above 1: " << count << std::endl;
}

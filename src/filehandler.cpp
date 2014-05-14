#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <filehandler.h>
#include <graph.h>

namespace facebook{

	Data readFile(std::ifstream &fin, Name &name){
		Data data;

		std::string line;
		std::getline(fin, line);

		//number of lines
		const int N = atoi(line.c_str());

		//resize vector
		data.resize(N);
		
		//retrieve name data
		for(int i=0; i<N; i++){
			std::getline(fin, line);
			std::string text = "";
			std::string id = "";
			int flag = 0;

			for(int j=0; j<=line.length(); j++){
				if(j == line.length()){
					flag = 0;
					id = "";
				}
				else if(line.compare(j, 1, ":") == 0){
					if(line.compare(j+1, 1, ":") == 0){
						flag = 1;
						name.push_back(text);
						text = "";
					}
				}
				else{
					if(flag == 0)
						text += line[j];
					else if(flag == 1)
						id += line[j];
				}
			}
		}
		
		//retrieve network data
		for(int i=0; i<N; i++){
			std::getline(fin, line);
			std::string text = "";

			for(int j=0; j<=line.length(); j++){
				if(line.compare(j, 1, " ") == 0 || j == line.length()){
					data[i].push_back(text);
					text = "";
				}
				else{
					text += line[j];
				}
			}
		}
		return data;
	}
	
	Map createMap(Data data){
		Map map;
		for(int i=0; i<data.size(); i++){
			
			map.insert(Map::value_type(i, data[i][0]));
		}
		return map;
	}

	socialgraph::Edges createEdges(Data data, Map map){
		socialgraph::Edges edges;

		for(int i=0; i<data.size(); i++){
			for(int j=1; j<data[i].size(); j++){
				edges.push_back(socialgraph::Edge(i, map.right.find(data[i][j])->second));
			}
		}

		return edges;
	}

	socialgraph::Graph createGraph(Data data, Map map){
		socialgraph::Graph g(data.size());

		socialgraph::Edges edges = createEdges(data, map);

		//add edges to graph
		for(int i=0; i<edges.size(); i++){
			boost::add_edge(edges[i].first, edges[i].second, g);
		}

		
		for(int i=0; i<data.size(); i++){
			put(vertex_name, g, vertex(i, g), data[i][0]);
		}
		

		return g;
	}

}

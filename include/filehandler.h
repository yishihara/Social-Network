#ifndef FILEHANDLER_H
#define FILEHANDLER_H

#include <vector>
#include <string>
#include <fstream>
#include <boost/bimap/bimap.hpp>
#include "graph.h"

namespace facebook{

	typedef std::vector<std::vector<std::string> > Data;

	typedef std::map<std::string, int> ID2INT;
	
	typedef boost::bimaps::bimap<int, std::string> Map;
	
	typedef std::vector<std::string> Name;

	Data readFile(std::ifstream &fin, Name &name);

	Map createMap(Data data);

	socialgraph::Edges createEdges(Data data, Map map);

	socialgraph::Graph createGraph(Data data, Map map);

}

#endif

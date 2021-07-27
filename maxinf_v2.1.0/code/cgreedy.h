#ifndef cgreedy_h__
#define cgreedy_h__

#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <string>
#include <cassert>
#include <vector>
#include "greedy.h"
#include "graph.h"

/// Topic-aware:
///    Greedy algorihtm with candidates
class CGreedy:
	public Greedy
{
public:
	CGreedy();
	void BuildWithCandidates(int num, std::set<int>& candidates, Graph& gf, ICascade& cascade);

	static void initialize_gwc(int argc, std::vector<std::string>&, std::set<int>& candidates, IGraph& gf);
};

#endif // cgreedy_h__

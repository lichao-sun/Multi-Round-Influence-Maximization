#ifndef RIM_CRGREEDY_H
#define RIM_CRGREEDY_H

#include <vector>
#include <set>
#include <string>
#include "algo_base.h"
#include "common.h"
#include "graph.h"
#include "cascade.h"

/// Greedy algorithm with lazy-forward optimization
class RIM_CRGreedy
	: public AlgoBase
{
protected:
	std::string file;

public:
	RIM_CRGreedy();

	void Build(IGraph& gf, int k, int rounds, ICascade& cascade);
	void BuildRanking(IGraph& gf, int k, ICascade& cascade);
	void BuildFromFile(IGraph& gf, const char* filename);
};

#endif


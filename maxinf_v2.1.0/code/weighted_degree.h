#ifndef WEIGHTED_DEGREE_H
#define WEIGHTED_DEGREE_H

#include <vector>
#include <string>
#include "graph.h"
#include "algo_base.h"
#include "common.h"

/// Algorithm to select k seeds based on their weighted degree.
class WeightedDegree
	: public AlgoBase
{
private:
	std::string file;

	void qsort_degree(int h, int t);

public:
	WeightedDegree();
	void Build(Graph& gf, int k = SET_SIZE);
	void BuildFromFile(Graph& gf);
};

#endif


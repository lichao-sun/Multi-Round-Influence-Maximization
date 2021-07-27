#ifndef DEGREE_H
#define DEGREE_H

#include <vector>
#include <string>
#include "common.h"
#include "graph.h"
#include "algo_base.h"

/// Algorithm selectes seeds with largest degree
class Degree
	: public AlgoBase
{
protected:
	std::string file;

	void qsort_degree(int h, int t);

public:
	Degree();

	void Build(Graph& gf, int k = SET_SIZE);
	void BuildFromFile(Graph& gf);
};

#endif


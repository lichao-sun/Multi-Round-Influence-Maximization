#ifndef DEGREE_DISCOUNT_IC_H
#define DEGREE_DISCOUNT_IC_H

#include <string>
#include <vector>
#include "common.h"
#include "graph.h"
#include "algo_base.h"
/// Algorithm selectes seeds by using largest degree with discount rate

class DegreeDiscount_IC
	: public AlgoBase
{
protected:
	std::string file;

//	static void qsort_degree(int h, int t);

public:
	DegreeDiscount_IC();

	void Build(Graph& gf, int k = SET_SIZE, double ratio = 0.01);
	void BuildFromFile(Graph& gf);
};

#endif


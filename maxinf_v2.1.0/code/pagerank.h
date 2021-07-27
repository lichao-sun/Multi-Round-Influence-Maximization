#ifndef PAGERANK_H
#define PAGERANK_H

#include <vector>
#include <string>

#include "common.h"
#include "graph.h"
#include "algo_base.h"

/// Pagerank algorithm
class Pagerank
	: public AlgoBase
{
protected:
	std::string file;
	std::vector<double> dp;  // pagerank values
	std::vector<int> dd;	 // neighbourCount

public:
	Pagerank();
	virtual double Build(Graph& gf, int num, double dampen = 0.15);
	void BuildFromFile(Graph& gf);
	int GetMax(int round, std::vector<double>& values);
};

/// Use pagerank to generate a graph
class PagerankRegen
	: protected Pagerank
{
public:
	double Build(Graph& gf, std::string& graph_file, int num, double dampen = 0.15);

protected:
	double NAdefault(double v, double defaultVal);
	double roundInInterval(double v, double left, double right);
	bool almostZero(double v, double threshold);
};


#endif

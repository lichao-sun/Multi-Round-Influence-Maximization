#ifndef MIA_H
#define MIA_H

#include <vector>
#include <string>
#include "common.h"
#include "graph.h"
#include "cascade.h"
#include "algo_base.h"

/// MIA algorithm.
/// See: Wei Chen, Chi Wang, and Yajun Wang. Scalable influence 
///    maximization for prevalent viral marketing in large-scale social networks.
///    In Proceedings of the 16th ACM SIGKDD Conference on Knowledge
///    Discovery and Data Mining (KDD'2010), Washington DC, U.S.A., July 2010.
class MIA
	: public AlgoBase
{
protected:
	char file[STR_LEN];
	std::vector<int> dd;
	double longest;
	std::vector<double> dp;
	std::vector<bool> used;
	std::vector<double*> self;
	std::vector<int> lastupdate;
	std::vector<double *>delta;
	std::vector<int *>children, path;
	int *S, *numchild, *queue;
	double *distance, *b;
	int *heap;
	int *childlist, *oldchildlist, *parent;
	std::vector<bool *> validlist;
	std::vector<int> *childnum;
	std::vector<double> *allb;

protected:
	void Initialize(int n, int k);

public:
	MIA();

	double Build(Graph& gf, int num, int bound);
	double Build(Graph& gf, int num, int k, int bound, 
		ICascade& simuCascade,
		ICascade& simuCascadeFast);
	void BuildFromFile(Graph& gf, int bound);
	int GetMax(int round);
	int GetMax0(int round);

	int generateMIAfrom(Graph& gf, int round, int node);
	int generateMIAto(Graph& gf, int node);
	int generateMIAto0(int node);
	int count(Graph& gf, int node);
	char* filename(int bound);
};

#endif


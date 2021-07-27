#ifndef SPT_new_H
#define SPT_new_H


#include "common.h"
#include "graph.h"
#include "cascade.h"
#include "algo_base.h"

/// PMIA algorithm.
/// See: Wei Chen, Chi Wang, and Yajun Wang. Scalable influence 
///    maximization for prevalent viral marketing in large-scale social networks.
///    In Proceedings of the 16th ACM SIGKDD Conference on Knowledge
///    Discovery and Data Mining (KDD'2010), Washington DC, U.S.A., July 2010.
class PMIA
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
	void Initialize(int totaln, int k);

public:
	PMIA();
	double Build(Graph& gf, int num, int bound);
	double Build(Graph& gf, int num, int k, int bound, 
			ICascade& simuCascade,
			ICascade& simuCascadeFast);
	void BuildFromFile(Graph& gf, int bound);

	int GetMax(int round);
	int GetMax0(int round);
	int generateSPT_newfrom(Graph& gf, int round, int node);
	int generateSPT_newto(Graph& gf, int node);
	int generateSPT_newto0(int node);
	int count(Graph& gf, int node);
	char* filename(int bound);
};

#endif


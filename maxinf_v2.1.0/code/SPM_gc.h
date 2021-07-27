#ifndef SPM_GC_H
#define SPM_GC_H

#include <vector>
#include "common.h"
#include "graph.h"
#include "cascade.h"

/// Estimate cascade spread, which is used for MIA and PMIA
class SPM_gc :
	public ICascade
{
private:
	int	n, m;
	int	targetSize;
	double	resultSize;
	std::vector<int> target;
	Graph* gf;

public:
	SPM_gc();
	void Build(Graph& gf);
	void SetTarget(int size, int set[]);
	double Run(int itr_num, int size, int set[]);
};

#endif

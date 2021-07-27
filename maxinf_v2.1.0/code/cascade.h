#ifndef __cascade_t_h__
#define __cascade_t_h__

#include <iostream>
#include<vector>
#include<set>
#include <cassert>
#include <algorithm>
#include "common.h"
#include "graph.h"
#include "algo_base.h"


/// Interface for cascade (abstract base class)
class ICascade 
{
public:
	virtual double Run(int num_iter, int size, int set[]) = 0;
	virtual double RRun(int num_iter, int size, int set[], bool act_nodes_arr[]) = 0;
	virtual double ERRun(int num_iter, int size, int set[], bool **act_nodes_arr) = 0;
	virtual double CRRun(int num_iter, int size, int set[], bool **act_nodes_arr, int T) = 0;
	virtual double ActRun(int top, std::set<int>& cur_act_nodes, std::vector<std::vector<int>>& seeds_list) = 0;
};

/// A general / template class of cascade
template<class TGraph=Graph>
class CascadeT : 
	public ICascade
{
public:
	typedef CascadeT<TGraph> self_type;
	typedef TGraph graph_type;

protected:
	int	n;
	int m;
	MIRandom random;
	TGraph* gf;

	virtual void _Build(TGraph& gf) {
        MI_STATIC_ASSERT(has_mem_GetN<graph_type>::value, "graph_type should has member GetN");
        MI_STATIC_ASSERT(has_mem_GetM<graph_type>::value, "graph_type should has member GetM");
        
		this->gf = &gf;
		this->n = gf.GetN();
		this->m = gf.GetN();
	}
public:
	CascadeT() : gf(NULL) {}
};



#endif // __cascade_t_h__

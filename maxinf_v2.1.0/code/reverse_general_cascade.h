#ifndef ReverseGeneralCascade_h__
#define ReverseGeneralCascade_h__


#include <random>
#include <vector>
#include <set>
#include <cassert>
#include "common.h"
#include "graph.h"
#include "mi_random.h"

/// Container for RR set.
typedef std::vector<int> RRVec;

/// Template class for Reverse General Cascade
template<class TGraph=Graph>
class ReverseGCascadeT
{
public:
	typedef TGraph graph_type;

protected:
	int	n, m;
	MIRandom random;
	graph_type* gf;

	

public:
	ReverseGCascadeT() : n(0), m(0), gf(NULL) {}

public:
	void Build(TGraph& gf)
	{
		MI_STATIC_ASSERT(has_mem_GetN<graph_type>::value, "graph_type should has member GetN");
		MI_STATIC_ASSERT(has_mem_GetM<graph_type>::value, "graph_type should has member GetM");

		this->gf = &gf;
		this->n = gf.GetN();
		this->m = gf.GetM();
	}

	int GenRandomNode()
	{
		int id = random.RandInt(0, n-1); // id: 0 ~ n-1
		return id;
	}

	double ReversePropagate(int num_iter, int target,
						std::vector< RRVec >& outRRSets,
                            int& outEdgeVisited)
	{
		MI_STATIC_ASSERT(has_mem_GetNeighborCount<graph_type>::value, "graph_type should has member GetNeighbor");
		MI_STATIC_ASSERT(has_mem_GetEdge<graph_type>::value, "graph_type should has member GetEdge");
		MI_STATIC_ASSERT(has_mem_edgeForm<graph_type>::value, "graph_type should has member edgeForm");

		if (gf == NULL) {
			throw NullPointerException("Please Build Graph first. (gf==NULL)");
		}

		int targetSize = 1;
		int resultSize = 0;
		outEdgeVisited = 0;
		ProbTransfom trans(gf->edgeForm);

		for (int it=0; it<num_iter; it++)
		{
            std::vector<bool> active(n, false);
			std::vector<int> RR;
			RR.push_back(target);
			active[target] = true;
			resultSize ++;

			int	h = 0;
			int t = targetSize;

			while (h<t) 
			{
				int k = gf->GetNeighborCount(RR[h]);
				for (int i=0; i<k; i++)
				{
					auto& e = gf->GetEdge(RR[h], i);
					// e = e(u, v), now RR[h] = u

					if (active[e.v]) continue;
					outEdgeVisited++;
					if (random.RandBernoulli(trans.Prob(e.w2)))
					{	
						RR.push_back(e.v);
						active[e.v] = true;
						t++;
						resultSize++;
					}
				}
				h++;
			}
			outRRSets.push_back(RR);
		}
		return (double)resultSize / (double)num_iter;
	}
};

typedef ReverseGCascadeT<Graph> ReverseGCascade;

#endif // ReverseGeneralCascade_h__

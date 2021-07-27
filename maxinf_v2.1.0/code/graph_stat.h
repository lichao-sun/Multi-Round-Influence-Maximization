#ifndef graph_stat_h__
#define graph_stat_h__

#include <iostream>
#include <sstream>
#include <vector>
#include <queue>
#include "graph.h"

/// Algorithm to generate basic statistics for graphs
template<class TGraph>
class GraphStatisticsT
{
public:
	typedef TGraph graph_type;

public:
	void Stats(TGraph& gf, std::ostream& out)
	{
		MI_STATIC_ASSERT(has_mem_n<graph_type>::value, "graph_type should has member n");
		MI_STATIC_ASSERT(has_mem_m<graph_type>::value, "graph_type should has member m");
		MI_STATIC_ASSERT(has_mem_degree<graph_type>::value, "graph_type should has member degree");
		MI_STATIC_ASSERT(has_mem_edges<graph_type>::value, "graph_type should has member edges");
		MI_STATIC_ASSERT(has_mem_index<graph_type>::value, "TGraph should has member index");
		MI_STATIC_ASSERT(has_mem_GetNeighborCount<graph_type>::value, "TGraph should has member GetNeighborCount");
		MI_STATIC_ASSERT(has_mem_GetEdge<graph_type>::value, "TGraph should has member GetEdge");
        
		out << "number of vertices:\t" << gf.n << std::endl;
		out << "number of edges:\t" << (gf.m / 2) << std::endl;
		out << "density:\t" << (double(gf.m) / gf.n / (gf.n - 1)) << std::endl;

		int maxdegree = 0;
		double tdegree = 0.0;
		int i, j, k;
		for (i = 0; i < gf.n; i++)
		{
			if (gf.GetDegree(i) > maxdegree) maxdegree = gf.GetDegree(i);
			tdegree += gf.GetDegree(i);
			//if (degree[i]%2) printf("%d\n", i);
		}
		out << "average degree:\t" << (tdegree / gf.n) << std::endl;
		out << "maximal degree:\t" << maxdegree << std::endl;
		int maxcmp = 0, ncmp = 0;

		std::vector<bool> used(gf.n, false);
		while (1)
		{
			std::queue<int> q;
			for (i = 0; i < gf.n; i++)
				if (!used[i]) break;
			if (i == gf.n) break;
			ncmp++;
			int cmpsize = 0;
			q.push(i);
			used[i] = true;
			while (!q.empty())
			{
				k = q.front();
				q.pop();
				cmpsize++;
				j = gf.GetNeighborCount(k);
				for (i = 0; i<j; i++)
				{
					auto& e = gf.GetEdge(k, i);
					if (used[e.v]) continue;
					q.push(e.v);
					used[e.v] = true;
				}
			}
			if (cmpsize>maxcmp) maxcmp = cmpsize;
		}
		out << "# of connected component:\t" << ncmp << std::endl;
		out << "largest component size:\t" << maxcmp << std::endl;
		out << "average component size:\t" << double(gf.n) / ncmp << std::endl;
	}
};

typedef GraphStatisticsT<Graph> GraphStatistics;

#endif // graph_stat_h__

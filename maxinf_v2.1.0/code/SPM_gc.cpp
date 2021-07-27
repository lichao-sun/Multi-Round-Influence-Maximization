#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include "graph.h"
#include "SPM_gc.h"
#include "cascade.h"


SPM_gc::SPM_gc()
{
	n = 0;
	m = 0;
	targetSize = 0;
	resultSize = 0;
	gf = NULL;
}


void SPM_gc::Build(Graph& gf)
{
	this->gf = &gf;
	this->n = gf.GetN();
	this->m = gf.GetM();
}

void SPM_gc::SetTarget(int size, int set[])
{
	targetSize = size;
	for (int i=0; i<size; i++)
		target[i] = set[i];
	resultSize = 0;
}

double SPM_gc::Run(int num_iter, int size, int set[])
{
	//printf("GreedyWC_SPM\n");
	if (gf == NULL) {
		throw NullPointerException("Please Build Graph first. (gf==NULL)");
	}

	SetTarget(size, set);
	int		h, t;
	int		*list=new int[n];
	int		*active=new int[n];
	double	*P=new double[n];
	ProbTransfom trans(gf->edgeForm);

	//for (int it=0; it<num_iter; it++)  //no simulation!
	{
		//printf("%d\n",Graph::GetN);
		memset(active, 0, sizeof(int)*n);
		for (int i=0; i<n; i++)
			P[i]=1.0;
		for (int i=0; i<targetSize; i++) 
		{
			list[i] = target[i];
			active[target[i]] = 1;
			P[target[i]] = 0;
		}
		//resultSize += targetSize;

		h = 0;
		t = targetSize;

		while (h<t) 
		{
			int k = gf->GetNeighborCount(list[h]);
			int ring = active[list[h]];
			P[list[h]]=1-P[list[h]];
			for (int i=0; i<k; i++)
			{
				auto e = gf->GetEdge(list[h], i);
				if (active[e.v]>0 && active[e.v]!=ring+1) continue;
				double puv = 1.0;
				puv = trans.Prob(e.w1);
				//double ratio=1.0/(double)Graph::GetDegree(e.v);
				//for (int j=0; j<e.c; j++)
				//	puv=puv*(1-ratio);
				//puv=1-puv;
					{
						P[e.v]=P[e.v]*(1-puv*P[e.u]);
						if (!active[e.v]) {
							list[t] = e.v;
							active[e.v] = ring+1;
							t++;
						}
					}
			}
			resultSize+=P[list[h]];
			h++;
		}
	}
	delete[] list;
	delete[] active;
	delete[] P;
	
	return (double)resultSize;
}

#include <cstdio>
#include <string>
#include <cstdlib>
#include <cmath>
//#include <time.h>
#include "SP1M_gc.h"
#include "graph.h"


SP1M_gc::SP1M_gc() 
{
	n = 0;
	m = 0;
	targetSize = 0;
	resultSize = 0;
	gf = NULL;
}

void SP1M_gc::Build(Graph& gf)
{
	this->gf = &gf;
	this->n = gf.GetN();
	this->m = gf.GetM();
}

void SP1M_gc::SetTarget(int size, int set[])
{
	targetSize = size;
	target.resize(targetSize);
	for (int i=0; i<size; i++)
		target[i] = set[i];
	resultSize = 0;
}

double SP1M_gc::Run(int num_iter, int size, int set[])
{
	//printf("GreedyIC_SPM\n");
	if (gf == NULL) {
		throw NullPointerException("Please Build Graph first. (gf==NULL)");
	}

	SetTarget(size, set);

	int		h, t;
	int		*list = new int[n];
	int		*active = new int[n];
	double	*P = new double[n];
	double  *P1 = new double[n];
	ProbTransfom trans(gf->edgeForm);

	//for (int it=0; it<num_iter; it++)
	{
		memset(active, 0, sizeof(int)*n);
		for (int i=0; i<n; i++){
			P[i]=1.0;
			P1[i]=1.0;
		}
		for (int i=0; i<targetSize; i++) 
		{
			list[i] = target[i];
			active[target[i]] = 1;
			P[target[i]] = 0;
		}
		//resultSize += targetSize;
		//printf("%d\n",n);

		h = 0;
		t = targetSize;

		int ring = 0;
		while (h<t) 
		{
			int u=list[h];
			if (active[u]>ring){
				ring=active[u];
				for (int h1=h;h1<t && active[list[h1]]==ring;h1++) {
					int k = gf->GetNeighborCount(list[h1]);
					P[list[h1]]=1-P[list[h1]];
					for (int i=0; i<k; i++)
					{
						auto& e = gf->GetEdge(list[h1], i);
						if (active[e.v]!=ring) continue;
						double puv = 1.0;
						double ratio=1.0 / (double) gf->GetDegree(e.v);
						for (int j=0; j<e.c; j++)
							puv=puv*(1-ratio);
							puv=1-puv;
							{
									P1[e.v]=P1[e.v]*(1-puv*P[e.u]);
							}
					}

				}

			}
			int k = gf->GetNeighborCount(u);
			P1[u]=(1-P[u])*(1-P1[u]);
			for (int i=0; i<k; i++)
			{
				auto& e = gf->GetEdge(u, i);
				if (active[e.v]>0 && active[e.v]!=ring+1) continue;
				double puv = trans.Prob(e.w1);
				//double ratio=1.0/(double)Graph::GetDegree(e.v);
				//for (int j=0; j<e.c; j++)
				//	puv=puv*(1-ratio);
				//puv=1-puv;
					{
						if (!active[e.v]) {
							list[t] = e.v;
							active[e.v] = ring+1;
							t++;
						}
						//if (active[e.v]==ring+1) {
							P[e.v]=P[e.v]*(1-puv*P[e.u]);
							P1[e.v]=P1[e.v]*(1-puv*P1[e.u]);
						//}
						//else if (active[e.v]==ring) {
						//	P1[e.v]=P1[e.v]*(1-puv*P[e.u]);
						//}
					}
			}
			resultSize+=P[u]+P1[u];
			h++;
		}
	}

	delete[] list;
	delete[] active;
	delete[] P;
	delete[] P1;
	
	return (double)resultSize;
}

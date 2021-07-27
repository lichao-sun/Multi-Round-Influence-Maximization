#include <vector>
#include "independ_cascade.h"

using namespace std;



void IndependCascade::Build(Graph& gf, double ratio)
{
	_Build(gf);
	this->ratio = ratio;
}

double IndependCascade::Run(int num_iter, int size, int set[])
{
	if (gf == NULL) {
		throw NullPointerException("Please Build Graph first. (gf==NULL)");
	}

	int targetSize = size;
	int resultSize = 0;

	int	h, t;
	int* list = new int[n];
	bool* active = new bool[n];

	for (int it = 0; it < num_iter; it++)
	{
		memset(active, 0, sizeof(bool)*n);
		for (int i = 0; i < targetSize; i++)
		{
			list[i] = set[i];
			active[list[i]] = true;
		}
		resultSize += targetSize;

		h = 0;
		t = targetSize;

		while (h < t)
		{
			int k = gf->GetNeighborCount(list[h]);
			//printf("%d \n",k);
			for (int i = 0; i < k; i++)
			{
				auto& e = gf->GetEdge(list[h], i);
				if (active[e.v]) continue;
				//printf("%d %d %g %g\n", e.u, e.v, exp(-e.w1), exp(-e.w2));
				//for (int j=0; j< e.c; j++)
				if (random.RandBernoulli(ratio))
				{
					list[t] = e.v;
					active[e.v] = true;
					t++;
					resultSize++;
					//break;
				}
			}
			h++;
		}
	}

	SAFE_DELETE_ARRAY(active);
	SAFE_DELETE_ARRAY(list);
	return (double)resultSize / (double)num_iter;
}


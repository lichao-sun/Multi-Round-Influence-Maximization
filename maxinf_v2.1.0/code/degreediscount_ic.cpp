#include <iostream>
#include <cmath>
#include <cstdio>
#include <string>
#include "degreediscount_ic.h"
#include "graph.h"
#include "common.h"

using namespace std;


DegreeDiscount_IC::DegreeDiscount_IC() {
	file = "degreediscount_ic.txt";
}

void DegreeDiscount_IC::Build(Graph& gf, int k, double ratio)
{
	n = gf.GetN();
	top = k;

	list.resize(k);
	d.resize(k);

	vector<double> weighted(n, 0.0);
	
	for (int i=0; i < n; i++)
		weighted[i] = gf.GetDegree(i);

	bool *used=new bool[n];
	memset(used, 0, sizeof(bool)*n);
	int *count=new int[n];
	memset(count, 0, sizeof(int)*n);
	
	
	for (int i=0; i < top; i++)
	{
		double max = -1000000.0;
		int mp = -1;
		for (int j=0; j<n; j++)
			if (!used[j])
			{
				double tmp = weighted[j] - 2 * count[j] - ratio*count[j] * (weighted[j] - count[j]);
				if (tmp >max)
				{
					max = tmp;
					mp = j;
				}
			}

		list[i] = mp;
		d[i] = max;
		used[mp] = true;

		for (int j=0; j < gf.GetNeighborCount(mp); j++)
		{
			auto& e = gf.GetEdge(mp, j);
			count[e.v]+=e.c;
		}
	}


	WriteToFile(file, gf);

	SAFE_DELETE_ARRAY(used);
	SAFE_DELETE_ARRAY(count);
}

void DegreeDiscount_IC::BuildFromFile(Graph& gf)
{
	ReadFromFile(file, gf);
}


#include <cstdio>
#include <cmath>
#include "weighted_degree.h"
#include "graph.h"
#include "common.h"

using namespace std;

WeightedDegree::WeightedDegree()
{
	file = "weighteddegree.txt";
}

void WeightedDegree::qsort_degree(int h, int t)
{
	if (h<t) 
	{
		int i = h, j = t;
		double midd = d[(i+j)/2];
		int midn = list[(i+j)/2];
		d[(i+j)/2] = d[i];
		list[(i+j)/2] = list[i];

		while (i<j) 
		{
			while ((i<j) && (d[j]<midd))
				j--;
			if (i<j) {
				d[i] = d[j];
				list[i] = list[j];
				i++;
			}
			while ((i<j) && (d[i]>midd))
				i++;
			if (i<j) {
				d[j] = d[i];
				list[j] = list[i];
				j--;
			}
		}

		d[i] = midd;
		list[i] = midn;
		qsort_degree(h, i-1);
		qsort_degree(i+1, t);
	}
}

void WeightedDegree::Build(Graph& gf, int k)
{
	n = gf.GetN();
	top = k;
	
	std::vector<double> weighted(n, 0.0);
	d.resize(top);
	list.resize(top);
	ProbTransfom trans(gf.edgeForm);

	for (int i=0; i<n; i++)
	{
		int degree = gf.GetNeighborCount(i);
		for (int j=0; j < degree; j++)
			weighted[i] += trans.Prob(gf.GetEdge(i, j).w1);
	}

	for (int i=0; i < top; i++)
	{	
		int max = 0;
		for (int j=0; j<n; j++)
			if (weighted[j]>weighted[max])
				max = j;
		//printf("%d %g\n",max,d[max]);
		
		list[i] = max;
		d[i] = weighted[max];
		weighted[max] = 0.0;
	}

	WriteToFile(file, gf);
}

void WeightedDegree::BuildFromFile(Graph& gf)
{
	ReadFromFile(file, gf);
}

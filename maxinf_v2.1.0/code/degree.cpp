#include <iostream>
#include <cmath>
#include <cstdio>
#include "degree.h"
#include "graph.h"
#include "common.h"

using namespace std;

Degree::Degree() {
	file = "degree.txt";
}

void Degree::qsort_degree(int h, int t)
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

void Degree::Build(Graph& gf, int k)
{
	n = gf.GetN();
	top = k;
	
	// d.resize(n, 0);
	list.resize(k, 0);
	d.resize(k, 0);

	vector<int> degrees(n, 0);

	for (int i=0; i<n; i++)
		degrees[i] = gf.GetDegree(i);

	for (int i=0; i < top; i++)
	{
		int max = 0;
		for (int j=0; j < n; j++)
			if (degrees[j] > degrees[max])
				max = j;

		list[i] = max;
		d[i] = degrees[max];

		degrees[max] = 0;
	}


	WriteToFile(file, gf);
}

void Degree::BuildFromFile(Graph& gf)
{
	ReadFromFile(file, gf);
}



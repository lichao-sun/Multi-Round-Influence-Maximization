#include <cstdio>
#include <string>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <fstream>
#include "pagerank.h"
#include "graph.h"
#include "event_timer.h"
#include "common.h"


using namespace std;

Pagerank::Pagerank()
{
	file = "pagerank.txt";
}



int Pagerank::GetMax(int round, std::vector<double>& values)
{
	double max = -1000000.0;
	int mp = -1;
	for (int j=0; j<n; j++)
	{
		double tmp = values[j];
		if (tmp >max)
		{
			max = tmp;
			mp = j;
		}
	}
	return mp;
}


double Pagerank::Build(Graph& gf, int num, double dampen)
{
	EventTimer timer;
    timer.SetTimeEvent("start");
    
	n = gf.GetN();
	top = num;
	d.resize(top);
	list.resize(top);

	dp.resize(n);
	dd.resize(n);

	vector<double> newdp(n,0);
	vector<double> self(n,0.0);

	vector<int> set(top);
	ProbTransfom trans(gf.edgeForm);

	int i=0;
	for (i=0; i<n; i++)
	{
		dp[i] = 1.0;
		dd[i] = gf.GetNeighborCount(i);
		for (int j=0;j<dd[i];j++)
			self[i] += trans.Prob(gf.GetEdge(i, j).w2);
	}
	

	double theta=1e-3;				
	double delta;
	long count=0;
	do
	{
		for (i=0;i<n;i++)
			newdp[i]=0;
		for (i=0;i<n;i++)
			for (int j=0;j<dd[i];j++)
			{
				auto& e = gf.GetEdge(i,j);
				newdp[e.v] += dp[i] * trans.Prob(e.w2) / self[i];
			}
		delta=0;
		for (i=0;i<n;i++)
		{
			delta+=fabs(dampen+(1-dampen)*newdp[i]-dp[i]);
			dp[i]=dampen+(1-dampen)*newdp[i];
		}
		//printf("%g\n",delta);
		
	} while (delta>theta &&++count<10000);
	
	vector<double> tmpDP = dp;

	i=0;
	double maxV = -1000000.0;
	int mp = 0;
	
	for (i=0; i<top; i++)
	{
		maxV = -1000000.0;
		int x = GetMax(i, tmpDP);
		set[i] = x;
		double improve = tmpDP[x];
		if (improve > maxV) {
			maxV = improve;
			mp = x;
		}
		tmpDP[mp] = 0;
		set[i] = mp;
		list[i] = mp;
		d[i] = maxV;
	}
	

    timer.SetTimeEvent("end");
    
	WriteToFile(file, gf);

	// QueryPerformanceCounter(&Eend);
	// double timer = (double)(Eend.QuadPart - start.QuadPart) / freq.QuadPart;
	FILE *timetmpfile;
	fopen_s(&timetmpfile, "time_pagerank.txt", "w");
	fprintf(timetmpfile,"%g\n", timer.TimeSpan("start", "end"));
	fclose(timetmpfile);

	return 1;
}


void Pagerank::BuildFromFile(Graph& gf)
{
	ReadFromFile(file, gf);
	dp.resize(n);
	dd.resize(n);
}

///////////////////////////////
double PagerankRegen::NAdefault(double v, double defaultVal) {
	if (v != v) {
		return defaultVal;
	}
	return v;
}
double PagerankRegen::roundInInterval(double v, double left, double right) {
	if (v > right) {
		return right;
	}
	if (v < left) {
		return left;
	}
	return v;
}
bool PagerankRegen::almostZero(double v, double threshold) {
	return abs(v) < threshold;
}

double PagerankRegen::Build(Graph& gf, std::string& graph_file, int num, double dampen /*= 0.15 */)
{
	double ret = Pagerank::Build(gf, num, dampen);

	printf("Begin Regenerate Graph\n");

	static double MIN_EPS = 1e-32;
	static double MIN_CLOSE = 1e-10;
	ProbTransfom trans(gf.edgeForm);

	ofstream fout(graph_file.c_str());

	int n = gf.GetN();
	int m = gf.GetM(); // m is the 2 times the number of undirected edges

	double ratio = double(n) / double(m);
	fout << n << " " << (m/2) << endl;  // in the file format, we use the count of undirected edges
	for (int u = 0; u < n; u++)
	{
		int neighborCount = gf.GetNeighborCount(u);
		for (int i = 0; i < neighborCount; i++)
		{
			Edge& e = gf.GetEdge(u, i);
			if (e.u >= e.v)
				continue;
			double p1 = almostZero(trans.Prob(e.w1), MIN_CLOSE) ? 0 : (dp[e.u] / max(dp[e.u] + dp[e.v], MIN_EPS)) * ratio;
			double p2 = almostZero(trans.Prob(e.w2), MIN_CLOSE) ? 0 : (dp[e.v] / max(dp[e.u] + dp[e.v], MIN_EPS)) * ratio;
			p1 = roundInInterval(NAdefault(p1, 0.0), 0.0, 1.0);
			p2 = roundInInterval(NAdefault(p2, 0.0), 0.0, 1.0);
			std::string& us = gf.MapIndexToNodeName(e.u);
			std::string& vs = gf.MapIndexToNodeName(e.v);
			fout << us << " " << vs << " " << p1 << " " << p2 << endl;
			fout << vs << " " << us << " " << p2 << " " << p1 << endl;
		}
	}

	return ret;
}


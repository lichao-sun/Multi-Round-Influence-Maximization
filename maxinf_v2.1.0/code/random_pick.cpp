#include <iostream>
#include <algorithm>
#include <cmath>
#include <random>
#include "random_pick.h"
#include "graph.h"
#include "common.h"

using namespace std;

RandomPick::RandomPick()
{
	n = 0;
	top = 0;
	file = "random_pick.txt";
}

void RandomPick::Build(Graph& gf, int k)
{
	n = gf.GetN();

	vector<int> nodes;
	int count = gf.GetRealNodeCount();
	nodes.resize(count);

	top = min(k, count);
	list.resize(top);
	d.resize(top);

	for (int i = 0; i < count; i++)
		nodes[i] = i;

	std::random_device rd;
	std::mt19937 g(rd());

	std::shuffle(nodes.begin(), nodes.end(), g);
	
	
	for (int i = 0; i < top; i++) {
		list[i] = nodes[i];
		d[i] = 0;
	}

	WriteToFile(file, gf);
}

void RandomPick::BuildFromFile(Graph& gf)
{
	ReadFromFile(file, gf);
}

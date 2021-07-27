#include <cstdio>
#include <string>
#include <cstdlib>
#include <vector>
#include <set>
#include "common.h"
#include "greedy_online.h"
#include "graph.h"


//#include <set>
using namespace std;

GreedyOnline::GreedyOnline()
{
	file = "greedy_online.txt";
	onlinebound_file = "greedy_online_bound_D.txt";
	isOnlineEstimated = true;
}

bool GreedyOnline::SortBestToTop(EstHeap& heap, int round, int* set, double old,
							Graph& gf, ICascade& cascade, int& ccc)
{
	if (heap.empty()) 
		return false;

	// until find next best
	while (heap.top()->last_update != round) 
	{
		EstNode* cur_top = heap.top();
		heap.pop();

		// to estimate
		set[round] = cur_top->id;
		cur_top->last_update = round;
		cur_top->d = cascade.Run(NUM_ITER, round + 1, set) - old;

		// re-insert
		heap.push(cur_top);
		ccc++;
	}
	return true;
}

void GreedyOnline::Build(Graph& gf, int num, ICascade& cascade)
{
	if (isOnlineEstimated) {
		onlineSeeds.clear();
		onlineD.clear();
	}
	n = gf.GetN();
	top = num;
	d.resize(top); 
	list.resize(top);

	bool *used= new bool[n];
	memset(used, 0, sizeof(bool)*n);
	
	// Change: Tian 2014/10/1
	// int set[SET_SIZE];
	int* set = new int[2*top];
	double old = 0.0;
	vector<double> old_history(SET_SIZE);

	EstHeap heap;
	vector<EstNode*> allNodes;
	for (int k = 0; k < n; k++) {
		//initialize: d = largest, last_update = -1
		EstNode* node = new EstNode(k, (double)(n+1), -1);
		allNodes.push_back(node);
		heap.push(node);
	}

	vector<int> firstSeeds;
	vector<double> firstD;
	vector<EstNode*> greedy_sequence;
	for (int i = 0; i < top; i++) {
		printf("Seed: %d", i);
		int round = i;
		int ccc = 0;

		// find the largest marginal
		bool succeed = false;
		
		// run on first time: estimate 50 nodes
		if (isOnlineEstimated && (round == 0)) {
			int estsize = (n >= top) ? top : n;
			vector<EstNode*> next_sequence;
			// estimate nexts
			for (int h = 0; h < estsize; ++h) {
				succeed = SortBestToTop(heap, 0, set, 0, gf, cascade, ccc);
				if (!succeed)
					break;

				// find the cur_top
				EstNode* cur_top = heap.top();
				next_sequence.push_back(cur_top);
				heap.pop();
			}

			// re-insert
			for (size_t h = 0; h < next_sequence.size(); ++h) {
				EstNode* node = next_sequence[h];
				firstSeeds.push_back(node->id);
				firstD.push_back(node->d);
				heap.push(node);
			}

			// write online bound: first time
			FILE *out_bound;
			fopen_s(&out_bound, "greedy_online_bound_D_first.txt", "w");
			double sumD = 0;
			for (size_t h = 0; h < firstD.size(); ++h) {
				sumD += firstD[h];
			}
			fprintf(out_bound, "Round[%d] Size[%d] SumD[%g] Total[%g]\n", 0, firstSeeds.size(), sumD, sumD);	//the nodes we want!
			for (size_t h = 0; h < firstSeeds.size(); ++h) {
				fprintf(out_bound, "\t%d\t%g\n", firstSeeds[h], firstD[h]);	//the nodes we want!	
			}
			fclose(out_bound);
		}

		succeed = SortBestToTop(heap, round, set, old, gf, cascade, ccc);
		if (!succeed) break;

		// already find currrent best
		EstNode* top = heap.top();
		set[round] = top->id;
		list[round] = top->id;
		d[round] = top->d;
		old += top->d;
		old_history[round] = old;
		// pop from the heap
		greedy_sequence.push_back(top);
		heap.pop();

		printf("\tcounters: search %d", ccc);
		ccc = 0;

		if (isOnlineEstimated) {
			// estimate the online bound
			int next_round = round + 1;
			int tmpsize = n-round-1;
			int estsize = (tmpsize >= round+1) ? round+1 : tmpsize;
			vector<EstNode*> next_sequence;
			// estimate nexts
			for (int h = 0; h < estsize; ++h) {
				succeed = SortBestToTop(heap, next_round, set, old, gf, cascade, ccc);
				if (!succeed)
					break;

				// find the cur_top
				EstNode* cur_top = heap.top();
				next_sequence.push_back(cur_top);
				heap.pop();
			}

			// re-insert
			vector<int> curSeeds;
			vector<double> curD;
			for (size_t h = 0; h < next_sequence.size(); ++h) {
				EstNode* node = next_sequence[h];
				curSeeds.push_back(node->id);
				curD.push_back(node->d);
				heap.push(node);
			}
			onlineSeeds.push_back(curSeeds);
			onlineD.push_back(curD);
		}
		printf("\testimate %d\n", ccc);
	}


	WriteToFile(file, gf);
	
	FILE* out_bound;
	// write online bound: one-by-one
	fopen_s(&out_bound, onlinebound_file.c_str(), "w");
	for (int i=0; i<top; i++) {
		vector<int>& tmpseeds = onlineSeeds[i];
		vector<double>& tmpD = onlineD[i];
		double sumD = 0;
		for (size_t h = 0; h < tmpseeds.size(); ++h) {
			sumD += tmpD[h];
		}
		fprintf(out_bound, "Round[%d] Size[%d] SumD[%g] Total[%g]\n", i+1, tmpseeds.size(), sumD, old_history[i]+sumD);	//the nodes we want!
		for (size_t h = 0; h < tmpseeds.size(); ++h) {
			fprintf(out_bound, "\t%d\t%g\n", tmpseeds[h], tmpD[h]);	//the nodes we want!	
		}
	}
	fclose(out_bound);

	for (int k = 0; k < n; k++) {
		delete allNodes[k];
	}

	SAFE_DELETE_ARRAY(used);
	SAFE_DELETE_ARRAY(set);
}

void GreedyOnline::BuildRanking(Graph& gf, int num, ICascade& cascade)
{
	n = gf.GetN();
	top = num;
	d.resize(top); 
	list.resize(top);

	//int* set = new int[top]; //alert 100
	//double* inf = new double[top];
	for(int i=0;i<top;i++)
	{
		list[i] =0;
		d[i] =0;
	}
	vector<double> improve(n, 0.0);
	int tmp[1];
	for (int i=0; i<n; i++)
	{
		tmp[0] = i;
		improve[i] = cascade.Run(NUM_ITER, 1, tmp);
		if(improve[i] > list[top-1])
		{
			d[top-1] = improve[i];
			list[top - 1] = i;
			int j = top - 2;
			while(j>=0)
			{
				if(improve[i] > d[j])
				{
					int int_tmp = list[j];
					double inf_tmp = d[j];
					list[j] = i;
					d[j] = improve[i];
					list[j + 1] = int_tmp;
					d[j+1] = inf_tmp;
				}
				else
					break;
				j--;
			}
		}
	}

	WriteToFile(file, gf);
}


void GreedyOnline::BuildFromFile(Graph& gf, const char* name)
{
	ReadFromFile(name, gf);
}

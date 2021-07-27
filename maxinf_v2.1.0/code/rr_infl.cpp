
#include <iostream>
#include <vector>
#include <set>
#include <random>
#include <functional>
#include <algorithm>
#include <cassert>

#include "rr_infl.h"
#include "reverse_general_cascade.h"
#include "graph.h"
#include "event_timer.h"
#include "common.h"

using namespace std;




// define a comparator for counts
struct CountComparator
{
public:
	vector<int>& counts;
	CountComparator(vector<int>& c) :counts(c){}
public:
	bool operator () (int a, int b) {
		return (counts[a] <= counts[b]);
	}
};



void RRInflBase::InitializeConcurrent()
{
	if (isConcurrent)
	{
#ifdef MI_USE_OMP
		/////////////////////////////////////////////
		// run concurrently
		const double DYNAMIC_RATIO = 0.25;
		int maxThreads = omp_get_max_threads();
		omp_set_num_threads(maxThreads);
		int dynamicThreads = (int)(maxThreads * DYNAMIC_RATIO);
		omp_set_dynamic(dynamicThreads);

		cout << "== Turn on omp optimization: == " << endl;
		cout << "#Max Threads = " << omp_get_max_threads() << "\t#Dynamic Threads = " << omp_get_dynamic() << endl;
#else
		cout << "== omp is not supported or enabled == " << endl;
#endif
	}
}


void RRInflBase::_AddRRSimulation(size_t num_iter,
	cascade_type& cascade,
	std::vector< RRVec >& refTable,
	std::vector<int>& refTargets)
{
	vector<int> edgeVisited; // discard
	_AddRRSimulation(num_iter, cascade, refTable, refTargets, edgeVisited);
}

void RRInflBase::_AddRRSimulation(size_t num_iter,
	cascade_type& cascade,
	std::vector< RRVec >& refTable,
	std::vector<int>& refTargets,
	std::vector<int>& refEdgeVisited)
{
#ifdef MI_USE_OMP
	if (!isConcurrent) {
#endif
		// run single thread

		for (size_t iter = 0; iter < num_iter; ++iter) {
			int id = cascade.GenRandomNode();
			int edgeVisited;
			cascade.ReversePropagate(1, id, refTable, edgeVisited);

			refTargets.push_back(id);
			refEdgeVisited.push_back(edgeVisited);
		}

#ifdef MI_USE_OMP
	}
	else {
		// run concurrently
		// #pragma omp parallel for private(iter, id, edge_visited_count, tmpTable) shared(refTable, refTargets, refEdgeVisited)

#pragma omp parallel for ordered
		for (int iter = 0; iter < num_iter; ++iter) {
			int id = cascade.GenRandomNode();
			int edgeVisited;
			vector<RRVec> tmpTable;
			cascade.ReversePropagate(1, id, tmpTable, edgeVisited);
#pragma omp critical
			{
				refTable.push_back(tmpTable[0]);
				refTargets.push_back(id);
				refEdgeVisited.push_back(edgeVisited);
			}
		}
	}
#endif

}


////////////////////////////////////////////////////////////////////
// CR_IMM RRInflBase: KDD 2018
void RRInflBase::_CRIMAddRRSimulation(size_t num_iter,
	cascade_type& cascade,
	std::vector<std::vector< RRVec >>& refTable,
	std::vector<int>& refTargets)
{
	vector<int> edgeVisited; // discard
	_CRIMAddRRSimulation(num_iter, cascade, refTable, refTargets, edgeVisited);
}

void RRInflBase::_CRIMAddRRSimulation(size_t num_iter,
	cascade_type& cascade,
	std::vector<std::vector< RRVec >>& refTable,
	std::vector<int>& refTargets,
	std::vector<int>& refEdgeVisited)
{
#ifdef MI_USE_OMP
	if (!isConcurrent) {
#endif
		// run single thread
		std::vector< RRVec > crTable;
		std::vector<int> RR;
		
		for (size_t iter = 0; iter < num_iter; ++iter) {
			crTable.clear();
			RR.clear();

			int id = cascade.GenRandomNode();
			int edgeVisited;

			cascade.ReversePropagate(T, id, crTable, edgeVisited);
			
			refTable.push_back(crTable);
			/*for (int i = 0; i < 30000; i++){
				act_nodes_arr[i] = 1;
			}*/
			/*memset(act_nodes_arr, 1, sizeof(bool)* 30000);
			for (int i = 0; i < crTable.size(); i++){
				const RRVec& RR_cr = crTable[i];
				for (int source : RR_cr) {
					if (act_nodes_arr[source]){
						RR.push_back(source);
						act_nodes_arr[source] = 0;
					}
				}
			}
			refTable.push_back(RR);*/
			//}
			//for (size_t i = 0; i < table.size(); ++i) {
			//	
			//	
			//		degrees[source]++;
			//		degreeRRIndices[source].push_back(i); // add index of table
			//	}
			//}

			refTargets.push_back(id);
			refEdgeVisited.push_back(edgeVisited);
			
		}

#ifdef MI_USE_OMP
	}
	else {
		// run concurrently
		// #pragma omp parallel for private(iter, id, edge_visited_count, tmpTable) shared(refTable, refTargets, refEdgeVisited)

#pragma omp parallel for ordered
		for (int iter = 0; iter < num_iter; ++iter) {
			int id = cascade.GenRandomNode();
			int edgeVisited;
			vector<RRVec> tmpTable;
			cascade.ReversePropagate(1, id, tmpTable, edgeVisited);
#pragma omp critical
			{
				refTable.push_back(tmpTable);
				refTargets.push_back(id);
				refEdgeVisited.push_back(edgeVisited);
			}
		}
	}
#endif

}

void RRInflBase::_CRIMRebuildRRIndices()
{
	degreesVec.clear();
	degreeRRIndicesVec.clear();
	sourceSetVec.clear();
	cout << tableVec.size() << endl;
	for (int j = 0; j < T; j++){
		degrees.clear();
		degrees.resize(n, 0);
		degreeRRIndices.clear();
		for (int i = 0; i < n; ++i) {
			degreeRRIndices.push_back(vector<int>());
		}
		// to count hyper edges:
		for (size_t i = 0; i < tableVec.size(); ++i) {
			const RRVec& RR = tableVec[i][j];
			for (int source : RR) {
				degrees[source]++;
				degreeRRIndices[source].push_back(i); // add index of table
			}
		}
		// add to sourceSet where node's degree > 0
		sourceSet.clear();
		for (size_t i = 0; i < degrees.size(); ++i) {
			if (degrees[i] > 0) {
				sourceSet.insert(i);
			}
		}
		//cout << degrees[559] << endl;
		degreesVec.push_back(degrees);
		degreeRRIndicesVec.push_back(degreeRRIndices);
		sourceSetVec.push_back(sourceSet);
	}
}


double RRInflBase::_CRIMRunGreedy(int seed_size,
	vector<vector<int>>& outSeedsVec,
	vector<double>& outEstSpread)
{
	//cout << degreesVec.size() << ' ' << degreeRRIndicesVec.size() << ' ' << sourceSetVec.size() << endl;
	outSeedsVec.clear();
	outEstSpread.clear();

	// set enables for table
	//vector<bool> enables;
	enables.clear();
	enables.resize(tableVec.size(), true);

	seeds_list.clear();
	for (int i = 0; i < T; ++i) {
		seeds_list.push_back(vector<int>());
	}

	vector<set<int>> candidatesVec(sourceSetVec);
	/*for (int i = 0; i < sourceSetVec.size(); i++){
		cout << "sourceS:" << candidatesVec[i].size() << endl;
	}
	for (int i = 0; i < candidatesVec.size(); i++){
		cout << "sourceC:"<< candidatesVec[i].size() << endl;
	}*/

	double spread = 0;
	vector<int> resRounds;

	for (int iter = 0; iter < seed_size * T; ++iter) {
		int max = 0;
		resRounds.clear();

		for (int j = 0; j < T; j++){
			if (seeds_list[j].size() < seed_size){
				resRounds.push_back(j);
			}
		}

		int seeds_count = 0;
		for (int i = 0; i < seeds_list.size(); i++){
			seeds_count += seeds_list[i].size();
		}
		//cout << "seeds:" << seeds_count << endl;
		int largestRound = -1;
		//cout << "resRound:" << resRounds.size() << endl;
		for (int j = 0; j < resRounds.size(); j++){
			CountComparator comp(degreesVec[j]);
			set<int>::const_iterator tmpMaxPt = max_element(candidatesVec[resRounds[j]].begin(), candidatesVec[resRounds[j]].end(), comp);
			int tmpMaxSource = *tmpMaxPt;
			//cout << resRounds[j] << ' ' << tmpMaxSource << endl;
			if (degreesVec[resRounds[j]][tmpMaxSource] > max){
				max = degreesVec[resRounds[j]][tmpMaxSource];
				largestRound = resRounds[j];
			}
		}
		
		CountComparator comp(degreesVec[largestRound]);
		set<int>::const_iterator maxPt = max_element(candidatesVec[resRounds[largestRound]].begin(), candidatesVec[resRounds[largestRound]].end(), comp);
		int maxSource = *maxPt;
		//cout << iter << ' ' << largestRound << ' ' << maxSource << endl;
		//cRound[largestRound]++;
		// selected one node
		//cout << maxSource << endl;
		seeds_list[largestRound].push_back(maxSource);
		//cout << iter << ' ' << outSeeds.size() << endl;

		//cout << "max:" << largestRound << ' ' << degreesVec[largestRound][maxSource] << ' ' << maxSource << endl;
		assert(degreesVec[largestRound][maxSource] >= 0);
		// estimate spread
		spread = spread + ((double)n * degreesVec[largestRound][maxSource] / tableVec.size());
		outEstSpread.push_back(spread);
		
		/*if (maxSource == 559){
			cout << degreesVec[largestRound][559] << endl;
		}*/
		//cout << "spread:" << ((double)n * degreesVec[largestRound][maxSource] / tableVec.size()) << endl;

		// clear values
		candidatesVec[largestRound].erase(maxPt);
		degreesVec[largestRound][maxSource] = -1;

		// deduct the counts from the rest nodes
		//cout << enables.size() << endl;

		
		if (!isConcurrent) {
			const vector<int>& idxList = degreeRRIndicesVec[largestRound][maxSource];
			//cout << idxList.size() << endl;
			for (int idx : idxList) {
				//if (enables[idx]) {
				//	const RRVec& RRset = tableVec[idx][largestRound];
				//	for (int rr : RRset) {
				//		if (rr == maxSource) continue;
				//		degreesVec[largestRound][rr]--; // deduct
				//		
				//	}
				//	enables[idx] = false;
				//}
				if (enables[idx]) {
					for (int j = 0; j < T; j++){
						const RRVec& RRset = tableVec[idx][j];
						for (int rr : RRset) {
							if (rr == maxSource && j==largestRound) continue;
							degreesVec[j][rr]--; // deduct
						}
					}
					enables[idx] = false;
				}
			}
			
			int cout_a = 0;
			int cout_b = 0;
			for (int i = 0; i < T; i++){
				if (i != largestRound){
					const vector<int>& idxList = degreeRRIndicesVec[i][maxSource];
					for (int idx : idxList) {
						/*if (i == 3 && maxSource == 559 && targets[idx] == maxSource)
							cout << iter << ' ' << largestRound << ' ' << maxSource << ' ' << targets[idx] << ' ' << enables[idx] << endl;*/
						if (!enables[idx]) {
							if (targets[idx] == maxSource){
								cout_b += 1;
								const RRVec& RRset = tableVec[idx][i];
								//cout << i << ' ' << targets[idx] << ' ' << maxSource << ' ' << RRset.size() << endl;
								for (int rr : RRset) {
									if (rr == maxSource) continue;
									degreesVec[i][rr]--; // deduct
								}
							}
							enables[idx] = false;
						}
					}

				}
				/*if (cout_b > 0){
					cout << cout_a << ' ' << cout_b << endl;
				}*/
			}
			
		}
		else {
			const vector<int>& idxList = degreeRRIndices[maxSource];
			// run concurrently
#pragma omp parallel for
			for (int idxIter = 0; idxIter < idxList.size(); ++idxIter) {
				int idx = idxList[idxIter];
				if (enables[idx]) {
					const RRVec& RRset = table[idx];
					for (int rr : RRset) {
						if (rr == maxSource) continue;

#pragma omp atomic
						degrees[rr]--; // deduct
					}
					enables[idx] = false;
				}
			}
		}
	}
	//cout << 11 << endl;
	int seeds_count = 0;
	for (int i = 0; i < seeds_list.size(); i++){
		seeds_count += seeds_list[i].size();
	}
	assert(seeds_count == seed_size * T);
	assert(outEstSpread.size() == seed_size * T);

	return spread;
}

void RRInflBase::_RebuildRRIndices()
{
	
	degrees.clear();
	degrees.resize(n, 0);
	degreeRRIndices.clear();
	for (int i = 0; i < n; ++i) {
		degreeRRIndices.push_back(vector<int>());
	}
	// to count hyper edges:
	for (size_t i = 0; i < table.size(); ++i) {
		const RRVec& RR = table[i];
		for (int source : RR) {
			degrees[source]++;
			degreeRRIndices[source].push_back(i); // add index of table
		}
	}
	// add to sourceSet where node's degree > 0
	sourceSet.clear();
	for (size_t i = 0; i < degrees.size(); ++i) {
		if (degrees[i] > 0) {
			sourceSet.insert(i);
		}
	}
}


// Apply Greedy to solve Max Cover
double RRInflBase::_RunGreedy(int seed_size,
	vector<int>& outSeeds,
	vector<double>& outEstSpread)
{
	outSeeds.clear();
	outEstSpread.clear();

	// set enables for table
	//vector<bool> enables;
	enables.clear();
	enables.resize(table.size(), true);

	set<int> candidates(sourceSet);
	CountComparator comp(degrees);

	double spread = 0;
	for (int iter = 0; iter < seed_size; ++iter) {
		set<int>::const_iterator maxPt = max_element(candidates.begin(), candidates.end(), comp);
		int maxSource = *maxPt;
		// cout << "" << maxSource << "\t";
		assert(degrees[maxSource] >= 0);

		// selected one node
		outSeeds.push_back(maxSource);

		// estimate spread
		spread = spread + ((double)n * degrees[maxSource] / table.size());
		outEstSpread.push_back(spread);

		// clear values
		candidates.erase(maxPt);
		degrees[maxSource] = -1;

		// deduct the counts from the rest nodes
		const vector<int>& idxList = degreeRRIndices[maxSource];
		if (!isConcurrent) {
			for (int idx : idxList) {
				if (enables[idx]) {
					const RRVec& RRset = table[idx];
					for (int rr : RRset) {
						if (rr == maxSource) continue;
						degrees[rr]--; // deduct
					}
					enables[idx] = false;
				}
			}
		}
		else {
			// run concurrently
#pragma omp parallel for
			for (int idxIter = 0; idxIter < idxList.size(); ++idxIter) {
				int idx = idxList[idxIter];
				if (enables[idx]) {
					const RRVec& RRset = table[idx];
					for (int rr : RRset) {
						if (rr == maxSource) continue;

#pragma omp atomic
						degrees[rr]--; // deduct
					}
					enables[idx] = false;
				}
			}
		}
	}

	assert(outSeeds.size() == seed_size);
	assert(outEstSpread.size() == seed_size);

	return spread;
}

void RRInflBase::_SetResults(const vector<int>& seeds, const vector<double>& cumu_spread)
{
	for (int i = 0; i < top; ++i) {
		list[i] = seeds[i];
		// d[i] = cumu_spread[i];  // use cumulative spread
		d[i] = (i > 0) ? (cumu_spread[i] - cumu_spread[i - 1]) : cumu_spread[i]; // use marginal spread
	}
}


////////////////////////////////////////////////////////////////////
// WR_IMM RRInflBase: KDD 2018
void RRInflBase::_WRIMAddRRSimulation(size_t num_iter,
	cascade_type& cascade,
	std::vector< RRVec >& refTable,
	std::vector<int>& refTargets,
	int* root_arr,
	int round)
{
	vector<int> edgeVisited; // discard
	_WRIMAddRRSimulation(num_iter, cascade, refTable, refTargets, edgeVisited, root_arr, round);
}

void RRInflBase::_WRIMAddRRSimulation(size_t num_iter,
	cascade_type& cascade,
	std::vector< RRVec >& refTable,
	std::vector<int>& refTargets,
	std::vector<int>& refEdgeVisited,
	int* root_arr,
	int round)
{
#ifdef MI_USE_OMP
	if (!isConcurrent) {
#endif
		// run single thread
		cout << "root size and R size:" << R << endl;
		if (round == 0){
			for (size_t iter = 0; iter < num_iter; ++iter) {
				int id = cascade.GenRandomNode();

				int edgeVisited;
				cascade.ReversePropagate(1, id, refTable, edgeVisited);

				refTargets.push_back(id);
				refEdgeVisited.push_back(edgeVisited);
			}
		}
		else{
			for (size_t iter = 0; iter < num_iter; ++iter) {

				int id = root_arr[rand() % R];

				int edgeVisited;
				cascade.ReversePropagate(1, id, refTable, edgeVisited);

				refTargets.push_back(id);
				refEdgeVisited.push_back(edgeVisited);
			}
			//cout << "size" << root_arr.size() << ' ' << targets.size() << endl;
		}

#ifdef MI_USE_OMP
	}
	else {
		// run concurrently
		// #pragma omp parallel for private(iter, id, edge_visited_count, tmpTable) shared(refTable, refTargets, refEdgeVisited)

#pragma omp parallel for ordered
		for (int iter = 0; iter < num_iter; ++iter) {
			int id = cascade.GenRandomNode();
			int edgeVisited;
			vector<RRVec> tmpTable;
			cascade.ReversePropagate(1, id, tmpTable, edgeVisited);
#pragma omp critical
			{
				refTable.push_back(tmpTable[0]);
				refTargets.push_back(id);
				refEdgeVisited.push_back(edgeVisited);
			}
		}
	}
#endif

}


////////////////////////////////////////////////////////////////////
// AdaRRInflBase: KDD 2018
void RRInflBase::_AdaAddRRSimulation(size_t num_iter,
	cascade_type& cascade,
	std::vector< RRVec >& refTable,
	std::vector<int>& refTargets)
{
	vector<int> edgeVisited; // discard
	_AdaAddRRSimulation(num_iter, cascade, refTable, refTargets, edgeVisited);
}

void RRInflBase::_AdaAddRRSimulation(size_t num_iter,
	cascade_type& cascade,
	std::vector< RRVec >& refTable,
	std::vector<int>& refTargets,
	std::vector<int>& refEdgeVisited)
{
#ifdef MI_USE_OMP
	if (!isConcurrent) {
#endif
		// run single thread

		for (size_t iter = 0; iter < num_iter; ++iter) {
			int id = cascade.GenRandomNode();
			int edgeVisited;

			if (act_nodes_arr[id] == 1){
				iter = iter - 1;
				continue;
			}
			cascade.ReversePropagate(1, id, refTable, edgeVisited);

			refTargets.push_back(id);
			refEdgeVisited.push_back(edgeVisited);
		}

#ifdef MI_USE_OMP
	}
	else {
		// run concurrently
		// #pragma omp parallel for private(iter, id, edge_visited_count, tmpTable) shared(refTable, refTargets, refEdgeVisited)

#pragma omp parallel for ordered
		for (int iter = 0; iter < num_iter; ++iter) {
			int id = cascade.GenRandomNode();
			int edgeVisited;
			vector<RRVec> tmpTable;
			cascade.ReversePropagate(1, id, tmpTable, edgeVisited);
#pragma omp critical
			{
				refTable.push_back(tmpTable[0]);
				refTargets.push_back(id);
				refEdgeVisited.push_back(edgeVisited);
			}
		}
	}
#endif

}

void RRInflBase::_AdaRebuildRRIndices()
{
	degrees.clear();
	degrees.resize(n, 0);
	degreeRRIndices.clear();

	for (int i = 0; i < n; ++i) {
		degreeRRIndices.push_back(vector<int>());
	}

	int count = 0;
	cout << "act nodes:" << act_nodes.size() << endl;
	// to count hyper edges:
	for (size_t i = 0; i < table.size(); ++i) {

		const RRVec& RR = table[i];
		if (act_nodes_arr[RR[0]] == 1){
			count++;
			continue;
		}
		for (int source : RR) {
			degrees[source]++;
			degreeRRIndices[source].push_back(i); // add index of table
		}
	}
	cout << "count:" << count << endl;
	// add to sourceSet where node's degree > 0
	sourceSet.clear();
	for (size_t i = 0; i < degrees.size(); ++i) {
		if (degrees[i] > 0) {
			sourceSet.insert(i);
		}
	}
}


// Same thing has been implemented as a part of _Greedy
double RRInflBase::_EstimateInfl(const vector<int>& seeds, vector<double>& out_cumu_spread)
{
	set<int> covered;
	double spd = 0;

	vector<bool> enables;
	enables.resize(table.size(), true);
	for (size_t i = 0; i < seeds.size(); ++i) {
		int sd = seeds[i];
		for (int idx : degreeRRIndices[sd]) {
			if (enables[idx]) {
				covered.insert(idx);
				enables[idx] = false;
			}
		}
		spd = (double)(n * covered.size()) / table.size();
		out_cumu_spread.push_back(spd);
	}
	return spd;
}



////////////////////////////////////////////////////////////////////
// RRInfl: paper [1]
double RRInfl::DefaultRounds(int n, int m, double epsilon)
{
	return max(144.0 * (n + m) / pow(epsilon, 3) * log(max(n, 1)), 1.0); // to make it positive
}

void RRInfl::BuildInError(graph_type& gf, int k, cascade_type& cascade, double epsilon/*=0.1*/)
{
	n = gf.GetN();
	m = gf.GetM();
	size_t num_iter = (size_t)(ceil(DefaultRounds(n, m, epsilon)));
	_Build(gf, k, cascade, num_iter);
}

void RRInfl::Build(graph_type& gf, int k, cascade_type& cascade, size_t num_iter/*=1000000*/)
{
	n = gf.GetN();
	m = gf.GetM();
	_Build(gf, k, cascade, num_iter);
}


void RRInfl::_Build(graph_type& gf, int k, cascade_type& cascade, size_t num_iter)
{
	InitializeConcurrent();

	n = gf.GetN();
	m = gf.GetM();

	top = k;
	d.resize(top, 0.0);
	list.resize(top, 0);
	cascade.Build(gf);

	cout << "#round = " << num_iter << endl;

	// table: target_id  RRset
	// id1  RR: rr1, rr2 ...
	// id2  RR: rr1, rr2 ...
	// ...
	table.clear();
	targets.clear();

	EventTimer pctimer;
	pctimer.SetTimeEvent("start");

	// Step 1:
	_AddRRSimulation(num_iter, cascade, table, targets);
	assert(targets.size() == num_iter);
	assert(table.size() == num_iter);
	pctimer.SetTimeEvent("step1");

	// Step 2:
	vector<int> seeds;
	vector<double> est_spread;
	_RebuildRRIndices();
	pctimer.SetTimeEvent("step2");

	// Step 3:
	double spd = _RunGreedy(k, seeds, est_spread);
	_SetResults(seeds, est_spread);
	pctimer.SetTimeEvent("end");

	cout << "  final (estimated) spread = " << spd << "\t round = " << num_iter << endl;

	// Write results to file:
	//FILE *out;
	//fopen_s(&out, file.c_str(), "w");
	//fprintf(out, "%d\n", top);
	//for (int i=0; i<top; i++)
	//	fprintf(out, "%d\t%g\n", list[i], d[i]);
	//fclose(out);
	WriteToFile(file, gf);

	FILE *timetmpfile;
	fopen_s(&timetmpfile, time_file.c_str(), "w");
	fprintf(timetmpfile, "%g\n", pctimer.TimeSpan("start", "end"));
	fprintf(timetmpfile, "Gen graph: %g\n", pctimer.TimeSpan("start", "step1"));
	fprintf(timetmpfile, "Build RR: %g\n", pctimer.TimeSpan("step1", "step2"));
	fprintf(timetmpfile, "Greedy: %g\n", pctimer.TimeSpan("step2", "end"));
	fclose(timetmpfile);
}




///////////////////////////////////////////////////////////////////////////
/// TimPlus: paper 2
void TimPlus::Build(graph_type& gf, int k, cascade_type& cascade, double eps/*=0.1*/, double ell/*=1.0*/)
{
	InitializeConcurrent();

	n = gf.GetN();
	m = gf.GetM();
	top = k;
	d.resize(top, 0.0);
	list.resize(top, 0);

	ell = ell + log(2) / log(n); // Original IMM has failure probability 2/n^ell, so use this transformation to 
	// make the failure probability 1/n^ell

	cascade.Build(gf);
	table.clear();
	targets.clear();

	EventTimer pctimer;
	EventTimer pctimer2;
	pctimer.SetTimeEvent("start");

	// Step1: 
	double est_opt1 = 0.0;
	size_t round1 = 0;
	{
		cout << "Step 1: estimate EPT" << endl;
		pctimer2.SetTimeEvent("est_start");
		int maxIth = (int)ceil(log(n) / log(2) - 1);
		for (int ith = 1; ith < maxIth; ++ith) {
			double lb = pow((double)0.5, ith);

			int loop = (int)(StepThreshold(n, lb, ell) + 1);

			vector<int> edgeVisited;
			_AddRRSimulation(loop, cascade, table, targets, edgeVisited);
			round1 += loop;

			double sumKappa = 0.0;
			for (int visit : edgeVisited) {
				double width = (double)visit / m;
				double kappa = 1.0 - pow(1.0 - width, k);
				sumKappa += kappa;
			}
			if (sumKappa / loop > lb) {
				// Algorithm 2 in [2] uses  (n * sumKappa) / (2.0 * loop), but their code does not have 2.0
				est_opt1 = (n * sumKappa) / (loop);
				break;
			}
		}
		est_opt1 = max(est_opt1, (double)k);
		pctimer2.SetTimeEvent("est_end");
	}
	cout << "  est_opt1 = " << est_opt1 << "\t round1 = " << round1 << endl;
	pctimer.SetTimeEvent("step1");
	double time_per_round = pctimer2.TimeSpan("est_start", "est_end") / ((double)round1);
	cout << "  (Time per round = " << time_per_round << ")" << endl;

	// Step2: Use greedy to estimate opt lower bound
	double eps_step2, spd2, est_opt2;
	size_t round2;
	{
		eps_step2 = EpsPrime(eps, k, ell);
		vector<int> seeds2;
		vector<double> est_spread2;
		double theta = RThreshold_0(eps, est_opt1, ell);
		round2 = (size_t)max(theta + 1.0, 1.0);
		cout << "Step 2: estimate opt by greedy. round = " << round2 << endl
			<< "  (Estimate time = " << time_per_round * round2 << ")" << endl;
		_AddRRSimulation(round2, cascade, table, targets);
		_RebuildRRIndices();
		spd2 = _RunGreedy(k, seeds2, est_spread2);
		est_opt2 = spd2 / (1 + eps_step2);
		est_opt2 = max(est_opt2, est_opt1);
	}
	cout << "  est_opt2 = " << est_opt2 << "\t round2 = " << round2 << endl;
	pctimer.SetTimeEvent("step2");

	// Step3: Set final runs
	vector<int> seeds3;
	vector<double> est_spread3;
	double spd3;
	size_t round3;
	{
		double theta = RThreshold(eps, est_opt2, k, ell);
		round3 = (size_t)max(theta + 1.0, 1.0);
		cout << "Step 3: final greedy. round = " << round3 << endl
			<< "  (Estimate time = " << time_per_round * round3 << ")" << endl;
		pctimer.SetTimeEvent("step3-1");
		_AddRRSimulation(round3, cascade, table, targets);
		pctimer.SetTimeEvent("step3-2");
		_RebuildRRIndices();
		pctimer.SetTimeEvent("step3-3");
		spd3 = _RunGreedy(top, seeds3, est_spread3);
	}
	cout << "  final spread = " << spd3 << "\t round3 = " << round3 << endl;
	_SetResults(seeds3, est_spread3);
	pctimer.SetTimeEvent("end");

	// Write results to file:
	//FILE *out;
	//fopen_s(&out, file.c_str(), "w");
	//fprintf(out, "%d\n", top);
	//for (int i=0; i<top; i++)
	//	fprintf(out, "%d\t%g\n", list[i], d[i]);
	//fclose(out);
	WriteToFile(file, gf);

	FILE *timetmpfile;
	fopen_s(&timetmpfile, time_file.c_str(), "w");
	fprintf(timetmpfile, "%g\n", pctimer.TimeSpan("start", "end"));
	fprintf(timetmpfile, "Step1: %g\n", pctimer.TimeSpan("start", "step1"));
	fprintf(timetmpfile, "Step2: %g\n", pctimer.TimeSpan("step1", "step2"));
	fprintf(timetmpfile, "Step3: %g\n", pctimer.TimeSpan("step2", "end"));
	fprintf(timetmpfile, "   Step3-1: %g\n", pctimer.TimeSpan("step3-1", "step3-2"));
	fprintf(timetmpfile, "   Step3-2: %g\n", pctimer.TimeSpan("step3-2", "step3-3"));
	fprintf(timetmpfile, "   Step3-3: %g\n", pctimer.TimeSpan("step3-3", "end"));
	fclose(timetmpfile);
}

double TimPlus::LogNChooseK(int n, int k){
	double ans = 0.0;
	assert(n >= k && k >= 0);
	for (int i = 0; i < k; i++){
		ans += log(n - i);
	}
	for (int i = 1; i <= k; i++){
		ans -= log(i);
	}
	return ans;
}

double TimPlus::RThreshold_0(double eps, double opt, double ell/*=1.0*/)
{
	double lambda = (double)(8 + 2 * eps) * n * (ell*log(n) + log(2)) / (eps*eps);
	double theta = lambda / (opt*4.0);
	return theta;
}

double TimPlus::RThreshold(double eps, double opt, int k, double ell/*=1.0*/)
{
	double lambda = (double)(8 + 2 * eps) * n * (ell*log(n) + LogNChooseK(n, k) + log(2)) / (eps*eps);
	double theta = lambda / opt;
	return theta;
}


double TimPlus::EpsPrime(double eps, int k, double ell/*=1.0*/)
{
	return 5.0 * pow(eps*eps*ell / ((double)ell + k), 1.0 / 3.0);
}

double TimPlus::StepThreshold(int n, double lb, double ell/*=1.0*/)
{
	double logValue = max(log(n) / log(2.0), 1.0);
	double loop = ((double)6 * ell * log(n) + 6 * log(logValue)) / lb;
	return loop;
}


///////////////////////////////////////////////////////////////////////////
/// IMM: paper 3
void IMM::Build(graph_type& gf, int k, cascade_type& cascade, double eps /* = 0.1 */, double ell /* = 1.0 */)
{
	InitializeConcurrent();

	EventTimer pctimer;
	vector<EventTimer> steptimers;
	pctimer.SetTimeEvent("start");

	n = gf.GetN();
	m = gf.GetM();
	top = k;
	d.resize(top, 0.0);
	list.resize(top, 0);

	ell = ell + log(2) / log(n); // Original IMM has failure probability 2/n^ell, so use this transformation to 
	// make the failure probability 1/n^ell

	cascade.Build(gf);

	table.clear();
	targets.clear();

	double epsprime = eps * sqrt(2.0); // eps'
	double LB = 1.0;   // lower bound
	double maxRounds = max(max(log2((double)n), 1.0) - 1.0, 1.0);

	// first a few passes of sampling
	double spread = 0.0;
	for (int r = 1; r < maxRounds; r++) {
		EventTimer stept;
		stept.SetTimeEvent("step_start");

		cout << "  Step" << r << ":";
		double x = max(double(n) / pow(2, r), 1.0);
		double theta = LambdaPrime(epsprime, k, ell, n) / x;

		// select nodes
		vector<int> seeds;
		vector<double> est_spread;
		LARGE_INT64 nNewSamples = LARGE_INT64(theta - table.size() + 1);
		if (table.size() < theta) {
			// generate samples
			_AddRRSimulation(nNewSamples, cascade, table, targets);
			_RebuildRRIndices();
			spread = _RunGreedy(top, seeds, est_spread);
			_SetResults(seeds, est_spread);
		}
		cout << " spread = " << spread << "\t round = " << ((nNewSamples > 0) ? nNewSamples : 0) << endl;

		stept.SetTimeEvent("step_end");
		steptimers.push_back(stept);

		// check influence and set lower bound
		if (spread >= ((1.0 + epsprime) * x)) {
			LB = spread / (1.0 + epsprime);
			break;
		}
	}

	// final pass of sampling
	{
		EventTimer stept;
		stept.SetTimeEvent("step_start");
		cout << "  Estimated Lower bound: " << LB << endl;
		cout << "  Final Step:";
		double theta = LambdaStar(eps, k, ell, n) / LB;

		vector<int> seeds;
		vector<double> est_spread;
		LARGE_INT64 nNewSamples = LARGE_INT64(theta - table.size() + 1);
		if (table.size() < theta) {
			// generate samples
			_AddRRSimulation(nNewSamples, cascade, table, targets);
			_RebuildRRIndices();
			spread = _RunGreedy(top, seeds, est_spread);
			_SetResults(seeds, est_spread);
		}
		cout << " spread = " << spread << "\t round = " << ((nNewSamples > 0) ? nNewSamples : 0) << endl;

		stept.SetTimeEvent("step_end");
		steptimers.push_back(stept);
	}
	pctimer.SetTimeEvent("end");

	// Write results to file:
	//FILE *out;
	//fopen_s(&out, file.c_str(), "w");
	//fprintf(out, "%d\n", top);
	//for (int i = 0; i<top; i++)
	//	fprintf(out, "%d\t%g\n", list[i], d[i]);
	//fclose(out);
	WriteToFile(file, gf);

	FILE *timetmpfile;
	fopen_s(&timetmpfile, time_file.c_str(), "w");
	fprintf(timetmpfile, "%g\n", pctimer.TimeSpan("start", "end"));
	for (size_t t = 0; t < steptimers.size(); t++) {
		if (t != steptimers.size() - 1) {
			fprintf(timetmpfile, "  Step%d: %g\n", t + 1, steptimers[t].TimeSpan("step_start", "step_end"));
		}
		else {
			fprintf(timetmpfile, "  Final Step: %g\n", steptimers[t].TimeSpan("step_start", "step_end"));
		}
	}
	fclose(timetmpfile);
}

double IMM::LambdaPrime(double epsprime, int k, double ell, int n)
{
	static double SMALL_DOUBLE = 1e-16;
	double cst = (2.0 + 2.0 / 3.0 * epsprime) / max(epsprime * epsprime, SMALL_DOUBLE);
	double part2 = LogNChooseK(n, k);
	part2 += ell * log(max((double)n, 1.0));
	part2 += log(max(log2(max((double)n, 1.0)), 1.0));
	// calc \lambda'
	double lambda = cst * part2 * n;
	return lambda;
}

double IMM::LambdaStar(double eps, int k, double ell, int n)
{
	static double SMALL_DOUBLE = 1e-16;
	static double APPRO_RATIO = (1.0 - 1.0 / exp(1));
	double logsum = ell * log(n) + log(2);
	// calc \alpha and \beta
	double alpha = sqrt(max(logsum, SMALL_DOUBLE));
	double beta = sqrt(APPRO_RATIO * (LogNChooseK(n, k) + logsum));
	// calc \lambda*
	double lambda = 2.0 * n / max(pow(eps, 2), SMALL_DOUBLE);
	lambda = lambda * pow(APPRO_RATIO * alpha + beta, 2);
	return lambda;
}


void ShapleyInfl::Shapley_Build(graph_type& gf, cascade_type& cascade, GeneralCascade & gc, double eps/*=0.1*/, double ell/*=1.0*/, int topk /*=50*/, bool isSingleInf /* = false */)
{

	InitializeConcurrent();

	EventTimer pctimer;
	EventTimer pctimer2;
	vector<EventTimer> steptimers;

	pctimer.SetTimeEvent("start");

	n = gf.GetN();
	m = gf.GetM();

	std::vector<double> shapleyV;
	std::vector<int> hitCount;
	LARGE_INT64 totalEdgeVisited = 0;
	LARGE_INT64 totalRRSetSoFar = 0; // the number of RR sets already generated

	double shpv = 0.0;

	shapleyV.assign(n, 0.0);
	hitCount.assign(n, 0);

	shapleyValueId nullSVI;
	nullSVI.sid = 0;
	nullSVI.value = 0.0;


	top = 1;  // for legacy reason

	cascade.Build(gf);
	table.clear();
	targets.clear();

	// when n is too small, need to use a smaller topk. A topk that is equal to n is also not good,
	// violating the assumption that \phi_k \ge 1. Thus use n/2. Only for very small graph.
	topk = min(n / 2, topk);

	// Phase 1: Estimate the number of RR sets needed
	cout << "Phase 1: Estimate the number of RR sets needed" << endl;
	pctimer2.SetTimeEvent("phase1_start");
	double LB = 1.0;   // lower bound

	// use the IMM sampling method to estimate lower bound LB
	cout << "lower bound LB estimated using the sampling method similar to IMM" << endl;
	double epsprime = eps * sqrt(2.0); // eps'
	double maxRounds = max(max(log2(n), 1.0) - 1.0, 1.0);

	// first a few passes of sampling, based on the IMM sampling method
	for (int r = 1; r < maxRounds; r++) {
		EventTimer stept;
		stept.SetTimeEvent("step_start");

		cout << "  Step" << r << ":";
		double x = max(double(n) / pow(2, r), 1.0);
		LARGE_INT64 theta = (LARGE_INT64)ceil(LambdaPrime(epsprime, top, ell + log(2) / log(n), n) / x);
		LARGE_INT64 nNewSamples = (LARGE_INT64)(theta - totalRRSetSoFar);

		if (totalRRSetSoFar < theta) {
			// generate samples
			Shapley_AddRRSimulation(nNewSamples, cascade, shapleyV, hitCount, table, totalEdgeVisited);
			totalRRSetSoFar = theta;

			// find the k-th largest Shapley value or SNI value
			shapleyVId.assign(n, nullSVI);

			if (!isSingleInf) {
				for (int i = 0; i < n; i++) {
					shapleyVId[i].sid = i;
					shapleyVId[i].value = shapleyV[i];
				}
			}
			else { /* if computing single node influence, use the shapleyVId array, but use hitCount collected from the simulation */
				for (int i = 0; i < n; i++) {
					shapleyVId[i].sid = i;
					shapleyVId[i].value = (double)hitCount[i];
				}
			}

			std::sort(shapleyVId.begin(), shapleyVId.end(), ShapleyComparator());

			shpv = 1.0 * shapleyVId[topk - 1].value / totalRRSetSoFar * n;
		}

		if (!isSingleInf) {
			cout << "No. " << topk << " largest Shapley value = " << shpv << "\t current round = " << ((nNewSamples > 0) ? nNewSamples : 0) << "\t total round = " << totalRRSetSoFar << endl;
		}
		else {
			cout << "No. " << topk << " largest SNI value = " << shpv << "\t current round = " << ((nNewSamples > 0) ? nNewSamples : 0) << "\t total round = " << totalRRSetSoFar << endl;
		}
		stept.SetTimeEvent("step_end");
		steptimers.push_back(stept);

		// check stopping condition and set lower bound
		if (shpv >= ((1.0 + epsprime) * x)) {
			LB = shpv / (1.0 + epsprime);
			break;
		}
	}


	LARGE_INT64 stheta = (LARGE_INT64)ceil(n*((ell + 1)* log(n) + log(4))* (2 + 2 / 3 * eps) / (eps * eps * LB));

	pctimer2.SetTimeEvent("phase1_end");
	if (!isSingleInf) {
		cout << "Phase 1 ends: Estimated lower bound of the k-th largest Shapley value: " << LB << endl;
	}
	else {
		cout << "Phase 1 ends: Estimated lower bound of the k-th largest SNI value: " << LB << endl;
	}
	cout << "              Number of RR sets needed for Phase 2: " << stheta << endl;

	// pctimer.SetTimeEvent("step1");

	// Phase 2: generate RR sets and estimate Shapley values of nodes
	// need to start from the scratch, because otherwise the estimates could be unbiased.
	// Martingale property cannot guarantee unbiasness

	if (!isSingleInf) {
		cout << "Phase 2: generate RR sets and estimate Shapley values of all nodes, theta = " << stheta << endl;
	}
	else {
		cout << "Phase 2: generate RR sets and estimate SNI values of all nodes, theta = " << stheta << endl;
	}

	shapleyV.assign(n, 0.0);
	hitCount.assign(n, 0);
	totalEdgeVisited = 0;

	Shapley_AddRRSimulation(stheta, cascade, shapleyV, hitCount, table, totalEdgeVisited);

	if (!isSingleInf) {
		for (int i = 0; i < n; i++) {
			shapleyV[i] = 1.0 * shapleyV[i] * n / stheta;
		}
	}
	else {
		for (int i = 0; i < n; i++) {
			shapleyV[i] = 1.0 * hitCount[i] * n / stheta;
		}
	}


	// pctimer.SetTimeEvent("step2");
	pctimer2.SetTimeEvent("phase2_end");

	double EPT = 1.0 * totalEdgeVisited / stheta;

	cout << "Phase 2 ends" << endl;


	// sort Shapley values. Since we need to sort node id with the value, we need another structure for this

	shapleyVId.assign(n, nullSVI);

	for (int i = 0; i < n; i++) {
		shapleyVId[i].sid = i;
		shapleyVId[i].value = shapleyV[i];
	}

	std::sort(shapleyVId.begin(), shapleyVId.end(), ShapleyComparator());

	if (!isSingleInf) {
		cout << "k-th largest Shapley value = " << shapleyVId[topk - 1].value << "; EPT = " << EPT << endl;
	}
	else {
		cout << "k-th largest SNI value = " << shapleyVId[topk - 1].value << "; EPT = " << EPT << endl;
	}

	// this is to indicate that we have found n seeds, from the Shapley value, to be used
	// later for simulation
	top = n;

	d.resize(top, 0.0);
	list.resize(top, 0);
	// fill list[] and d[] arrays
	for (int i = 0; i < n; i++) {
		list[i] = shapleyVId[i].sid;
		d[i] = shapleyVId[i].value;
	}



	// Write results to file:
	FILE *out;
	int nshapley = n;
	fopen_s(&out, file.c_str(), "w");
	fprintf(out, "%d\n", nshapley);
	if (!isSingleInf) {
		fprintf(out, "# nodeid\t Shapley value\t  epsilon = %Lg; ell = %Lg; EPT = %Lg\n", eps, ell, EPT);
	}
	else {
		fprintf(out, "# nodeid\t SNI value\t  epsilon = %Lg; ell = %Lg; EPT = %Lg\n", eps, ell, EPT);
	}
	for (int i = 0; i<nshapley; i++)
		fprintf(out, "%s\t%Lg\n", gf.MapIndexToNodeName(list[i]).c_str(), d[i]);
	fclose(out);


	FILE *timetmpfile;
	fopen_s(&timetmpfile, time_file.c_str(), "w");
	fprintf(timetmpfile, "%lg\n", pctimer2.TimeSpan("phase1_start", "phase2_end"));
	fclose(timetmpfile);
}


void ShapleyInfl::Shapley_AddRRSimulation(LARGE_INT64 num_iter,
	cascade_type& cascade,
	std::vector<double>& shapleyV,
	std::vector<int>& hitCount,  // hitCount[i] is the number of times the RR sets hit node i
	std::vector< std::vector<int> >& refTable,
	LARGE_INT64 & totalEdgeVisited)
	// generate one RR set, and getting the statistics and the estimate of the Shapley values.
{
	for (LARGE_INT64 iter = 0; iter < num_iter; ++iter) {
		int id = cascade.GenRandomNode();
		int edge_visited_count;

		// we do not need the keep the RR sets, so we can disgard them before generating the next one. 
		// I do not want to redefine ReverseGCascade::ReversePropagate, so use clear() here.
		refTable.clear();

		cascade.ReversePropagate(1, id, refTable, edge_visited_count);

		vector<int>& lastRRset = *(refTable.rbegin());
		for (int i = 0; i < lastRRset.size(); i++) {
			shapleyV[lastRRset[i]] += 1.0 / lastRRset.size();
			hitCount[lastRRset[i]] += 1;
		}

		totalEdgeVisited += edge_visited_count;

	}
}


bool ShapleyInfl::ShapleyComparator::operator() (shapleyValueId a, shapleyValueId b) {
	return (a.value > b.value);
}

///////////////////////////////////////////////////////////////////////////
/// AdaIMM: KDD 2018
void RIM_IMM::Build(graph_type& gf, int k, cascade_type& cascade, double eps /* = 0.1 */, double ell /* = 1.0 */, int rounds)
{
	InitializeConcurrent();

	EventTimer pctimer;
	vector<EventTimer> steptimers;
	pctimer.SetTimeEvent("start");
	vector<int> cur_act_nodes;
	vector<double> res_time;

	n = gf.GetN();
	m = gf.GetM();
	top = k;
	T = rounds;
	d.resize(top, 0.0);
	list.resize(top, 0);
	int previous_act_nodes = 0;
	ell = ell + log(2*T) / log(n); // Original IMM has failure probability 2/n^ell, so use this transformation to 
	// make the failure probability 1/n^ell

	cascade.Build(gf);

	table.clear();
	targets.clear();
	act_nodes.clear();
	seeds_list.clear();
	cur_act_nodes.clear();

	double epsprime = eps * sqrt(2.0); // eps'
	double LB = 1.0;   // lower bound
	double maxRounds = max(max(log2((double)(n - act_nodes.size())), 1.0) - 1.0, 1.0);


	double xprime = n;
	for (int round = 0; round <= T; round++){
		table.clear();
		targets.clear();
		memset(act_nodes_arr, 0, sizeof(bool)* 30000);
		EventTimer timer;
		timer.SetTimeEvent("start");
		cout << "round:" << round << endl;
		// first a few passes of sampling
		if (round > 0){

			int targetSize = k;
			ProbTransfom trans(gf.edgeForm);

			// single thread
			int* tmplist = new int[n];
			bool* active = new bool[n];
			for (int it = 0; it < 1; it++)
			{
				memset(active, 0, sizeof(bool)*n);
				for (int i = 0; i < targetSize; i++){
					tmplist[i] = seeds_list[round - 1][i];
					active[tmplist[i]] = true;
					act_nodes.insert(tmplist[i]);
					act_nodes_arr[tmplist[i]] = 1;
				}

				int h = 0;
				int t = targetSize;
				while (h < t)
				{
					int kl = gf.GetNeighborCount(tmplist[h]);
					//printf("%d \n",k);

					for (int i = 0; i < kl; i++)
					{
						auto& e = gf.GetEdge(tmplist[h], i);
						if (active[e.v]) continue;
						if (random.RandBernoulli(trans.Prob(e.w1)))
						{
							tmplist[t] = e.v;
							active[e.v] = true;
							t++;
							act_nodes.insert(e.v);
							act_nodes_arr[e.v] = 1;
							//break;
						}
					}
					h++;
				}
			}
			SAFE_DELETE_ARRAY(active);
			SAFE_DELETE_ARRAY(tmplist);
		}

		if (round > 0)
			cur_act_nodes.push_back(act_nodes.size());
		if (round == rounds)
			break;

		double spread = 0.0;

		for (int r = 1; r < maxRounds; r++) {
			EventTimer stept;
			stept.SetTimeEvent("step_start");

			cout << "  Step" << r << ":";
			double x;
			int marginal_act_nodes = act_nodes.size() - previous_act_nodes;
			if (round == 0)
				x = max((double(xprime) / pow(2, r)), 1.0);
			else
				//x = max(double(xprime - marginal_act_nodes) / pow(2, (r - 1)), 1.0);
				x = max(min(double(xprime) / pow(2, (r - 1)), (n - act_nodes.size())), 1.0);

			double theta = LambdaPrime(epsprime, k, ell, n - act_nodes.size()) / x;

			// select nodes
			vector<int> seeds;
			vector<double> est_spread;
			LARGE_INT64 nNewSamples = LARGE_INT64(theta - table.size() + 1);
			cout << "  theta and table size:" << theta << " " << table.size() << endl;
			if (table.size() < theta) {
				// generate samples
				_AdaAddRRSimulation(nNewSamples, cascade, table, targets);
				_AdaRebuildRRIndices();
				spread = _RunGreedy(top, seeds, est_spread);
			}
			else if (round > 0){
				_AdaRebuildRRIndices();
				spread = _RunGreedy(top, seeds, est_spread);
			}

			cout << " spread = " << spread << "\t round = " << ((nNewSamples > 0) ? nNewSamples : 0) << endl;
			stept.SetTimeEvent("step_end");
			steptimers.push_back(stept);

			// check influence and set lower bound
			if (spread >= ((1.0 + epsprime) * x)) {
				LB = spread / (1.0 + epsprime);
				xprime = x;
				break;
			}
		}

		// final pass of sampling
		{
			EventTimer stept;
			stept.SetTimeEvent("step_start");
			cout << "  Estimated Lower bound: " << LB << endl;
			cout << "  Final Step:";
			double theta = LambdaStar(eps, k, ell, n) / LB;

			vector<int> seeds;
			vector<double> est_spread;
			LARGE_INT64 nNewSamples = LARGE_INT64(theta - table.size() + 1);
			cout << "  theta and table size:" << theta << " " << table.size() << endl;
			if (table.size() < theta) {
				// generate samples
				_AdaAddRRSimulation(nNewSamples, cascade, table, targets);
				_AdaRebuildRRIndices();
				//_SetResults(seeds, est_spread);
			}
			spread = _RunGreedy(top, seeds, est_spread);
			cout << " spread = " << spread << "\t round = " << ((nNewSamples > 0) ? nNewSamples : 0) << endl;

			stept.SetTimeEvent("step_end");
			steptimers.push_back(stept);
			seeds_list.push_back(seeds);
		}
		previous_act_nodes = act_nodes.size();
		timer.SetTimeEvent("end");
		res_time.push_back(timer.TimeSpan("start", "end"));
		cout << "time: " << timer.TimeSpan("start", "end") << endl;

	}

	for (const auto &s : cur_act_nodes) std::cout << s << '\t';
	std::cout << std::endl;

	for (const auto &s : res_time) std::cout << s << '\t';
	std::cout << std::endl;

	for (const auto &row : seeds_list)
	{
		for (const auto &s : row) std::cout << s << '\t';
		std::cout << std::endl;
	}

	{
		ofstream myfile;
		std::stringstream ss;
		ss << "rim_imm" << ".csv";
		std::string fullName = ss.str();

		myfile.open(fullName);
		for (const auto &s : cur_act_nodes) myfile << s << ',';
		myfile << '\n';

		for (const auto &s : res_time) myfile << s << ',';
		myfile << '\n';

		for (const auto &row : seeds_list)
		{
			for (const auto &s : row) myfile << s << ',';
			myfile << '\n';
		}
		myfile.close();
	}
	pctimer.SetTimeEvent("end");
}



///////////////////////////////////////////////////////////////////////////
/// WR_IMM: KDD 2018
void WR_IMM::Build(graph_type& gf, int k, cascade_type& cascade, double eps /* = 0.1 */, double ell /* = 1.0 */, int rounds)
{
	InitializeConcurrent();

	EventTimer pctimer;
	vector<EventTimer> steptimers;
	vector<double> res_time;
	pctimer.SetTimeEvent("start");

	n = gf.GetN();
	m = gf.GetM();
	top = k;
	T = rounds;
	d.resize(top, 0.0);
	list.resize(top, 0);

	ell = ell + log(2 * T) / log(n); // Original IMM has failure probability 2/n^ell, so use this transformation to 
	// make the failure probability 1/n^ell

	cascade.Build(gf);
	enables.clear();
	seeds_list.clear();

	double epsprime = eps * sqrt(2.0); // eps'
	double LB = 1.0;   // lower bound
	double maxRounds = max(max(log2((double)n), 1.0) - 1.0, 1.0);

	// first a few passes of sampling
	double spread = 0.0;
	double alpha = 1.0;

	for (int t = 0; t < T; t++){

		EventTimer timer;
		timer.SetTimeEvent("start");

		LARGE_INT64 length = targets.size();
		cout << length << endl;

		int *root_arr = new int[length];
		R = 0;
		for (int i = 0; i < enables.size(); i++){
			if (enables[i]){
				root_arr[R] = targets[i];
				R++;
			}
		}
		cout << "root and target size:" << R << ' ' << targets.size() << ' ' << endl;
		if (t > 0){
			alpha = alpha * R / targets.size();
		}


		table.clear();
		targets.clear();
		cout << "round:" << t << endl;

		for (int r = 1; r < maxRounds; r++) {
			EventTimer stept;
			stept.SetTimeEvent("step_start");

			cout << "  Step" << r << ":";
			double x = max(double(n) / pow(2, r), 1.0);
			double theta = LambdaPrime(epsprime, k, ell, n) / x;

			// select nodes
			vector<int> seeds;
			vector<double> est_spread;
			LARGE_INT64 nNewSamples = LARGE_INT64(theta - table.size() + 1);

			if (table.size() < theta) {
				// generate samples
				_WRIMAddRRSimulation(nNewSamples, cascade, table, targets, root_arr, t);
				_RebuildRRIndices();
				spread = _RunGreedy(top, seeds, est_spread);
				//_SetResults(seeds, est_spread);
			}

			cout << " spread = " << spread << "\t round = " << ((nNewSamples > 0) ? nNewSamples : 0) << endl;

			stept.SetTimeEvent("step_end");
			steptimers.push_back(stept);

			// check influence and set lower bound
			if (alpha * spread >= ((1.0 + epsprime) * x)) {
				LB = alpha * spread / (1.0 + epsprime);
				break;
			}
		}
		// final pass of sampling
		{
			EventTimer stept;
			stept.SetTimeEvent("step_start");
			cout << "  Estimated Lower bound: " << LB << endl;
			cout << "  Final Step:";
			double theta = LambdaStar(eps, k, ell, n) / LB;

			vector<int> seeds;
			vector<double> est_spread;
			LARGE_INT64 nNewSamples = LARGE_INT64(theta - table.size() + 1);
			if (table.size() < theta) {
				// generate samples
				_WRIMAddRRSimulation(nNewSamples, cascade, table, targets, root_arr, t);
				_RebuildRRIndices();
				spread = _RunGreedy(top, seeds, est_spread);
				//_SetResults(seeds, est_spread);
			}
			seeds_list.push_back(seeds);
			cout << " spread = " << spread << "\t round = " << ((nNewSamples > 0) ? nNewSamples : 0) << endl;

			stept.SetTimeEvent("step_end");
			steptimers.push_back(stept);
		}
		pctimer.SetTimeEvent("end");

		delete[] root_arr;

		/*root_arr.clear();
		for (int i = 0; i < enables.size(); i++){
		if (enables[i]){
		root_arr.push_back(targets[i]);
		R++;
		}
		}*/

		/*for (int i = 0; i < R; i++){
		cout << "root arr:" << i << ' ' << root_arr[i] << endl;
		}*/
		timer.SetTimeEvent("end");
		res_time.push_back(timer.TimeSpan("start", "end"));
		cout << "time: " << timer.TimeSpan("start", "end") << endl;


	}

	for (const auto &row : seeds_list)
	{
		for (const auto &s : row) std::cout << s << '\t';
		std::cout << std::endl;
	}

	{
		ofstream myfile;
		std::stringstream ss;
		ss << "rr_wr_imm" << ".csv";
		std::string fullName = ss.str();
		myfile.open(fullName);

		for (const auto &s : res_time) myfile << s << ',';
		myfile << '\n';

		for (const auto &row : seeds_list)
		{
			for (const auto &s : row) myfile << s << ',';
			myfile << '\n';
		}
		myfile.close();
	}

	{
		ofstream myfile;
		myfile.open("rr_wr_imm.txt");

		for (const auto &row : seeds_list)
		{
			for (const auto &s : row) myfile << s << '\t';
			myfile << '\n';
		}
		myfile.close();
	}
	// Write results to file:
	//FILE *out;
	//fopen_s(&out, file.c_str(), "w");
	//fprintf(out, "%d\n", top);
	//for (int i = 0; i<top; i++)
	//	fprintf(out, "%d\t%g\n", list[i], d[i]);
	//fclose(out);
	WriteToFile(file, gf);

	FILE *timetmpfile;
	fopen_s(&timetmpfile, time_file.c_str(), "w");
	fprintf(timetmpfile, "%g\n", pctimer.TimeSpan("start", "end"));
	for (size_t t = 0; t < steptimers.size(); t++) {
		if (t != steptimers.size() - 1) {
			fprintf(timetmpfile, "  Step%d: %g\n", t + 1, steptimers[t].TimeSpan("step_start", "step_end"));
		}
		else {
			fprintf(timetmpfile, "  Final Step: %g\n", steptimers[t].TimeSpan("step_start", "step_end"));
		}
	}
	fclose(timetmpfile);
}


///////////////////////////////////////////////////////////////////////////
/// IMM: paper 3
void CR_IMM::Build(graph_type& gf, int k, cascade_type& cascade, double eps /* = 0.1 */, double ell /* = 1.0 */, int rounds)
{
	InitializeConcurrent();

	EventTimer pctimer;
	vector<EventTimer> steptimers;
	pctimer.SetTimeEvent("start");

	n = gf.GetN();
	m = gf.GetM();
	top = k;
	T = rounds;
	//cRound = new int[T];
	d.resize(top * T, 0.0);
	list.resize(top * T, 0);

	ell = ell + log(2*T) / log(n); // Original IMM has failure probability 2/n^ell, so use this transformation to 
	// make the failure probability 1/n^ell

	cascade.Build(gf);

	tableVec.clear();
	targets.clear();
	seeds_list.clear();

	double epsprime = eps * sqrt(2.0); // eps'
	double LB = 1.0;   // lower bound
	double maxRounds = max(max(log2((double)n), 1.0) - 1.0, 1.0);
	
	EventTimer timer;
	timer.SetTimeEvent("start");

	// first a few passes of sampling
	double spread = 0.0;
	for (int r = 1; r < maxRounds; r++) {
		EventTimer stept;
		stept.SetTimeEvent("step_start");

		cout << "  Step" << r << ":";
		double x = max(double(n) / pow(2, r), 1.0);
		double theta =  LambdaPrime_CR(epsprime, k, ell, n) / x;

		// select nodes
		vector<int> seeds;
		vector<double> est_spread;
		LARGE_INT64 nNewSamples = LARGE_INT64(theta - tableVec.size() + 1);
		if (tableVec.size() < theta) {
			// generate samples
			_CRIMAddRRSimulation(nNewSamples, cascade, tableVec, targets);
			_CRIMRebuildRRIndices();
			spread = _CRIMRunGreedy(top, seeds_list, est_spread);
		}
		cout << " spread = " << spread << "\t round = " << ((nNewSamples > 0) ? nNewSamples : 0) << endl;

		stept.SetTimeEvent("step_end");
		steptimers.push_back(stept);

		// check influence and set lower bound
		if (spread >= ((1.0 + epsprime) * x)) {
			LB = spread / (1.0 + epsprime);
			break;
		}
	}

	// final pass of sampling
	{
		EventTimer stept;
		stept.SetTimeEvent("step_start");
		cout << "  Estimated Lower bound: " << LB << endl;
		cout << "  Final Step:";
		double theta = LambdaStar_CR(eps, k, ell, n) / LB;

		vector<int> seeds;
		vector<double> est_spread;
		LARGE_INT64 nNewSamples = LARGE_INT64(theta - tableVec.size() + 1);
		if (tableVec.size() < theta) {
			// generate samples
			_CRIMAddRRSimulation(nNewSamples, cascade, tableVec, targets);
			_CRIMRebuildRRIndices();
			spread = _CRIMRunGreedy(top, seeds_list, est_spread);
		}
		cout << " spread = " << spread << "\t round = " << ((nNewSamples > 0) ? nNewSamples : 0) << endl;

		stept.SetTimeEvent("step_end");
		steptimers.push_back(stept);
	}
	pctimer.SetTimeEvent("end");

	timer.SetTimeEvent("end");
	cout << "time: " << timer.TimeSpan("start", "end") << endl;

	for (const auto &row : seeds_list)
	{
		for (const auto &s : row) std::cout << s << '\t';
		std::cout << std::endl;
	}

	{
		ofstream myfile;
		std::stringstream ss;
		ss << "rr_cr_imm" << ".txt";
		std::string fullName = ss.str();

		myfile.open(fullName);

		//for (const auto &s : res_time) myfile << s << ',';
		//myfile << '\n';

		for (const auto &row : seeds_list)
		{
			for (const auto &s : row) myfile << s << '\t';
			myfile << '\n';
		}
		myfile.close();
	}

	// Write results to file:
	//FILE *out;
	//fopen_s(&out, file.c_str(), "w");
	//fprintf(out, "%d\n", top);
	//for (int i = 0; i<top; i++)
	//	fprintf(out, "%d\t%g\n", list[i], d[i]);
	//fclose(out);
	WriteToFile(file, gf);

	FILE *timetmpfile;
	fopen_s(&timetmpfile, time_file.c_str(), "w");
	fprintf(timetmpfile, "%g\n", pctimer.TimeSpan("start", "end"));
	for (size_t t = 0; t < steptimers.size(); t++) {
		if (t != steptimers.size() - 1) {
			fprintf(timetmpfile, "  Step%d: %g\n", t + 1, steptimers[t].TimeSpan("step_start", "step_end"));
		}
		else {
			fprintf(timetmpfile, "  Final Step: %g\n", steptimers[t].TimeSpan("step_start", "step_end"));
		}
	}
	fclose(timetmpfile);
}

double CR_IMM::LambdaPrime_CR(double epsprime, int k, double ell, int n)
{
	static double SMALL_DOUBLE = 1e-16;
	double cst = (2.0 + 2.0 / 3.0 * epsprime) / max(epsprime * epsprime, SMALL_DOUBLE);
	double part2 = T * LogNChooseK(n, k);
	part2 += ell * log(max((double)n, 1.0));
	part2 += log(max(log2(max((double)n, 1.0)), 1.0));
	// calc \lambda'
	double lambda = cst * part2 * n;
	return lambda;
}

double CR_IMM::LambdaStar_CR(double eps, int k, double ell, int n)
{
	static double SMALL_DOUBLE = 1e-16;
	static double APPRO_RATIO = (1.0 - 1.0 / exp(1));
	double logsum = ell * log(n) + log(2);
	// calc \alpha and \beta
	double alpha = sqrt(max(logsum, SMALL_DOUBLE));
	double beta = sqrt(APPRO_RATIO * (LogNChooseK(n, k) + logsum));
	// calc \lambda*
	double lambda = 2.0 * n * T / max(pow(eps, 2), SMALL_DOUBLE);
	lambda = lambda * pow(APPRO_RATIO * alpha + beta, 2);
	return lambda;
}
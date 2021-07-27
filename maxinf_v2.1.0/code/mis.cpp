#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <sstream>
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>

#include <vector>
#include <string>
#include <map>
#include "mis.h"
#include "graph.h"

// #include <windows.h>
#include "event_timer.h"
#include "common.h"

using namespace std;

MIS::MIS()
{
	file = "mis_estimate.txt";
	timerfile = "time_mis.txt";
}

bool MIS::MINodeGreater::operator() (const MINode& left, const MINode& right) const
{
	if(left.influence > right.influence)
		return true;
	else if(left.influence == right.influence && left.num > right.num)
		return true;
	return false;
}

int MIS::FindSeedByNum(vector<MINode>& candidates, int num)
{
	for (size_t i = 0; i < candidates.size(); ++i) {
		if (candidates[i].num == num) {
			return i;
		}
	}
	return -1;
} 

void MIS::EstimateSeeds(std::vector<double>& distribution, Graph& gf)
{
	vector<int> each_indices;
	vector<MINode> candidates;

	// LARGE_INTEGER start, tend, freq;
	// QueryPerformanceFrequency(&freq);
	// QueryPerformanceCounter(&start);
    EventTimer timer;
    timer.SetTimeEvent("start");

	// begin process:
	// find valid dimensions:
	for (int i = 0; i < TopicAwareBase::ntopics; ++i) {
		double dist = distribution[i];
		if (dist >= EPS) {
			vector<TASeeds>& ss = TopicAwareBase::topicVector[i];
			// cout <<  "Topic[" <<  i+1 <<  "] " << dist << " ";
			int index = TopicAwareBase::RoundSeedsToIndex(ss, dist);
			if (index >= 0) {
				each_indices.push_back(index);
				// cout << "\t" << ss.vlambda << " " << ss.seeds.size() << endl; 
				continue;
			}
		}
		each_indices.push_back(-1); // invalid
	}

	// sum up marginal influence: heuristic
	for (int t = 0; t < TopicAwareBase::ntopics; ++t) {
		int index = each_indices[t];
		if (index < 0) { continue; } // skip invalid 
		
		TASeeds& ss = TopicAwareBase::topicVector[t][index];
		vector<MINode>& possible = ss.seeds;
		
		// try possible seeds, and add its marginal influence 
		for (size_t i = 0; i < possible.size(); ++i) {
			int num = possible[i].num;
			double influence = possible[i].influence;
			int pos = FindSeedByNum(candidates, num);
			if (pos == -1) {
				MINode temp;
				temp.num = num;
				temp.influence = 0;
				candidates.push_back(temp);
				pos = candidates.size() - 1;
			}
			candidates[pos].influence += influence;
		}
	}
	
	// sort by influence
	sort(candidates.begin(), candidates.end(), MINodeGreater());
	
	// end process
	// QueryPerformanceCounter(&tend);
    timer.SetTimeEvent("end");

	// print rounding results:
	cout << "Pick candidates after rounding to:" << endl;
	for (int t = 0; t < TopicAwareBase::ntopics; ++t) {
		int index = each_indices[t];
		if (index < 0) { cout << "/ "; } // skip invalid
		else {
			TASeeds& ss = TopicAwareBase::topicVector[t][index];
			cout << ss.vlambda << " ";
		}
	}
	cout << endl;


	// double timer = (double)(tend.QuadPart - start.QuadPart) / freq.QuadPart;
	
	FILE* out;
	fopen_s(&out, timerfile.c_str(), "w");
	fprintf(out, "%g", timer.TimeSpan("start", "end"));
	fclose(out);

	// prepare write to file
	size_t max_len = (candidates.size() <= (size_t)top) ? candidates.size() : (size_t)top;
	d.resize(max_len);
	list.resize(max_len);
	for (size_t i = 0; i < max_len; ++i) {
		list[i] = candidates[i].num;
		d[i] = candidates[i].influence;
	}
	top = max_len;


	//FILE* dout;
	//fopen_s(&dout, file.c_str(), "w");
	//fprintf(dout, "%d\n", max_len);
	//for (size_t i = 0; i < max_len; ++i) {
	//	fprintf(dout, "%d\t%g\n", list[i], d[i]);
	//}
	//fclose(dout);

	WriteToFile(file, gf);
}

void MIS::Build(Graph& gf, int set_size, std::vector<double>& distribution)
{
	n = gf.GetN();
	top = set_size;
	EstimateSeeds(distribution, gf);
}


#include <vector>
// #include "windows.h"
#include <iostream>
#include "common.h"
#include "event_timer.h"
#include "top_selection.h"
#include "graph.h"

using namespace std;

TopSelection::TopSelection()
{
	file = "top_selection_estimate.txt";
	timerfile = "time_top_selection.txt";
}

void TopSelection::SelectTop(std::vector<double>& distribution, Graph& gf)
{
	// LARGE_INTEGER start, tend, freq;
    EventTimer timer;
	vector<int> each_indices;
	vector<MINode> candidates;
	
    timer.SetTimeEvent("start");
	// QueryPerformanceFrequency(&freq);
    // QueryPerformanceCounter(&start);

	// begin process:
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

	double max_influence = -99999999;
	int best_topic = 0;
	int select_index = 0;
	for (int t = 0; t < TopicAwareBase::ntopics; ++t) {
		int index = each_indices[t];
		if (index < 0) { continue; } // skip invalid 
		
		TASeeds& ss = TopicAwareBase::topicVector[t][index];
		vector<MINode>& nodes = ss.seeds;
		// valid topic seed
		double curr = 0;
		for (size_t i = 0; i < nodes.size(); ++i) {
			curr += nodes[i].influence;
		}
		if ( max_influence < curr ) {
			best_topic = t;
			select_index = index;
			max_influence = curr;
		}
	}
	candidates = TopicAwareBase::topicVector[best_topic][select_index].seeds;

    timer.SetTimeEvent("end");
	// QueryPerformanceCounter(&tend);
    double duration = timer.TimeSpan("start", "end");
    // (double)(tend.QuadPart - start.QuadPart) / freq.QuadPart;

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
	cout << "Choose seeds from topic: " << (best_topic+1) << endl;

	FILE* out;
	fopen_s(&out, timerfile.c_str(), "w");
	fprintf(out, "%g", duration);
	fclose(out);

	// write to file

	size_t max_len = (candidates.size() <= (size_t)top) ? candidates.size() : (size_t)top;
	d.resize(max_len);
	list.resize(max_len);
	for (size_t i = 0; i < max_len; ++i) {
		list[i] = candidates[i].num;
		d[i] = candidates[i].influence;
	}
	top = max_len;

	WriteToFile(file, gf);
}

void TopSelection::Build(Graph& gf, int set_size, std::vector<double>& distribution)
{
	n = gf.GetN();
	top = set_size;

	SelectTop(distribution, gf);

}

#ifndef GREEDY_ONLINE_H
#define GREEDY_ONLINE_H

#include <set>
#include <vector>
#include <queue>

#include "common.h"
#include "graph.h"
#include "cascade.h"
#include "algo_base.h"



/// Greedy with online bound
class GreedyOnline
	: public AlgoBase
{
private:
	class EstNode{
	public:
		int id;
		double d;
		int last_update;

	public:
		EstNode(int id = -1, double d = -1, int last_update = -1)
			: id(id), d(d), last_update(last_update) {}
	};

	/// in descending order
	struct EstGreater {
	public:
		bool operator() (const EstNode* a, const EstNode* b) {
			return a->d < b->d;
		}
	};

	typedef std::priority_queue<EstNode*, std::vector<EstNode*>, EstGreater> EstHeap;
	

protected:
	std::string file;
	std::string onlinebound_file;
	bool isOnlineEstimated; // switch to estimate online bound
	std::vector< std::vector<int> >  onlineSeeds;
	std::vector< std::vector<double> > onlineD;

public:
	GreedyOnline();

	void Build(Graph& gf, int k, ICascade& cascade);
	void BuildRanking(Graph& gf, int k, ICascade& cascade);
	void BuildFromFile(Graph& gf, const char*);

protected:
	bool SortBestToTop(EstHeap& heap, int round, int* set, double old,
		Graph& gf, ICascade& cascade, int& ccc);
};

#endif ///:~GREEDY_ONLINE_H


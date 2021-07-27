#ifndef top_selection_h__
#define top_selection_h__

#include <vector>
#include <string>
#include "topic_aware.h"
#include "graph.h"
#include "algo_base.h"

/// Topic-Aware:
/// Top Selection Algorithm
class TopSelection
	: public AlgoBase
{
private:
	std::string file;
	std::string timerfile;
	void SelectTop(std::vector<double>& distribution, Graph& gf);
	

public:
	TopSelection();

	void Build(Graph& gf, int set_size, std::vector<double>& distribution);
};

#endif // top_selection_h__

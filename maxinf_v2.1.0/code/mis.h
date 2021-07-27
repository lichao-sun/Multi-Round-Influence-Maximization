#ifndef MIS_H
#define MIS_H

#include <vector>
#include <string>
#include "topic_aware.h"
#include "graph.h"
#include "algo_base.h"

/// Topic-Aware:
/// Marginal Influence Sort (MIS) Algorithm
class MIS
	: public AlgoBase
{
protected:
	std::string file;
	std::string timerfile;
	void EstimateSeeds(std::vector<double>& distribution, Graph& gf);
	int FindSeedByNum(std::vector<MINode>& candidates, int num);

	struct MINodeGreater
	{
		bool operator () (const MINode& left, const MINode& right) const;
	};

public:
	MIS();
	void Build(Graph& gf, int set_size, std::vector<double>& distribution);
};


#endif ///:~ MIS_H


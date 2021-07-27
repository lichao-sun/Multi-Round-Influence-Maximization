#ifndef mi_command_line_h__
#define mi_command_line_h__

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <set>
#include <fstream>
#include <sstream>
#include <random>

#include "common.h"
#include "graph.h"
#include "graph_stat.h"
#include "random_pick.h"
#include "degree.h"
#include "greedy.h"
#include "rim_greedy.h"
#include "rim_crgreedy.h"
#include "rim_agreedy.h"
#include "greedy_online.h"
#include "degreediscount_ic.h"
#include "weighted_degree.h"
#include "independ_cascade.h"

#include "SPM_gc.h"
#include "SP1M_gc.h"
#include "pmia.h"
#include "pagerank.h"

#include "mia.h"

#include "topic_aware.h"
#include "top_selection.h"
#include "cgreedy.h"
#include "mis.h"
#include "event_timer.h"
#include "simulate.h"
#include "rim_simulate.h"

#include "rr_infl.h"
#include "cascade.h"
#include "contcascade.h"
#include "general_cascade.h"


/// Command line for a set of max_influence algorithms
class MICommandLine
{
public:
	int Main(int argc, char* argv[]);
	int Main(int argc, std::vector<std::string>& argv);
	std::string Help();
	void BuildRanking(int argc, std::vector<std::string>& argv);
	void TopicAwareCGreedy(int argc, std::vector<std::string>& argv);
	void TopicAwareCGreedyFromList(int argc, std::vector<std::string>& argv);
	void TopicAwareMIS(int argc, std::vector<std::string>& argv);
	void TopicAwareTopSelection(int argc, std::vector<std::string>& argv);
	void TestSeeds(int argc, std::vector<std::string>& argv);
	void RIM_TestSeeds(int argc, std::vector<std::string>& argv);
	void GraphStat(int argc, std::vector<std::string>& argv);
	void BaselineAlg(int argc, std::vector<std::string>& argv);
	void GreedyAlg(int argc, std::vector<std::string>& argv);
	void RIMGreedyAlg(int argc, std::vector<std::string>& argv);
	void RIMCRGreedyAlg(int argc, std::vector<std::string>& argv);
	void RIMAGreedyAlg(int argc, std::vector<std::string>& argv);
	void GreedyOnlineAlg(int argc, std::vector<std::string>& argv);
	void PMIAAlg(int argc, std::vector<std::string>& argv);
	void MIAAlg(int argc, std::vector<std::string>& argv);
	void RRAlg(int argc, std::vector<std::string>& argv);
	void ContinousICAlg(int argc, std::vector<std::string>& argv);
	void PagerankGraphGen(int argc, std::vector<std::string>& argv);
};


#endif // mi_command_line_h__

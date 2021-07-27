#ifndef topic_aware_h__
#define topic_aware_h__

#include <string>
#include <vector>
#include <map>
#include <set>
#include "graph.h"
#include "common.h"


/// Marginal influence node
struct MINode
{
	int num;
	double influence;
};

/// Topic Aware Seeds with its lambda
struct TASeeds
{
	double vlambda;
	std::vector<MINode> seeds;
};

/// TopicAwareBase used for topic-aware influence maximization algorithms
class TopicAwareBase {
private:
	struct TASeedLambdaSmaller
	{
		bool operator() (const TASeeds& m1, const TASeeds& m2)
		{
			return (m1.vlambda < m2.vlambda);
		}
	};

	/// Read seeds from file <br/>
	/// [Format] <br/>
	/// <seeds_num> <br/>
	/// <seed1>  <maringal_influence1> <br/>
	/// <seed2>  <maringal_influence2> <br/>
	/// ... <br/>
	static void ReadSeeds(std::string file, TASeeds& ta, IGraph& gf);
	
public:
	typedef std::map< int, std::vector<TASeeds> > TASeedsMap;

	static int ntopics;
	static int nlandmarks;
	static TASeedsMap topicMap;
	static std::vector< std::vector<TASeeds> > topicVector; // another representation started with 0

	/// Round down seeds to the closest vlambda, and return its index
	/// (If only it is only one vlambda, then round to this vlambda, especially for 1.0.)
	static int RoundSeedsToIndex(std::vector<TASeeds>& seeds, double vlambda);
	
	/// Read seeds list: (topic's index starts with 1) <br/>
	/// <br/>
	/// [Format] <br/>
	/// <number_topic> <number_landmarks>  <br/>
	/// <topic1> <vlambda1>  <seed_file1>  <br/>
	/// <topic2> <vlambda2>  <seed_file2>  <br/>
	/// ... <br/>
	/// [Example] <br/>
	/// 2 2 <br/>
	/// 1   0.3  topic1_0.3-seeds.txt <br/>
	/// 1   1.0  topic1_1.0-seeds.txt <br/>
	/// 2   0.1  topic1_0.1-seeds.txt <br/>
	/// 2   1.0  topic1_0.1-seeds.txt <br/>
	static void ReadTopicFile(std::string seeds_list_file, IGraph& gf);

	static std::set<int> GetAllSeedsIndices();
};

#endif // topic_aware_h__

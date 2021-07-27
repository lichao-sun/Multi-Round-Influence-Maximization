#ifndef algo_base_h__
#define algo_base_h__

#include <vector>
#include <set>
#include <iostream>
#include "graph.h"

class AlgoBase
{
protected:
	/// node num
	int n;
	/// seed num
	int top;
	/// rounds num
	int T;
	/// root num;
	int R;
	/// seeds
	std::vector<int> list;
	/// all nodes before i-th round
	std::set<int> act_nodes;
	/// seeds list in first i-th round
	bool act_nodes_arr[30000];
	/// seeds list in first i-th round
	std::vector<std::vector<int>> seeds_list;
	/// influence spread
	std::vector<double> d; 

public:
	AlgoBase();
	virtual ~AlgoBase();

public:
	/// Get seed by index
	virtual int GetSeed(int i);
	/// Get seed list
	virtual std::vector<int>& GetSeedList();

protected:
	/// Write seeds and their influence to file
	virtual void WriteToFile(const std::string& filename, IGraph& gf);
	/// Read seeds and their influence from file
	virtual void ReadFromFile(const std::string& filename, IGraph& gf);
};


#endif ///:~ algo_base_h__
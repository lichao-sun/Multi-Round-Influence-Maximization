#ifndef RANDOM_PICK_H
#define RANDOM_PICK_H

#include <vector>
#include <string>
#include "common.h"
#include "graph.h"
#include "algo_base.h"

/// Randomly pick k nodes as seeds
class RandomPick
	: public AlgoBase
{
private:
	std::string file;

public:
	RandomPick();

	void Build(Graph& gf, int k = SET_SIZE);
	void BuildFromFile(Graph& gf);
};

#endif


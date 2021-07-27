#ifndef __SIMULATE__H__
#define __SIMULATE__H__

#include <iostream>
#include <functional>
#include <vector>
#include <fstream>
#include "common.h"
#include "graph.h"
#include "cascade.h"


/// simulator
class Simu
{
public:
	std::string file;
	std::vector<int> seeds;
	int simuSize;
	int simuIterNum;

	Simu();
	Simu(std::vector<int>& seeds, 
			std::string filename = "", 
			int simulateSize = SET_SIZE, 
			int simuIterNum = NUM_ITER);
	Simu(int* seeds, 
			std::string filename = "",
			int simulateSize = SET_SIZE, 
			int simuIterNum = NUM_ITER);
	Simu(int(*GetNode)(int i), 
			std::string filename = "", 
			int simulateSize = SET_SIZE, 
			int simuIterNum = NUM_ITER);

public:
	double toSimulate(ICascade& cascade);
	double toSimulateOnce(ICascade& cacade);

protected:
	inline bool isWriteFile()
	{
		return !file.empty();
	}
};

#endif // __SIMULATE__H__

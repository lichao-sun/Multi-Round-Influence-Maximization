#ifndef __RIM_Simu__H__
#define __RIM_Simu__H__

#include <iostream>
#include <functional>
#include <vector>
#include <fstream>
#include "common.h"
#include "graph.h"
#include "cascade.h"


/// simulator
class RIM_Simu
{
public:
	std::string file;
	std::vector<std::vector<int>> seeds_list;
	int simuSize;
	int simuIterNum;

	RIM_Simu();
	RIM_Simu(std::vector<std::vector<int>>& seeds_list,
			std::string filename = "", 
			int simulateSize = SET_SIZE, 
			int simuIterNum = NUM_ITER);
	RIM_Simu(int* seeds,
			std::string filename = "",
			int simulateSize = SET_SIZE, 
			int simuIterNum = NUM_ITER);
	RIM_Simu(int(*GetNode)(int i),
			std::string filename = "", 
			int simulateSize = SET_SIZE, 
			int simuIterNum = NUM_ITER);

public:
	double toSimulate(ICascade& cascade, int T);
	double toSimulateOnce(ICascade& cacade);

protected:
	inline bool isWriteFile()
	{
		return !file.empty();
	}
};

#endif // __SIMULATE__H__

#include "rim_simulate.h"

RIM_Simu::RIM_Simu()
{
	file = "";
	simuSize = 0;
	simuIterNum = NUM_ITER;
}

RIM_Simu::RIM_Simu(std::vector<std::vector<int>>& seeds_list,
			std::string filename /*= ""*/, 
			int simulateSize /*= SET_SIZE*/, 
			int simuIterNum /*= NUM_ITER*/)
{
	this->file = filename;
	this->seeds_list = seeds_list;
	this->simuSize = ((int)seeds_list.size() < simulateSize) ? (int)seeds_list.size() : simulateSize;
	this->simuIterNum = simuIterNum;
}

//RIM_Simu::RIM_Simu(int* seeds,
//			std::string filename /*= ""*/, 
//			int simulateSize /*= SET_SIZE*/, 
//			int simuIterNum /*= NUM_ITER*/)
//{
//	this->file = filename;
//	for (int i = 0; i < simulateSize; i++) {
//		this->seeds_list.push_back(seeds[i]);
//	}
//	this->simuSize = simulateSize;
//	this->simuIterNum = simuIterNum;
//}
//
//
//RIM_Simu::RIM_Simu(int(*GetNode)(int i),
//			std::string filename /*= ""*/, 
//			int simulateSize /*= SET_SIZE*/, 
//			int simuIterNum /*= NUM_ITER*/)
//{
//	this->file = filename;
//	for (int i = 0; i < simulateSize; i++) {
//		seeds.push_back(GetNode(i));
//	}
//	this->simuSize = simulateSize;
//	this->simuIterNum = simuIterNum;
//}

double RIM_Simu::toSimulate(ICascade& cascade, int T)
{
	/*std::ofstream out;
	if (isWriteFile()) {
		out.open(file.c_str());
	}*/
	std::vector<double> result;
	for (const auto &row : seeds_list)
	{
		for (const auto &s : row) std::cout << s << '\t';
		std::cout << std::endl;
	}

	std::set<int> act_nodes;
	std::vector<std::vector<int>> tmp_seeds_list;
	double spread = 0.0;
	int rounds = T;
	for (int i = 0; i < NUM_ITER; i++){
		for (int round = 0; round < rounds; round++){
			tmp_seeds_list.push_back(seeds_list[round]);
			cascade.ActRun(seeds_list[round].size(), act_nodes, tmp_seeds_list);
			tmp_seeds_list.clear();
		}
		//std::cout << i << ' ' << act_nodes.size() << std::endl;
		spread = spread + act_nodes.size();
		result.push_back(act_nodes.size());
		act_nodes.clear();
	}
	{
		std::ofstream myfile;
		std::stringstream ss;
		ss << "rim_tna_" << rounds << ".csv";
		std::string fullName = ss.str();

		myfile.open(fullName);
		for (const auto &s : result) myfile << s << ',';
		myfile << '\n';

		myfile.close();
	}
	std::cout << "ESpread:" << double(spread) / NUM_ITER;
	return spread;
}



double RIM_Simu::toSimulateOnce(ICascade& cacade)
{
	int* set = new int[simuSize];
	for (int t = 0; t < simuSize; t++)
	{
		set[t] = seeds_list[1][t];
	}
	double spread = cacade.Run(simuIterNum, simuSize, set);
	SAFE_DELETE_ARRAY(set);
	return spread;
}
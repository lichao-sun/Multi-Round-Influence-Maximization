#include <cstdio>
#include <string>
#include <cstdlib>
#include <vector>
#include <set>
#include <sstream>

#include "rim_agreedy.h"
#include "graph.h"
#include "common.h"
#include "cascade.h"
#include "event_timer.h"

//#include <set>
using namespace std;


RIM_AGreedy::RIM_AGreedy()
{
	file = "rim_agreedy.txt";
}

void RIM_AGreedy::Build(IGraph& gf, int num, int rounds, ICascade& cascade)
{
	int flag = 1;
	list.clear();
	act_nodes.clear();
	seeds_list.clear();

	n = gf.GetN();
	top = num;
	T = rounds;
	d.resize(num, 0);
	list.resize(num, 0);
	memset(act_nodes_arr, 0, sizeof(bool)*n);
	vector<vector<int>> empty_list;
	bool *act_nodes_arr = new bool[n];
	vector<int> cur_act_nodes;
	vector<double> res_time;
	vector<vector<set<int>>> previous_act_nodes;

	for (int round = 0; round <= T; round++){
		EventTimer timer;
		timer.SetTimeEvent("start");

		cout << "round: " << round << ' ' << act_nodes.size() << ' ' << seeds_list.size() << endl;
		if (round > 0){
			cascade.ActRun(top, act_nodes, seeds_list);
		}
		for (auto const& value : act_nodes) {
			act_nodes_arr[value] = 1;
		}

		cur_act_nodes.push_back(act_nodes.size());

		cout << "round: " << round << ' ' << act_nodes.size() << ' ' << seeds_list.size() << endl;
		if (round == T)
			continue;

		bool *used = new bool[n];
		memset(used, 0, sizeof(bool)*n);
		int* set = new int[num];

		double old = 0.0;

		double *improve = new double[n];
		int *lastupdate = new int[n];
		int *heap = new int[n];
		for (int i = 0; i<n; i++)
		{
			heap[i] = i;
			lastupdate[i] = -1;
			improve[i] = (double)(n + 1);//initialize largest
		}

		for (int i = 0; i<top; i++)
		{
			int ccc = 0;
			//printf("%d\n", i);
			while (lastupdate[heap[0]] != i)
			{
				/*printf("%d %d %d\n", i, heap[0], ccc);*/
				ccc++;
				lastupdate[heap[0]] = i;
				set[i] = heap[0];
				//printf("GreedyGC_SPM %d %d\n", heap[0], improve[heap[0]]);

				improve[heap[0]] = cascade.RRun(NUM_ITER, i + 1, set, act_nodes_arr) - old;
				char tmpfilename[200];
				sprintf_s(tmpfilename, "tmp/%02d%05d.txt", i, heap[0]);
				//FILE *tmpfile;
				//fopen_s(&tmpfile, tmpfilename,"w");
				//fprintf(tmpfile, "%g\n", improve[heap[0]]); 
				//fclose(tmpfile);

				int x = 0;
				while (x * 2 + 2 <= n - i)
				{
					int newx = x * 2 + 1;
					if ((newx + 1<n - i) && (improve[heap[newx]]<improve[heap[newx + 1]]))
						newx++;
					if (improve[heap[x]]<improve[heap[newx]])
					{
						int t = heap[x];
						heap[x] = heap[newx];
						heap[newx] = t;
						x = newx;
					}
					else
						break;
				}
			}

			used[heap[0]] = true;
			set[i] = heap[0];
			list[i] = heap[0];
			act_nodes.insert(heap[0]);
			d[i] = improve[heap[0]];
			old += d[i];
			cout << list[i] << ' ' << improve[heap[0]] << endl;
			//char bakname[200];
			//sprintf(bakname, "greedychoice%02d.txt", i+1);
			//FILE *bak = fopen(bakname, "w");
			//fprintf(bak, "%6d\t%d\t%g\n", i+1, heap[0], improve[heap[0]]);
			//fclose(bak);
			//printf("%d\t%g\n", i+1, improve[131]);

			heap[0] = heap[n - i - 1];
			int x = 0;
			while (x * 2 + 2 <= n - i)//bug should-1
			{
				int newx = x * 2 + 1;
				if ((newx + 1<n - i) && (improve[heap[newx]]<improve[heap[newx + 1]]))	//bug should-1
					newx++;
				if (improve[heap[x]]<improve[heap[newx]])
				{
					int t = heap[x];
					heap[x] = heap[newx];
					heap[newx] = t;
					x = newx;
				}
				else
					break;
			}

		}

		timer.SetTimeEvent("end");
		res_time.push_back(timer.TimeSpan("start", "end"));
		cout << "time: " << timer.TimeSpan("start", "end") << endl;
		WriteToFile(file, gf);
		seeds_list.push_back(list);
		cout << ' ' << seeds_list.size() << endl;
		SAFE_DELETE_ARRAY(set);
		SAFE_DELETE_ARRAY(heap);
		SAFE_DELETE_ARRAY(lastupdate);
		SAFE_DELETE_ARRAY(improve);
		SAFE_DELETE_ARRAY(used);
	}

	for (const auto &s : cur_act_nodes) std::cout << s << '\t';
	std::cout << std::endl;

	for (const auto &s : res_time) std::cout << s << '\t';
	std::cout << std::endl;

	for (const auto &row : seeds_list)
	{
		for (const auto &s : row) std::cout << s << '\t';
		std::cout << std::endl;
	}

	{
		ofstream myfile;
		std::stringstream ss;
		ss << "rim_agreedy" << ".csv";
		std::string fullName = ss.str();

		myfile.open(fullName);
		for (const auto &s : cur_act_nodes) myfile << s << ',';
		myfile << '\n';

		for (const auto &s : res_time) myfile << s << ',';
		myfile << '\n';

		for (const auto &row : seeds_list)
		{
			for (const auto &s : row) myfile << s << ',';
			myfile << '\n';
		}
		myfile.close();
	}
}


void RIM_AGreedy::BuildFromFile(IGraph& gf, const char* name)
{
	ReadFromFile(name, gf);
}

#include <cstdio>
#include <string>
#include <cstdlib>
#include <vector>
#include <set>
#include <random>
#include "mi_random.h"
#include "rim_greedy.h"
#include "graph.h"
#include "common.h"
#include "cascade.h"
#include "event_timer.h"

//#include <set>
using namespace std;


RIM_Greedy::RIM_Greedy() 
{
	file = "rim_greedy.txt";
}

void RIM_Greedy::Build(IGraph& gf, int num, int rounds, ICascade& cascade)
{
	n = gf.GetN();
	top = num;
	T = rounds;
	d.resize(num, 0);
	list.resize(num, 0);
	vector<double> res_time;

	vector<vector<set<int>>> previous_act_nodes;
	MIRandom random;
	double total_spread = 0;

	for (int round = 0; round < T; round++){
		EventTimer timer;
		timer.SetTimeEvent("start");
		cout << "round: " << round << ' ' << act_nodes.size() << ' ' << seeds_list.size() << endl;
		if (round > 0){
			vector<set<int>> tmp_vector_set;
			for (int iter = 0; iter < 100; iter++){
				act_nodes.clear();
				vector<vector<int>> cur_seeds_list;
				cur_seeds_list.push_back(seeds_list[round - 1]);
				cascade.ActRun(top, act_nodes, cur_seeds_list);
				tmp_vector_set.push_back(act_nodes);
			}
			previous_act_nodes.push_back(tmp_vector_set);
			tmp_vector_set.clear();
			act_nodes.clear();
		}
		bool** previous_act_nodes_arr = new bool*[NUM_ITER];
		for (int i = 0; i < NUM_ITER; ++i)
			previous_act_nodes_arr[i] = new bool[n];

		if (round > 0){
			for (int iter = 0; iter < NUM_ITER; iter++){
				for (int cur_round = 0; cur_round < round; cur_round++){
					int id = random.RandInt(0, 99);
					for (auto const& value : previous_act_nodes[round - 1][id]) {
						previous_act_nodes_arr[iter][value] = 1;
					}
				}
			}
		}

		cout << "round: " << round << ' ' << act_nodes.size() << ' ' << seeds_list.size() << endl;

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
				improve[heap[0]] = cascade.ERRun(NUM_ITER, i + 1, set, previous_act_nodes_arr) - old;
				//improve[heap[0]] = cascade.Run(NUM_ITER, i + 1, set) - old;
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
			total_spread = total_spread + improve[heap[0]];
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

		SAFE_DELETE_ARRAY(set);
		SAFE_DELETE_ARRAY(heap);
		SAFE_DELETE_ARRAY(lastupdate);
		SAFE_DELETE_ARRAY(improve);
		SAFE_DELETE_ARRAY(used);
	}

	cout << "total spread: " << total_spread << endl;

	for (const auto &s : res_time) std::cout << s << '\t';
	std::cout << std::endl;

	for (const auto &row : seeds_list)
	{
		for (const auto &s : row) std::cout << s << '\t';
		std::cout << std::endl;
	}

	{
		ofstream myfile;
		myfile.open("rim_greedy.csv");

		for (const auto &s : res_time) myfile << s << ',';
		myfile << '\n';

		for (const auto &row : seeds_list)
		{
			for (const auto &s : row) myfile << s << ',';
			myfile << '\n';
		}
		myfile.close();
	}

	{
		ofstream myfile;
		myfile.open("rim_greedy.txt");

		for (const auto &row : seeds_list)
		{
			for (const auto &s : row) myfile << s << '\t';
			myfile << '\n';
		}
		myfile.close();
	}

}


void RIM_Greedy::BuildFromFile(IGraph& gf, const char* name)
{
	ReadFromFile(name, gf);
}

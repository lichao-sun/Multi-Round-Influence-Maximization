#include <cstdio>
#include <string>
#include <cstdlib>
#include <vector>
#include <set>
#include <random>
#include "mi_random.h"
#include "rim_crgreedy.h"
#include "graph.h"
#include "common.h"
#include "cascade.h"
#include "event_timer.h"

//#include <set>
using namespace std;


RIM_CRGreedy::RIM_CRGreedy()
{
	file = "rim_greedy.txt";
}

void RIM_CRGreedy::Build(IGraph& gf, int num, int rounds, ICascade& cascade)
{
	n = gf.GetN();
	top = num;
	T = rounds;
	d.resize(num*T, 0);
	list.resize(num*T, 0);
	vector<double> res_time;

	vector<vector<set<int>>> previous_act_nodes;
	MIRandom random;
	double total_spread = 0;

	bool *used = new bool[n*T];
	memset(used, 0, sizeof(bool)*n*T);
	int* set = new int[num*T];

	double old = 0.0;

	double *improve = new double[n*T];
	int *lastupdate = new int[n*T];
	int *heap = new int[n*T];
	bool*resRounds = new bool[T];

	for (int j = 0; j < T; j++){
		resRounds[j] = 0;
	}

	seeds_list.clear();
	for (int i = 0; i < T; ++i) {
		seeds_list.push_back(vector<int>());
	}
	for (int i = 0; i<n*T; i++)
	{
		heap[i] = i;
		lastupdate[i] = -1;
		improve[i] = (double)(n + 1);//initialize largest
	}

	std::vector<std::set<int>> tmp_vector_set;
	vector<vector<int>> cur_seeds_list;
	vector<int> last_seed;
	bool** previous_act_nodes_arr = new bool*[NUM_ITER];
	for (int i = 0; i < NUM_ITER; ++i)
		previous_act_nodes_arr[i] = new bool[n];
	
	EventTimer timer;
	timer.SetTimeEvent("start");
	for (int i = 0; i<top*T; i++)
	{	
		//timer.SetTimeEvent("start");
		//if (i > 0){
		//	tmp_vector_set.clear();
		//	for (int iter = 0; iter < NUM_ITER; iter++){
		//		act_nodes.clear();
		//		for (int j = 0; j < T; j++){
		//			cur_seeds_list.clear();
		//			cur_seeds_list.push_back(seeds_list[j]);
		//			cascade.ActRun(seeds_list[j].size(), act_nodes, cur_seeds_list);
		//		}
		//		tmp_vector_set.push_back(act_nodes);
		//		/*act_nodes.clear();
		//		cur_seeds_list.clear();
		//		last_seed.clear();
		//		last_seed.push_back(set[i - 1] % n);
		//		cur_seeds_list.push_back(last_seed);
		//		cascade.ActRun(1, act_nodes, cur_seeds_list);
		//		tmp_vector_set.push_back(act_nodes);*/
		//	}
		//	//previous_act_nodes.push_back(tmp_vector_set);
		//	act_nodes.clear();
		//}

		//if (i > 0){
		//	for (int iter = 0; iter < NUM_ITER; iter++){
		//		//for (int cur_round = 0; cur_round < T; cur_round++){
		//		for (auto const& value : tmp_vector_set[iter]) {
		//			previous_act_nodes_arr[iter][value] = 1;
		//		}
		//		//}
		//	}
		//	tmp_vector_set.clear();
		//}

		for (int j = 0; j < T; j++){
			cout << j << ' ' << seeds_list[j].size() << ' ' << resRounds[j]  << endl;
			if (seeds_list[j].size() >= top){
				resRounds[j] = 1;
				cout << j << ' ' << seeds_list[j].size() << endl;
			}
		}

		int ccc = 0;
		//printf("%d\n", i);
		while (lastupdate[heap[0]] != i)
		{
			//cout << i << ' ' << heap[0] << endl;
			/*printf("%d %d %d\n", i, heap[0], ccc);*/
			ccc++;
			lastupdate[heap[0]] = i;
			set[i] = heap[0];
			if (resRounds[int(set[i] / n)])
				improve[heap[0]] = -10000;
			else{
				improve[heap[0]] = cascade.CRRun(NUM_ITER, i+1, set, previous_act_nodes_arr, T) - old;
				//improve[heap[0]] = cascade.Run(NUM_ITER, i + 1, set) - old;
			}
			if (heap[0] % n == 28)
				cout << heap[0] / n << ' ' << improve[heap[0]] << endl;
			//improve[heap[0]] = cascade.Run(NUM_ITER, i + 1, set) - old;
			char tmpfilename[200];
			sprintf_s(tmpfilename, "tmp/%02d%05d.txt", i, heap[0]);
			//FILE *tmpfile;
			//fopen_s(&tmpfile, tmpfilename,"w");
			//fprintf(tmpfile, "%g\n", improve[heap[0]]); 
			//fclose(tmpfile);

			int x = 0;
			while (x * 2 + 2 <= n * T - i)
			{
				int newx = x * 2 + 1;
				if ((newx + 1<n * T - i) && (improve[heap[newx]]<improve[heap[newx + 1]]))
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
		std::cout << list[i] << ' ' << int(list[i] / n) << ' ' << list[i] % n << ' ' << improve[heap[0]] << endl;
		seeds_list[int(list[i] / n)].push_back(list[i] % n);
		total_spread = total_spread + improve[heap[0]];
		//char bakname[200];
		//sprintf(bakname, "greedychoice%02d.txt", i+1);
		//FILE *bak = fopen(bakname, "w");
		//fprintf(bak, "%6d\t%d\t%g\n", i+1, heap[0], improve[heap[0]]);
		//fclose(bak);
		//printf("%d\t%g\n", i+1, improve[131]);

		heap[0] = heap[n - i - 1];
		int x = 0;
		while (x * 2 + 2 <= n * T - i)//bug should-1
		{
			int newx = x * 2 + 1;
			if ((newx + 1<n * T - i) && (improve[heap[newx]]<improve[heap[newx + 1]]))	//bug should-1
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
	//res_time.push_back(timer.TimeSpan("start", "end"));
	std::cout << "time: " << timer.TimeSpan("start", "end") << endl;
	//WriteToFile(file, gf);

	SAFE_DELETE_ARRAY(set);
	SAFE_DELETE_ARRAY(heap);
	SAFE_DELETE_ARRAY(lastupdate);
	SAFE_DELETE_ARRAY(improve);
	SAFE_DELETE_ARRAY(used);
	

	cout << "total spread: " << total_spread << endl;

	/*for (const auto &s : res_time) std::cout << s << '\t';
	std::cout << std::endl;*/

	for (const auto &row : seeds_list)
	{
		for (const auto &s : row) std::cout << s << '\t';
		std::cout << std::endl;
	}

	{
		ofstream myfile;
		myfile.open("rim_crgreedy.txt");

		/*for (const auto &s : res_time) myfile << s << '\t';
		myfile << '\n';*/

		for (const auto &row : seeds_list)
		{
			for (const auto &s : row) myfile << s << '\t';
			myfile << '\n';
		}
		myfile.close();
	}

}


void RIM_CRGreedy::BuildFromFile(IGraph& gf, const char* name)
{
	ReadFromFile(name, gf);
}

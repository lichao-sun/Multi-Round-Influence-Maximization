#ifndef GENERAL_CASCADE_H
#define GENERAL_CASCADE_H


#include "cascade.h"
#include <vector>

/// Template class that implements general cascade diffusion
template<class TGraph>
class GeneralCascadeT:
	public CascadeT<TGraph>
{
public:
	int nthreads;

public:
	GeneralCascadeT() : nthreads(1) {}

	virtual void Build(TGraph& gf)
	{
		_Build(gf);
	}

	double Run(int num_iter, int size, int set[])
	{
		InitializeConcurrent();
		return _Run(num_iter, size, set);
	}

	double CRRun(int num_iter, int size, int set[], bool** previous_act_nodes_arr, int T)
	{
		InitializeConcurrent();
		return _CRRun(num_iter, size, set, T);
	}

	double RRun(int num_iter, int size, int set[], bool act_nodes_arr[])
	{
		return _RRun(num_iter, size, set, act_nodes_arr);
	}

	double ERRun(int num_iter, int size, int set[], bool** previous_act_nodes_arr)
	{
		return _ERRun(num_iter, size, set, previous_act_nodes_arr);
	}

	double ActRun(int top, std::set<int>& act_nodes, std::vector<std::vector<int>>& seeds_list)
	{
		return _ActRun(top, act_nodes, seeds_list);
	}

protected:
	void InitializeConcurrent() 
	{
		if (IsConcurrent())
		{
#ifdef MI_USE_OMP
			/////////////////////////////////////////////
			// run concurrently
			const double DYNAMIC_RATIO = 0.25;
			omp_set_num_threads(nthreads);
			int dynamicThreads = (int)(nthreads * DYNAMIC_RATIO);
			omp_set_dynamic(dynamicThreads);

			std::cout << "== Turn on omp optimization: == " << std::endl;
			std::cout << "#Max Threads = " << omp_get_max_threads() << "\t#Dynamic Threads = " << omp_get_dynamic() << std::endl;
#else
			std::cout << "== omp is not supported or enabled == " << std::endl;
#endif
		}
	}

	inline bool IsConcurrent()
	{
		return (nthreads > 1);
	}

	double _ActRun(int top, std::set<int>& act_nodes, std::vector<std::vector<int>>& seeds_list){
		MI_STATIC_ASSERT(has_mem_GetNeighborCount<graph_type>::value, "TGraph should has member GetNeighborCount");
		MI_STATIC_ASSERT(has_mem_GetEdge<graph_type>::value, "TGraph should has member GetEdge");
		MI_STATIC_ASSERT(has_mem_edgeForm<graph_type>::value, "graph_type should has member edgeForm");

		if (gf == NULL) {
			throw NullPointerException("Please Build Graph first. (gf==NULL)");
		}

		int targetSize = top;
		int resultSize = 0;
		ProbTransfom trans(gf->edgeForm);

		int count = 0;
		int* list = new int[n];
		bool* active = new bool[n];
		int cur_round = seeds_list.size();

		memset(active, 0, sizeof(bool)*n);
		for (int i = 0; i < targetSize; i++)
		{
			list[i] = seeds_list[cur_round - 1][i];
			active[list[i]] = true;
			//cur_act_nodes.push_back(list[i]);
			act_nodes.insert(list[i]);
		}
		int h = 0;
		int t = targetSize;

		while (h < t)
		{
			int k = gf->GetNeighborCount(list[h]);
			//printf("%d \n",k);
			for (int i = 0; i < k; i++)
			{
				auto& e = gf->GetEdge(list[h], i);
				if (active[e.v]) continue;

				if (random.RandBernoulli(trans.Prob(e.w1)))
				{
					list[t] = e.v;
					active[e.v] = true;
					t++;
					act_nodes.insert(e.v);
				}
			}
			h++;
		}
		SAFE_DELETE_ARRAY(active);
		SAFE_DELETE_ARRAY(list);
	return 0.0;
	}

	double _RRun(int num_iter, int size, int set[], bool act_nodes_arr[])
	{
		MI_STATIC_ASSERT(has_mem_GetNeighborCount<graph_type>::value, "TGraph should has member GetNeighborCount");
		MI_STATIC_ASSERT(has_mem_GetEdge<graph_type>::value, "TGraph should has member GetEdge");
		MI_STATIC_ASSERT(has_mem_edgeForm<graph_type>::value, "graph_type should has member edgeForm");

		if (gf == NULL) {
			throw NullPointerException("Please Build Graph first. (gf==NULL)");
		}

		int targetSize = size;
		int resultSize = 0;
		ProbTransfom trans(gf->edgeForm);

			int count = 0;
			int* list = new int[n];
			bool* active = new bool[n];
			for (int it = 0; it < num_iter; it++)
			{
				memset(active, 0, sizeof(bool)*n);
				for (int i = 0; i < targetSize; i++)
				{
					list[i] = set[i];
					active[list[i]] = true;
				}
				resultSize += targetSize;

				int h = 0;
				int t = targetSize;

				while (h < t)
				{
					int k = gf->GetNeighborCount(list[h]);
					//printf("%d \n",k);
					for (int i = 0; i < k; i++)
					{
						auto& e = gf->GetEdge(list[h], i);
						if (active[e.v]) continue;

						if (random.RandBernoulli(trans.Prob(e.w1)))
						{
							list[t] = e.v;
							active[e.v] = true;
							t++;
							//if (tmp_act_nodes[e.v] == 1){
							if (act_nodes_arr[e.v] == 1){
								continue;
								count++;
							}	
							else{
								resultSize++;
							}
						}
					}
					h++;
				}
			}
			SAFE_DELETE_ARRAY(active);
			SAFE_DELETE_ARRAY(list);

		return (double)resultSize / (double)num_iter;
	}

	double _ERRun(int num_iter, int size, int set[], bool** previous_act_nodes_arr)
	{
		MI_STATIC_ASSERT(has_mem_GetNeighborCount<graph_type>::value, "TGraph should has member GetNeighborCount");
		MI_STATIC_ASSERT(has_mem_GetEdge<graph_type>::value, "TGraph should has member GetEdge");
		MI_STATIC_ASSERT(has_mem_edgeForm<graph_type>::value, "graph_type should has member edgeForm");

		if (gf == NULL) {
			throw NullPointerException("Please Build Graph first. (gf==NULL)");
		}

		int targetSize = size;
		int resultSize = 0;
		ProbTransfom trans(gf->edgeForm);

		int count = 0;
		int* list = new int[n];
		bool* active = new bool[n];
		for (int it = 0; it < num_iter; it++)
		{
			memset(active, 0, sizeof(bool)*n);
			for (int i = 0; i < targetSize; i++)
			{
				list[i] = set[i];
				active[list[i]] = true;
			}
			resultSize += targetSize;

			int h = 0;
			int t = targetSize;

			while (h < t)
			{
				int k = gf->GetNeighborCount(list[h]);
				//printf("%d \n",k);
				for (int i = 0; i < k; i++)
				{
					auto& e = gf->GetEdge(list[h], i);
					if (active[e.v]) continue;

					if (random.RandBernoulli(trans.Prob(e.w1)))
					{
						list[t] = e.v;
						active[e.v] = true;
						t++;
						if (previous_act_nodes_arr[it][e.v] == 1){
							continue;
							count++;
						}
						else{
							resultSize++;
						}
					}
				}
				h++;
			}
		}

		SAFE_DELETE_ARRAY(active);
		SAFE_DELETE_ARRAY(list);

		return (double)resultSize / (double)num_iter;
	}

	double _CRRun(int num_iter, int size, int set[], int T)
	{
		MI_STATIC_ASSERT(has_mem_GetNeighborCount<graph_type>::value, "TGraph should has member GetNeighborCount");
		MI_STATIC_ASSERT(has_mem_GetEdge<graph_type>::value, "TGraph should has member GetEdge");
		MI_STATIC_ASSERT(has_mem_edgeForm<graph_type>::value, "graph_type should has member edgeForm");

		if (gf == NULL) {
			throw NullPointerException("Please Build Graph first. (gf==NULL)");
		}


		int resultSize = 0;
		vector<vector<int>> tmp_seed_list;
		//cout << cur_round << endl;
		int* tmp_list = new int[n];
		int targetSize = 0;
		int cround = int(set[size - 1] / n);
		//vector<int> orderT;

		for (int i = 0; i < T; ++i) {
			tmp_seed_list.push_back(vector<int>());
			//if (cround != i)
				//orderT.push_back(i);
		}
		//orderT.push_back(cround);
		for (int j = 0; j < size; j++){
			tmp_seed_list[int(set[j] / n)].push_back(set[j] % n);
		}

		ProbTransfom trans(gf->edgeForm);


		// single thread
		//int T = = sizeof(previous_act_nodes_arr) / sizeof(previous_act_nodes_arr[0]);
		
		//bool* active_T = new bool[n];
		//vector<vector<int>> cr_seeds_list;

		//cout << "rounds: " << T << endl;
		//for (int i = 0; i < T; ++i) {
		//	cr_seeds_list.push_back(vector<int>());
		//}
		bool* active_T = new bool[n];
		int* list = new int[n];
		bool* active = new bool[n];
		memset(active, false, sizeof(bool)*n);

		for (int it = 0; it < num_iter; it++)
		{
			memset(active_T, false, sizeof(bool)*n);
			/*for (int j = 0; j < n; j++)
				active_T[j] = 0; */
			memset(active, false, sizeof(bool)*n);
			for (int cur_round = 0; cur_round < T; cur_round++){//orderT.size()
				//cout << orderT.size() << endl;
				/*for (int j = 0; j < n; j++){
					if (active[j]){
						active_T[j] = true;
					}
				}*/
				int targetSize = tmp_seed_list[cur_round].size();
				/*if (it == 1){
					cout << cur_round <<' '<< targetSize << endl;
				}*/
				memset(list, -1, sizeof(int)*n);
				memset(active, false, sizeof(bool)*n);
				for (int i = 0; i < targetSize; i++)
				{
					list[i] = tmp_seed_list[cur_round][i];
					/*if (it == 1 && cur_round == 1)
						cout << list[i] << ' ' << targetSize << endl;*/
					active[list[i]] = true;
					//cout << "set[i]:" << set[i] << ' ' << int(set[i] / n) << ' ' << set[i] % n << endl;
					//cr_seeds_list[int(set[i] / n)] = set[i] % n;
				}
				resultSize += targetSize;

				int h = 0;
				int t = targetSize;

				while (h < t)
				{
					int k = gf->GetNeighborCount(list[h]);
					//printf("%d \n",k);
					for (int i = 0; i < k; i++)
					{
						auto& e = gf->GetEdge(list[h], i);
						if (active[e.v]) continue;

						if (random.RandBernoulli(trans.Prob(e.w1)))
						{
							list[t] = e.v;
							active[e.v] = true;
							t++;
							if (active_T[e.v]){
								//cout << 1 << endl;
								continue;
								//count++;
							}
							else{
								active_T[e.v] = true;
								resultSize++;
							}
							//break;
						}
					}
					h++;
				}
			}
		}
		//cout << set[size-1] << ' ' <<resultSize << endl;
		SAFE_DELETE_ARRAY(active_T);
		SAFE_DELETE_ARRAY(active);
		SAFE_DELETE_ARRAY(list);
		return (double)resultSize / (double)num_iter;
	}

	//double _CRRun(int num_iter, int size, int set[], bool** previous_act_nodes_arr, int T)
	//{
	//	MI_STATIC_ASSERT(has_mem_GetNeighborCount<graph_type>::value, "TGraph should has member GetNeighborCount");
	//	MI_STATIC_ASSERT(has_mem_GetEdge<graph_type>::value, "TGraph should has member GetEdge");
	//	MI_STATIC_ASSERT(has_mem_edgeForm<graph_type>::value, "graph_type should has member edgeForm");

	//	if (gf == NULL) {
	//		throw NullPointerException("Please Build Graph first. (gf==NULL)");
	//	}

	//	
	//	int resultSize = 0;
	//	int cur_round = int(set[size - 1] / n);
	//	//cout << cur_round << endl;
	//	int* tmp_list = new int[n];
	//	int targetSize = 0;
	//	for (int j = 0; j < size; j++){
	//		if (int(set[j] / n) == cur_round){
	//			tmp_list[targetSize] = set[j] % n;
	//			targetSize++;
	//		}
	//	}

	//	ProbTransfom trans(gf->edgeForm);


	//	// single thread
	//	//int T = = sizeof(previous_act_nodes_arr) / sizeof(previous_act_nodes_arr[0]);
	//	int* list = new int[n];
	//	bool* active = new bool[n];
	//	//bool* active_T = new bool[n];
	//	//vector<vector<int>> cr_seeds_list;

	//	//cout << "rounds: " << T << endl;
	//	//for (int i = 0; i < T; ++i) {
	//	//	cr_seeds_list.push_back(vector<int>());
	//	//}

	//	for (int it = 0; it < num_iter; it++)
	//	{
	//		memset(active, 0, sizeof(bool)*n);
	//		for (int i = 0; i < targetSize; i++)
	//		{
	//			list[i] = tmp_list[i];
	//			active[list[i]] = true;
	//			//cout << "set[i]:" << set[i] << ' ' << int(set[i] / n) << ' ' << set[i] % n << endl;
	//			//cr_seeds_list[int(set[i] / n)] = set[i] % n;
	//		}
	//		resultSize += targetSize;

	//		int h = 0;
	//		int t = targetSize;

	//		while (h < t)
	//		{
	//			int k = gf->GetNeighborCount(list[h]);
	//			//printf("%d \n",k);
	//			for (int i = 0; i < k; i++)
	//			{
	//				auto& e = gf->GetEdge(list[h], i);
	//				if (active[e.v]) continue;

	//				if (random.RandBernoulli(trans.Prob(e.w1)))
	//				{
	//					list[t] = e.v;
	//					active[e.v] = true;
	//					t++;
	//					if (previous_act_nodes_arr[cur_round][e.v] == 1){
	//						continue;
	//						//count++;
	//					}
	//					else{
	//						resultSize++;
	//					}
	//					//break;
	//				}
	//			}
	//			h++;
	//		}
	//	}

	//	SAFE_DELETE_ARRAY(active);
	//	SAFE_DELETE_ARRAY(list);
	//	return (double)resultSize / (double)num_iter;
	//}


	double _Run(int num_iter, int size, int set[])
	{
		MI_STATIC_ASSERT(has_mem_GetNeighborCount<graph_type>::value, "TGraph should has member GetNeighborCount");
		MI_STATIC_ASSERT(has_mem_GetEdge<graph_type>::value, "TGraph should has member GetEdge");
		MI_STATIC_ASSERT(has_mem_edgeForm<graph_type>::value, "graph_type should has member edgeForm");

		if (gf == NULL) {
			throw NullPointerException("Please Build Graph first. (gf==NULL)");
		}

		int targetSize = size;
		int resultSize = 0;
		ProbTransfom trans(gf->edgeForm);

#ifdef MI_USE_OMP
		if (!IsConcurrent()) {
#endif
			// single thread
			int* list = new int[n];
			bool* active = new bool[n];
			for (int it = 0; it < num_iter; it++)
			{
				memset(active, 0, sizeof(bool)*n);
				for (int i = 0; i < targetSize; i++)
				{
					list[i] = set[i];
					active[list[i]] = true;
				}
				resultSize += targetSize;

				int h = 0;
				int t = targetSize;

				while (h < t)
				{
					int k = gf->GetNeighborCount(list[h]);
					//printf("%d \n",k);
					for (int i = 0; i < k; i++)
					{
						auto& e = gf->GetEdge(list[h], i);
						if (active[e.v]) continue;

						if (random.RandBernoulli(trans.Prob(e.w1)))
						{
							list[t] = e.v;
							active[e.v] = true;
							t++;
							resultSize++;
							//break;
						}
					}
					h++;
				}
			}

			SAFE_DELETE_ARRAY(active);
			SAFE_DELETE_ARRAY(list);

#ifdef MI_USE_OMP
		}
		else {
			// concurrent

#pragma omp parallel for
			for (int it = 0; it < num_iter; it++) {
				std::vector<int> list(n);
				std::vector<bool> active(n, false);
				for (int i = 0; i < targetSize; i++)
				{
					list[i] = set[i];
					active[list[i]] = true;
				}
				int curResult = targetSize;
				int h = 0;
				int t = targetSize;

				while (h < t)
				{
					int k = gf->GetNeighborCount(list[h]);
					//printf("%d \n",k);
					for (int i = 0; i < k; i++)
					{
						auto& e = gf->GetEdge(list[h], i);
						if (active[e.v]) continue;

						if (random.RandBernoulli(trans.Prob(e.w1)))
						{
							list[t] = e.v;
							active[e.v] = true;
							t++;
							curResult++;
						}
					}
					h++;
				}

#pragma omp atomic
				resultSize += curResult;
			}
		}
#endif

		return (double)resultSize / (double)num_iter;
	}
};


typedef GeneralCascadeT<Graph> GeneralCascade;


#endif

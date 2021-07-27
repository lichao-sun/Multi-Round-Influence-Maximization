///
/// rr_infl.h
/// 
/// In this file, we implement three algorithms as follows.
///
/// [1] Borgs, Christian, et al. "Maximizing social influence in nearly optimal time." 
///     Proceedings of the Twenty-Fifth Annual ACM-SIAM Symposium on Discrete Algorithms. SIAM, 2014.
/// [2] Tang, Youze, Xiaokui Xiao, and Yanchen Shi. "Influence maximization: Near-optimal time complexity meets practical efficiency." 
///     Proceedings of the 2014 ACM SIGMOD international conference on Management of data. ACM, 2014.
/// [3] Youze Tang, Yanchen Shi, and Xiaokui Xiao. "Influence maximization in near-linear time: a martingale approach".
///     Proceedings of the 2015 ACM SIGMOD International Conference on Management of Data. ACM, 2015.


#ifndef __rr_infl_h__
#define __rr_infl_h__


#include <set>
#include <vector>
#include "graph.h"
#include "common.h"
#include "reverse_general_cascade.h"
#include "algo_base.h"
#include "general_cascade.h"
#include "mi_random.h"


/// Base class for Reverse Influence Maximization 
class RRInflBase
	: public AlgoBase
{
public:
	typedef Graph graph_type;
	typedef ReverseGCascade cascade_type;

	RRInflBase() : m(0),
			isConcurrent(false) {
	}

	// for concurrent optimization: using omp
	bool isConcurrent; // turn on to use openmp

protected:
	int m;
	int* cRound;
	
	std::vector< RRVec > table;
	std::vector<int> targets;
	std::vector<bool> enables;
	//std::vector<int> root_arr;
	// degree of hyper-edges v, where e(u, v) in the hyper graph
	// source id --> degrees
	std::vector<int> degrees;
	std::vector<std::vector<RRVec >> tableVec;
	std::vector<std::vector<int>> degreesVec;
	std::vector<std::vector< std::vector<int> >> degreeRRIndicesVec;
	std::vector<std::set<int>> sourceSetVec; // all the source node ids
	std::vector< std::vector<int> > degreeRRIndices;
	std::set<int> sourceSet; // all the source node ids
	

	void InitializeConcurrent();
	
	void _AddRRSimulation(size_t num_iter, 
		cascade_type& cascade, 
		std::vector< RRVec >& refTable, 
		std::vector<int>& refTargets);
	void _AddRRSimulation(size_t num_iter,
		cascade_type& cascade, 
		std::vector< RRVec >& refTable,
		std::vector<int>& refTargets,
		std::vector<int>& refEdgeVisited);
	void _AdaAddRRSimulation(size_t num_iter,
		cascade_type& cascade,
		std::vector< RRVec >& refTable,
		std::vector<int>& refTargets);
	void _AdaAddRRSimulation(size_t num_iter,
		cascade_type& cascade,
		std::vector< RRVec >& refTable,
		std::vector<int>& refTargets,
		std::vector<int>& refEdgeVisited);
	void _WRIMAddRRSimulation(size_t num_iter,
		cascade_type& cascade,
		std::vector< RRVec >& refTable,
		std::vector<int>& refTargets,
		int* root_arr,
		int round);
	void _WRIMAddRRSimulation(size_t num_iter,
		cascade_type& cascade,
		std::vector< RRVec >& refTable,
		std::vector<int>& refTargets,
		std::vector<int>& refEdgeVisited,
		int* root_arr,
		int round);
	void _CRIMAddRRSimulation(size_t num_iter,
		cascade_type& cascade,
		std::vector<std::vector< RRVec >>& refTable,
		std::vector<int>& refTargets);
	void _CRIMAddRRSimulation(size_t num_iter,
		cascade_type& cascade,
		std::vector<std::vector< RRVec >>& refTable,
		std::vector<int>& refTargets,
		std::vector<int>& refEdgeVisited);
	double _RunGreedy(int seed_size, 
		std::vector<int>& outSeeds, 
		std::vector<double>& outMarginalCounts);
	double _WRIMRunGreedy(int seed_size,
		std::vector<int>& outSeeds,
		std::vector<double>& outMarginalCounts);
	double _CRIMRunGreedy(int seed_size,
		std::vector<std::vector<int>>& outSeedsVec,
		std::vector<double>& outMarginalCounts);
	double _RunRimGreedy(int seed_size,
		std::vector<int>& outSeeds,
		std::vector<double>& outMarginalCounts);
	void _RebuildRRIndices();
	void _CRIMRebuildRRIndices();
	void _AdaRebuildRRIndices();
	/*void _TableRemove(vector<int>& activated_seeds,
		std::vector< RRVec >& refTable,
		std::vector<int>& refTargets);*/
	double _EstimateInfl(const std::vector<int>& seeds, std::vector<double>& out_cumuInfl);
	void _SetResults(const std::vector<int>& seeds, const std::vector<double>& cumu_spread);

};


/// Reverse Influence Maximization
/// Implementation of Paper [1] (same idea, see: rr_infl.h)
///
/// * Use concurrent optimization for multi-cores, turn on switch /openmp
class RRInfl : 
	public RRInflBase
{
protected:
	std::string file;
	std::string time_file;
	
public:
	RRInfl() : RRInflBase(), 
		file("rr_infl.txt"),
		time_file("time_rr_infl.txt")	
	{}
	

public:
	virtual void Build(graph_type& gf, int k, cascade_type& cascade, size_t num_iter = 1000000); // [1]
	virtual void BuildInError(graph_type& gf, int k, cascade_type& cascade, double epsilon = 0.1); // [1] 0 < epsilon < 1
	
protected:
	void _Build(graph_type& gf, int k, cascade_type& cascade, size_t num_iter = 0); // internal
	// methods for finding the solution
	double DefaultRounds(int n, int m, double epsilon = 0.2); // [1]

};

/// TimPlus Algorithm
/// Implementation of Paper [2] (see: rr_infl.h)
///
class TimPlus : 
	public RRInflBase
{
protected:
	std::string file;
	std::string time_file;

public:
	TimPlus() : RRInflBase(), 
		file("rr_timplus_infl.txt"), 
		time_file("time_rr_timplus_infl.txt") {
	} 

public:
	virtual void Build(graph_type& gf, int k, cascade_type& cascade, double eps = 0.1, double ell = 1.0); // [2]

protected:
	double StepThreshold(int n, double lb, double ell=1.0);
	double RThreshold_0(double eps, double opt, double ell=1.0);
	double LogNChooseK(int n, int k);
	double RThreshold(double eps, double opt, int k, double ell=1.0); // Lemma 3 in [2]
	double EpsPrime(double eps, int k, double ell=1.0); // Last equation in Section 4.1 [2]
};

/// IMM algorithm
/// Implementation of Paper [3] (see: rr_infl.h)
class IMM :
	public TimPlus
{
public:
	IMM() : TimPlus() {
		file = "rr_imm_infl.txt";
		time_file = "time_rr_imm_infl.txt";
	}

public:
	/// override Build
	void Build(graph_type& gf, int k, cascade_type& cascade, double eps = 0.1, double ell = 1.0); // [3]

	double LambdaPrime(double epsprime, int k, double ell, int n); // Equation(9) in [3]
	double LambdaStar(double eps, int k, double ell, int n); // Equation(6) in [3]
};

class RIM_IMM :
	public IMM
{
public:
	MIRandom random;
	RIM_IMM() : IMM() {
		file = "rr_rimm_infl.txt";
		time_file = "time_rr_rimm_infl.txt";
	}

public:
	/// override Build
	void Build(graph_type& gf, int k, cascade_type& cascade, double eps = 0.1, double ell = 1.0, int rounds = 5); // KDD 18
};

class WR_IMM :
	public IMM
{
public:
	MIRandom random;
	WR_IMM() : IMM() {
		file = "rr_wrimm_infl.txt";
		time_file = "time_rr_wrimm_infl.txt";
	}

public:
	/// override Build
	void Build(graph_type& gf, int k, cascade_type& cascade, double eps = 0.1, double ell = 1.0, int rounds = 5); // KDD 18
};

class CR_IMM :
	public IMM
{
public:
	MIRandom random;
	CR_IMM() : IMM() {
		file = "rr_crimm_infl.txt";
		time_file = "time_rr_crimm_infl.txt";
	}

public:
	/// override Build
	void Build(graph_type& gf, int k, cascade_type& cascade, double eps = 0.1, double ell = 1.0, int rounds = 5); // KDD 18
	double LambdaPrime_CR(double epsprime, int k, double ell, int n);
	double LambdaStar_CR(double eps, int k, double ell, int n);
};

/// ASV-RR algorithm <br/>
/// using the RR set method to compute Shapley values of nodes
class ShapleyInfl :  
	public IMM
{
public:
	struct shapleyValueId {
		int sid;
		double value;
	};

	struct ShapleyComparator
	{
		bool operator () (shapleyValueId a, shapleyValueId b);
	};

public:
	ShapleyInfl() : IMM() {
		file = "rrs_ASVRR_infl.txt";
		time_file = "time_rrs_ASVRR_infl.txt";
	}

public:

	std::vector<shapleyValueId> shapleyVId;

public:
	/// Build for Shapley value computation
	void Shapley_Build(graph_type& gf, cascade_type& cascade, GeneralCascade & gc, double eps = 0.1, double ell = 1.0, int topk=50, bool isSingleInf = false);

	/// RR set generation together with Shapley value computation
	void Shapley_AddRRSimulation(LARGE_INT64 num_iter,
		cascade_type& cascade,
		std::vector<double>& shapleyV,
		std::vector<int>& hitCount,  // hitCount[i] is the number of times the RR sets hit node i
		std::vector< std::vector<int> >& refTable,
		LARGE_INT64 & totalEdgeVisited);

};

/// Compute single node influence, using adapted Shapley computation
class SNIInfl :
	public ShapleyInfl   /* single node influence actually inherit from Shapley computation, since it computes
							only single node influence, and no need to store RR sets, and thus it is closer to
							Shapley rather than IMM. Shapley_AddRRSimulation() is reused here.
							*/
{
public:
	SNIInfl() : ShapleyInfl() {
		file = "rr_sni_infl.txt";
		time_file = "time_rr_sni_infl.txt";
	}

};



#endif ///:~ __rr_infl_h__


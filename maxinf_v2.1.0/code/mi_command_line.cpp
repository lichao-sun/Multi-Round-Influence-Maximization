#include "mi_command_line.h"

using namespace std;



std::string MICommandLine::Help()
{
	std::string help =
		"[Help]\n"
		MI_SOFTWARE_META "\n"
		"\n"
		"-h: print the help \n"
		"-st : statistics of the weighted cascade graph\n"
		"-stg : statistics of the general independent cascade graph\n"
		"-b : baseline: random (-br), degree (-bd), degreediscount (-bdd), weighteddegree (-bw), pagerank (-bp) for general ic\n"
		"-g : greedy, SPM and SP1M for general ic\n"
		"-go : greedy with online bound for general ic\n"
		"-p <bound1=10> <bound2=2000> : PMIA with 1/theta from bound1 to bound2 \n"
		"-m <bound1=10> <bound2=2000> : MIA with 1/theta from bound1 to bound2 \n"
		"-gwc filename1 filename2... : c-greedy algorithm\n"
		"-gwcl seeds_list_file : c-greedy algorithm\n"
		"-mis seeds_list_file dist1 dist2 ... distd : MIS(marginal influence sort) algorithm \n"
		"-ts seeds_list_file dist1 dist2 ... distd : BTS(top-selection) algorithm \n"
		"-t seeds_file <num_iter=10000> <seed_set_size = 50> <output_file=GC_spread.txt> <nthreads=1>: test influence spread with seeds \n"
		"-rr : reverse influence maximization algorithms.\n"
		"      -rr1 <sample_num=1000000> (SODA'14) \n"
		"      -rr2 <eps=0.1> <ell=1.0>  (SIGMOD'14, TIM-Plus) \n"
		"      -rr3 <eps=0.1> <ell=1.0>  (SIGMOD'15, IMM) \n"
		"      -rrs <eps=0.1> <ell=1.0>	<k = 50>	compute node Shapley values, with top k nodes ensuring relative error of eps \n"
		"      -rrsn <eps=0.1> <ell=1.0> <k = 50>   compute single node influence ranking, with top k nodes ensuring relative error of eps, adapted from Shapley computation \n"
		"   -rro / -rr1o / -rr2o / -rr3o : parallel optimization of the above \n"
		"-cic <max_time=10.0>: continuos independent cascade model. -cics : for graph statistics \n"
		"-prgen <dampen=0.15> <out_graph=pagerank_regen.txt>: use pagerank to re-generate graph \n"
		"\n"
		"example: max_influence -p 20 2000 < hep.inf \n"
	;

	// read from file
	// Example input:          // Nodes starts from 1, and <edge_size>*2 lines follows
	//	15233 58891 ---------> #nodes #edges(undirected)
	//	1 2 2.380952e-002 2.500000e-001
	//	2 1 2.500000e-001 2.380952e-002
	//	3 4 7.692308e-002 5.000000e-001
	//	4 3 5.000000e-001 7.692308e-002
	//	бн
	//	15232 15233 3.127006e-001 3.127006e-001
	//	15233 15232 3.127006e-001 3.127006e-001
	//
	//  printf("example: python gen_graph.py <parm1> <parm2> ... | max_influence -p 20 2000 \n"); // read from other program
	
	return help;
}

int MICommandLine::Main(int argc, char * argv[])
{
	// create an arguments vector to include all strings
	std::vector<std::string> argVec;
	for (int i = 0; i < argc; i++) {
		std::string param = argv[i];
		argVec.push_back(param);
	}

	return Main(argc, argVec);
}

int MICommandLine::Main(int argc, std::vector<std::string>& argv)
{
	srand((unsigned)time(NULL));
	
	if (argc <= 1) {
		std::cout << Help() << std::endl;
		return 0;
	}

	std::string arg1 = argv[1]; // the switch string like "-abc"
	std::string s;

	// for -h, print help
	s = "-h";
	if (s.compare(arg1) == 0) {
		std::cout << Help() << std::endl;
		return 0;
	}

	// the following contains switches for algorithms
	system("mkdir tmp");
	system("cd tmp");

	// create empty _running_.log to indicate running
	system("del /Q _finished_.log");
	system("echo. 2> _running_.log");

	s = "-r";
	if (s.compare(arg1) == 0){
		BuildRanking(argc, argv);
	}

	s = "-gwc";
	if (s.compare(arg1) == 0){
		TopicAwareCGreedy(argc, argv);
	}

	s = "-gwcl";
	if (s.compare(arg1) == 0){
		TopicAwareCGreedyFromList(argc, argv);
	}

	s = "-mis";
	if (s.compare(arg1) == 0){
		TopicAwareMIS(argc, argv);
	}

	s = "-ts";
	if (s.compare(arg1) == 0){
		TopicAwareTopSelection(argc, argv);
	}

	s = "-t";
	if (s.compare(arg1) == 0){
		TestSeeds(argc, argv);
	}

	s = "-rimt";
	if (s.compare(arg1) == 0){
		RIM_TestSeeds(argc, argv);
	}

	s = "-st";
	if (s.compare(arg1) == 0) {
		GraphStat(argc, argv);
	}

	s = "-b";
	if (s.compare(arg1.substr(0, 2)) == 0) {
		BaselineAlg(argc, argv);
	}

	s = "-g";
	if (s.compare(arg1) == 0) {
		GreedyAlg(argc, argv);
	}

	s = "-rimg";
	if (s.compare(arg1) == 0) {
		RIMGreedyAlg(argc, argv);
	}

	s = "-rimgc";
	if (s.compare(arg1) == 0) {
		RIMCRGreedyAlg(argc, argv);
	}

	s = "-rimga";
	if (s.compare(arg1) == 0) {
		RIMAGreedyAlg(argc, argv);
	}

	s = "-go";
	if (s.compare(arg1) == 0) {
		GreedyOnlineAlg(argc, argv);
	}


	s = "-p";
	if (s.compare(arg1) == 0) {
		PMIAAlg(argc, argv);
	}

	s = "-m";
	if (s.compare(arg1) == 0) {
		MIAAlg(argc, argv);
	}

	s = "-rr";
	if (s.compare(arg1.substr(0, 3)) == 0) {
		RRAlg(argc, argv);
	}

	/*s = "-cic";
	if (s.compare(arg1.substr(0, 4)) == 0) {
		ContinousICAlg(argc, argv);
	}*/

	s = "-prgen";
	if (s.compare(arg1) == 0) {
		PagerankGraphGen(argc, argv);
	}

	// delete _running_.log to indicate finish
	system("del /Q _running_.log");
	system("echo. 2> _finished_.log");

	return 0;
}

void MICommandLine::TestSeeds(int argc, std::vector<std::string>& argv)
{
	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);

	string seeds_file = "Test.txt";
	if (argc >= 3) {
		seeds_file = argv[2];
	}
	int num_iter = NUM_ITER;
	if (argc >= 4) {
		num_iter = std::stoi(argv[3]);
	}
	int seed_set_size = SET_SIZE;
	if (argc >= 5) {
		seed_set_size = std::stoi(argv[4]);
	}
	string outfile = "GC_spread.txt";
	if (argc >= 6) {
		outfile = argv[5];
	}
	int nthreads = 1;
	if (argc >= 7) {
		nthreads = std::stoi(argv[6]);
	}
	cascade.nthreads = (nthreads <= 1) ? 1 : nthreads;

	SeedIO io;
	std::vector<int> seeds = io.Read(seeds_file, gf);
	int size = min(seeds.size(), seed_set_size);

	EventTimer timer;
	timer.SetTimeEvent("start");
	Simu(seeds, outfile, size).toSimulate(cascade);
	timer.SetTimeEvent("end");
	char timefilename[] = "time_test.txt";
	FILE *out;
	fopen_s(&out, timefilename, "w");
	fprintf(out, "%g\n", timer.TimeSpan("start", "end"));
	fclose(out);
}

void MICommandLine::RIM_TestSeeds(int argc, std::vector<std::string>& argv)
{
	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);

	string seeds_file = "Seeds.txt";
	if (argc >= 3) {
		seeds_file = argv[2];
	}
	int rounds = 5;
	if (argc >= 4) {
		rounds = std::stoi(argv[3]);
	}
	int num_iter = NUM_ITER;
	if (argc >= 5) {
		num_iter = std::stoi(argv[4]);
	}
	int seed_set_size = SET_SIZE;
	if (argc >= 6) {
		seed_set_size = std::stoi(argv[5]);
	}
	string outfile = "RIM_GC_spread.txt";
	if (argc >= 7) {
		outfile = argv[6];
	}
	int nthreads = 1;
	if (argc >= 8) {
		nthreads = std::stoi(argv[7]);
	}
	cascade.nthreads = (nthreads <= 1) ? 1 : nthreads;

	/*SeedIO io;
	std::vector<int> seeds_list = io.Read(seeds_file, gf);
	int size = min(seeds_list.size(), seed_set_size);*/
	std::ifstream infile(seeds_file);
	
	//int a, b, c, d, e, f, g, h, i, j;
	//while (infile >> a >> b >> c >> d >> e >> f >> g >> h >> i >> j)
	//{
	//	vector<int> tmp;
	//	tmp.push_back(a);
	//	tmp.push_back(b);
	//	tmp.push_back(c);
	//	tmp.push_back(d);
	//	tmp.push_back(e);
	//	tmp.push_back(f);
	//	tmp.push_back(g);
	//	tmp.push_back(h);
	//	tmp.push_back(i);
	//	tmp.push_back(j);
	//	seeds_list.push_back(tmp);
	//	// process pair (a,b)
	//}

	std::string line;
	vector<vector<int>> seeds_list;
	while (getline(infile, line)) {
		std::istringstream is(line);
		seeds_list.push_back(
			std::vector<int>(std::istream_iterator<int>(is),
			std::istream_iterator<int>()));
	}


	EventTimer timer;
	timer.SetTimeEvent("start");
	RIM_Simu(seeds_list, outfile, 10).toSimulate(cascade, rounds);
	timer.SetTimeEvent("end");
	char timefilename[] = "time_test.txt";
	FILE *out;
	fopen_s(&out, timefilename, "w");
	fprintf(out, "%g\n", timer.TimeSpan("start", "end"));
	fclose(out);
}

void MICommandLine::BuildRanking(int argc, std::vector<std::string>& argv)
{
	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);
	Greedy alg;
	alg.BuildRanking(gf, 100, cascade);
}

void MICommandLine::TopicAwareCGreedy(int argc, std::vector<std::string>& argv)
{
	//int* seed;
	//printf("-gwc filename1 filename2... : c-greedy algorithm\n");
	
	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);

	// read seeds
	set<int> candidates;
	CGreedy::initialize_gwc(argc, argv, candidates, gf);

	EventTimer timer;
	timer.SetTimeEvent("start");
	// CGreedy: candidate greeedy
	CGreedy alg;
	alg.BuildWithCandidates(SET_SIZE, candidates, gf, cascade);
	timer.SetTimeEvent("end");
	char SPTfilename[] = "time_cgreedy.txt";
	FILE *out;
	fopen_s(&out, SPTfilename, "w");
	fprintf(out, "%g\n", timer.TimeSpan("start", "end"));
	fclose(out);

	Simu(alg.GetSeedList(), "GC_cgreedy.txt").toSimulate(cascade);
}

void MICommandLine::TopicAwareCGreedyFromList(int argc, std::vector<std::string>& argv)
{
	set<int> candidates;
	//int* seed;
	printf("-gwcl seeds_list_file : c-greedy algorithm\n");

	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);

	string seeds_file_name = argv[2];
	TopicAwareBase::ReadTopicFile(seeds_file_name, gf);
	candidates = TopicAwareBase::GetAllSeedsIndices();

	EventTimer timer;
	timer.SetTimeEvent("start");
	// CGreedy: candidate greedy
	CGreedy alg;
	alg.BuildWithCandidates(SET_SIZE, candidates, gf, cascade);
	timer.SetTimeEvent("end");
	char SPTfilename[] = "time_cgreedy.txt";
	FILE *out;
	fopen_s(&out, SPTfilename, "w");
	fprintf(out, "%g\n", timer.TimeSpan("start", "end"));
	fclose(out);

	Simu(alg.GetSeedList(), "GC_cgreedy.txt").toSimulate(cascade);
}

void MICommandLine::TopicAwareMIS(int argc, std::vector<std::string>& argv)
{
	string seeds_file_name = argv[2];

	std::vector<double> distribution;
	printf("distribution: ( ");
	for (int i = 3; i < argc; ++i) {
		double dis = std::stof(argv[i]);
		distribution.push_back(dis);
		printf(" %g ", dis);
	}
	printf(")\n");
	// = { 1.0, 0.3, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);

	TopicAwareBase::ReadTopicFile(seeds_file_name, gf);
	// Performance test is in Build
	MIS mis;
	mis.Build(gf, min(SET_SIZE, gf.GetN()), distribution);

	Simu(mis.GetSeedList(), "GC_mis.txt").toSimulate(cascade);
}

void MICommandLine::TopicAwareTopSelection(int argc, std::vector<std::string>& argv)
{
	string seeds_file_name = argv[2];

	std::vector<double> distribution;
	printf("distribution: ( ");
	for (int i = 3; i < argc; ++i) {
		double dis = std::stof(argv[i]);
		distribution.push_back(dis);
		printf(" %g ", dis);
	}
	printf(")\n");
	// double distribution[] = { 1.0, 0.3, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);
	TopicAwareBase::ReadTopicFile(seeds_file_name, gf);
	// Performance test is in Build
	TopSelection ts;
	ts.Build(gf, min(SET_SIZE, gf.GetN()), distribution);

	Simu(ts.GetSeedList(), "GC_top_selection.txt").toSimulate(cascade);
}


void MICommandLine::GraphStat(int argc, std::vector<std::string>& argv)
{
	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GraphStatistics sta;
	sta.Stats(gf, std::cout);
}

void MICommandLine::BaselineAlg(int argc, std::vector<std::string>& argv)
{
	// ("-b : baseline: random (-br), degree (-bd), degreediscount (-bdd), 
	//   weighteddegree (-bw), pagerank (-bp) for general ic\n")

	string arg1(argv[1]);
	bool isRandom = false, isDegree = false, isDegreediscount = false, isWeighteddegree = false, isPagerank = false;
	if (arg1.compare("-b") == 0) {
		// test all baseline
		isRandom = isDegree = isDegreediscount = isWeighteddegree = isPagerank = true;
	}
	else if (arg1.compare("-br") == 0) {
		isRandom = true;
	}
	else if (arg1.compare("-bd") == 0) {
		isDegree = true;
	}
	else if (arg1.compare("-bdd") == 0) {
		isDegreediscount = true;
	}
	else if (arg1.compare("-bw") == 0) {
		isWeighteddegree = true;
	}
	else if (arg1.compare("-bp") == 0) {
		isPagerank = true;
	}

	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);

	// Random
	if (isRandom) {
		cout << "Start: Random" << endl;
		EventTimer timer;
		timer.SetTimeEvent("start");
		RandomPick alg;
		alg.Build(gf, min(SET_SIZE, gf.GetN()));
		timer.SetTimeEvent("end");

		FILE* timetmpfile;
		fopen_s(&timetmpfile, "time_random.txt", "w");
		fprintf(timetmpfile, "%g\n", timer.TimeSpan("start", "end"));
		fclose(timetmpfile);
		//Simu(alg.GetSeedList(), "GC_Random.txt").toSimulate(cascade);
	}
	// Weighted Degree
	if (isWeighteddegree) {
		cout << "Start: Weighted Degree" << endl;
		EventTimer timer;
		timer.SetTimeEvent("start");
		WeightedDegree alg;
		alg.Build(gf, min(SET_SIZE, gf.GetN()));
		timer.SetTimeEvent("end");
		FILE* timetmpfile;
		fopen_s(&timetmpfile, "time_weighteddegree.txt", "w");
		fprintf(timetmpfile, "%g\n", timer.TimeSpan("start", "end"));
		fclose(timetmpfile);
		//Simu(alg.GetSeedList(), "GC_WeightedDegree.txt").toSimulate(cascade);
	}
	// Highest Degree
	if (isDegree) {
		cout << "Start: Highest Degree" << endl;
		EventTimer timer;
		timer.SetTimeEvent("start");
		Degree alg;
		alg.Build(gf, min(SET_SIZE, gf.GetN()));
		timer.SetTimeEvent("end");
		FILE* timetmpfile;
		fopen_s(&timetmpfile, "time_degree.txt", "w");
		fprintf(timetmpfile, "%g\n", timer.TimeSpan("start", "end"));
		fclose(timetmpfile);
		//Simu(alg.GetSeedList(), "GC_Degree.txt").toSimulate(cascade);
	}
	// DegreeDiscountIC
	if (isDegreediscount) {
		cout << "Start: DegreeDiscountIC" << endl;
		EventTimer timer;
		timer.SetTimeEvent("start");
		DegreeDiscount_IC alg;
		alg.Build(gf, min(SET_SIZE, gf.GetN()), 0.01);
		timer.SetTimeEvent("end");

		FILE* timetmpfile;
		fopen_s(&timetmpfile, "time_degreediscount_ic.txt", "w");
		fprintf(timetmpfile, "%g\n", timer.TimeSpan("start", "end"));
		fclose(timetmpfile);
		//Simu(alg.GetSeedList(), "GC_DiscountIC.txt").toSimulate(cascade);
	}
	// Pagerank
	if (isPagerank) {
		cout << "Start: Pagerank" << endl;
		EventTimer timer;
		timer.SetTimeEvent("start");
		Pagerank alg;
		alg.Build(gf, min(SET_SIZE, gf.GetN()));
		timer.SetTimeEvent("end");

		FILE* timetmpfile;
		fopen_s(&timetmpfile, "time_pagerank.txt", "w");
		fprintf(timetmpfile, "%g\n", timer.TimeSpan("start", "end"));
		fclose(timetmpfile);
		//Simu(alg.GetSeedList(), "GC_pagerank.txt").toSimulate(cascade);
	}
}

void MICommandLine::GreedyAlg(int argc, std::vector<std::string>& argv)
{
	// GreedyGC (improved by CELF)
	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);

	EventTimer timer;
	timer.SetTimeEvent("start");
	Greedy alg;
	alg.Build(gf, min(SET_SIZE, gf.GetN()), cascade);
	timer.SetTimeEvent("end");

	FILE* timetmpfile;
	fopen_s(&timetmpfile, "time_greedy_gc.txt", "w");
	fprintf(timetmpfile, "%g\n", timer.TimeSpan("start", "end"));
	fclose(timetmpfile);
	system("copy greedy.txt greedy_gc.txt");
	//Simu(alg.GetSeedList(), "GC_Greedy.txt").toSimulate(cascade);
	system("del /Q tmp\\*");

	//system("pause");

	/*
	// GreedyGC_SPM (improved by CELF)
	QueryPerformanceCounter(&start);
	Greedy::Build(SET_SIZE,SPM_gc::Run);
	QueryPerformanceCounter(&Eend);
	timer = (double)(Eend.QuadPart - start.QuadPart) / freq.QuadPart;
	fopen_s(&timetmpfile, "time_greedy_gc_spm.txt", "w");
	fprintf(timetmpfile,"%g\n", timer);
	fclose(timetmpfile);
	system("copy greedy.txt greedy_gc_spm.txt");
	system("del /Q tmp\\*");
	toSimulate("GC_SPM.txt", Greedy::GetNode, GeneralCascade::Run);

	// GreedyGC_SP1M (improved by CELF)
	QueryPerformanceCounter(&start);
	Greedy::Build(SET_SIZE,SP1M_gc::Run);
	QueryPerformanceCounter(&Eend);
	timer = (double)(Eend.QuadPart - start.QuadPart) / freq.QuadPart;
	fopen_s(&timetmpfile,"time_greedy_gc_sp1m.txt", "w");
	fprintf(timetmpfile,"%g\n", timer);
	fclose(timetmpfile);
	system("copy greedy.txt greedy_gc_sp1m.txt");
	system("del /Q tmp\\*");
	toSimulate("GC_SP1M.txt", Greedy::GetNode, GeneralCascade::Run);*/
}

void MICommandLine::RIMGreedyAlg(int argc, std::vector<std::string>& argv)
{
	int seeds = 10;
	int rounds = 5;
	if (argc >= 3) seeds = std::stoi(argv[3]);
	if (argc >= 4) rounds = std::stoi(argv[4]);
	// GreedyGC (improved by CELF)
	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);

	EventTimer timer;
	timer.SetTimeEvent("start");
	RIM_Greedy alg;
	alg.Build(gf, min(seeds, gf.GetN()), rounds, cascade);
	timer.SetTimeEvent("end");

	FILE* timetmpfile;
	fopen_s(&timetmpfile, "time_greedy_gc.txt", "w");
	fprintf(timetmpfile, "%g\n", timer.TimeSpan("start", "end"));
	fclose(timetmpfile);
	system("copy within-round greedy.txt greedy_gc.txt");
	//Simu(alg.GetSeedList(), "GC_Greedy.txt").toSimulate(cascade);
	system("del /Q tmp\\*");
}

void MICommandLine::RIMCRGreedyAlg(int argc, std::vector<std::string>& argv)
{
	int seeds = 10;
	int rounds = 5;
	//cout << argc << ' ' << stoi(argv[3]) << ' ' << stoi(argv[4]) << endl;
	if (argc >= 3) seeds = std::stoi(argv[2]);
	if (argc >= 4) rounds = std::stoi(argv[3]);
	// GreedyGC (improved by CELF)
	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);

	EventTimer timer;
	timer.SetTimeEvent("start");
	RIM_CRGreedy alg;
	alg.Build(gf, min(seeds, gf.GetN()), rounds, cascade);
	timer.SetTimeEvent("end");

	FILE* timetmpfile;
	fopen_s(&timetmpfile, "time_crgreedy_gc.txt", "w");
	fprintf(timetmpfile, "%g\n", timer.TimeSpan("start", "end"));
	fclose(timetmpfile);
	system("copy cross-round greedy.txt greedy_gc.txt");
	//Simu(alg.GetSeedList(), "GC_Greedy.txt").toSimulate(cascade);
	system("del /Q tmp\\*");
}

void MICommandLine::RIMAGreedyAlg(int argc, std::vector<std::string>& argv)
{
	int seeds = 10;
	int rounds = 5;
	
	if (argc >= 3) seeds = std::stoi(argv[3]);
	if (argc >= 4) rounds = std::stoi(argv[4]);
	// GreedyGC (improved by CELF)
	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);

	EventTimer timer;
	timer.SetTimeEvent("start");
	RIM_AGreedy alg;
	alg.Build(gf, min(seeds, gf.GetN()), rounds, cascade);
	timer.SetTimeEvent("end");

	FILE* timetmpfile;
	fopen_s(&timetmpfile, "time_agreedy_gc.txt", "w");
	fprintf(timetmpfile, "%g\n", timer.TimeSpan("start", "end"));
	fclose(timetmpfile);
	system("copy adaptive greedy.txt greedy_gc.txt");
	//Simu(alg.GetSeedList(), "GC_Greedy.txt").toSimulate(cascade);
	system("del /Q tmp\\*");
}


void MICommandLine::GreedyOnlineAlg(int argc, std::vector<std::string>& argv)
{
	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);
	EventTimer timer;
	timer.SetTimeEvent("start");
	GreedyOnline alg;
	alg.Build(gf, min(SET_SIZE, gf.GetN()), cascade);
	timer.SetTimeEvent("end");

	FILE* timetmpfile;
	fopen_s(&timetmpfile, "time_greedy_online_gc.txt", "w");
	fprintf(timetmpfile, "%g\n", timer.TimeSpan("start", "end"));
	fclose(timetmpfile);
	system("copy greedy_online.txt greedy_online_gc.txt");
	//Simu(alg.GetSeedList(), "GC_Greedy_online.txt").toSimulate(cascade);
	system("del /Q tmp\\*");
}

void MICommandLine::PMIAAlg(int argc, std::vector<std::string>& argv)
{
	//control bound to test PMIA_GC
	int bound1 = 10, bound2 = 2000;
	if (argc >= 3) bound1 = std::stoi(argv[2]);
	if (argc >= 4) bound2 = std::stoi(argv[3]);

	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);

	char SPTfilename[] = "PMIA_control.txt";
	FILE *out;
	fopen_s(&out, SPTfilename, "w");
	char timefilename[] = "time_PMIA_0000.txt";
	char SPT_new_WC[] = "GC_PMIA_0000.txt";
	for (int bound = bound1; bound < bound2; bound += bound){
		printf("%d ", bound);
		double spread, treesize = 0;
#ifdef COUNT
		{
			spread = PMIA::Build(SET_SIZE, bound);
			printf("%g\n", spread);
			continue;
		}
#endif
		sprintf_s(timefilename, "time_PMIA_%04d.txt", bound);
		sprintf_s(SPT_new_WC, "GC_PMIA_%04d.txt", bound);
		EventTimer timer;
		{
			timer.SetTimeEvent("start");
			PMIA alg;
			treesize = alg.Build(gf, min(SET_SIZE, gf.GetN()), bound);
			timer.SetTimeEvent("end");
			printf("\n");
			FILE* timetmpfile;
			fopen_s(&timetmpfile, timefilename, "w");
			fprintf(timetmpfile, "%g\n", timer.TimeSpan("start", "end"));
			fclose(timetmpfile);
			spread = Simu(alg.GetSeedList()).toSimulateOnce(cascade);
			//Simu(alg.GetSeedList(), SPT_new_WC).toSimulate(cascade);
		}
		fprintf(out, "%g %g %d %g\n", timer.TimeSpan("start", "end"), spread, bound, treesize);
	}
	fclose(out);
}

void MICommandLine::MIAAlg(int argc, std::vector<std::string>& argv)
{
	//control bound to test MIA_GC
	int bound1 = 10, bound2 = 2000;
	if (argc >= 3) bound1 = std::stoi(argv[2]);
	if (argc >= 4) bound2 = std::stoi(argv[3]);

	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);

	char SPTfilename[] = "MIA_control.txt";
	FILE *out;
	fopen_s(&out, SPTfilename, "w");
	char timefilename[] = "time_MIA_0000.txt";
	char SPT_new_WC[] = "GC_MIA_0000.txt";
	for (int bound = bound1; bound < bound2; bound += bound){
		printf("%d ", bound);
		double spread, treesize = 0;
#ifdef COUNT
		{
			spread = PMIA::Build(SET_SIZE, bound);
			printf("%g\n", spread);
			continue;
		}
#endif
		sprintf_s(timefilename, "time_MIA_%04d.txt", bound);
		sprintf_s(SPT_new_WC, "GC_MIA_%04d.txt", bound);
		EventTimer timer;
		{
			timer.SetTimeEvent("start");
			MIA alg;
			treesize = alg.Build(gf, min(SET_SIZE, gf.GetN()), bound);
			timer.SetTimeEvent("end");
			printf("\n");

			FILE* timetmpfile;
			fopen_s(&timetmpfile, timefilename, "w");
			fprintf(timetmpfile, "%g\n", timer.TimeSpan("start", "end"));
			fclose(timetmpfile);
			spread = Simu(alg.GetSeedList()).toSimulateOnce(cascade);
			//Simu(alg.GetSeedList(), SPT_new_WC).toSimulate(cascade);
		}
		fprintf(out, "%g %g %d %g\n", timer.TimeSpan("start", "end"), spread, bound, treesize);
	}
	fclose(out);
}

void MICommandLine::RRAlg(int argc, std::vector<std::string>& argv)
{
	string arg1(argv[1]);

	bool isSODA14 = false, isSIGMOD14 = false, isSIGMOD15 = false, isWWW16 = false, isKDD18wr = false, isKDD18cr = false;
	int num_iter = 1000000;
	double eps = 0.1;
	double ell = 1.0;
	int topk = 10;
	int nout = 10;
	int rounds = 5;
	bool isConcurrent = false;
	bool isSingleInf = false;
	bool isShapley = false;
	bool isRepeated = false;
	// switches:
	// -rr  -rro
	// -rr1  -rr1o
	// -rr2  -rr2o
	// -rr3  -rr3o 
	// -rrs  -rrsn -rrso
	if (arg1.compare("-rr") == 0) {
		isSODA14 = isSIGMOD14 = isSIGMOD15 = true; // test algorithms
	}
	else if (arg1.compare("-rro") == 0) {
		isSODA14 = isSIGMOD14 = isSIGMOD15 = true; // test algorithms with parallel optimization
		isConcurrent = true;
	}
	else if (arg1.substr(0, 4).compare("-rr1") == 0) {
		isSODA14 = true;
		if (argc >= 3) num_iter = std::stoi(argv[2]);
		if (arg1.compare("-rr1o") == 0)
			isConcurrent = true;
	}
	else if (arg1.substr(0, 4).compare("-rr2") == 0) {
		isSIGMOD14 = true;
		if (argc >= 3) eps = std::stod(argv[2]);
		if (argc >= 4) ell = std::stod(argv[3]);
		if (arg1.compare("-rr2o") == 0)
			isConcurrent = true;
	}
	else if (arg1.substr(0, 4).compare("-rr3") == 0) {
		isSIGMOD15 = true;
		if (argc >= 3) eps = std::stod(argv[2]);
		if (argc >= 4) ell = std::stod(argv[3]);
		if (argc >= 5) nout = std::stoi(argv[4]);
		if (arg1.compare("-rr3o") == 0)
			isConcurrent = true;
	}
	else if (arg1.substr(0, 4).compare("-rrr") == 0) {
		isWWW16 = true;
		if (argc >= 3) rounds = std::stoi(argv[2]);
		if (argc >= 4) topk = std::stoi(argv[3]);
		if (argc >= 5) eps = std::stod(argv[4]);
		if (argc >= 6) ell = std::stod(argv[5]);
		if (argc >= 7) nout = std::stoi(argv[6]);
		
		if (arg1.compare("-rrro") == 0)
			isConcurrent = true;
	}
	else if (arg1.substr(0, 4).compare("-rrw") == 0) {
		isKDD18wr = true;
		if (argc >= 3) rounds = std::stoi(argv[2]);
		if (argc >= 4) topk = std::stoi(argv[3]);
		if (argc >= 5) eps = std::stod(argv[4]);
		if (argc >= 6) ell = std::stod(argv[5]);
		if (argc >= 7) nout = std::stoi(argv[6]);

		if (arg1.compare("-rrwo") == 0)
			isConcurrent = true;
	}
	else if (arg1.substr(0, 4).compare("-rrc") == 0) {
		isKDD18cr = true;
		if (argc >= 3) rounds = std::stoi(argv[2]);
		if (argc >= 4) topk = std::stoi(argv[3]);
		if (argc >= 5) eps = std::stod(argv[4]);
		if (argc >= 6) ell = std::stod(argv[5]);
		if (argc >= 7) nout = std::stoi(argv[6]);

		if (arg1.compare("-rrco") == 0)
			isConcurrent = true;
	}
	else if (arg1.substr(0, 4).compare("-rrs") == 0) {
		isShapley = true;
		if (argc >= 3) eps = std::stod(argv[2]);
		if (argc >= 4) ell = std::stod(argv[3]);
		if (argc >= 5) topk = std::stoi(argv[4]);
		else if (arg1.compare("-rrsn") == 0){
			isSingleInf = true;
		}
		if (arg1.compare("-rrso") == 0)
			isConcurrent = true;
	}

	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	ReverseGCascade cascade;
	cascade.Build(gf);
	

	if (isSODA14) {
		int maxK = min(SET_SIZE, gf.GetN());
		cout << "=== Algorithm 1: SODA'14 ===" << endl;
		cout << "#seeds = " << maxK << endl;
		RRInfl infl;
		infl.isConcurrent = isConcurrent;
		infl.Build(gf, maxK, cascade, num_iter);
		char rrinfl_simu_file[] = "GC_rr_infl.txt";
		// toSimulate(rrinfl_simu_file, RRInfl::GetNode, GeneralCascade::Run);
	}

	if (isSIGMOD14) {
		int maxK = min(SET_SIZE, gf.GetN());
		cout << "=== Algorithm 2: TimPlus, SIGMOD'14 ===" << endl;
		cout << "#seeds = " << maxK << endl;
		cout << "eps = " << eps << endl;
		cout << "ell = " << ell << endl;
		TimPlus infl;
		infl.isConcurrent = isConcurrent;
		infl.Build(gf, maxK, cascade, eps, ell);
		//char rrinfl_simu_file[] = "GC_rr_timplus_infl.txt";
		// toSimulate(rrinfl_simu_file, TimPlus::GetNode, GeneralCascade::Run);
	}

	if (isSIGMOD15) {
		int maxK = min(SET_SIZE, gf.GetN());
		cout << "=== Algorithm 3: IMM, SIGMOD'15 ===" << endl;
		cout << "#seeds = " << maxK << endl;
		cout << "eps = " << eps << endl;
		cout << "ell = " << ell << endl;
		cout << "isConcurrent = " << isConcurrent << endl;
		IMM infl;
		infl.isConcurrent = isConcurrent;
		infl.Build(gf, maxK, cascade, eps, ell);
		// char rrinfl_simu_file[] = "GC_rr_imm_infl.txt";
		// toSimulate(rrinfl_simu_file, IMM::GetNode, GeneralCascade::Run);
	}

	if (isWWW16) {
		int maxK = min(topk, gf.GetN());
		cout << "=== Algorithm Repeated: IMM ===" << endl;
		cout << "#seeds = " << maxK << endl;
		cout << "eps = " << eps << endl;
		cout << "ell = " << ell << endl;
		cout << "rounds = " << rounds << endl;
		cout << "isConcurrent = " << isConcurrent << endl;
		RIM_IMM infl;
		infl.isConcurrent = isConcurrent;
		infl.Build(gf, maxK, cascade, eps, ell, rounds);
		// char rrinfl_simu_file[] = "GC_rr_imm_infl.txt";
		// toSimulate(rrinfl_simu_file, IMM::GetNode, GeneralCascade::Run);
	}

	if (isKDD18wr) {
		int maxK = min(topk, gf.GetN());
		cout << "=== Algorithm Repeated: WRIMM ===" << endl;
		cout << "#seeds = " << maxK << endl;
		cout << "eps = " << eps << endl;
		cout << "ell = " << ell << endl;
		cout << "rounds = " << rounds << endl;
		cout << "isConcurrent = " << isConcurrent << endl;
		WR_IMM infl;
		infl.isConcurrent = isConcurrent;
		infl.Build(gf, maxK, cascade, eps, ell, rounds);
		// char rrinfl_simu_file[] = "GC_rr_imm_infl.txt";
		// toSimulate(rrinfl_simu_file, IMM::GetNode, GeneralCascade::Run);
	}

	if (isKDD18cr) {
		int maxK = min(topk, gf.GetN());
		cout << "=== Algorithm Repeated: CRIMM ===" << endl;
		cout << "#seeds = " << maxK << endl;
		cout << "eps = " << eps << endl;
		cout << "ell = " << ell << endl;
		cout << "rounds = " << rounds << endl;
		cout << "isConcurrent = " << isConcurrent << endl;
		CR_IMM infl;
		infl.isConcurrent = isConcurrent;
		infl.Build(gf, maxK, cascade, eps, ell, rounds);
		// char rrinfl_simu_file[] = "GC_rr_imm_infl.txt";
		// toSimulate(rrinfl_simu_file, IMM::GetNode, GeneralCascade::Run);
	}


	if (isShapley && !isSingleInf) {
		cout << "=== Algorithm ASV-RR: Approximated Shapley Value using RR sets ===" << endl;
		cout << "eps = " << eps << endl;
		cout << "ell = " << ell << endl;
		cout << "k = " << topk << endl;
		GeneralCascade gc;
		gc.Build(gf);

		ShapleyInfl infl;
		infl.isConcurrent = isConcurrent;
		infl.Shapley_Build(gf, cascade, gc, eps, ell, topk, isSingleInf);
	}

	if (isShapley && isSingleInf) {
		cout << "=== Adapt ASV-RR to compute single node influence ===" << endl;
		cout << "eps = " << eps << endl;
		cout << "ell = " << ell << endl;
		cout << "k = " << topk << endl;
		GeneralCascade gc;
		gc.Build(gf);

		SNIInfl infl;
		infl.isConcurrent = isConcurrent;
		infl.Shapley_Build(gf, cascade, gc, eps, ell, topk, isSingleInf); // use the Shapley computation, but adapt it for compute single node influence
	}

}

//void MICommandLine::ContinousICAlg(int argc, std::vector<std::string>& argv)
//{
//	// Run continuous-time independent cascade model.
//	cout << "=== Continous IC model ===" << endl;
//
//	ContGraphFactory fact;
//	ContGraph gf = fact.Build(std::cin);
//
//	string s = "-cics"; // statistics
//	string arg1 = argv[1]; 
//
//	if (s.compare(arg1) == 0) {
//		GraphStatisticsT<ContGraph> sta;
//		sta.Stats(gf, std::cout);
//	}
//	else {
//		double maxTime = 10.0;
//		if (argc >= 3) maxTime = std::stod(argv[2]);
//		ContGeneralCascade contcascade;
//		contcascade.Build(gf, maxTime);
//
//		EventTimer timer;
//
//		cout << "Run Greedy" << endl;
//		timer.SetTimeEvent("greedy_start");
//		Greedy alg;
//		alg.Build(gf, min(SET_SIZE, gf.GetN()), contcascade);
//		timer.SetTimeEvent("greedy_end");
//		cout << "  Time: " << timer.TimeSpan("greedy_start", "greedy_end") << endl;
//
//		cout << "Run simulation (CIC) for Greedy. Max time=" << maxTime << endl;
//		char f4[] = "GC_Greedy_cont.txt";
//		timer.SetTimeEvent("g2_start");
//		Simu(alg.GetSeedList(), f4).toSimulate(contcascade);
//		timer.SetTimeEvent("g2_end");
//		cout << "  Time: " << timer.TimeSpan("g2_start", "g2_end") << endl;
//	}
//}

void MICommandLine::PagerankGraphGen(int argc, std::vector<std::string>& argv)
{
	double dampen = 0.15;
	std::string out_graph = "pagerank_regen.txt";
	if (argc >= 3) dampen = std::stod(argv[2]);
	if (argc >= 4) out_graph = argv[3];

	GraphFactory fact;
	Graph gf = fact.Build(std::cin);

	cout << "Start: Pagerank" << endl;
	EventTimer timer;
	timer.SetTimeEvent("start");
	PagerankRegen alg;
	alg.Build(gf, out_graph, min(SET_SIZE, gf.GetN()), dampen);
	timer.SetTimeEvent("end");
	cout << "  time: " << timer.TimeSpan("start", "end") << endl;
}


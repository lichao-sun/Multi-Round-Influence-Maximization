import sys
import unittest
from inspect import getmembers, isfunction, isclass
import max_influence as MI


TEST_GRAPH = "test_graph.inf"

class TestMIModule(unittest.TestCase):
    def test_module(self):
        functions_list = [o for o in getmembers(MI) if isfunction(o[1])]
        classes_list = [o for o in getmembers(MI) if isclass(o[1])]
        print("=== Module functions ===")
        for o in functions_list:
            print(o)

        print("=== Module classes ===")
        for o in classes_list:
            print(o)

    def test_version(self):
        print("=== Build Version ===")
        MI.PrintVersion()

    def test_basic(self):
        print("--- General Cascade ---")
        pb = MI.PdfCdfConverter()
        # print(pb)
        print("prob_pos: %s" % str(pb.ExpPDF(0.3, 0.4)))

        pb = MI.MIRandom()
        for i in xrange(10):
            print(pb.RandUnit())

        timer = MI.EventTimer()
        timer.SetTimeEvent("start")
        timer.SetTimeEvent("end")
        print(timer.TimeSpan("start", "end"))

    def test_edge(self):
        e = MI.Edge()
        e.u = 1
        e.v = 2
        e.c = 100
        e.w1 = 0.8
        # e.w2 = 0.3
        print("%d %d %d %f %f" % (e.u, e.v, e.c, e.w1, e.w2))


class TestMIGraphFactory(unittest.TestCase):
    def test_graph_factory_stdin(self):
        # this is static
        print("--- Read Graph from stdin ---")
        fact = MI.GraphFactory()
        stdin = MI.PyHelper().GetStdin()
        gf = fact.Build(stdin)
        print("Node:%d" % gf.GetN())
        print("Edge:%d" % gf.GetM())

    def test_graph_factory_input_filestream(self):
        print("--- Read Graph from File ---")
        fact = MI.GraphFactory()
        fin = MI.PyHelper().GetInputFileStream(TEST_GRAPH)
        gf = fact.Build(fin)
        print("Node:%d" % gf.GetN())
        print("Edge:%d" % gf.GetM())

    def test_graph_general_cascade(self):
        print("--- General Cascade ---")
        fact = MI.GraphFactory()
        fin = MI.PyHelper().GetInputFileStream(TEST_GRAPH)
        gf = fact.Build(fin)
        gc = MI.GeneralCascade()
        gc.Build(gf)


class TestMI_HeuristicAlgos(unittest.TestCase):
    def setUp(self):
        fact = MI.GraphFactory()
        fin = MI.PyHelper().GetInputFileStream(TEST_GRAPH)
        self.gf = fact.Build(fin)
        self.gc = MI.GeneralCascade()
        self.gc.Build(self.gf)

    def test_pmia(self):
        print("--- PMIA ---")
        alg = MI.PMIA()
        alg.Build(self.gf, 50, 250)
        print alg.filename(230)
        ls = alg.GetSeedList()
        print ls
        print [ls[i] for i in xrange(len(ls))]

    def test_mia(self):
        print("--- MIA ---")
        alg = MI.MIA()
        alg.Build(self.gf, 50, 250)
        print alg.filename(230)
        ls = alg.GetSeedList()
        print [ls[i] for i in xrange(len(ls))]

    def test_graph_stat(self):
        print("--- GraphStatistics ---")
        alg = MI.GraphStatistics()
        stdout = MI.PyHelper().GetStdout()
        alg.Stats(self.gf, stdout)

    def test_random_pick(self):
        print("--- RandomPick ---")
        alg = MI.RandomPick()
        alg.Build(self.gf, 50)
        ls = alg.GetSeedList()
        print [ls[i] for i in xrange(len(ls))]

    def test_degree(self):
        print("--- Degree ---")
        alg = MI.Degree()
        alg.Build(self.gf, 50)
        ls = alg.GetSeedList()
        print [ls[i] for i in xrange(len(ls))]

    def test_degreediscount_ic(self):
        print("--- DegreeDiscount_IC ---")
        alg = MI.DegreeDiscount_IC()
        alg.Build(self.gf, 50)
        ls = alg.GetSeedList()
        print [ls[i] for i in xrange(len(ls))]

    def test_weighted_degree(self):
        print("--- WeightedDegree ---")
        alg = MI.WeightedDegree()
        alg.Build(self.gf, 50)
        ls = alg.GetSeedList()
        print [ls[i] for i in xrange(len(ls))]

    def test_pagerank(self):
        print("--- Pagerank ---")
        alg = MI.Pagerank()
        alg.Build(self.gf, 50)
        ls = alg.GetSeedList()
        print [ls[i] for i in xrange(len(ls))]


class TestMI_GreedyAlgos(unittest.TestCase):
    def setUp(self):
        fact = MI.GraphFactory()
        fin = MI.PyHelper().GetInputFileStream(TEST_GRAPH)
        self.gf = fact.Build(fin)
        self.gc = MI.GeneralCascade()
        self.gc.Build(self.gf)

    def test_Greedy(self):
        print("--- Greedy ---")
        alg = MI.Greedy()
        cascade = MI.GeneralCascade()
        cascade.Build(self.gf)
        alg.Build(self.gf, 50, cascade)
        ls = alg.GetSeedList()
        print [ls[i] for i in xrange(len(ls))]


class TestMI_RRAlgos(unittest.TestCase):
    def setUp(self):
        fact = MI.GraphFactory()
        fin = MI.PyHelper().GetInputFileStream(TEST_GRAPH)
        self.gf = fact.Build(fin)
        self.gc = MI.GeneralCascade()
        self.gc.Build(self.gf)

    def test_RRInfl(self):
        print("--- RRInfl ---")
        alg = MI.RRInfl()
        cascade = MI.ReverseGCascade()
        cascade.Build(self.gf)
        alg.Build(self.gf, 50, cascade, 100000)
        ls = alg.GetSeedList()
        print [ls[i] for i in xrange(len(ls))]

    def test_RRInfl_InError(self):
        print("--- RRInfl ---")
        alg = MI.RRInfl()
        cascade = MI.ReverseGCascade()
        cascade.Build(self.gf)
        alg.BuildInError(self.gf, 50, cascade, 5.0)
        ls = alg.GetSeedList()
        print [ls[i] for i in xrange(len(ls))]

    def test_TimPlus(self):
        print("--- TimPlus ---")
        alg = MI.TimPlus()
        cascade = MI.ReverseGCascade()
        cascade.Build(self.gf)
        alg.Build(self.gf, 50, cascade, eps=0.1, ell=1.0)
        ls = alg.GetSeedList()
        print [ls[i] for i in xrange(len(ls))]

    def test_IMM(self):
        print("--- IMM ---")
        alg = MI.IMM()
        cascade = MI.ReverseGCascade()
        cascade.Build(self.gf)
        alg.Build(self.gf, 50, cascade, eps=0.1, ell=1.0)
        ls = alg.GetSeedList()
        print [ls[i] for i in xrange(len(ls))]
        return ls

    def test_Simulation(self):
        print("--- Simulation ---")
        ls = self.test_IMM()
        simu = MI.Simu(ls, "test_simu.txt", 50)
        cs = MI.GeneralCascade()
        cs.Build(self.gf)
        simu.toSimulate(cs)


if __name__ == "__main__":
    print("============ In test ==============")
    unittest.main()

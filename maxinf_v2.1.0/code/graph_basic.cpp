#include "graph_basic.h"
#include "graph.h"




void Node::Serialize(std::ostream& sout)
{
	sout << label;
}

void Node::Deserialize(std::istream& sin)
{
	sin >> label;
}

bool Node::operator == (const Node& n) const
{
	return label == n.label;
}


/// Edge adapter for downward compatibility.
///  (1) Node id in input starts from 1; and the storage index starts from 0. (thus minus 1)
///  (2) Probability is stored in a logrithmic form. (deal with values that are close to 0)
class _EdgeDownwardAdapter
{
public:
	/// input starts from 1 and take log
	static void ApplyTransform(Edge& e);
	static inline double logStar(double v, double INF_POS_MIN = 1e-32) {
		return log(v > INF_POS_MIN ? v : INF_POS_MIN);
	}
};


void _EdgeDownwardAdapter::ApplyTransform(Edge& e)
{
	//e.u--;
	//e.v--;
	//e.w1 = -logStar(e.w1);
	//e.w2 = -logStar(e.w2);
}

//////////////////////////////////////////
// Edge

void Edge::Serialize(std::ostream& sout, IGraph& gf) {
	sout << u << "\t" << v << "\t"
		<< w1 << "\t" << w2;
}
void Edge::Deserialize(std::istream& sin, IGraph& gf) {
	Node su;
	Node sv;
	su.Deserialize(sin);
	sv.Deserialize(sin);
	u = gf.InsertNode(su);
	v = gf.InsertNode(sv);

	// sin >> u >> v;
	sin >> w1 >> w2;
	// _EdgeDownwardAdapter::ApplyTransform(*this);
}

//////////////////////////////////////////
// ContEdge
void ContEdge::Serialize(std::ostream& sout, IGraph& gf)  {
	sout << u << "\t" << v << "\t"
		<< u_a << "\t" << u_b << "\t"
		<< v_a << "\t" << v_b << "\t"
		<< w1 << "\t" << w2;
}

void ContEdge::Deserialize(std::istream& sin, IGraph& gf) {
	Node su;
	Node sv;
	su.Deserialize(sin);
	sv.Deserialize(sin);
	u = gf.InsertNode(su);
	v = gf.InsertNode(sv);

	// sin >> u >> v;
	sin >> u_a >> u_b
		>> v_a >> v_b
		>> w1 >> w2;
	// _EdgeDownwardAdapter::ApplyTransform(*this);
}

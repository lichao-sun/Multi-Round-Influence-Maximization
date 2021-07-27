#include "contcascade.h"


using namespace std;

///////////////////////////////////
// ContGeneralCascade

void ContGeneralCascade::Build(ContGraph& gf, double maxTime)
{
	_Build(gf);
	this->maxTime = maxTime;
}

double ContGeneralCascade::Run(int num_iter, int size, int set[])
{
	if (gf == NULL) {
		throw NullPointerException("Please Build Graph first. (gf==NULL)");
	}

	int resultSize = 0;
	for (int it = 0; it < num_iter; it++)
	{
		// In every iteration: it
		TimeQueue Q;
		vector<bool> active(n, false);

		for (int i = 0; i < size; i++)
		{
			// old: active[set[i]] = true;
			Q.push(TimeItem(set[i], 0.0));
			resultSize++;
		}

		// start simulation:
		while (!Q.empty())
		{
			TimeItem cur = Q.top();
			// not timed out
			if (cur.time > maxTime) {
				break;
			}
			Q.pop();

			// Q may have redundant node (activated in different time)
			// And we use counted to sift out the earliest node.
			if (active[cur.id]) {
				continue;
			}
			active[cur.id] = true;
			resultSize++;

			//if (it==0) {
			//	cout << "  " << cur.id << "\t" << cur.time << endl;
			//}
			int k = gf->GetNeighborCount(cur.id);
			for (int i = 0; i < k; i++)
			{
				graph_type::edge_type& e = gf->GetEdge(cur.id, i);
				if (active[e.v]) continue;

				double delta = random.RandWeibull(e.u_a, e.u_b);
				double newTime = cur.time + delta;
				if (newTime < maxTime)
				{
					// add new node e.v
					Q.push(TimeItem(e.v, newTime));
					// old: active[e.v] = true;
					// resultSize++;
					// break;
				}
			}
		}
	}
	return (double)resultSize / (double)num_iter;
}


bool TimeItem::operator<(const TimeItem& b) const
{
	return (time > b.time);
}

TimeItem::TimeItem(int id/*=0*/, double time/*=0*/) : id(id), time(time)
{
}



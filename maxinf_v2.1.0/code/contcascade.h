#ifndef __contcascade_h__
#define __contcascade_h__


#include <iostream>
#include <vector>
#include <queue>
#include "common.h"
#include "graph.h"
#include "cascade.h"


/// Activated time (used in Continous-time Cascade-related model)
class TimeItem
{
public:
	int id;
	double time;
	TimeItem(int id = 0, double time = 0);
public:
	bool operator < (const TimeItem& b) const;
};


/// Continous-time General Cascade
class ContGeneralCascade :
	public CascadeT<ContGraph>
{
public:
	typedef ContGeneralCascade self_type;
	typedef CascadeT<ContGraph> base_type;

protected:
	typedef std::priority_queue<TimeItem> TimeQueue; // smaller to larger

	double maxTime;

public:
	void Build(ContGraph& gf, double maxTime);
	double Run(int num_iter, int size, int set[]);
};



#endif // __contcascade_h__
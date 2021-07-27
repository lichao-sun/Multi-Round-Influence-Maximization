#ifndef INDEPEND_CASCADE_H
#define INDEPEND_CASCADE_H

#include "cascade.h"

/// Implement Independent Cascade
class IndependCascade :
	public CascadeT<Graph>
{
protected:
	double ratio;

public:
	void Build(Graph& gf, double ratio);
	double Run(int num_iter, int size, int set[]);
};



#endif

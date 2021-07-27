#include <algorithm>
#include "cgreedy.h"
#include "graph.h"
#include "common.h"


using namespace std;

CGreedy::CGreedy() 
{
	file = "cgreedy.txt";
}

void CGreedy::BuildWithCandidates(int num, set<int>& candidates, Graph& gf, ICascade& cascade)
{
	n = gf.GetN();
	top = num;
	d.resize(num, 0);
	list.resize(num, 0);

	int* set = new int[num];

	double old = 0.0;

	double *improve=new double[n];
	int *lastupdate=new int[n];
	int *heap=new int[candidates.size()];
	std::set<int>::iterator iter=candidates.begin();
	for(size_t i=0;i<candidates.size();i++)
	{
		heap[i] = *iter;
		lastupdate[*iter] = -1;
		improve[*iter] = (double)(n+1);//initialize largest
		iter++;
	}
	//printf("%d\n", 666);

	for (int i=0; i<top; i++)
	{
		int ccc = 0;
		//printf("%d\n",i);
		while (lastupdate[heap[0]] != i)
		{
			//printf("%d %d %d\n", i, heap[0], ccc);
			ccc++;
			lastupdate[heap[0]] = i;
			set[i] = heap[0];
			//printf("GreedyGC_SPM %d %d\n",heap[0],improve[heap[0]]);
			improve[heap[0]] = cascade.Run(NUM_ITER, i + 1, set) - old;

			char tmpfilename[200];
			sprintf_s(tmpfilename, "tmp/%02d%05d.txt", i, heap[0]);
			//FILE *tmpfile;
			//fopen_s(&tmpfile, tmpfilename,"w");
			//fprintf(tmpfile, "%g\n", improve[heap[0]]); 
			//fclose(tmpfile);

			int x = 0;
			while (x*2+2<= ((int)candidates.size()) - i)
			{
				int newx=x*2+1;
				if ((newx+1< ((int) candidates.size())-i) && (improve[heap[newx]]<improve[heap[newx+1]]))
					newx++;
				if (improve[heap[x]]<improve[heap[newx]])
				{
					int t=heap[x];
					heap[x] = heap[newx];
					heap[newx] = t;
					x = newx;
				}
				else
					break;
			}
		}
		//printf("i:%d, heap[0]: %d\n",i,heap[0]);
		//printf("4020: %g\n",improve[4020]);
		//printf("1771: %g\n",improve[1771]);

		//used[heap[0]] = true;
		set[i] = heap[0];
		list[i] = heap[0];
		d[i] = improve[heap[0]];
		old+=d[i];

		//char bakname[200];
		//sprintf(bakname, "greedychoice%02d.txt", i+1);
		//FILE *bak = fopen(bakname, "w");
		//fprintf(bak, "%6d\t%d\t%g\n", i+1, heap[0], improve[heap[0]]);
		//fclose(bak);
		//printf("%d\t%g\n", i+1, improve[131]);

		heap[0] = heap[candidates.size()-i-1];//suibian nayige bushang
		int x = 0;
		while (x*2+2<= ((int)candidates.size())-i-1)
		{
			int newx=x*2+1;
			if ((newx+1<((int)candidates.size())-i-1) && (improve[heap[newx]]<improve[heap[newx+1]]))
				newx++;
			if (improve[heap[x]]<improve[heap[newx]])
			{
				int t=heap[x];
				heap[x] = heap[newx];
				heap[newx] = t;
				x = newx;
			}
			else
				break;
		}
	}

	//FILE *out;
	//fopen_s(&out, file.c_str(), "w");
	//fprintf(out, "%d\n", top);
	//for (int i=0; i<top; i++)
	//	fprintf(out, "%d\t%g\n", list[i], d[i]);	//the nodes we want!
	//fclose(out);
	WriteToFile(file, gf);

	SAFE_DELETE_ARRAY(heap);
	SAFE_DELETE_ARRAY(lastupdate);
	SAFE_DELETE_ARRAY(improve);
	SAFE_DELETE_ARRAY(set);
	//	delete[] used;
}


void CGreedy::initialize_gwc(int argc, std::vector<std::string>& argv, std::set<int>& candidates, IGraph& gf)
{
	for (int i = 2; i < argc; i++)
	{
		SeedIO io;
		std::vector<int>& seeds = io.Read(argv[i], gf);
		candidates.insert(seeds.begin(), seeds.end());
	}
}

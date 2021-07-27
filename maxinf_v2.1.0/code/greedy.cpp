#include <cstdio>
#include <string>
#include <cstdlib>
#include "event_timer.h"
#include "greedy.h"
#include "graph.h"
#include "common.h"
#include "cascade.h"

//#include <set>
using namespace std;


Greedy::Greedy() 
{
	file = "greedy.txt";
}

void Greedy::Build(IGraph& gf, int num, ICascade& cascade)
{
	n = gf.GetN();
	top = num;
	d.resize(num, 0);
	list.resize(num, 0);

	bool *used= new bool[n];
	memset(used, 0, sizeof(bool)*n);
	int* set = new int[num];

	double old = 0.0;

	EventTimer timer;
	timer.SetTimeEvent("start");

	double *improve = new double[n];
	int *lastupdate = new int[n];
	int *heap = new int[n];
	for (int i=0; i<n; i++)
	{
		heap[i] = i;
		lastupdate[i] = -1;
		improve[i] = (double)(n+1);//initialize largest
	}

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
			improve[heap[0]] = cascade.Run(NUM_ITER, i+1, set) - old;
			char tmpfilename[200];
			sprintf_s(tmpfilename, "tmp/%02d%05d.txt", i, heap[0]);
			//FILE *tmpfile;
			//fopen_s(&tmpfile, tmpfilename,"w");
			//fprintf(tmpfile, "%g\n", improve[heap[0]]); 
			//fclose(tmpfile);

			int x = 0;
			while (x*2+2<=n-i)
			{
				int newx=x*2+1;
				if ((newx+1<n-i) && (improve[heap[newx]]<improve[heap[newx+1]]))
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

		used[heap[0]] = true;
		set[i] = heap[0];
		list[i] = heap[0];
		d[i] = improve[heap[0]];
		old+=d[i];
		//cout << heap[0] << ' ' << improve[heap[0]] << endl;
		cout << heap[0] << endl;
		//char bakname[200];
		//sprintf(bakname, "greedychoice%02d.txt", i+1);
		//FILE *bak = fopen(bakname, "w");
		//fprintf(bak, "%6d\t%d\t%g\n", i+1, heap[0], improve[heap[0]]);
		//fclose(bak);
		//printf("%d\t%g\n", i+1, improve[131]);

		heap[0] = heap[n-i-1];
		int x = 0;
		while (x*2+2<=n-i)//bug should-1
		{
			int newx=x*2+1;
			if ((newx+1<n-i) && (improve[heap[newx]]<improve[heap[newx+1]]))	//bug should-1
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
		//cout << lastupdate[heap[0]] << ' ' << i << ' ' << ccc << endl;
	}

	//FILE *out;
	//fopen_s(&out, file.c_str(), "w");
	//fprintf(out, "%d\n", top);
	//for (int i=0; i<top; i++)
	//	fprintf(out, "%d\t%g\n", list[i], d[i]);	//the nodes we want!
	//fclose(out);

	//cout << 1 << endl;
	timer.SetTimeEvent("end");
	cout << "time: " << timer.TimeSpan("start", "end") << endl;

	for (const auto &s : list) std::cout << s << '\t';
	std::cout << std::endl;

	WriteToFile(file, gf);

	SAFE_DELETE_ARRAY(set);
	SAFE_DELETE_ARRAY(heap);
	SAFE_DELETE_ARRAY(lastupdate);
	SAFE_DELETE_ARRAY(improve);
	SAFE_DELETE_ARRAY(used);
}

void Greedy::BuildRanking(IGraph& gf, int num, ICascade& cascade)
{
	this->n = gf.GetN();
	this->top = num;
	d.resize(num, 0);
	list.resize(num, 0);

	/*int* set = new int[num];
	double* inf = new double[top];*/
	//for(int i=0;i < top;i++)
	//{
	//	set[i] =0;
	//	inf[i] =0;
	//}

	vector<double> improve(n, 0.0);
	int tmp[1];
	for (int i=0; i<n; i++)
	{
		tmp[0] = i;
		improve[i] = cascade.Run(NUM_ITER, 1, tmp);
		if (improve[i] > list[top - 1])
		{
			d[top - 1] = improve[i];
			list[top - 1] = i;
			int j = top - 2;
			while(j>=0)
			{
				if (improve[i] > d[j])
				{
					int int_tmp = list[j];
					double inf_tmp = d[j];
					list[j] = i;
					d[j] = improve[i];
					list[j + 1] = int_tmp;
					d[j + 1] = inf_tmp;
				}
				else
					break;
				j--;
			}
		}
	}

	cout << 2 << endl;
	WriteToFile(file, gf);
	//SAFE_DELETE_ARRAY(set);
	//SAFE_DELETE_ARRAY(improve);
	//SAFE_DELETE_ARRAY(inf);
	//delete[] used;
}


void Greedy::BuildFromFile(IGraph& gf, const char* name)
{
	ReadFromFile(name, gf);
}

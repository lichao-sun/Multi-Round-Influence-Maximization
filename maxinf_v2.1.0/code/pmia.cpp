#include <cstdio>
#include <string>
#include <cstdlib>
#include <cmath>
#include <ctime>
// #include <windows.h>
#include <vector>
#include <cassert>
#include <algorithm>
#include <fstream>
#include "event_timer.h"
#include "pmia.h"
#include "graph.h"
#include "common.h"


using namespace std;

PMIA::PMIA() {
	Initialize(0, 0);
}

void PMIA::Initialize(int totaln, int topk) {
	n = totaln;
	top = topk;
	sprintf_s(file, "SPT_new_0000.txt");
	longest = log(100.0);
	
	if (n > 0) 
	{
		dd.resize(n, 0);
		dp.resize(n, 1.0);
		used.resize(n);
		self.resize(n);
		lastupdate.resize(n, -1);
		delta.resize(n);
		children.resize(n, NULL);
		path.resize(n, NULL);
	}

	S = NULL;
	distance = NULL;
	b = NULL;
	numchild = NULL;
	queue = NULL;
	heap = NULL;
	childlist = NULL;
	oldchildlist = NULL;
	parent = NULL;
	childnum = NULL;
	allb = NULL;

	if (top > 0) {
		d.resize(top);
		list.resize(top);
		validlist.resize(top, NULL);
	}
}



bool largepair(pair<int,double> a,pair<int,double> b)
{
	if(a.second<b.second)
		return false;
	else 
		return true;
}

int PMIA::GetMax(int round)
{
		vector< pair<int,double> > person_influence;
		double max = -1000000.0;
		int mp = -1;
		for (int j=0; j<n; j++)
			if (!used[j] && lastupdate[j]!=round)
			{
				person_influence.push_back(pair<int,double>(j,dp[j]));
				double tmp = dp[j];
				if (tmp >max)
				{
					max = tmp;
					mp = j;
				}
			}
		/*	if(round == 0){
		sort(person_influence.begin(),person_influence.end(),largepair);//should pass a function '<',but i pass a function '>'
		ofstream fout("first.txt");
		vector<pair<int,double>>::iterator iter = person_influence.begin();
		int cnt=0;
		while(iter!=person_influence.end() && cnt!=1000)
		{
			fout<<iter->first<<" "<<iter->second<<endl;
			cnt++;
			iter++;
		}
		printf("%d",round);
		fout.close();
		person_influence.clear();
		system("pause");
		}*/
		return mp;


}

int PMIA::GetMax0(int round)
{
		int mp = -1;
		return mp;
}

int PMIA::generateSPT_newfrom(Graph& gf, int round, int node){
	ProbTransfom trans(gf.edgeForm);
	int top=0, bottom=0;
	distance[node]=0;
	heap[0]=node;
	top++;
	while (true){
		//pop out of heap
		if (distance[heap[0]]<longest) S[heap[0]]=-1; 
		else break;
		childlist[bottom++]=heap[0];
		for (int i=0;i < gf.GetNeighborCount(heap[0]);i++){
			auto& e = gf.GetEdge(heap[0],i);
			if (used[e.v] || S[e.v]<0) continue;
			if (distance[e.v]>distance[heap[0]] + trans.LogNegProb(e.w1) + EPS) {
				if (S[e.v]>=n){
					distance[e.v] = distance[heap[0]] + trans.LogNegProb(e.w1);
					heap[top]=e.v;
					int j=top++, x=(j-1)/2;
					double temp=distance[heap[j]];
					while (j>0) {
						if (distance[heap[x]]>temp){
							heap[j]=heap[x];
							if (S[heap[j]]<n) S[heap[j]]=j;
							j=x;
							x=(j-1)/2;
						}
						else break;
					}
					heap[j]=e.v;
					S[e.v]=j;
				}
				else{
					distance[e.v] = distance[heap[0]] + trans.LogNegProb(e.w1);
					int j=S[e.v], x=(j-1)/2;
					double temp=distance[heap[j]];
					while (j>0) {
						if (distance[heap[x]]>temp){
							heap[j]=heap[x];
							if (S[heap[j]]<n) S[heap[j]]=j;
							j=x;
							x=(j-1)/2;
						}
						else break;
					}
					heap[j]=e.v;
					S[e.v]=j;
				}
			}//endif
			
		}// end for
		heap[0]=heap[--top];
		if (!top) break;
		//siftdown
		int j=0, x=j*2+1;
		double temp=distance[heap[j]];
		while (x<top){
			if (x+1<top && distance[heap[x+1]]<distance[heap[x]]) x++;
			if (distance[heap[x]]<temp){
				heap[j]=heap[x];
				S[heap[j]]=j;
				j=x; x=j*2+1;
			}
			else break;
		}
		heap[j]=heap[top];
		if (S[heap[j]]<n) S[heap[j]]=j;
	}
	for (int i=0;i<bottom;i++)
	{
		int child=childlist[i];
		distance[child]=longest;
		S[child]=n;
		parent[child]=-1;
	}
	//update tree node set and heuristic
	validlist[round]=new bool[n];
	memset(validlist[round],false,n);
	for (int i=0;i<bottom;i++){
		int child=childlist[i];
		oldchildlist[i]=child;
		if (dd[child]>0) {
			for (int j=0;j<dd[child];j++) dp[children[child][j]]-=delta[child][j]*(1-self[child][j]);

			//block some nodes in child's inverse tree
			for (int j=0;j<round;j++) {
				if (validlist[j][child]) {
					int k;
					for (k=0;k<dd[child] && children[child][k]!=list[j];k++) ;
					if (k==dd[child]) {
						continue;
					}

					for (;k && children[child][k]!=node;k=path[child][k]) {
					}
					validlist[j][child]=!k; //k==0 means ok
				}
			}
			validlist[round][child]=true;
		}
	}
	for (int i=1;i<bottom;i++) 
		generateSPT_newto(gf, oldchildlist[i]);
	return bottom;
}

int PMIA::generateSPT_newto(Graph& gf, int node){
	ProbTransfom trans(gf.edgeForm);
	int top=0, bottom=0;
	if (used[node]) {
		dd[node]=1;
		path[node][0]=0;
		self[node][0]=1;
		delta[node][0]=0;
		return 1;
	}
	distance[node]=0;
	heap[0]=node;
	top++;
	parent[node]=node;
	b[node]=1;
	while (true){
		//stack out of heap
		if (distance[heap[0]]<longest) S[heap[0]]=-bottom-1; 
		else break;
		childlist[bottom++]=heap[0];
		if (parent[heap[0]]!=heap[0]) numchild[parent[heap[0]]]++;
		if (!used[heap[0]])
		for (int i=0;i < gf.GetNeighborCount(heap[0]);i++){
			auto& e = gf.GetEdge(heap[0], i);
			if ((used[e.v] && !validlist[lastupdate[e.v]][node]) || S[e.v]<0) continue;
			if (distance[e.v]>distance[heap[0]] + trans.LogNegProb(e.w2) + EPS) {
				parent[e.v]=heap[0];
				b[e.v] = trans.Prob(e.w2);
				if (S[e.v]>=n){
					distance[e.v] = distance[heap[0]] + trans.LogNegProb(e.w2);
					heap[top]=e.v;
					int j=top++, x=(j-1)/2;
					double temp=distance[heap[j]];
					while (j>0) {
						if (distance[heap[x]]>temp){
							heap[j]=heap[x];
							if (S[heap[j]]<n) S[heap[j]]=j;
							j=x;
							x=(j-1)/2;
						}
						else break;
					}
					heap[j]=e.v;
					S[heap[j]]=j;
				}
				else{
					distance[e.v] = distance[heap[0]] + trans.LogNegProb(e.w2);
					int j=S[e.v], x=(j-1)/2;
					double temp=distance[heap[j]];
					while (j>0) {
						if (distance[heap[x]]>temp){
							heap[j]=heap[x];
							if (S[heap[j]]<n) S[heap[j]]=j;
							j=x;
							x=(j-1)/2;
						}
						else break;
					}
					heap[j]=e.v;
					S[e.v]=j;
				}
			}//endif
			
		}// end for
		heap[0]=heap[--top];
		if (!top) break;
		//siftdown
		int j=0, x=j*2+1;
		double temp=distance[heap[j]];
		while (x<top){
			if (x+1<top && distance[heap[x+1]]<distance[heap[x]]) x++;
			if (distance[heap[x]]<temp){
				heap[j]=heap[x];
				S[heap[j]]=j;
				j=x; x=j*2+1;
			}
			else break;
		}
		heap[j]=heap[top];
		if (S[heap[j]]<n) S[heap[j]]=j;
	}
	//update tree node set and heuristic
	if (!dd[node]){
		children[node]=new int[bottom];
		delta[node]=new double[bottom];
		self[node]=new double[bottom];
		path[node]=new int[bottom];
	}
	dd[node]=bottom;
	int head=0, tail=0;
	for (int i=0;i<bottom;i++){
		children[node][i]=childlist[i];
		if (numchild[childlist[i]])	self[node][i]=1;
		else {
			self[node][i]=used[childlist[i]]?1:0;
			queue[tail++]=i;
		}
		path[node][i]=-S[parent[childlist[i]]]-1;
	}
	for (int i=0;i<bottom;i++)
	{
		int child=childlist[i];
		distance[child]=longest;
		S[child]=n;
		parent[child]=-1;
	}
	
	int x = 0,u;
	while (head<tail) {
		x=queue[head++];
		u=path[node][x];
		self[node][u]*=(1-self[node][x]*b[childlist[x]]);
		if (!--numchild[childlist[u]]) {
			self[node][u]=1-self[node][u];
			queue[tail++]=u;
		}		
	}
	numchild[node]=0;
	delta[node][queue[--head]]=1;
	dp[node]+=1-self[node][x];
	for (head--;head>=0;head--) {
		x=queue[head], u=path[node][x];
		delta[node][x]=(1-self[node][u])/(1-self[node][x]*b[childlist[x]])*b[childlist[x]]*delta[node][u];
		dp[childlist[x]]+=delta[node][x]*(1-self[node][x]);
	}

	return bottom;
}

int PMIA::generateSPT_newto0(int node){
	int top=0, bottom=0;
	distance[node]=0;
	heap[0]=node;
	top++;
	bottom=dd[node];
	int head=0, tail=0;
	for (int i=0;i<bottom;i++){
		childlist[i]=children[node][i];
		numchild[i]=childnum[node][i];
		b[i]=allb[node][i];
		
		if (numchild[i])	self[node][i]=1;
		else {
			self[node][i]=used[childlist[i]]?1:0;
			queue[tail++]=i;
		}
	}
	for (int i=0;i<bottom;i++)
	{
		int child=childlist[i];
		distance[child]=longest;
		S[child]=n;
		parent[child]=-1;
	}
	int x = 0,u;
	while (head<tail) {
		x=queue[head++];
		u=path[node][x];
		self[node][u]*=(1-self[node][x]*b[x]);
		if (!--numchild[u]) {
			self[node][u]=1-self[node][u];
			queue[tail++]=u;
		}		
	}
	delta[node][queue[--head]]=1;
	dp[node]+=1-self[node][x];
	for (head--;head>=0;head--) {
		x=queue[head], u=path[node][x];
		delta[node][x]=(1-self[node][u])/(1-self[node][x]*b[x])*b[x]*delta[node][u];
		dp[childlist[x]]+=delta[node][x]*(1-self[node][x]);
	}

	return bottom;
}
int PMIA::count(Graph& gf, int node){
	ProbTransfom trans(gf.edgeForm);
	int top=0, bottom=0;
	distance[node]=0;
	heap[0]=node;
	top++;
	parent[node]=node;
	b[node]=1;
	int count=0;
	while (true){
		//stack out of heap
		if (distance[heap[0]]<longest) S[heap[0]]=-bottom-1; 
		else break;
		childlist[bottom++]=heap[0];
		if (parent[heap[0]]!=heap[0]) numchild[parent[heap[0]]]++;
		if (used[heap[0]])
			if (count++>0) return count;
		for (int i=0;i < gf.GetNeighborCount(heap[0]);i++){
			auto& e = gf.GetEdge(heap[0], i);
			if (S[e.v]<0) continue;
			if (distance[e.v]>distance[heap[0]] + trans.LogNegProb(e.w2) + EPS) {
				parent[e.v]=heap[0];
				b[e.v] = trans.Prob(e.w2);
				if (S[e.v]>=n){
					distance[e.v] = distance[heap[0]] + trans.LogNegProb(e.w2);
					heap[top]=e.v;
					int j=top++, x=(j-1)/2;
					double temp=distance[heap[j]];
					while (j>0) {
						if (distance[heap[x]]>temp){
							heap[j]=heap[x];
							if (S[heap[j]]<n) S[heap[j]]=j;
							j=x;
							x=(j-1)/2;
						}
						else break;
					}
					heap[j]=e.v;
					S[heap[j]]=j;
				}
				else{
					distance[e.v] = distance[heap[0]] + trans.LogNegProb(e.w2);
					int j=S[e.v], x=(j-1)/2;
					double temp=distance[heap[j]];
					while (j>0) {
						if (distance[heap[x]]>temp){
							heap[j]=heap[x];
							if (S[heap[j]]<n) S[heap[j]]=j;
							j=x;
							x=(j-1)/2;
						}
						else break;
					}
					heap[j]=e.v;
					S[e.v]=j;
				}
			}//endif
			
		}// end for
		heap[0]=heap[--top];
		if (!top) break;
		//siftdown
		int j=0, x=j*2+1;
		double temp=distance[heap[j]];
		while (x<top){
			if (x+1<top && distance[heap[x+1]]<distance[heap[x]]) x++;
			if (distance[heap[x]]<temp){
				heap[j]=heap[x];
				S[heap[j]]=j;
				j=x; x=j*2+1;
			}
			else break;
		}
		heap[j]=heap[top];
		if (S[heap[j]]<n) S[heap[j]]=j;
	}
	//update tree node set and heuristic
	return count;
    /*
	if (!dd[node]){
		children[node]=new int[bottom];
		delta[node]=new double[bottom];
		self[node]=new double[bottom];
		path[node]=new int[bottom];
	}
	dd[node]=bottom;
	int head=0, tail=0;
	for (int i=0;i<bottom;i++){
		children[node][i]=childlist[i];
		if (numchild[childlist[i]])	self[node][i]=1;
		else {
			self[node][i]=used[childlist[i]]?1:0;
			queue[tail++]=i;
		}
		path[node][i]=-S[parent[childlist[i]]]-1;
	}
	for (int i=0;i<bottom;i++)
	{
		int child=childlist[i];
		distance[child]=longest;
		S[child]=n;
		parent[child]=-1;
	}
	
	int x,u;
	while (head<tail) {
		x=queue[head++];
		u=path[node][x];
		self[node][u]*=(1-self[node][x]*b[childlist[x]]);
		if (!--numchild[childlist[u]]) {
			self[node][u]=1-self[node][u];
			queue[tail++]=u;
		}		
	}
	numchild[node]=0;
	delta[node][queue[--head]]=1;
	dp[node]+=1-self[node][x];
	for (head--;head>=0;head--) {
		x=queue[head], u=path[node][x];
		delta[node][x]=(1-self[node][u])/(1-self[node][x]*b[childlist[x]])*b[childlist[x]]*delta[node][u];
		dp[childlist[x]]+=delta[node][x]*(1-self[node][x]);
	}

	return bottom;
    */
}

double PMIA::Build(Graph& gf, int num, int bound)
{
	// LARGE_INTEGER start, end, freq;
    // QueryPerformanceFrequency(&freq);
	// QueryPerformanceCounter(&start);
    EventTimer timer;
    timer.SetTimeEvent("start");
    
	n = gf.GetN();
	top = num;
	Initialize(n, top);
	longest = log(double(bound));
	double treesize=0;
	S = new int[n];
	distance = new double[n];
	b = new double[n];
	heap = new int[n];
	childlist = new int[n];
	oldchildlist = new int[n];
	parent = new int[n];
	numchild = new int[n];
	queue = new int[n];
	childnum=new vector<int>[n];
	allb=new vector<double>[n];

	used.resize(n);
	lastupdate.resize(n);
	children.resize(n);
	dp.resize(n);
	self.resize(n);
	dd.resize(n);
	delta.resize(n);
	path.resize(n);
	int set[SET_SIZE];

	double old = 0.0;

	int i=0;
	for (i=0; i<n; i++)
	{
		lastupdate[i] = -1;
		children[i]=NULL;
		dp[i]=0.0;
		self[i]=NULL;
		used[i]=false;
	}
	for (i=0; i<n; i++)
		dd[i] = 0;
	for (int i=0;i<n;i++) distance[i]=longest;
	for (int i=0;i<n;i++) S[i]=n;
	for (int i=0;i<n;i++) parent[i]=-1;
	for (int i=0;i<n;i++) numchild[i]=0;

	for (i=0;i<n;i++)
	{
		double size=generateSPT_newto(gf, i);
		treesize+=size*size;
	}
    // QueryPerformanceCounter(&end);
    timer.SetTimeEvent("end");
	// double	timer = (double)(end.QuadPart - start.QuadPart) / freq.QuadPart;
    printf("%g ", timer.TimeSpan("start", "end"));


	i=0;
	double max = -1000000.0;
	int mp = 0;
	{
		int x=GetMax(i);
		set[i] = x;
		lastupdate[x] = i;			
		double improve=dp[x];
		if (improve > max) {
			max=improve;
			mp=x;
		}
	}
	used[mp] = true;
	set[i] = mp;
	list[i] = mp;
	d[i] = max;
	old+=d[i];
	generateSPT_newfrom(gf, i, mp);

	for (i=1; i<top; i++)
	{
		max = -1000000.0;
			int x=GetMax(i);
			set[i] = x;
			lastupdate[x] = i;			
			double improve=dp[x];
			if (improve > max) {
				max=improve;
				mp=x;
			}
		used[mp] = true;
		set[i] = mp;
		list[i] = mp;
		d[i] = max;
		old+=d[i];
		generateSPT_newfrom(gf, i, mp);
	}
	
    delete[] childlist;
	delete[] oldchildlist;
	delete[] distance;
	delete[] S;
	delete[] heap;
	delete[] b;
	delete[] parent;
	delete[] numchild;
	delete[] queue;

	for (i=0;i<n;i++)
		if (dd[i]) {
			delete[] children[i];
			delete[] delta[i];
			delete[] self[i];
			delete[] path[i];
			dd[i]=0;
		}
	for (i=0;i<top;i++)
		delete[] validlist[i];

	sprintf_s(file, "SPT_new_%04d.txt", bound);
	WriteToFile(file, gf);

	printf("%g ",treesize/n);
	return treesize/n;
}

double PMIA::Build(Graph& gf, int num, int k0, int bound, 
	ICascade& simuCascade,
	ICascade& simuCascadeFast)
{
	n = gf.GetN();
	top = num;
	Initialize(n, top);

	longest=log(double(bound));
	int k=k0;
	double treesize=0;
	S = new int[n];
	distance = new double[n];
	b = new double[n];
	heap = new int[n];
	childlist = new int[n];
	oldchildlist = new int[n];
	parent = new int[n];

	used.resize(n);
	lastupdate.resize(n);
	children.resize(n);
	dp.resize(n);
	self.resize(n);
	dd.resize(n);
	int set[SET_SIZE];

	double old = 0.0;

	int i=0;
	for (i=0; i<n; i++)
	{
		lastupdate[i] = -1;
		children[i]=NULL;
		dp[i]=0.0;
		self[i]=NULL;
	}
	for (i=0; i<n; i++)
		dd[i] = 0;
	for (i=0;i<n;i++)
	{
		treesize+=generateSPT_newto(gf, i);
	}


	double max = -1000000.0;
	int mp = 0;
	for (int j=0;j<k;j++){
		int x=GetMax(i);
		set[i] = x;
		lastupdate[x] = i;			
		double improve = simuCascadeFast.Run(NUM_ITER, i + 1, set);
		if (improve > max) {
			max=improve;
			mp=x;
		}
	}
	used[mp] = true;
	set[i] = mp;
	list[i] = mp;
	d[i] = max;
	old+=d[i];
	generateSPT_newfrom(gf, i, mp);

	for (i=1; i<top; i++)
	{
		max = -1000000.0;
		for (int j=0;j<k;j++){
			int x=GetMax(i);
			set[i] = x;
			lastupdate[x] = i;
			double improve = (i>top) ?
				simuCascade.Run(NUM_ITER/100, i+1, set) - old : 
				simuCascadeFast.Run(NUM_ITER/100, i+1, set) - old;
			if (improve > max) {
				max=improve;
				mp=x;
			}
		}
		used[mp] = true;
		set[i] = mp;
		list[i] = mp;
		d[i] = max;
		old+=d[i];
		generateSPT_newfrom(gf, i, mp);
		// update trees.
	}
	delete[] childlist;
	delete[] oldchildlist;
	delete[] distance;
	delete[] S;
	delete[] heap;
	delete[] b;
	delete[] parent;

	for (i=0;i<n;i++)
		if (dd[i]) {
			delete[] children[i];
			delete[] delta[i];
			delete[] self[i];
			delete[] path[i];
		}
	for (i=0;i<top;i++)
		delete[] validlist[i];

	sprintf_s(file,"SPT_new_%04d.txt", bound);
	WriteToFile(file, gf);

	return treesize/n;
}

void PMIA::BuildFromFile(Graph& gf, int bound)
{
	sprintf_s(file, "SPT_new_%04d.txt", bound);
	ReadFromFile(file, gf);
	Initialize(n, top);
}

char* PMIA::filename(int bound)
{
	sprintf_s(file,"SPT_new_%04d.txt", bound);
	return file;
}

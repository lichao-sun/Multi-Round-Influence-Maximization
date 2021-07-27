#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include "topic_aware.h"
#include "graph.h"

using namespace std;

int TopicAwareBase::ntopics = -1;
int TopicAwareBase::nlandmarks = -1;
TopicAwareBase::TASeedsMap TopicAwareBase::topicMap;
vector< vector<TASeeds> > TopicAwareBase::topicVector;

int TopicAwareBase::RoundSeedsToIndex(std::vector<TASeeds>& seeds, double vlambda)
{
	// find the right seeds:
	int sz = (int)seeds.size();
	if (sz <= 0) { return -1; }

	if (vlambda <= EPS) {
		return -1;
	} else if ( vlambda >= seeds[sz-1].vlambda - EPS ) {
		return sz-1;
	} 

	for (int i = 1; i < sz; ++i) {
		if (vlambda < seeds[i].vlambda  - EPS) {
			return i-1;
		}
	}
	// [Special] if it is too small, and less than [0], we return index 0.
	return 0;
}

void TopicAwareBase::ReadSeeds(std::string file, TASeeds& ta, IGraph& gf)
{
	SeedIO io;
	std::vector<double> infl;
	std::vector<int> seeds = io.Read(file, infl, gf);
	for (size_t i = 0; i < seeds.size(); i++) {
		MINode temp;
		temp.num = seeds[i];
		temp.influence = infl[i];
		ta.seeds.push_back(temp);
	}
}

void TopicAwareBase::ReadTopicFile(std::string seeds_list_file, IGraph& gf)
{
	vector<int> topics;
	vector<double> vlambdas;
	vector<string> files;

	ifstream flist(seeds_list_file.c_str());
	string line;

	getline(flist, line);
	stringstream header(line);
	header >> TopicAwareBase::ntopics >> TopicAwareBase::nlandmarks;

	int counter = 0;
	while ( getline(flist, line) )
	{
		stringstream stream(line);
		int topic;
		double vlambda;
		string file;
		stream >> topic >> vlambda >> file;
		topics.push_back(topic);
		vlambdas.push_back(vlambda);
		files.push_back(file);
		counter++;
	}

	TASeedsMap& tmap = TopicAwareBase::topicMap;
	// topic index start from 1; initialize
	for (int t = 1; t <= TopicAwareBase::ntopics; ++t) 
	{
		tmap[t] = vector<TASeeds>();
	}

	for (int i = 0; i < counter; ++i) 
	{
		vector<TASeeds>& vec = tmap[topics[i]];

		TASeeds s;
		s.vlambda = vlambdas[i];
		TopicAwareBase::ReadSeeds(files[i], s, gf);
		
		vec.push_back(s);
		/*cout << topics[i] << "\t" 
			<< vlambdas[i] << "\t"
			<< files[i] << endl;
		for (int s = 0; s < seeds.size(); ++s) {
			cout << seeds[s].num << " " << seeds[s].influence << endl;
		}*/
	}

	// sort by vlambda:
	for (int t = 1; t <= TopicAwareBase::ntopics; ++t)
	{
		vector<TASeeds>& vec = tmap[t];
		sort(vec.begin(), vec.end(), TASeedLambdaSmaller());
		
		//cout << t << "\t" << vec.size() << endl;
		//for (int i = 0; i < vec.size(); i++) {
		//	cout <<  "  " << vec[i].vlambda << endl;
		//}
	}
	for (int t = 1; t <= TopicAwareBase::ntopics; ++t) {
		topicVector.push_back( tmap[ t ] );
	}
}

std::set<int> TopicAwareBase::GetAllSeedsIndices()
{
	set<int> allseeds;
	for (size_t j = 1; j < topicVector.size(); j++) {
		vector<TASeeds>& s = topicVector[j];
		for (size_t i = 0; i < s.size(); i++) {
			TASeeds& sss = s[i];
			for (size_t k = 0; k < sss.seeds.size(); k++) {
				allseeds.insert(sss.seeds[k].num);
			}
		}
	}
	return allseeds;
}

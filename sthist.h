
#ifndef _STHIST_H
#define _STHIST_H

#include "histogram.h"
#include "glibwrap.h"

#include<vector>
#include<forward_list>

using namespace std;

struct coordenada{
	double x;
	double y;
};

struct{
        Envelope mbr;
	unsigned int F;
        forward_list<coordenada> O;
        vector<int> filhos;
}hotSpot;

typedef struct hotSpot STBucket;

typedef vector<STBucket> hotSpotTree;
typedef hotSpotTree STHist;

STHist * STHist_generate(dataset *ds, int buckets_num);

double STHist_search(STHist *hist, Envelope query);

void STHist_print(STHist *hist);

#endif

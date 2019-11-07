
#ifndef _STHIST_H
#define _STHIST_H

#include "histogram.h"
#include "glibwrap.h"

#include<vector>
#include<forward_list>

struct coordenada{
	double x;
	double y;
}

struct{
        Envelope mbr;
	usingned int F;
        forward_list<coordenada> O;
        vector<int> filhos;
}hotSpot;
hotSpot typedef STBucket;

vector<STHist_bucket>> typedef hotSpotTree;
hostSpotTree typedef STHist;

STHist * STHist_generate(dataset *ds, int buckets_num);

double STHist_search(STHist *hist, Envelope query);

void STHist_print(STHist *hist);

#endif

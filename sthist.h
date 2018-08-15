
#ifndef _STHIST_H
#define _STHIST_H

#include "histogram.h"
#include "glibwrap.h"

typedef struct {
	Envelope mbr;
	double cardin;
} st_hist_bucket;

st_hist_bucket *st_generate_hist(dataset *ds, int buckets_num/*****/);

double st_search_hist(st_hist_bucket *hist, Envelope query);

void st_print_hist(dataset *ds, st_hist_bucket *hist);

#endif


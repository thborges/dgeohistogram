
#ifndef _MINSKEW_H
#define _MINSKEW_H

#include "histogram.h"
#include "glibwrap.h"

typedef struct {
	double skew;
	char split_axis;
	short split_point;
	double skew_reduction;
	Envelope mbr;
	double cardin;
} minskew_bucket;

typedef struct {
	double mean;
	double variance;
	int n;
} variance_result;

GList *minskew_generate_hist(dataset *ds, int buckets_num);

double minskew_search_hist(GList *hist, Envelope query);

void minskew_print_hist(dataset *ds, GList *hist);

variance_result calculate_skew_row_col(dataset_histogram *dh, 
	int xini, int xfim, int yini, int yfim);

#endif


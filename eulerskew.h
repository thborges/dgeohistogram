
#ifndef _EULERKEW_H
#define _EULERKEW_H

#include "histogram.h"
#include "glibwrap.h"


typedef struct {
} eulerskew_edge;

typedef struct {
} eulerskew_vertex;

typedef struct {
	double skew;
	char split_axis;
	short split_point;
	double skew_reduction;

	Envelope mbr;
	double cardin;
} eulerskew_face;


typedef struct {
	GList *faces;
	GList *edges;
	GList *vertex;
} eulerskew_hist;


//GList *eulerskew_generate_hist(eh, int buckets_num);

#endif

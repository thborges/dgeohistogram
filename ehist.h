
#ifndef _STHIST_H
#define _STHIST_H

#include "histogram.h"
#include "glibwrap.h"

typedef struct {
	double cardin;
} euler_face;

typedef struct {
	Envelope mbr;
	double cardin;
} euler_edge;

typedef struct {
	double x;
	double y;
	double cardin;
} euler_vertex;

typedef	struct {
	Envelope mbr;
	int xqtd;
	int yqtd;
	double xsize;
	double ysize;
	double *xtics;
	double *ytics;
	euler_face *faces;
	euler_edge *edges;
	euler_vertex *vertexes;
} euler_histogram;

euler_histogram *eh_generate_hist(dataset *ds, HistogramGenerateSpec spec, enum JoinPredicateCheck pcheck); 

euler_histogram *eh_generate_intermed(dataset *ds,euler_histogram *ehA,euler_histogram *ehB);

int euler_search_hist(euler_histogram *eh, Envelope query);
void euler_print_hist(char *name, euler_histogram *eh);


#endif



#ifndef _STHIST_H
#define _STHIST_H

#include "histogram.h"
#include "glibwrap.h"

typedef struct {
	double cardin;
    double avg_heigth;
    double avg_width;
    double avg_area;
} euler_face;

typedef struct {
	Envelope mbr;
	double cardin;
    double avg_projection;
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

int euler_search_hist(euler_histogram *eh, Envelope query);
void euler_print_hist(dataset *ds, euler_histogram *eh);

int euler_spatial_join(euler_histogram* ehr, euler_histogram* ehs);
#endif


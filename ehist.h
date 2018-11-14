
#ifndef _STHIST_H
#define _STHIST_H

#include "histogram.h"
#include "glibwrap.h"
#include "rtree.h"
#include "rtree-star.h"


#define GET_VERT_EDGE(x, y) ((x == eh->xqtd) ? (x * (2*eh->yqtd+1) + y) : (x * (2*eh->yqtd+1) + 2*y + 1))
#define GET_HORZ_EDGE(x, y) (x * (2*eh->yqtd+1) + 2*y)

#define GET_VERT_EDGE_EHR(x, y) ((x == ehr->xqtd) ? (x * (2*ehr->yqtd+1) + y) : (x * (2*ehr->yqtd+1) + 2*y + 1))
#define GET_HORZ_EDGE_EHR(x, y) (x * (2*ehr->yqtd+1) + 2*y)
#define GET_VERT_EDGE_EHS(x, y) ((x == ehs->xqtd) ? (x * (2*ehs->yqtd+1) + y) : (x * (2*ehs->yqtd+1) + 2*y + 1))
#define GET_HORZ_EDGE_EHS(x, y) (x * (2*ehs->yqtd+1) + 2*y)

typedef struct {
	double cardin;
    double avg_height;
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
int euler_join_cardinality(dataset *dr, 
        dataset *ds,
        euler_histogram* ehr, 
        euler_histogram* ehs,
        rtree_root* rtree_r,
        rtree_root* rtree_s,
        double* stddev);
#endif


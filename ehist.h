
#ifndef _STHIST_H
#define _STHIST_H

#include "histogram.h"
#include "glibwrap.h"
#include "rtree.h"
#include "rtree-star.h"

typedef struct {
	double cardin;
	double cardin_prop;
    double avg_height;
    double avg_width;
    double avg_area;
    double areasum;
    Envelope usedarea;
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

euler_histogram *eh_generate_hreal(dataset *ds,dataset *dsb,euler_histogram *ehA,enum JoinPredicateCheck pcheck);

euler_histogram *eh_generate_intermed(dataset *ds,dataset *dsb,euler_histogram *ehA,euler_histogram *ehB);
euler_histogram *eh_generate_intermed_real(dataset *ds,dataset *dsb,dataset *dsc,euler_histogram *ehA,euler_histogram *ehB,euler_histogram *ehC);

double estimate_intersections_mp_edges_vert(Envelope el, Envelope er, Envelope inters,
        euler_edge *ehr_face, euler_edge *ehs_face);
double estimate_intersections_mp_edges_horz(Envelope el, Envelope er, Envelope inters,
		euler_edge *ehr_face, euler_edge *ehs_face);
double estimate_intersections_mamoulis_papadias(Envelope el, Envelope er, Envelope inters,
		euler_face *ehr_face, euler_face *ehs_face);

double estimate_intersections_zu(Envelope el, Envelope er,
Envelope inters, euler_face *ehl_face, euler_face *ehr_face,
dataset *dh_l, dataset *dh_r);


int euler_search_hist(euler_histogram *eh, Envelope query);
void euler_print_hist(char *name, euler_histogram *eh);
int euler_spatial_join(euler_histogram* ehr, euler_histogram* ehs);

#endif


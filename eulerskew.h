
#ifndef _STEULERSKEW_H
#define _STEULERSKEW_H

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
    Envelope mbr;
	double cardin;
    double avg_projection;
} eulerskew_edge;

typedef struct {
    double x;
	double y;
	double cardin;
} eulerskew_vertex;

typedef struct {
	double skew;
	char split_axis;
	short split_point;
	double skew_reduction;
    double avg_height;
    double avg_width;
    double avg_area;
	Envelope mbr;
	double cardin;
} eulerskew_face;

//eulerskew_face é o bucket

typedef	struct {
	Envelope mbr;
	int xqtd;
	int yqtd;
	double xsize;
	double ysize;
	double *xtics;
	double *ytics;
	eulerskew_face *faces;
	eulerskew_edge *edges;
	eulerskew_vertex *vertexes;
} eulerskew_histogram;

typedef struct {
	double mean;
	double variance;
	int n;
} variance_result_eulerskew;

void eulerskew_get_ini_fim(dataset_histogram *dh, Envelope ev, int *xini, int *xfim, int *yini, int *yfim);

variance_result_eulerskew eulerskew_calculate_skew_row_col(dataset_histogram *dh, int xini, int xfim, int yini, int yfim);

void eulerskew_calculate_bucket_with_mbr(dataset_histogram *dh, eulerskew_face *b);

void eulerskew_calculate_skew_reduction(dataset_histogram *dh, eulerskew_face *bucket);

GList *eulerskew_generate_hist(dataset *ds, int buckets_num);

//void eulerskew_alloc(dataset *ds, eulerskew_histogram *eh, int xqtd, int yqtd, double psizex, double psizey);
void eulerskew_alloc(dataset *ds, eulerskew_histogram *eh, int xqtd, int yqtd, double psizex, double psizey, GList *minskewhist);

void eulerskew_hash_ds_objects(dataset *ds, eulerskew_histogram *eh, enum JoinPredicateCheck pcheck);

//void eulerskew_generate_hw(dataset *ds, eulerskew_histogram *eh, double x, double y, enum JoinPredicateCheck pcheck);
void eulerskew_generate_hw(dataset *ds, eulerskew_histogram *eh, double x, double y, enum JoinPredicateCheck pcheck, GList *minskewhist) ;

//void eulerskew_generate_fix(dataset *ds, eulerskew_histogram *eh, int fsizex, int fsizey, enum JoinPredicateCheck pcheck);
void eulerskew_generate_fix(dataset *ds, eulerskew_histogram *eh, int fsizex, int fsizey, enum JoinPredicateCheck pcheck, GList *minskewhist);

eulerskew_histogram *eulerskew_generate_hist_with_euler(dataset *ds, HistogramGenerateSpec spec, enum JoinPredicateCheck pcheck, GList *minskewhist);
//eulerskew_histogram *eulerskew_generate_hist_with_euler(dataset *ds, HistogramGenerateSpec spec, enum JoinPredicateCheck pcheck);

int eulerskew_search_hist(eulerskew_histogram *eh, Envelope query2);

void eulerskew_print_hist(dataset *ds, eulerskew_histogram *eh);

double eulerskew_estimate_intersections_mp_edges_vert(Envelope el, Envelope er, Envelope inters,
        eulerskew_edge *ehr_face, eulerskew_edge *ehs_face);

double eulerskew_estimate_intersections_mp_edges_horz(Envelope el, Envelope er, Envelope inters,
        eulerskew_edge *ehr_face, eulerskew_edge *ehs_face);

double eulerskew_estimate_intersections_mamoulis_papadias(Envelope el, Envelope er, Envelope inters,
        eulerskew_face *ehr_face, eulerskew_face *ehs_face);

int eulerskew_join_cardinality(dataset *dr,
        dataset *ds,
        eulerskew_histogram* ehr,
        eulerskew_histogram* ehs,
        rtree_root* rtree_r,
        rtree_root* rtree_s,
        double* stddev);


int eulerskew_spatial_join(eulerskew_histogram* ehr, eulerskew_histogram* ehs);

int eulerskew_cardinality_per_face(dataset *dr,
        dataset *ds,
        eulerskew_histogram* ehr,
        eulerskew_histogram* ehs,
        rtree_root* rtree_r,
        rtree_root* rtree_s);

int real_cardin_eulerskew_histogram_cell(rtree_root* rtree_r, rtree_root* rtree_s, Envelope inters);




//double eulerskew_search_hist(GList *hist, Envelope query);
//qual?
//no main.c a função search é chamada pra estimar a cardinalidade das intersecções dos poligonos com os buckets, a contagem multipla será feita atraves do metodo de euler
//possibilidade: gerar o histograma por minskew, utilizar as funções do histograma de euler substituindo as edges pelo bucket, certificando que o mbr e as variaveis necessarias possuem os mesmos valores.
//É necessario analisar a função de contagem de objetos interceptados, do euler, e as funções necessarias do minskew




#endif

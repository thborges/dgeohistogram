
#ifndef JOINCOMMON_H
#define JOINCOMMON_H

enum JoinPredicateCheck { CHECKR, CHECKS };
#define get_join_pair_leaf(ir, pcheck) pcheck == CHECKR ? &ir[0] : &ir[1]

typedef struct {
	double to_pnts;
	double io_objs;
	double io_pnts;
	double io_save;
	double inters;
} multiway_histogram_estimate;

void multiway_print_estimate(char *name, multiway_histogram_estimate *estimate, int servers);

void multiway_totalize_estimate(multiway_histogram_estimate *estimate, int servers,
	double *totalpnts, double *totalcomm, 
	double *mkspan, double *max_comm,
	double *stdev_mkspan, double *stdev_comm,
	double *mkspan_gap);

#endif


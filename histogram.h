/*
 * histogram.h
 *
 *  Created on: 09/08/2014
 *      Author: thborges
 */

#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "glibwrap.h"
#include "dataset.h"
#include "geosext.h"

#define GET_HISTOGRAM_CELL(h, x, y) &(h)->hcells[(x)*(h)->yqtd + (y)]
#define SET_IN_PLACE(var, place) (var = var | 1<<(place-1))
#define IS_IN_PLACE(var, place) (var & 1<<(place-1))

typedef struct {
	double to_pnts;
	double io_objs;
	double io_pnts;
	double io_save;
	double inters;
} multiway_histogram_estimate;

void histogram_generate(dataset *ds, enum JoinPredicateCheck pcheck);
void histogram_distribute(dataset *ds);
void histogram_print(dataset *ds, histogram_type type);
void histogram_print_geojson(dataset *ds);
int histogram_join_cardinality(dataset *dr, dataset *ds);
void histogram_alloc(dataset_histogram *dh, int xqtd, int yqtd);
void histogram_print_estimative(char *name, multiway_histogram_estimate *estimate, int servers);

#ifdef __cplusplus
}
#endif


#endif /* HISTOGRAM_H_ */

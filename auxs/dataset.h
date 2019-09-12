/*
 * dataset.h
 *
 *  Created on: 26/07/2014
 *      Author: thborges
 */

#ifndef DATASET_H_
#define DATASET_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "uthash.h"
#include <ogrext.h>
#include <stdbool.h>
#include <pthread.h>

#define DATASET_NAME_MAX 20
#define DATASET_PAGE_SIZE 4096
#define DATASET_CACHE_MAX 10
#define DATASET_HASH_FULL(dh, seg) !dh->metadata.memory_dataset && HASH_COUNT(seg->page_cache) > DATASET_CACHE_MAX
#define DATASET_PAGE_CAPACITY(dh) (DATASET_PAGE_SIZE / (dh->metadata.geocount * sizeof(dataset_leaf)) -1)

enum JoinPredicateCheck { CHECKR, CHECKS };
#define get_join_pair_leaf(ir, pcheck) pcheck == CHECKR ? &ir[0] : &ir[1]


#define AVGLENGTH_FIND_BUCKET(olength,dim_length) ((olength) / ((dim_length) / AVGL_HISTO_DIV))

typedef GEOSGeometry* GEOSGeometryH;

struct leaf {
	bool cloned;
	long long gid;
	unsigned points;
	Envelope mbr;
	//GEOSGeometryH chull;
	GEOSGeometryH geo;
	const GEOSPreparedGeometry *pgeo;
};

typedef struct leaf dataset_leaf;

typedef enum {CLEAN, DIRTY} dataset_page_wstate;
typedef enum {LOCKED, UNLOCKED} dataset_page_state;

typedef struct {
	union {
		unsigned used;
		void *__; // align
		// check DATASET_PAGE_CAPACITY after adding fields here!
	};
	dataset_leaf leaves[1];
} dataset_page_data;

typedef struct {
	unsigned id;
	dataset_page_wstate wstate;
	dataset_page_state state;
	UT_hash_handle hh;
	union {
		dataset_page_data *data;
		unsigned char *memory;
	};
} dataset_page;

typedef struct _dataset_segment_item_ {
	unsigned pid;
	struct _dataset_segment_item_ *next;
} dataset_segment_item;

typedef struct {
	unsigned sid;
	unsigned count;
	unsigned itemcount;

	dataset_segment_item *fillingitem;
	dataset_page *fillingpage;

	dataset_segment_item *items;
	dataset_page *page_cache;

	UT_hash_handle hh;
} dataset_segment;

typedef enum { CARDIN, POINTS, PLACES } histogram_type;

typedef struct {
	double cardin;
	double points;
	double avgwidth;
	double avgheight;
	unsigned place;
	unsigned copies;
    int objcount;
	float areasum;
    //void *extra_data;
} histogram_cell;


double OLD_AVG_LENGTH_X(histogram_cell *cell);
double OLD_AVG_LENGTH_Y(histogram_cell *cell);


typedef	struct dataset_histogram{
	Envelope mbr;
	int xqtd;
	int yqtd;
	double xsize;
	double ysize;
	double *xtics;
	double *ytics;
	//OGRwkbGeometryType geom_type;
    //void *extra_data;
	histogram_cell *hcells;
} dataset_histogram;

#define dataset_histogram_size sizeof(dataset_histogram)

typedef struct dataset_head_ {
	struct {
		char *name;
		double x_average;
		double y_average;
		double x_psa; //x power sum average for stddev
		double y_psa; //y power sum average for stddev
		unsigned count;
		unsigned geocount;
		unsigned pagecount;
		char memory_dataset;
		char cloned;
		int servers;
		dataset_histogram hist;
	} metadata;

	// temp data
	OGRLayerH temp_ogr_layer;
	int fpages;
	dataset_segment *segments;
	pthread_mutex_t *mem_get_lock;
} dataset;

typedef struct {
	dataset *dh;
    dataset_segment *seg;
	dataset_segment_item *si;
	dataset_page *page;
	unsigned position;
	dataset_leaf *item;
    int fpages;
} dataset_iter_seg, dataset_iter;

dataset*		dataset_create(const char *name, unsigned short geocount);
dataset*		dataset_create_mem(const char *name, unsigned short geocount);
void			dataset_destroy(dataset *dh);
void			dataset_destroy_full(dataset *dh, char destroy_geos);
dataset_leaf*	dataset_add(dataset *dh);
void			dataset_concat(dataset *dh, dataset *di);
dataset_iter	dataset_first(dataset *dh);
void			dataset_fetch(dataset_iter *i);
void			dataset_next(dataset_iter *i);
dataset*		dataset_clone(dataset *dh);
void                    dataset_seg_sync(dataset *dh, dataset_segment *seg);
dataset*		dataset_clone_seg(dataset *dh, dataset_segment *seg, dataset_segment **clone_seg);
void			dataset_fill_leaf_ext(dataset_leaf *leaf, int index, long long gid, bool cloned, GEOSGeometryH ggeo, Envelope *e, GEOSPreparedGeometry *pgeo);
void			dataset_fill_leaf_id(dataset_leaf *leaf, int index, long long gid, Envelope *e);
void			dataset_fill_leaf(dataset_leaf *leaf, int index, GEOSGeometryH ggeo, long long gid, Envelope *e, GEOSPreparedGeometry *pgeo);
void			dataset_set_histogram(dataset *dh, dataset_histogram *hist);
dataset_histogram *dataset_get_histogram(dataset *ds);

dataset_segment* dataset_new_seg(dataset *dh, unsigned sid);
void             dataset_destroy_seg(dataset *dh, dataset_segment *seg);
dataset_segment* dataset_get_seg(dataset *dh, unsigned sid);
dataset_leaf*    dataset_add_seg(dataset *dh, dataset_segment *seg);
dataset_iter_seg dataset_first_seg(dataset *dh, dataset_segment *seg);
void             dataset_fetch_seg(dataset_iter_seg *i);
void             dataset_next_seg(dataset_iter_seg *i);
void             dataset_release_iter(dataset_iter_seg *i);

GEOSGeometryH    dataset_get_leaf_geo(dataset *dh, dataset_leaf *leaf);

void             dataset_write_csv(dataset *dh, FILE *f);

#define	dataset_foreach(var, dataset) dataset_foreach_seg(var, dataset, dataset_get_seg(dataset, 0))

#define	dataset_foreach_seg(var, dh, segment)		\
	for((var) = dataset_first_seg((dh), (segment));	\
	   (var).item;								\
	   dataset_next_seg(&(var)))

#define dataset_is_last(i) ((i.si->pid+1 == i.dh->metadata.pagecount) && (i.page->data->used == i.position+1))
#define dataset_iter_is_same(i, j) (i.dh == j.dh && i.si->pid == j.si->pid && i.position == j.position)

#define dataset_meta_stddev(m, axis) (sqrt(m.axis##_psa / (m.count)))

#ifdef __cplusplus
}
#endif

#endif /* DATASET_H_ */

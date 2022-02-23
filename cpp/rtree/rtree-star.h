
#ifndef RTREE_STAR_H
#define RTREE_STAR_H

#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>
#include "rtree.h"
#include "rtree-reinsert.h"

// R0-tree definitions used on reinsert
#define ENV_WIDTH(r) (r.MaxX - r.MinX)
#define ENV_HEIGHT(r) (r.MaxY - r.MinY)
#define R0_QUALITY(r) ((1.0/(ENV_WIDTH(r)*ENV_HEIGHT(r))) * powf( (MIN(ENV_WIDTH(r),ENV_HEIGHT(r))/MAX(ENV_WIDTH(r),ENV_HEIGHT(r))), 0.5))
#define R0_GAIN(r1, r2) (1-(R0_QUALITY(r1)/R0_QUALITY(r2)))
#define R0_LOSS(r1, r2) R0_GAIN(r2,r1)

void planesweep_split(rtree_leaf *node, const int m);

int sortleafx_min(const void *x, const void *y);

int rtree_choose_subtree_star(rtree_root *root, rtree_node *n, const EnvelopeC e);

int rtree_split_star_choose_index_leaf(rtree_leaf *leaves, int m, int usedtotal);
int rtree_split_star_choose_index(rtree_node **dirs, int m, int usedtotal);

rtree_node **rtree_split_star_choose_axis(rtree_node **dirsx, rtree_node **dirsy, int m, int usedtotal, int *minperimeter);
rtree_leaf *rtree_split_star_choose_axis_leaf(rtree_leaf *leavesx, rtree_leaf *leavesy, int m, int usedtotal, int *minperimeter);

rtree_node *rtree_split_star_leaf(rtree_root *root, rtree_node *n, rtree_node *parent, const GEOSGeometryH ngeo, const EnvelopeC e);
rtree_node *rtree_split_star(rtree_root *root, rtree_node *n, rtree_node *newn, rtree_node *parent, int m);

rtree_root *rtree_new_rstar(const int m, const int servers);

#ifdef __cplusplus
}
#endif

#endif

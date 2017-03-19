
#ifndef RTREE_STAR_H
#define RTREE_STAR_H

#ifdef __cplusplus
extern "C" {
#endif

#include "rtree.h"
#include "ogrext.h"
#include <assert.h>

void planesweep_split(rtree_leaf *node, const int m);

int sortleafx_min(const void *x, const void *y);

int rtree_choose_subtree_star(rtree_root *root, rtree_node *n, const Envelope e);

int rtree_split_star_choose_index_leaf(rtree_leaf *leaves, int m, int usedtotal);
int rtree_split_star_choose_index(rtree_node **dirs, int m, int usedtotal);

rtree_node **rtree_split_star_choose_axis(rtree_node **dirsx, rtree_node **dirsy, int m, int usedtotal, int *minperimeter);
rtree_leaf *rtree_split_star_choose_axis_leaf(rtree_leaf *leavesx, rtree_leaf *leavesy, int m, int usedtotal, int *minperimeter);

rtree_node *rtree_split_star_leaf(rtree_root *root, rtree_node *n, rtree_node *parent, const GEOSGeometryH ngeo, const Envelope e);
rtree_node *rtree_split_star(rtree_root *root, rtree_node *n, rtree_node *newn, rtree_node *parent, int m);

rtree_root *rtree_new_rstar(const int m, const int servers);

#ifdef __cplusplus
}
#endif

#endif

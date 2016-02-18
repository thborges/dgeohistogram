
#ifndef RTREE_GUT_H
#define RTREE_GUT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "rtree.h"
#include "ogrext.h"
#include <assert.h>

int rtree_choose_subtree_gut(rtree_root *root, rtree_node *n, const Envelope e);

rtree_node *rtree_split_gut_quad(rtree_root *root, rtree_node *n, rtree_node *newn, rtree_node *parent, int m);

rtree_root *rtree_new_rtree_gut_quad(const int m, const int servers);


int rtree_choose_subtree_linear(rtree_root *root, rtree_node *n, const Envelope e);

rtree_node *rtree_split_gut_linear(rtree_root *root, rtree_node *n, rtree_node *newn, rtree_node *parent, int m);

rtree_root *rtree_new_rtree_gut_linear(const int m, const int servers);





int rtree_choose_subtree_gut_corner(rtree_root *root, rtree_node *n, const Envelope e);
rtree_node *rtree_split_corner(rtree_root *root, rtree_node *n, rtree_node *newn, rtree_node *parent, int m);

rtree_node *rtree_split_corner_leaf(rtree_root *root, rtree_node *n, rtree_node *parent, const GEOSGeometryH ngeo, const Envelope e);


rtree_root *rtree_new_rtree_corner(const int m, const int servers);




#ifdef __cplusplus
}
#endif

#endif

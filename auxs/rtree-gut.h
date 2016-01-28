
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

rtree_root *rtree_new_rtree_gut_quad(const int m, const int servers);

#ifdef __cplusplus
}
#endif

#endif

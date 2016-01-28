#include "rtree.h"

#ifndef RTREE_TEST_H
#define RTREE_TEST_H

#ifdef __cplusplus
extern "C" {
#endif

count_node_struct count_nodes(rtree_root *rtree);

void check_mbrs(rtree_root *rtree);

void print_overlap(rtree_root *rtree, OGRSpatialReferenceH refspatial);

void write_dirs_shapefile(rtree_root *rtree);

#ifdef __cplusplus
}
#endif

#endif


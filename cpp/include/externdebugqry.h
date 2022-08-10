
#pragma once
#include <rtree.h>

#define DEBUG_QUERYNO
#ifdef DEBUG_QUERYNO
extern int global_qryno;
extern double global_query_size;
extern int global_debug_qryno;
extern double global_debug_qrysize;
extern rtree_root *global_search_rtree;
#endif

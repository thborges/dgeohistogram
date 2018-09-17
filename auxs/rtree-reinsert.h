#ifndef RTREE_REINSERT_H
#define RTREE_REINSERT_H

void rstar_reinsert_on_parent(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g, const Envelope gmbr);
void border_reinsert_on_parent(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g, const Envelope gmbr);
void greedy_reinsert_on_parent(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g, const Envelope gmbr);
void exponen_reinsert_on_parent(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g, const Envelope gmbr);
void rlazy_reinsert_on_parent(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g, const Envelope gmbr);

#endif


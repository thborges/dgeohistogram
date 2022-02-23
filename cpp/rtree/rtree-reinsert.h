#ifndef RTREE_REINSERT_H
#define RTREE_REINSERT_H

void rstar_reinsert_on_parent(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g, const EnvelopeC gmbr);
void border_reinsert_on_parent(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g, const EnvelopeC gmbr);
void greedy_reinsert_on_parent(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g, const EnvelopeC gmbr);
void exponen_reinsert_on_parent(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g, const EnvelopeC gmbr);
void rlazy_reinsert_on_parent(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g, const EnvelopeC gmbr);

#endif


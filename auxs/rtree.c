#include <float.h>
#include <math.h>
#include <assert.h>
#include "glibwrap.h"
#include "rtree.h"
#include "rtree-distributor.h"
#include "ogrext.h"

int nodeindex = 0;
lrubuffer *lru = NULL;

rtree_node *rtree_new_node(rtree_root *root, const enum NodeTypeEnum type) {
	rtree_node *n = g_new0(rtree_node, 1);
	ENVELOPE_INIT(n->mbr);
	n->index = nodeindex++;
	n->server = 0; 
	n->type = type;
	n->used = 0;
	if (type == LEAF)
		n->leaves = g_new0(rtree_leaf, root->m);
	else
		n->dirs = g_new0(rtree_node*, root->m);

	//Where to put the node?
	rtree_distributor_rr(root, n);

	return n;
}

rtree_root *rtree_new(
	const int m,
	const int servers, 
	int         (*choose_subtree_func)   (rtree_root *, rtree_node *, const Envelope e),
	rtree_node *(*rtree_split_func)      (rtree_root *, rtree_node *, rtree_node *, rtree_node *, int),
	rtree_node *(*rtree_split_leaf_func) (struct rtree_root_ *, rtree_node *, rtree_node *, const GEOSGeometryH, const Envelope e),
	void        (*reinsert_on_parent)    (rtree_node *, struct rtree_root_ *, GList **, const GEOSGeometryH, const Envelope),
	void        (*rtree_redivide)        (rtree_node *, rtree_node *, const int),
	void        (*rtree_redivide_leaf)   (rtree_node *, rtree_node *, const int)

) {
	rtree_root *r = g_new(rtree_root, 1);
	r->m = m;
	r->servers = servers;
	r->rr_next_server = 0;
	r->root = rtree_new_node(r, LEAF);
	r->choose_subtree_func = choose_subtree_func;
	r->rtree_split_func = rtree_split_func;
	r->rtree_split_leaf_func = rtree_split_leaf_func;
	r->reinsert_on_parent = reinsert_on_parent;
	r->rtree_redivide = rtree_redivide;
	r->rtree_redivide_leaf = rtree_redivide_leaf;
	return r;
}

Envelope rtree_compute_mbr(const rtree_node *n) {
	Envelope e;
	int i;
	if (n->type == DIRECTORY) {
		e = n->dirs[0]->mbr;
		for(i=1; i < n->used; i++) {
			ENVELOPE_MERGE(e, n->dirs[i]->mbr);
		}
	}
	else {
		e = n->leaves[0].mbr;
		for(i=1; i < n->used; i++) {
			ENVELOPE_MERGE(e, n->leaves[i].mbr);
		}
	}
	return e;
}

rtree_node *rtree_split(rtree_root *root, rtree_node *n, rtree_node *newn, rtree_node *parent, int m) {
	return root->rtree_split_func(root, n, newn, parent, m);
}

int sort_node_used(const void *x, const void *y) {
	rtree_node *iy = *(rtree_node**)y;
	rtree_node *ix = *(rtree_node**)x;
	if (ix->used == iy->used) return 0;
	if (ix->used > iy->used) return 1;
	return -1;
}

rtree_node *rtree_append_aux(rtree_root *root, rtree_node *n, GList **parents, const GEOSGeometryH g, const Envelope gmbr, rtree_node *at_reinsert) {

	if (n->type == LEAF) {
		if (n->used < root->m) {
			//if (n->used < minm(root->m))
			//	printf("Outlier %04d insert: %02d used\n", n->index, n->used);

			n->leaves[n->used].mbr = gmbr;
			//n->leaves[n->used].chull = NULL;//OGR_G_ConvexHull(g);
			n->leaves[n->used].geo = g;
			n->leaves[n->used].pgeo = NULL;
			n->used++;
			n->mbr = rtree_compute_mbr(n);

			//if (n->used == minm(root->m))
			//	printf("Deixou de ser outlier\n");
		}
		else {
			if (n == root->root || !root->reinsert_on_parent || at_reinsert != NULL) {
				//split leaf
				GList *glp = g_list_first(*parents);
				rtree_node *parent = glp ? (rtree_node*)(glp)->data : NULL;
				rtree_node *new_leaf = root->rtree_split_leaf_func(root, n, parent, g, gmbr);
				return new_leaf;
			}
			else
				root->reinsert_on_parent(n, root, parents, g, gmbr);
		}
	}
	else {
		//verify the best bound box 
		int index = root->choose_subtree_func(root, n, gmbr);
		assert(index != -1); // choose_subtree inserted the item
		//	return NULL;

		*parents = g_list_prepend(*parents, n);
		rtree_node *new_children = rtree_append_aux(root, n->dirs[index], parents, g, gmbr, at_reinsert);
		*parents = g_list_remove(*parents, n);
		if (new_children) {
			if (n->used < root->m) {
				n->dirs[n->used++] = new_children;
				n->mbr = rtree_compute_mbr(n);

				/*if (root->rtree_redivide) {
					qsort(VPOINTER(n->dirs), n->used, sizeof(struct rtree_node*), sort_node_used);
					int totalused = 0;
					for(int i = 0; i < n->used; i++)
						totalused += n->dirs[i]->used;
					int minimum_used = totalused / n->used;
					int i = 0;
					while (n->dirs[i]->used < minimum_used) {
						//printf("E menor que minimo %d %d\n", n->dirs[i]->used, minimum_used);
						double maximun_overlap = 0.0;
						int mergewith = -1;
						for(int j = 0; j < n->used; j++) {
							if (j == i || n->dirs[i]->used > (root->m - n->dirs[j]->used)) continue;
							double overlap = ENVELOPE_AREA(EnvelopeIntersection(n->dirs[i]->mbr, n->dirs[j]->mbr));
							if (overlap > maximun_overlap) {
								mergewith = j;
								maximun_overlap = overlap;
							}
						}
						if (mergewith != -1) {
							printf("Will merge %d with %d, at %d\n", i, mergewith, n->index);
							rtree_node *nd = n->dirs[mergewith];
							rtree_node *no = n->dirs[i];
							if (n->dirs[0]->type == DIRECTORY) {
								for (int x = 0; x < no->used; x++)
									nd->dirs[nd->used++] = no->dirs[x];
							}
							else {
								for (int x = 0; x < no->used; x++)
									nd->leaves[nd->used++] = no->leaves[x];
							}

							for(int x = i; x < n->used-1; x++)
								n->dirs[x] = n->dirs[x+1];

							n->dirs[n->used] = NULL;
							n->used--;
							no->used = 0; // will be deleted latter

							nd->mbr = rtree_compute_mbr(nd);
						}
						i++;
						mergesort(VPOINTER(n->dirs), n->used, sizeof(struct rtree_node*), sort_node_used);
					}

					double maximun_overlap = 0.0;
					int merge = -1;
					int mergewith = -1;
					for(int i = 0; i < n->used; i++) {
						for(int j = i; j < n->used; j++) {
							if (j == i) continue;
							double overlap = ENVELOPE_AREA(EnvelopeIntersection(n->dirs[i]->mbr, n->dirs[j]->mbr));
							if (overlap > maximun_overlap) {
								merge = i;
								mergewith = j;
								maximun_overlap = overlap;
							}
						}
					}
					if (merge != -1) {
						if (n->dirs[0]->type == DIRECTORY)
							root->rtree_redivide(n->dirs[merge], n->dirs[mergewith], root->m);
						else
							root->rtree_redivide_leaf(n->dirs[merge], n->dirs[mergewith], root->m);
					}

				}*/

				/*if (root->rtree_redivide) {
					for(int i = 0; i < n->used-1; i++) {
						for(int j = i+1; j < n->used-1; j++) {
							//if (j == i) continue;
							Envelope oti = n->dirs[i]->mbr;
							Envelope otj = n->dirs[j]->mbr;
							const float buff = 0.20;
							ENVELOPE_BUFFER(oti, buff * (oti.MaxX - oti.MinX), buff * (oti.MaxY - oti.MinY));
							ENVELOPE_BUFFER(otj, buff * (otj.MaxX - otj.MinX), buff * (otj.MaxY - otj.MinY));
							if (ENVELOPE_CONTAINS(oti, n->dirs[j]->mbr) ||
								ENVELOPE_CONTAINS(otj, n->dirs[i]->mbr)) {
									printf("Dirs %d(%d) contains %d(%d)\n", i, n->dirs[i]->used, j, n->dirs[j]->used);
									if (n->dirs[0]->type == DIRECTORY)
										root->rtree_redivide(n->dirs[i], n->dirs[j], root->m);
									else
										root->rtree_redivide_leaf(n->dirs[i], n->dirs[j], root->m);
								//}
							}
						}
					}
				}*/

				if (root->rtree_redivide) {
					for(int i = 0; i < n->used; i++) {
						for(int j = 0; j < n->used; j++) {
							if (i == j) continue;
							double overlap = ENVELOPE_AREA(EnvelopeIntersection(n->dirs[i]->mbr, n->dirs[j]->mbr));
							double total_area = ENVELOPE_AREA(n->dirs[i]->mbr) + ENVELOPE_AREA(n->dirs[j]->mbr);
							//printf("Fraction: %lf\n", overlap / total_area);
							if (overlap / total_area > 0.02) {
								if (n->dirs[0]->type == DIRECTORY)
									root->rtree_redivide(n->dirs[i], n->dirs[j], root->m);
								else
									root->rtree_redivide_leaf(n->dirs[i], n->dirs[j], root->m);
							}
						}
					}
				}
			}
			else {
				GList *glp = g_list_first(*parents);
				rtree_node *parent = glp ? (rtree_node*)(glp)->data : NULL;
				rtree_node *node = rtree_split(root, n, new_children, parent, root->m);
				return node;
			}
		}
		else {
			n->mbr = rtree_compute_mbr(n);
		}
	}
	
	return NULL;
}

void rtree_append(rtree_root *r, const GEOSGeometryH g) {
	Envelope gmbr;
	GEOSGetEnvelope(g, &gmbr);
	rtree_append_internal(r, g, gmbr, NULL);
}

void rtree_append_internal(rtree_root *r, const GEOSGeometryH g, const Envelope gmbr, rtree_node *at_reinsert) {
	GList *parents = NULL;
	rtree_node *new_leaf = rtree_append_aux(r, r->root, &parents, g, gmbr, at_reinsert);
	if (new_leaf) { 
		//root split 
		rtree_node *new_root = rtree_new_node(r, DIRECTORY);
		new_root->dirs[0] = r->root;
		new_root->dirs[1] = new_leaf;
		new_root->used = 2;
		r->root = new_root;
		r->root->mbr = rtree_compute_mbr(r->root);

		rtree_node *n = new_root;
		if (r->rtree_redivide) {
			for(int i = 0; i < n->used; i++) {
				for(int j = 0; j < n->used; j++) {
					if (i == j) continue;
					double overlap = 0.0;
					if (ENVELOPE_INTERSECTS(n->dirs[i]->mbr, n->dirs[j]->mbr))
						overlap = ENVELOPE_AREA(EnvelopeIntersection2(n->dirs[i]->mbr, n->dirs[j]->mbr));
					double total_area = ENVELOPE_AREA(n->dirs[i]->mbr) + ENVELOPE_AREA(n->dirs[j]->mbr);
					if (overlap / total_area > 0.01) {
						if (n->dirs[0]->type == DIRECTORY)
							r->rtree_redivide(n->dirs[i], n->dirs[j], r->m);
						else
							r->rtree_redivide_leaf(n->dirs[i], n->dirs[j], r->m);
					}
				}
			}
		}
	}
}

int rtree_height_recursive(const rtree_node *n) {
	return n->type == LEAF ? 1 : 1 + rtree_height_recursive(n->dirs[0]);
}

int rtree_height(const rtree_root *r) {
	return rtree_height_recursive(r->root);
}

void rtree_destroy_recursive(rtree_node *n) {
	if (n->type == DIRECTORY) {
		for(int i=0; i<n->used; i++)
			rtree_destroy_recursive(n->dirs[i]);
	}
	else { //LEAF
		for(int i=0; i<n->used; i++) {
			GEOSPreparedGeom_destroy(n->leaves[i].pgeo);
			// GEOM is released with the dataset
			// GEOSGeom_destroy(n->leaves[i].geo);
		}
	}
	g_free(n->dirs); // union with leaves
	g_free(n);
}

void rtree_destroy(rtree_root *r) {
	rtree_destroy_recursive(r->root);
	g_free(r);
}



#include <float.h>
#include "rtree-lazy.h"
#include "rtree-star.h"
#include "rtree-reinsert.h"

struct auxsplit {
	Envelope mbr;
	double value;
	int index;
};

struct bestchoice {
	int index;
	int border;
};

int rlazy_sortcomparer(const double i) {
	if (i < 0.0) return -1;
	if (i > 0.0) return 1;
	return 0;
}

int rlazy_sortx_min(const void *x, const void *y) {
	return rlazy_sortcomparer((*(rtree_node**)x)->mbr.MinX - (*(rtree_node**)y)->mbr.MinX);
}

int rtree_choose_subtree_lazy_new(rtree_root *root, rtree_node *n, const Envelope e) {
	
	int index = -1;
	int index2 = -1;
	int usedcurrent = 0;
	//int outlier_count = 0;
	//char no_intersects = TRUE;
	double minorloss = DBL_MAX;
	double minorloss2 = DBL_MAX;

	//mergesort(VPOINTER(n->dirs), n->used, sizeof(rtree_node*), rlazy_sortx_min);

	char check2 = (n->type == DIRECTORY && n->dirs[0]->type == LEAF);

	for(int i = 0; i < n->used; i++) {
		Envelope eaux = n->dirs[i]->mbr;
		ENVELOPE_MERGE(eaux, e);

		double oldarea = ENVELOPE_AREA(n->dirs[i]->mbr);
		double newarea = ENVELOPE_AREA(eaux);
		double loss;

		if (newarea - oldarea > 0.0) {

			double oldoverlap = 0.0;
			double overlap = 0.0;
			for(int j = 0; j < n->used; j++) {
				if (j == i) continue;
				//if ((eaux.MaxX >= n->dirs[j]->mbr.MinX) && (eaux.MaxY >= n->dirs[j]->mbr.MinY))
				if (ENVELOPE_INTERSECTS(eaux, n->dirs[j]->mbr)) {
					overlap += ENVELOPE_AREA(EnvelopeIntersection2(eaux, n->dirs[j]->mbr)); 
					
					if (ENVELOPE_INTERSECTS(n->dirs[i]->mbr, n->dirs[j]->mbr)) {
						oldoverlap += ENVELOPE_AREA(EnvelopeIntersection2(n->dirs[i]->mbr, n->dirs[j]->mbr));
					}
				}
			}

			loss = R1_LOSS(n->dirs[i]->mbr, oldoverlap, eaux, overlap);
			//loss = overlap; //MAX(oldoverlap, FLT_MIN);

		} else {
			loss = 0.0;	
		}

		if (loss < minorloss || (loss == minorloss && usedcurrent > n->dirs[i]->used)) {
			if (index != -1 && n->dirs[index]->used < root->m) {
				index2 = index;
				minorloss2 = minorloss;
			}
			index = i;
			minorloss = loss;
			usedcurrent = n->dirs[i]->used;
		}
	}

	if (check2 &&
		index2 != -1 &&
		n->dirs[index]->used == root->m) {
		index = index2;
	}

	return index;
}

int rtree_choose_subtree_lazy_new2(rtree_root *root, rtree_node *n, const Envelope e) {

	int index = -1;
	int index2 = -1;
	double minorloss = DBL_MAX;
	double minorloss2 = DBL_MAX;
	int usedcurrent = 0;

	char check2 = (n->type == DIRECTORY && n->dirs[0]->type == LEAF);

	for(int i = 0; i < n->used; i++) {
		Envelope eaux = n->dirs[i]->mbr;
		ENVELOPE_MERGE(eaux, e);
		double loss = R0_LOSS(n->dirs[i]->mbr, eaux);

		if (loss < minorloss || (loss == minorloss && usedcurrent > n->dirs[i]->used)) {
			if (index != -1 && n->dirs[index]->used < root->m) {
				index2 = index;
				minorloss2 = minorloss;
			}
			index = i;
			minorloss = loss;
			usedcurrent = n->dirs[i]->used;
		}
	}

	if (check2 &&
		index2 != -1 &&
		n->dirs[index]->used == root->m) {
		index = index2;
	}

	return index;
}

int rtree_choose_subtree_lazy(rtree_root *root, rtree_node *n, const Envelope e) {
	
	int index = -1;
	double minorloss = DBL_MAX;
	int usedcurrent = 0;

	for(int i = 0; i < n->used; i++) {
		Envelope eaux = n->dirs[i]->mbr;
		ENVELOPE_MERGE(eaux, e);
		double loss = R0_LOSS(n->dirs[i]->mbr, eaux);

		//printf("Loss: %f\n", loss);
		//			if (quality > bestquality || (quality == bestquality && area > currentarea)) {

		if (loss < minorloss || (loss == minorloss && usedcurrent > n->dirs[i]->used)) {
			index = i;
			minorloss = loss;
			usedcurrent = n->dirs[i]->used;
		}
	}

	return index;
}

int sort_split_aux(const void *x, const void *y) {
	double iy = *(double*)y;
	double ix = *(double*)x;
	if (ix == iy) return 0;
	if (ix > iy) return 1;
	return -1;
}

int sort_split_aux_rev(const void *x, const void *y) {
	double iy = *(double*)y;
	double ix = *(double*)x;
	if (ix == iy) return 0;
	if (ix < iy) return 1;
	return -1;
}

struct bestchoice rtree_split_lazy_choose_index(struct auxsplit **borders, const Envelope totalmbr, const int usedtotal, const int m) {

	qsort(VPOINTER(borders[0]), usedtotal, sizeof(struct auxsplit), sort_split_aux);
	qsort(VPOINTER(borders[1]), usedtotal, sizeof(struct auxsplit), sort_split_aux);
	qsort(VPOINTER(borders[2]), usedtotal, sizeof(struct auxsplit), sort_split_aux);
	qsort(VPOINTER(borders[3]), usedtotal, sizeof(struct auxsplit), sort_split_aux);

	//double totalarea = ENVELOPE_AREA(totalmbr);
	double mbrquality;
	double split_major_quality = DBL_MIN;
	double current_mbr_quality = DBL_MAX;
	struct bestchoice bc;
	bc.border = 0;
	bc.index = 0;

	int min = MAX(minm(m)-1, usedtotal-m-1);
	int pmax = usedtotal - min - 1;

	for(int b = 0; b < 4; b++) {

		Envelope newmbr1 = borders[b][0].mbr;
		for(int i = 1; i < min; i++) {
			ENVELOPE_MERGE(newmbr1, borders[b][i].mbr);
		}

		for(int i = min; i < pmax; i++) {
			ENVELOPE_MERGE(newmbr1, borders[b][i].mbr);

			Envelope newmbr2 = borders[b][i+1].mbr;
			for(int j = i+2; j < usedtotal; j++) {
				ENVELOPE_MERGE(newmbr2, borders[b][j].mbr);
			}

			double overlap = ENVELOPE_AREA(EnvelopeIntersection(newmbr1, newmbr2));
			double sarea = ENVELOPE_AREA(newmbr1) + ENVELOPE_AREA(newmbr2);
			double squality = SPLIT_QUALITY(overlap, sarea, newmbr1, newmbr2);

			//printf("Area: %f, Overlap: %f, Split Quality: %f : %f\n", sarea,
			//    overlap, squality, SPLIT_QUALITY(overlap, 1, newmbr1, newmbr2));
			//
		
			if (squality > split_major_quality ||
				(squality == split_major_quality && (mbrquality = (R0_QUALITY(newmbr1) + R0_QUALITY(newmbr2)) / 2.0) > current_mbr_quality)) {
				//squality == split_major_quality ? printf("Desempate %d: %f %f\n", desempate++, squality, split_major_quality) : 0;
				bc.border = b;
				bc.index = i;
				split_major_quality = squality;
				current_mbr_quality = mbrquality;
			}
		}
	}
	//print_geojson_footer();

	//fprintf(stderr, "Choosed: %d.%d", bc.border, bc.index);
	return bc;
}

rtree_node *rtree_split_lazy_leaf(rtree_root *root, rtree_node *n, rtree_node *parent, const GEOSGeometryH ngeo, const Envelope e) {

	/*if (parent && root->rtree_redivide) {
		char reorganized = FALSE;
		for(int i = 0; i < parent->used; i++) {
			if (parent->dirs[i] == n) continue;
			double overlap = ENVELOPE_AREA(EnvelopeIntersection(parent->dirs[i]->mbr, n->mbr));
			double total_area = ENVELOPE_AREA(parent->dirs[i]->mbr) + ENVELOPE_AREA(n->mbr);
			if (overlap / total_area > 0.01) {
				root->rtree_redivide_leaf(parent->dirs[i], n, root->m);
				reorganized = TRUE;
			}
		}
		if (reorganized) {
			int index = root->choose_subtree_func(root, parent, e);
			n = parent->dirs[index];
		}
		if (reorganized) {
			int index = root->choose_subtree_func(root, parent, e);
			n = parent->dirs[index];
			if (n->used < root->m) {
				n->leaves[n->used].mbr = e;
				n->leaves[n->used].chull = NULL;//OGR_G_ConvexHull(g);
				n->leaves[n->used].geo = ngeo;
				n->leaves[n->used].pgeo = NULL;
				n->used++;
				n->mbr = rtree_compute_mbr(n);
				return NULL;
			}
		}
	}*/

	rtree_node *new_leaf = rtree_new_node(root, LEAF);
	int usedtotal = n->used+1;

	//GEOSGeometryH ngeo_chull = NULL;//OGR_G_ConvexHull(ngeo);
	
	//add new entry (being added to node)
	rtree_leaf *leaves = g_new0(rtree_leaf, usedtotal);
	memcpy(leaves, n->leaves, sizeof(rtree_leaf) * n->used);
	leaves[n->used].geo = ngeo;
	//leaves[n->used].chull = ngeo_chull;
	leaves[n->used].mbr = e;

	struct auxsplit borderl[usedtotal];
	struct auxsplit borderr[usedtotal];
	struct auxsplit borderb[usedtotal];
	struct auxsplit bordert[usedtotal];

	for(int i = 0; i < usedtotal; i++) {
		borderl[i].value = leaves[i].mbr.MinX;
		borderr[i].value = leaves[i].mbr.MaxX;
		borderb[i].value = leaves[i].mbr.MinY;
		bordert[i].value = leaves[i].mbr.MaxY;
		borderl[i].index = borderr[i].index = borderb[i].index = bordert[i].index = i;
		borderl[i].mbr = borderr[i].mbr = borderb[i].mbr = bordert[i].mbr = leaves[i].mbr;
	}

	struct auxsplit *borders[4];
	borders[0] = borderl;
	borders[1] = borderr;
	borders[2] = bordert;
	borders[3] = borderb;

	Envelope totalmbr = n->mbr;
	ENVELOPE_MERGE(totalmbr, leaves[n->used].mbr);

	struct bestchoice best = rtree_split_lazy_choose_index(borders, totalmbr, usedtotal, usedtotal-1);

	//distribute
	n->used = 0;
	for(int i=0; i<usedtotal; i++) {
		if (i <= best.index)
			n->leaves[n->used++] = leaves[borders[best.border][i].index];
		else
			new_leaf->leaves[new_leaf->used++] = leaves[borders[best.border][i].index];
	}
	
	//printf("1: %d, 2: %d\n", n->used, new_leaf->used);
	assert(n->used + new_leaf->used == root->m+1);
	
	n->mbr = rtree_compute_mbr(n);
	new_leaf->mbr = rtree_compute_mbr(new_leaf);
	
#ifdef DEBUG
	double finaloverlap = ENVELOPE_AREA(EnvelopeIntersection(n->mbr, new_leaf->mbr));
	printf("Final overlap %f\n\n", finaloverlap);
#endif
	
	g_free(leaves);
	return new_leaf;


}

void rtree_redivide_lazy_leaf(rtree_node *n, rtree_node *s, int m) {

	assert(n != s);
	assert(n->type == s->type == LEAF);
#ifdef DEBUG
	double initialoverlap = ENVELOPE_AREA(EnvelopeIntersection(n->mbr, s->mbr));
#endif

	int usedtotal = n->used + s->used;

	rtree_leaf *leaves = g_new0(rtree_leaf, usedtotal);
	memcpy(leaves, n->leaves, sizeof(rtree_leaf) * n->used);
	memcpy(&leaves[n->used], s->leaves, sizeof(rtree_leaf) * s->used);

	struct auxsplit borderl[usedtotal];
	struct auxsplit borderr[usedtotal];
	struct auxsplit borderb[usedtotal];
	struct auxsplit bordert[usedtotal];

	for(int i = 0; i < usedtotal; i++) {
		borderl[i].value = leaves[i].mbr.MinX;
		borderr[i].value = leaves[i].mbr.MaxX;
		borderb[i].value = leaves[i].mbr.MinY;
		bordert[i].value = leaves[i].mbr.MaxY;
		borderl[i].index = borderr[i].index = borderb[i].index = bordert[i].index = i;
		borderl[i].mbr = borderr[i].mbr = borderb[i].mbr = bordert[i].mbr = leaves[i].mbr;
	}

	struct auxsplit *borders[4];
	borders[0] = borderl;
	borders[1] = borderr;
	borders[2] = bordert;
	borders[3] = borderb;

	Envelope totalmbr = n->mbr;
	ENVELOPE_MERGE(totalmbr, s->mbr);

	struct bestchoice best = rtree_split_lazy_choose_index(borders, totalmbr, usedtotal, m);

	//distribute
	n->used = 0;
	s->used = 0;
	for(int i=0; i<usedtotal; i++) {
		if (i <= best.index)
			n->leaves[n->used++] = leaves[borders[best.border][i].index];
		else
			s->leaves[s->used++] = leaves[borders[best.border][i].index];
	}

	//printf("1: %d, 2: %d\n", n->used, new_leaf->used);
	assert(n->used + s->used == usedtotal);

	n->mbr = rtree_compute_mbr(n);
	s->mbr = rtree_compute_mbr(s);

#ifdef DEBUG
	double finaloverlap = ENVELOPE_AREA(EnvelopeIntersection(n->mbr, s->mbr));
	printf("Leaf overlap: initial %f final %f. %d %d\n", initialoverlap, finaloverlap, n->used, s->used);
#endif

	g_free(leaves);
}

rtree_node *rtree_split_lazy(rtree_root *root, rtree_node *n, rtree_node *newn, rtree_node *parent, int m) {

	/*if (parent && root->rtree_redivide) {
		char reorganized = FALSE;
		for(int i = 0; i < parent->used; i++) {
			if (parent->dirs[i] == n) continue;
			double overlap = ENVELOPE_AREA(EnvelopeIntersection(parent->dirs[i]->mbr, n->mbr));
			double total_area = ENVELOPE_AREA(parent->dirs[i]->mbr) + ENVELOPE_AREA(n->mbr);
			if (overlap / total_area > 0.01) {
				root->rtree_redivide(parent->dirs[i], n, root->m);
				reorganized = TRUE;
			}
		}
		if (reorganized) {
			int index = root->choose_subtree_func(root, parent, newn->mbr);
			n = parent->dirs[index];
			if (n->used < root->m) {
				n->dirs[n->used++] = newn;
				return NULL;
			}
		}
	}*/

	rtree_node *new_dir = rtree_new_node(root, DIRECTORY);
	int usedtotal = n->used+1;

	//add new entry (being added to node)
	rtree_node **dirs = g_new(rtree_node*, usedtotal);
	memcpy(dirs, n->dirs, sizeof(rtree_node*) * n->used);
	dirs[n->used] = newn;

	struct auxsplit borderl[usedtotal];
	struct auxsplit borderr[usedtotal];
	struct auxsplit borderb[usedtotal];
	struct auxsplit bordert[usedtotal];

	for(int i = 0; i < usedtotal; i++) {
		rtree_node *li = dirs[i];
		borderl[i].value = li->mbr.MinX;
		borderr[i].value = li->mbr.MaxX;
		borderb[i].value = li->mbr.MinY;
		bordert[i].value = li->mbr.MaxY;
		borderl[i].index = borderr[i].index = borderb[i].index = bordert[i].index = i;
		borderl[i].mbr = borderr[i].mbr = borderb[i].mbr = bordert[i].mbr = li->mbr;
	}

	struct auxsplit *borders[4];
	borders[0] = borderl;
	borders[1] = borderr;
	borders[2] = bordert;
	borders[3] = borderb;

	Envelope totalmbr = n->mbr;
	ENVELOPE_MERGE(totalmbr, newn->mbr);

	struct bestchoice best = rtree_split_lazy_choose_index(borders, totalmbr, usedtotal, m);

	//distribute
	n->used = 0;
	for(int i=0; i<usedtotal; i++) {
		if (i <= best.index)
			n->dirs[n->used++] = dirs[borders[best.border][i].index];
		else
			new_dir->dirs[new_dir->used++] = dirs[borders[best.border][i].index];
	}

	assert(n->used + new_dir->used == m+1);

	n->mbr = rtree_compute_mbr(n);
	new_dir->mbr = rtree_compute_mbr(new_dir);

#ifdef DEBUG
	double finaloverlap = ENVELOPE_AREA(EnvelopeIntersection(n->mbr, new_dir->mbr));
	printf("Final overlap %f\n\n", finaloverlap);
#endif

	g_free(dirs);
	return new_dir;

}

void rtree_redivide_lazy(rtree_node *n, rtree_node *s, int m) {

	assert(n != s);
	assert(n->type == s->type);
	assert(n->type == DIRECTORY);
#ifdef DEBUG
	double initialoverlap = ENVELOPE_AREA(EnvelopeIntersection(n->mbr, s->mbr));
#endif

	int usedtotal = n->used + s->used;
	rtree_node **dirs = g_new(rtree_node*, usedtotal);
	memcpy(dirs, n->dirs, sizeof(rtree_node*) * n->used);
	memcpy(&dirs[n->used], s->dirs, sizeof(rtree_node*) * s->used);

	struct auxsplit borderl[usedtotal];
	struct auxsplit borderr[usedtotal];
	struct auxsplit borderb[usedtotal];
	struct auxsplit bordert[usedtotal];

	for(int i = 0; i < usedtotal; i++) {
		rtree_node *li = dirs[i];
		borderl[i].value = li->mbr.MinX;
		borderr[i].value = li->mbr.MaxX;
		borderb[i].value = li->mbr.MinY;
		bordert[i].value = li->mbr.MaxY;
		borderl[i].index = borderr[i].index = borderb[i].index = bordert[i].index = i;
		borderl[i].mbr = borderr[i].mbr = borderb[i].mbr = bordert[i].mbr = li->mbr;
	}

	struct auxsplit *borders[4];
	borders[0] = borderl;
	borders[1] = borderr;
	borders[2] = bordert;
	borders[3] = borderb;

	Envelope totalmbr = n->mbr;
	ENVELOPE_MERGE(totalmbr, s->mbr);

	struct bestchoice best = rtree_split_lazy_choose_index(borders, totalmbr, usedtotal, m);

	//distribute
	n->used = 0;
	s->used = 0;
	for(int i=0; i<usedtotal; i++) {
		if (i <= best.index)
			n->dirs[n->used++] = dirs[borders[best.border][i].index];
		else
			s->dirs[s->used++] = dirs[borders[best.border][i].index];
	}

	assert(n->used + s->used == usedtotal);

	n->mbr = rtree_compute_mbr(n);
	s->mbr = rtree_compute_mbr(s);

#ifdef DEBUG
	double finaloverlap = ENVELOPE_AREA(EnvelopeIntersection(n->mbr, s->mbr));
	printf("Dir overlap: initial %f final %f. %d %d\n", initialoverlap, finaloverlap, n->used, s->used);
#endif

	g_free(dirs);
}

rtree_root *rtree_new_r0(const int m, const int servers) {
	return rtree_new(
		m,
		servers,
		rtree_choose_subtree_lazy,
		rtree_split_star,
		rtree_split_star_leaf,
		greedy_reinsert_on_parent,
		NULL, NULL
	);
}

rtree_root *rtree_new_rlazy(const int m, const int servers) {
	return rtree_new(
		m,
		servers,
		rtree_choose_subtree_lazy_new2,
		rtree_split_lazy,
		rtree_split_lazy_leaf,
		NULL,
		rtree_redivide_lazy,
		rtree_redivide_lazy_leaf
	);
}



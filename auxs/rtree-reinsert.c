
#include <float.h>
#include "rtree.h"
#include "ogrext.h"
#include <ogr_api.h>
#include <ogr_core.h>
#include "rtree-lazy.h"

struct idis {
	double distance; //need to stay first. check sort_idis.
	GEOSGeometryH geo;
	Envelope mbr;
	double gain;
	double mygain;
	int bestidx;
	int index;
};

struct maxchoice {
	double gain;
	Envelope mbr;
	int *idxs;
	struct idis *d[4];
	int count;
};


struct maxgain {
	double gain;
	struct idis *distances;
	double oborder;
	double (*get_border_start)(const Envelope e);
	int bestidx;
	int remaining;
};

double getMinX(const Envelope e) {
	return e.MinX;
}
double getMinY(const Envelope e) {
	return e.MinY;
}
double getMaxX(const Envelope e) {
	return e.MaxX;
}
double getMaxY(const Envelope e) {
	return e.MaxY;
}

int sort_idis(const void *x, const void *y) {
	double iy = *(double*)y;
	double ix = *(double*)x;
	if (ix == iy) return 0;
	if (ix < iy) return 1;
	return -1;
}

int sort_maxgain(const void *x, const void *y) {
	double iy = *(double*)y;
	double ix = *(double*)x;
	if (ix == iy) return 0;
	if (ix < iy) return 1;
	return -1;
}

int sort_border(const void *x, const void *y) {
	double iy = *(double*)y;
	double ix = *(double*)x;
	if (ix == iy) return 0;
	if (ix > iy) return 1;
	return -1;
}

int sort_border_rev(const void *x, const void *y) {
	double iy = *(double*)y;
	double ix = *(double*)x;
	if (ix == iy) return 0;
	if (ix < iy) return 1;
	return -1;
}

void reinsert_items(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g,
	const Envelope gmbr, struct idis *distances, const int items_to_reinsert) {
	
	// reorganize the geos sequencially
	rtree_leaf *copy = (rtree_leaf *)g_memdup(VPOINTER(n->leaves), sizeof(rtree_leaf)*n->used);
	int last = n->used;
	n->used = 0;
	for(int i = 0; i < last; i++) {
		if (copy[i].geo != NULL)
			n->leaves[n->used++] = copy[i];
	}
	g_free(copy);

	n->mbr = rtree_compute_mbr(n);

	// recompute the parent mbr
	for(GList *parent = *parents; parent; parent = parent->next) {
		rtree_node *myparent = (rtree_node*)parent->data;
		myparent->mbr = rtree_compute_mbr(myparent);
	}

	for(int i = 0; i < items_to_reinsert; i++) {
		rtree_append_internal(root, distances[i].geo, distances[i].mbr, n);
	}
	rtree_append_internal(root, g, gmbr, n);
}

void rstar_reinsert_on_parent(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g, const Envelope gmbr) {

	struct idis distances[n->used];
	double x, y;

	// new strategy: inferior corner leads to better results!
	x = n->mbr.MinX;
	y = n->mbr.MinY;

	// original method: distance to the center
	//x = (n->mbr.MinX + n->mbr.MaxX) / 2.0;
	//y = (n->mbr.MinY + n->mbr.MaxY) / 2.0;

	for(int i = 0; i < n->used; i++) {
		rtree_leaf li = n->leaves[i];
		int xi = (li.mbr.MinX + li.mbr.MaxX) / 2.0;
		int yi = (li.mbr.MinY + li.mbr.MaxY) / 2.0;
		distances[i].index = i;
		distances[i].distance = sqrt(pow(xi-x, 2) + pow(yi-y, 2));
		distances[i].geo = li.geo;
		distances[i].mbr = li.mbr;
	}

	qsort(VPOINTER(distances), n->used, sizeof(struct idis), sort_idis);

#ifdef DEBUG
	for(int i = 0; i < n->used; i++) {
		printf("idis %d - %f\n", distances[i].index, distances[i].distance);
	}
	printf("\n\n");
#endif

	//remove p outiliers objects from n
	int p = (int)(reinsertthreshold * n->used);
	for(int i = 0; i < p; i++) // set null on removed geometries
		n->leaves[distances[i].index].geo = NULL;

	reinsert_items(n, root, parents, g, gmbr, distances, p);

}

int get_best_p(struct idis *distances, double oborder, double first, const int pmax, double *maxgain) {

	for(int i = 0; i < pmax; i++) {
		distances[i].gain = (oborder * fabs(distances[i+1].distance - first)) / (i+1);
		//printf("gain %d: %f ", i, gain_per_object[i]);
	}

	*maxgain = distances[0].gain;
	int maxgain_index = 0;
	for(int i = 1; i < pmax; i++) {
		if (distances[i].gain > *maxgain) {
			*maxgain = distances[i].gain;
			maxgain_index = i;
		}
	}
	//printf("max gain choosed: %d\n", maxgain_index);
	return maxgain_index+1 ;
}

void border_reinsert_on_parent(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g, const Envelope gmbr) {

	struct idis distancesl[n->used];
	struct idis distancesr[n->used];
	struct idis distancesb[n->used];
	struct idis distancest[n->used];

	for(int i = 0; i < n->used; i++) {
		rtree_leaf li = n->leaves[i];
		distancesl[i].distance = li.mbr.MinX;
		distancesr[i].distance = li.mbr.MaxX;
		distancesb[i].distance = li.mbr.MinY;
		distancest[i].distance = li.mbr.MaxY;
		distancesl[i].index = distancesr[i].index = distancesb[i].index = distancest[i].index = i;
		distancesl[i].geo = distancesr[i].geo = distancesb[i].geo = distancest[i].geo = li.geo;
		distancesl[i].mbr = distancesr[i].mbr = distancesb[i].mbr = distancest[i].mbr = li.mbr;
	}

	int pmax = (int)(reinsertthreshold * n->used);

	qsort(VPOINTER(distancesl), n->used, sizeof(struct idis), sort_border);
	qsort(VPOINTER(distancesr), n->used, sizeof(struct idis), sort_border_rev);
	qsort(VPOINTER(distancesb), n->used, sizeof(struct idis), sort_border);
	qsort(VPOINTER(distancest), n->used, sizeof(struct idis), sort_border_rev);

	double gainl, gainr, gainb, gaint;
	int lef = get_best_p(distancesl, (n->mbr.MaxY - n->mbr.MinY), n->mbr.MinX, pmax, &gainl);
	int rig = get_best_p(distancesr, (n->mbr.MaxY - n->mbr.MinY), n->mbr.MaxX, pmax, &gainr);
	int bot = get_best_p(distancesb, (n->mbr.MaxX - n->mbr.MinX), n->mbr.MinY, pmax, &gainb);
	int top = get_best_p(distancest, (n->mbr.MaxX - n->mbr.MinX), n->mbr.MaxY, pmax, &gaint);

	// choose best border
	struct maxgain mg[4];
	mg[0].gain = gainl; mg[0].bestidx = lef; mg[0].distances = distancesl;
	mg[1].gain = gainr; mg[1].bestidx = rig; mg[1].distances = distancesr;
	mg[2].gain = gaint; mg[2].bestidx = top; mg[2].distances = distancest;
	mg[3].gain = gainb; mg[3].bestidx = bot; mg[3].distances = distancesb;
	qsort(VPOINTER(mg), 4, sizeof(struct maxgain), sort_maxgain);
	int p = mg[0].bestidx;
	int p2 = mg[1].bestidx;
	struct idis *distances = mg[0].distances;
	struct idis *distances2 = mg[1].distances;

	/*int p = lef;
	int p2 = -1;
	double gain = gainl;
	struct idis *distances = distancesl;
	struct idis *distances2;
	if (gain < gainr) {
		p2 = p;
		distances2 = distances;
		p = rig;
		distances = distancesr;
		gain = gainr;
	}
	if (gain < gainb) {
		p2 = p;
		distances2 = distances;
		p = bot;
		distances = distancesb;
		gain = gainb;
	}
	if (gain < gaint) {
		p2 = p;
		distances2 = distances;
		p = top;
		distances = distancest;
		gain = gaint;
	}*/

	// remove p objects
	for(int i = 0; i < p; i++) // set null on removed geometries
		n->leaves[distances[i].index].geo = NULL;
	if (p2 != -1 && p2 <= (pmax-p)) {
		//printf("usou p2\n");
		int max = MIN(pmax-p, p2);
		for(int i = 0; i < max; i++) {
			if (n->leaves[distances2[i].index].geo != NULL) {
				n->leaves[distances2[i].index].geo = NULL;
				distances[p++] = distances2[i];
			}
		}
	}
	//else
	//	printf("nao usou p2\n");

	reinsert_items(n, root, parents, g, gmbr, distances, p);
}

void get_best_p_lazy_new(struct idis *shrink, struct idis *dt, struct idis *db, struct idis *dl, struct idis *dr, const Envelope mbr, const int pmax, struct maxgain *m) {
	int ib = 0;
	int il = 0;
	int ir = 0;
	int it = 0;
	Envelope newmbr = mbr;
	double startarea = ENVELOPE_AREA(mbr);
	double prior_gain = 0.0;
	m->gain = 0.0;

	for(int i = 0; i < pmax; i++) {
		if (shrink[i].index == dt[it].index)
			newmbr.MaxY = dt[++it].mbr.MaxY;
		if (shrink[i].index == dl[il].index)
			newmbr.MinX = dl[++il].mbr.MinX;
		if (shrink[i].index == dr[ir].index)
			newmbr.MaxX = dr[++ir].mbr.MaxX;
		if (shrink[i].index == db[ib].index)
			newmbr.MinY = db[++ib].mbr.MinY;
		double newarea = ENVELOPE_AREA(newmbr);
		double gain = (startarea - newarea) / (i+1);
		
		if (gain > m->gain) {
			m->gain = gain;
			m->bestidx = i;
		}

		shrink[i].gain = m->gain;
		shrink[i].bestidx = m->bestidx;
		shrink[i].mygain = (startarea - newarea) - prior_gain;
		prior_gain += shrink[i].mygain;
	}

	m->distances = shrink;
	m->remaining = pmax;
	m->oborder = 0;
	m->get_border_start = NULL;
}


void get_best_p_lazy(struct idis *distances, double oborder, double (*getfirst)(const Envelope e), const Envelope mbr, const int pmax, struct maxgain *m) {

	for(int i = 0; i < pmax; i++) {
		distances[i].gain = (oborder * fabs(distances[i+1].distance - getfirst(mbr))) / (i+1);
		//printf("gain %d: %f ", i, distances[i].gain);
	}

	m->gain = distances[0].gain;
	int maxgain_index = 0;
	for(int i = 1; i < pmax; i++) {
		if (distances[i].gain > m->gain) {
			m->gain = distances[i].gain;
			maxgain_index = i;
		}
		distances[i].bestidx = maxgain_index;
	}
	//printf("max gain choosed: %d\n", maxgain_index);
	m->distances = distances;
	m->bestidx = maxgain_index;
	m->remaining = pmax;
	m->oborder = oborder;
	m->get_border_start = getfirst;
}

void greedy_2bestb_reinsert_on_parent(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g, const Envelope gmbr) {

	struct idis distancesl[n->used];
	struct idis distancesr[n->used];
	struct idis distancesb[n->used];
	struct idis distancest[n->used];

	for(int i = 0; i < n->used; i++) {
		rtree_leaf li = n->leaves[i];
		distancesl[i].distance = li.mbr.MinX;
		distancesr[i].distance = li.mbr.MaxX;
		distancesb[i].distance = li.mbr.MinY;
		distancest[i].distance = li.mbr.MaxY;
		distancesl[i].index = distancesr[i].index = distancesb[i].index = distancest[i].index = i;
		distancesl[i].geo = distancesr[i].geo = distancesb[i].geo = distancest[i].geo = li.geo;
		distancesl[i].mbr = distancesr[i].mbr = distancesb[i].mbr = distancest[i].mbr = li.mbr;
	}

	int pmax = (int)(reinsertthreshold * n->used);

	qsort(VPOINTER(distancesl), n->used, sizeof(struct idis), sort_border);
	qsort(VPOINTER(distancesr), n->used, sizeof(struct idis), sort_border_rev);
	qsort(VPOINTER(distancesb), n->used, sizeof(struct idis), sort_border);
	qsort(VPOINTER(distancest), n->used, sizeof(struct idis), sort_border_rev);

	struct maxgain mg[4];
	get_best_p_lazy_new(distancest, distancest, distancesb, distancesl, distancesr, n->mbr, pmax, &mg[0]);
	get_best_p_lazy_new(distancesb, distancest, distancesb, distancesl, distancesr, n->mbr, pmax, &mg[1]);
	get_best_p_lazy_new(distancesl, distancest, distancesb, distancesl, distancesr, n->mbr, pmax, &mg[2]);
	get_best_p_lazy_new(distancesr, distancest, distancesb, distancesl, distancesr, n->mbr, pmax, &mg[3]);
	qsort(VPOINTER(mg), 4, sizeof(struct maxgain), sort_maxgain);

	struct idis distances[pmax];
	int items_to_reinsert = 0;
	int i = 0;
	int mg_atu = 0;
	while (items_to_reinsert < pmax && mg_atu < 2) {
		int leaf_index_to_remove = mg[mg_atu].distances[i].index;
		if (n->leaves[leaf_index_to_remove].geo != NULL) {
			n->leaves[leaf_index_to_remove].geo = NULL;
			distances[items_to_reinsert++] = mg[mg_atu].distances[i];
		}
		i++;

		if (i >= mg[mg_atu].bestidx) { // all itens from mg[mg_atu] was taken
			mg_atu++;
			
			int space = (pmax-items_to_reinsert);
			while (mg[mg_atu].bestidx > space && mg_atu < 3) {
				double g1 = mg[mg_atu].distances[space].gain;
				double g2 = mg[mg_atu+1].distances[MIN(space, mg[mg_atu+1].bestidx)].gain;
				if (g2 > g1) {
					mg_atu++;
					//printf("SKIPED %d\n", yy++);
				}
				else break;
			}
				
			i = 0;
			//printf("mg_atu %d\n", mg_atu);
		}
		
	}

	reinsert_items(n, root, parents, g, gmbr, distances, items_to_reinsert);

}

void greedy_2bestborder_reinsert_on_parent(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g, const Envelope gmbr) {

	/*print_geojson_header();	
	print_geojson_mbr(n->mbr, 99999999);
	for(int i = 0; i<n->used; i++) {
		print_geojson_mbr(n->leaves[i].mbr, i);
	}
	print_geojson_footer();*/

	struct idis distancesl[n->used];
	struct idis distancesr[n->used];
	struct idis distancesb[n->used];
	struct idis distancest[n->used];

	for(int i = 0; i < n->used; i++) {
		rtree_leaf li = n->leaves[i];
		distancesl[i].distance = li.mbr.MinX;
		distancesr[i].distance = li.mbr.MaxX;
		distancesb[i].distance = li.mbr.MinY;
		distancest[i].distance = li.mbr.MaxY;
		distancesl[i].index = distancesr[i].index = distancesb[i].index = distancest[i].index = i;
		distancesl[i].geo = distancesr[i].geo = distancesb[i].geo = distancest[i].geo = li.geo;
		distancesl[i].mbr = distancesr[i].mbr = distancesb[i].mbr = distancest[i].mbr = li.mbr;
	}

	int pmax = (int)(reinsertthreshold * n->used);
	//int pmax = n->used - minm(n->used);

	qsort(VPOINTER(distancesl), n->used, sizeof(struct idis), sort_border);
	qsort(VPOINTER(distancesr), n->used, sizeof(struct idis), sort_border_rev);
	qsort(VPOINTER(distancesb), n->used, sizeof(struct idis), sort_border);
	qsort(VPOINTER(distancest), n->used, sizeof(struct idis), sort_border_rev);

	struct maxgain mg[4];
	get_best_p_lazy_new(distancest, distancest, distancesb, distancesl, distancesr, n->mbr, pmax, &mg[0]);
	get_best_p_lazy_new(distancesb, distancest, distancesb, distancesl, distancesr, n->mbr, pmax, &mg[1]);
	get_best_p_lazy_new(distancesl, distancest, distancesb, distancesl, distancesr, n->mbr, pmax, &mg[2]);
	get_best_p_lazy_new(distancesr, distancest, distancesb, distancesl, distancesr, n->mbr, pmax, &mg[3]);

	qsort(VPOINTER(mg), 4, sizeof(struct maxgain), sort_maxgain);
	struct maxgain *mga = &mg[0];

	struct idis distances[pmax];
	int items_to_reinsert = 0;
	int i = 0;
	int space = pmax;
	int bordersused = 0;

	while (space > 0) {

		int leaf_index_to_remove = mga->distances[i].index;
		if (n->leaves[leaf_index_to_remove].geo != NULL) {
			n->leaves[leaf_index_to_remove].geo = NULL;
			distances[items_to_reinsert++] = mga->distances[i];
			space--;
			//printf("Removed item at index: %d\n", leaf_index_to_remove);
		}
		i++;
		
		if (bordersused < 1 && space > 0 && i > mga->bestidx) { // all itens from mg[mg_atu] was taken
			bordersused++;

			// remove used items
			mga->distances = &mga->distances[i];
			mga->remaining -= i;
			mga->gain = 0.0;
			double cumulative_gain = 0.0;
			for(int j=0; j < space; j++) {
				double gain = (mga->distances[j].mygain + cumulative_gain) / (j+1);	
				cumulative_gain += mga->distances[j].mygain;
				if (gain > mga->gain) {
					mga->gain = gain;
					mga->bestidx = j;
				}
			}
			
			// best idx again, for pmax-items_to_reinsert
			for(int j=1; j < 4; j++) {
				mg[j].bestidx = mg[j].distances[space-1].bestidx;
				mg[j].gain = mg[j].distances[space-1].gain;
			}
		
			qsort(VPOINTER(mg), 4, sizeof(struct maxgain), sort_maxgain);
			mga = &mg[0];
			i = 0;
		}
		else break;
	}

	reinsert_items(n, root, parents, g, gmbr, distances, items_to_reinsert);
}


void produced(const int k, const int *p, void *data) {
	struct maxchoice *mc = (struct maxchoice*)data;	
	
	int borders[4];
	borders[0] = borders[1] = borders[2] = borders[3] = 0;
	for(int i = 0; i < k; i++) {
		borders[p[i]]++;
	}

	Envelope newmbr;
	newmbr.MinX = mc->d[0][borders[0]].mbr.MinX;
	newmbr.MaxX = mc->d[1][borders[1]].mbr.MaxX;
	newmbr.MaxY = mc->d[2][borders[2]].mbr.MaxY;
	newmbr.MinY = mc->d[3][borders[3]].mbr.MinY;
	double gain = R0_GAIN(mc->mbr, newmbr) / k;
	if (gain > mc->gain) {
		mc->gain = gain;
		mc->count = k;
		memcpy(mc->idxs, p, sizeof(int)*k);
	}
}

void enumerate(void *data, const int n, const int k, const int start, const int level, int *current, void(*produced)(const int k, const int *p, void *data)) {
			
	if (level == k) {
		produced(k, current, data);
		return;
	}

	for(int i = start; i < n; i++) {
		current[level] = i;
		enumerate(data, n, k, start > i ? start : i, level+1, current, produced);
	}
}

void exponen_reinsert_on_parent(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g, const Envelope gmbr) {

	/*print_geojson_header();	
	print_geojson_mbr(n->mbr, 99999999);
	for(int i = 0; i<n->used; i++) {
		print_geojson_mbr(n->leaves[i].mbr, i);
	}
	print_geojson_footer();
	*/

	struct idis distancesl[n->used];
	struct idis distancesr[n->used];
	struct idis distancesb[n->used];
	struct idis distancest[n->used];

	for(int i = 0; i < n->used; i++) {
		rtree_leaf li = n->leaves[i];
		distancesl[i].distance = li.mbr.MinX;
		distancesr[i].distance = li.mbr.MaxX;
		distancesb[i].distance = li.mbr.MinY;
		distancest[i].distance = li.mbr.MaxY;
		distancesl[i].index = distancesr[i].index = distancesb[i].index = distancest[i].index = i;
		distancesl[i].geo = distancesr[i].geo = distancesb[i].geo = distancest[i].geo = li.geo;
		distancesl[i].mbr = distancesr[i].mbr = distancesb[i].mbr = distancest[i].mbr = li.mbr;
	}

	int pmax = (int)(reinsertthreshold * n->used);
	qsort(VPOINTER(distancesl), n->used, sizeof(struct idis), sort_border);
	qsort(VPOINTER(distancesr), n->used, sizeof(struct idis), sort_border_rev);
	qsort(VPOINTER(distancesb), n->used, sizeof(struct idis), sort_border);
	qsort(VPOINTER(distancest), n->used, sizeof(struct idis), sort_border_rev);

	int current[pmax];
	int max[pmax];
	struct maxchoice mc;
	struct idis distances[pmax];
	mc.gain = 0.0;
	mc.mbr = n->mbr;
	mc.idxs = max;
	mc.count = 0;
	mc.d[0] = distancesl;
	mc.d[1] = distancesr;
	mc.d[2] = distancest;
	mc.d[3] = distancesb;
	int k = 1;
	while (k <= pmax) {
		enumerate(&mc, 4 /*borders*/, k, 0, 0, current, produced);
		k++;
	}

	int borders[4];
	borders[0] = borders[1] = borders[2] = borders[3] = 0;
	for(int i = 0; i < mc.count; i++) {
		borders[mc.idxs[i]]++;
	}

	int d = 0;
	for(int i = 0; i < 4; i++) {
		for(int j = 0; j < borders[i]; j++) {
			int leaf_index_to_remove = mc.d[i][j].index;
			if (n->leaves[leaf_index_to_remove].geo != NULL) {
				n->leaves[leaf_index_to_remove].geo = NULL;
				distances[d++] = mc.d[i][j];	
			}
		}
	}

	//if (d < pmax) printf("D < PMAX: %d < %d\n", d, pmax);
	
	reinsert_items(n, root, parents, g, gmbr, distances, d);
}

void greedy_reinsert_on_parent(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g, const Envelope gmbr) {

	/*print_geojson_header();	
	print_geojson_mbr(n->mbr, 99999999);
	for(int i = 0; i<n->used; i++) {
		print_geojson_mbr(n->leaves[i].mbr, i);
	}
	print_geojson_footer();
	*/

	struct idis distancesl[n->used];
	struct idis distancesr[n->used];
	struct idis distancesb[n->used];
	struct idis distancest[n->used];

	for(int i = 0; i < n->used; i++) {
		rtree_leaf li = n->leaves[i];
		distancesl[i].distance = li.mbr.MinX;
		distancesr[i].distance = li.mbr.MaxX;
		distancesb[i].distance = li.mbr.MinY;
		distancest[i].distance = li.mbr.MaxY;
		distancesl[i].index = distancesr[i].index = distancesb[i].index = distancest[i].index = i;
		distancesl[i].geo = distancesr[i].geo = distancesb[i].geo = distancest[i].geo = li.geo;
		distancesl[i].mbr = distancesr[i].mbr = distancesb[i].mbr = distancest[i].mbr = li.mbr;
	}

	int pmax = (int)(reinsertthreshold * n->used);
	int remaining = pmax;
	qsort(VPOINTER(distancesl), n->used, sizeof(struct idis), sort_border);
	qsort(VPOINTER(distancesr), n->used, sizeof(struct idis), sort_border_rev);
	qsort(VPOINTER(distancesb), n->used, sizeof(struct idis), sort_border);
	qsort(VPOINTER(distancest), n->used, sizeof(struct idis), sort_border_rev);

	struct maxchoice mc;
	struct idis distances[pmax];
	int dist = 0;
	int remaining_sizes[4] = {pmax, pmax, pmax, pmax};

	mc.d[0] = distancesl;
	mc.d[1] = distancesr;
	mc.d[2] = distancest;
	mc.d[3] = distancesb;
	
	// greed_min_p
	int m = 5;
	int bestborder = 0;
	Envelope bestmbr;
	Envelope currentmbr = n->mbr;
	mc.gain = 0.0;
	do {
		mc.count = 0;

		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < MIN(remaining_sizes[i], MIN(m, remaining)); j++) {
				Envelope newmbr = currentmbr;
				switch (i) {
					case 0: newmbr.MinX = mc.d[i][j+1].mbr.MinX; break;
					case 1: newmbr.MaxX = mc.d[i][j+1].mbr.MaxX; break;
					case 2: newmbr.MaxY = mc.d[i][j+1].mbr.MaxY; break;
					case 3: newmbr.MinY = mc.d[i][j+1].mbr.MinY; break;
				}
				double gain = R0_GAIN(n->mbr, newmbr) / (dist+j+1);
				if (gain > mc.gain) {
					mc.gain = gain;
					mc.count = j+1;
					bestborder = i;
					bestmbr = newmbr;
				}
			}
		}

		if (mc.count == 0.0)
			break;
		
		for(int i = 0; i < mc.count; i++) {
			distances[dist++] = mc.d[bestborder][i];
			for(int b = 0; b < 4; b++) {
				if (b == bestborder) continue;
				for(int j = 0; j < remaining_sizes[b]; j++) {
					if (mc.d[b][j].index == mc.d[bestborder][i].index) {
						for(int k = j; k < remaining_sizes[b]-1; k++)
							mc.d[b][k] = mc.d[b][k+1];
						remaining_sizes[b]--;
						break;
					}
				}
			}
		}
		mc.d[bestborder] = &mc.d[bestborder][mc.count];
		remaining_sizes[bestborder] -= mc.count;

		currentmbr = bestmbr;
		remaining -= mc.count;

	} while (1);

	for(int i = 0; i < dist; i++) {
		n->leaves[distances[i].index].geo = NULL;
	}

	//if (dist < pmax) printf("DIST < PMAX: %d < %d\n", dist, pmax);
	
	reinsert_items(n, root, parents, g, gmbr, distances, dist);
}


void rlazy_reinsert_on_parent(rtree_node *n, rtree_root *root, GList **parents, const GEOSGeometryH g, const Envelope gmbr) {
	GList *parent = *parents;
	rtree_node *myparent = (rtree_node*)parent->data;

	double min_area_inc = DBL_MAX;
	int choosed_to_insert = -1;
	for(int i = 0; i < myparent->used; i++) {
		if (myparent->dirs[i]->used < root->m) {
			Envelope eaux = gmbr;
			ENVELOPE_MERGE(eaux, myparent->dirs[i]->mbr);
			double area_inc = ENVELOPE_AREA(eaux) - ENVELOPE_AREA(myparent->dirs[i]->mbr);
			if (area_inc < min_area_inc) {
				choosed_to_insert = i;
				min_area_inc = area_inc;
			}
		}
	}

	if (choosed_to_insert != -1) {
		rtree_node *nc = myparent->dirs[choosed_to_insert];
		nc->leaves[nc->used].mbr = gmbr;
		//nc->leaves[nc->used].chull = NULL;//OGR_G_ConvexHull(g);
		nc->leaves[nc->used].geo = g;
		nc->leaves[nc->used].pgeo = NULL;
		nc->used++;
		nc->mbr = rtree_compute_mbr(nc);

		/*for(int i = 0; i < myparent->used; i++) {
			if (nc == myparent->dirs[i]) continue;
			if (ENVELOPE_CONTAINS(nc->mbr, myparent->dirs[i]->mbr) ||
			    ENVELOPE_CONTAINS(myparent->dirs[i]->mbr, nc->mbr)) {
				root->rtree_redivide_leaf(nc, myparent->dirs[i], root->m);
				return;
			}
		}*/

		int t=0;
		int lastoverlap = -1;
		int choosed = -1;
		do {
			double max_overlap = 0.0;
			lastoverlap = choosed;
			choosed = -1;

			for(int i = 0; i < myparent->used; i++) {
				if (i == choosed_to_insert) continue;
				double overlap = ENVELOPE_AREA(EnvelopeIntersection(myparent->dirs[i]->mbr, nc->mbr));

				if (overlap > max_overlap) {
					max_overlap = overlap;
					choosed = i;
				}
			}

			if (choosed != -1) {
				//printf("Run times; %d", t++);
				root->rtree_redivide_leaf(nc, myparent->dirs[choosed], root->m);
			}
		} while (choosed != -1 && choosed != lastoverlap);
	}
	else
		greedy_reinsert_on_parent(n, root, parents, g, gmbr);
}

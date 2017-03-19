
#include <float.h>
#include "rtree-star.h"
#include "rtree-reinsert.h"

int sortcomparer(const double i) {
	if (i < 0.0) return -1;
	if (i > 0.0) return 1;
	return 0;
}

int sortx_min(const void *x, const void *y) {
	return sortcomparer((*(rtree_node**)x)->mbr.MinX - (*(rtree_node**)y)->mbr.MinX);
}

int sorty_min(const void *x, const void *y) {
	return sortcomparer((*(rtree_node**)x)->mbr.MinY - (*(rtree_node**)y)->mbr.MinY);
}

int sortx_max(const void *x, const void *y) {
	return sortcomparer((*(rtree_node**)x)->mbr.MaxX - (*(rtree_node**)y)->mbr.MaxX);
}

int sorty_max(const void *x, const void *y) {
	return sortcomparer((*(rtree_node**)x)->mbr.MaxY - (*(rtree_node**)y)->mbr.MaxY);
}

int sortleafx_min(const void *x, const void *y) {
	return sortcomparer(((rtree_leaf*)x)->mbr.MinX - ((rtree_leaf*)y)->mbr.MinX);
}

int sortleafy_min(const void *x, const void *y) {
	return sortcomparer(((rtree_leaf*)x)->mbr.MinY - ((rtree_leaf*)y)->mbr.MinY);
}

int sortleafx_max(const void *x, const void *y) {
	return sortcomparer(((rtree_leaf*)x)->mbr.MaxX - ((rtree_leaf*)y)->mbr.MaxX);
}

int sortleafy_max(const void *x, const void *y) {
	return sortcomparer(((rtree_leaf*)x)->mbr.MaxY - ((rtree_leaf*)y)->mbr.MaxY);
}

int rtree_choose_subtree_star(rtree_root *root, rtree_node *n, const Envelope e) {
	int index = -1;

	if (n->dirs[0]->type == LEAF) { //if the node points to leaves. r* trick
		double minoroverlap = DBL_MAX;
		double currentarea = DBL_MAX;
		int minoridx = 0;
	
		for(int i = 0; i < n->used; i++) {
			double overlap = 0.0;
			Envelope mbrnew = e;
			ENVELOPE_MERGE(mbrnew, n->dirs[i]->mbr);

			for(int j = 0; j < n->used; j++) {
				if (j == i) continue;
				if (ENVELOPE_INTERSECTS(mbrnew, n->dirs[j]->mbr))
					overlap += ENVELOPE_AREA(EnvelopeIntersection2(mbrnew, n->dirs[j]->mbr));
			}

			double area = ENVELOPE_AREA(mbrnew);
			if (overlap < minoroverlap || (overlap == minoroverlap && area < currentarea)) {
				minoroverlap = overlap;
				currentarea = area; 
				minoridx = i;
			}
		}
		index = minoridx;
	}
	else {
		double minor_enlargement = DBL_MAX;
		double current_area = DBL_MAX;
		for(int i = 0; i < n->used; i++) {
			Envelope eaux = n->dirs[i]->mbr;
			double initarea = ENVELOPE_AREA(eaux);
			ENVELOPE_MERGE(eaux, e);
			double newarea = ENVELOPE_AREA(eaux);
			double increase = newarea - initarea;
			if (increase < minor_enlargement || (increase == minor_enlargement && newarea < current_area)) {
				index = i;
				minor_enlargement = increase;
				current_area = newarea;
			}
		}
	}
	return index;
}

rtree_leaf *rtree_split_star_choose_axis_leaf(rtree_leaf *leavesx, rtree_leaf *leavesy, int m, int usedtotal, int *minperimeter) {
	
	int totaldistrib = m - 2*minm(m)+2;
	double min_perimeter_x = DBL_MAX;
	double min_perimeter_y = DBL_MAX;

	for(int distrib = 0; distrib < totaldistrib; distrib++) {
		int split_index = minm(m) + distrib;

		Envelope group1x = leavesx[0].mbr; // first mbr of group 1 to begin with
		Envelope group1y = leavesy[0].mbr;
		for(int i=1; i < split_index; i++) {
			ENVELOPE_MERGE(group1x, leavesx[i].mbr);
			ENVELOPE_MERGE(group1y, leavesy[i].mbr);
		}

		Envelope group2x = leavesx[split_index].mbr; // first mbr of group 2 to begin with
		Envelope group2y = leavesy[split_index].mbr;
		for(int i=split_index+1; i < usedtotal; i++) {
			ENVELOPE_MERGE(group2x, leavesx[i].mbr);
			ENVELOPE_MERGE(group2y, leavesy[i].mbr);
		}

		int perim_x = ENVELOPE_PERIMETER(group1x) + ENVELOPE_PERIMETER(group2x);
		int perim_y = ENVELOPE_PERIMETER(group1y) + ENVELOPE_PERIMETER(group2y);
		min_perimeter_x = perim_x < min_perimeter_x ? perim_x : min_perimeter_x;
		min_perimeter_y = perim_y < min_perimeter_y ? perim_y : min_perimeter_y;
	}
	
	*minperimeter = min_perimeter_x < min_perimeter_y ? min_perimeter_x : min_perimeter_y;
	if (min_perimeter_x < min_perimeter_y)
		return leavesx;
	else
		return leavesy;
}

int rtree_split_star_choose_index_leaf(rtree_leaf *leaves, int m, int usedtotal) {
	int result = 0;
	double minoroverlap = DBL_MAX;
	double currentarea = DBL_MAX;
	int totaldistrib = m - 2*minm(m)+2;

	for(int distrib = 0; distrib < totaldistrib; distrib++) {
		int split_index = minm(m) + distrib;

		Envelope group1 = leaves[0].mbr;
		for(int i=1; i < split_index; i++) {
			ENVELOPE_MERGE(group1, leaves[i].mbr);
		}

		Envelope group2 = leaves[split_index].mbr;
		for(int i=split_index+1; i < usedtotal; i++) {
			ENVELOPE_MERGE(group2, leaves[i].mbr);
		}

		double overlap = ENVELOPE_AREA(EnvelopeIntersection(group1, group2));
		double area = ENVELOPE_AREA(group1) + ENVELOPE_AREA(group2);

		if (overlap < minoroverlap || (overlap == minoroverlap && area < currentarea)) {
			minoroverlap = overlap;
			currentarea = area; 
			result = split_index;
		}
	}

	assert(result > 0 && result < usedtotal);
	return result;
}

rtree_node *rtree_split_star_leaf(rtree_root *root, rtree_node *n, rtree_node *parent, const GEOSGeometryH ngeo, const Envelope e) {
	rtree_node *new_leaf = rtree_new_node(root, LEAF);
	int i;
	int usedtotal = n->used+1;
	
	GEOSGeometryH ngeo_chull = NULL;//OGR_G_ConvexHull(ngeo);

	//add new entry
	rtree_leaf *leavesx = g_new0(rtree_leaf, usedtotal);
	memcpy(leavesx, n->leaves, sizeof(rtree_leaf) * n->used);
	leavesx[n->used].geo = ngeo;
	//leavesx[n->used].chull = ngeo_chull;
	leavesx[n->used].mbr = e;

	rtree_leaf *leavesy = g_new0(rtree_leaf, usedtotal);
	memcpy(leavesy, n->leaves, sizeof(rtree_leaf) * n->used);
	leavesy[n->used].geo = ngeo;
	//leavesy[n->used].chull = ngeo_chull;
	leavesy[n->used].mbr = e;

	//choose split axis
	rtree_leaf *leaves;
	{
		//choose split axis - by max
		qsort(VPOINTER(leavesx), usedtotal, sizeof(rtree_leaf), sortleafx_max);		
		qsort(VPOINTER(leavesy), usedtotal, sizeof(rtree_leaf), sortleafy_max);
		int perimeter_max;
		rtree_leaf *leaves_max = rtree_split_star_choose_axis_leaf(leavesx, leavesy, root->m, usedtotal, &perimeter_max);
		rtree_leaf *leaves_max_copy = g_new(rtree_leaf, usedtotal);
		memcpy(leaves_max_copy, leaves_max, sizeof(rtree_leaf) * usedtotal);

		//choose split axis - by min
		qsort(VPOINTER(leavesx), usedtotal, sizeof(rtree_leaf), sortleafx_min);		
		qsort(VPOINTER(leavesy), usedtotal, sizeof(rtree_leaf), sortleafy_min);
		int perimeter_min;
		rtree_leaf *leaves_min = rtree_split_star_choose_axis_leaf(leavesx, leavesy, root->m, usedtotal, &perimeter_min);

        if (perimeter_min > perimeter_max)
            leaves = leaves_max_copy;
        else {
            leaves = leaves_min;
        }

        // free the unused leaves array
        if (leaves != leavesx) g_free(leavesx);
        if (leaves != leavesy) g_free(leavesy);
        if (leaves != leaves_max_copy) g_free(leaves_max_copy);
	}

    //choose split index
	int splitindex = rtree_split_star_choose_index_leaf(leaves, root->m, usedtotal);

	//TODO: Test planesweep_split(leaves, usedtotal);

	//distribute
	n->used = 0;
	for(i=0; i<usedtotal; i++) {
		if (i < splitindex)
			n->leaves[n->used++] = leaves[i];
		else
			new_leaf->leaves[new_leaf->used++] = leaves[i];
	}
	
	assert(n->used + new_leaf->used == root->m+1);
	
	n->mbr = rtree_compute_mbr(n);
	new_leaf->mbr = rtree_compute_mbr(new_leaf);
   
	g_free(leaves);
	return new_leaf;
}

int rtree_split_star_choose_index(rtree_node **dirs, int m, int usedtotal) {

	int result = 0;
	double minoroverlap = DBL_MAX;
	double currentarea = DBL_MAX;
	int totaldistrib = m - 2*minm(m)+2;

	for(int distrib = 0; distrib < totaldistrib; distrib++) {
		int split_index = minm(m) + distrib;
		//printf("Split index: %d\n", split_index-1);

		Envelope group1 = dirs[0]->mbr;
		for(int i=1; i < split_index; i++) {
			ENVELOPE_MERGE(group1, dirs[i]->mbr);
		}

		Envelope group2 = dirs[split_index]->mbr;
		for(int i=split_index+1; i < usedtotal; i++) {
			ENVELOPE_MERGE(group2, dirs[i]->mbr);
		}

		double overlap = ENVELOPE_AREA(EnvelopeIntersection(group1, group2));
		double area = ENVELOPE_AREA(group1) + ENVELOPE_AREA(group2);

		if (overlap < minoroverlap || (overlap == minoroverlap && area < currentarea)) {
			minoroverlap = overlap;
			currentarea = area; 
			result = split_index;
		}
	}

	assert(result > 0 && result < usedtotal);
	return result;
}

rtree_node **rtree_split_star_choose_axis(rtree_node **dirsx, rtree_node **dirsy, int m, int usedtotal, int *minperimeter) {
	
	int totaldistrib = m - 2*minm(m)+2;
	double min_perimeter_x = DBL_MAX;
	double min_perimeter_y = DBL_MAX;

	for(int distrib = 0; distrib < totaldistrib; distrib++) {
		int split_index = minm(m) + distrib;

		Envelope group1x = dirsx[0]->mbr; // first mbr of group 1 to begin with
		Envelope group1y = dirsy[0]->mbr;
		for(int i=1; i < split_index; i++) {
			ENVELOPE_MERGE(group1x, dirsx[i]->mbr);
			ENVELOPE_MERGE(group1y, dirsy[i]->mbr);
		}

		Envelope group2x = dirsx[split_index]->mbr; // first mbr of group 2 to begin with
		Envelope group2y = dirsy[split_index]->mbr;
		for(int i=split_index+1; i < usedtotal; i++) {
			ENVELOPE_MERGE(group2x, dirsx[i]->mbr);
			ENVELOPE_MERGE(group2y, dirsy[i]->mbr);
		}

		int perim_x = ENVELOPE_PERIMETER(group1x) + ENVELOPE_PERIMETER(group2x);
		int perim_y = ENVELOPE_PERIMETER(group1y) + ENVELOPE_PERIMETER(group2y);
		min_perimeter_x = perim_x < min_perimeter_x ? perim_x : min_perimeter_x;
		min_perimeter_y = perim_y < min_perimeter_y ? perim_y : min_perimeter_y;
	}
	
	*minperimeter = min_perimeter_x < min_perimeter_y ? min_perimeter_x : min_perimeter_y;
	if (min_perimeter_x < min_perimeter_y)
		return dirsx;
	else
		return dirsy;
}

rtree_node *rtree_split_star(rtree_root *root, rtree_node *n, rtree_node *newn, rtree_node *parent, int m) {
	rtree_node *new_dir = rtree_new_node(root, DIRECTORY);
	int i;
	int usedtotal = n->used+1;
	
#ifdef DEBUG
	for(i=0; i<n->used; i++) printf("%p ", n->dirs[i]);
	printf("\n");
	for(i=0; i<n->used; i++) printf("%f ", n->dirs[i]->mbr.MinX);
	printf("\n");
#endif

	//add new entry (being added to node)
	rtree_node **dirsx = g_new(rtree_node*, usedtotal);
	memcpy(dirsx, n->dirs, sizeof(rtree_node*) * n->used);
	dirsx[n->used] = newn;

	rtree_node **dirsy = g_new(rtree_node*, usedtotal);
	memcpy(dirsy, n->dirs, sizeof(rtree_node*) * n->used);
	dirsy[n->used] = newn;
		
	//choose split axis
	rtree_node **dirs;
	{
		//choose split axis - by max
		qsort(VPOINTER(dirsx), usedtotal, sizeof(rtree_node*), sortx_max);
		qsort(VPOINTER(dirsy), usedtotal, sizeof(rtree_node*), sorty_max);
		int perimeter_max;
		rtree_node **dirs_max = rtree_split_star_choose_axis(dirsx, dirsy, m, usedtotal, &perimeter_max);
		rtree_node **dirs_max_copy = g_new(rtree_node*, sizeof(rtree_node*) * usedtotal);
		memcpy(dirs_max_copy, dirs_max, sizeof(rtree_node*) * usedtotal);

		//choose split axis - by min 
		qsort(VPOINTER(dirsx), usedtotal, sizeof(rtree_node*), sortx_min);
		qsort(VPOINTER(dirsy), usedtotal, sizeof(rtree_node*), sorty_min);
		int perimeter_min;
		rtree_node **dirs_min = rtree_split_star_choose_axis(dirsx, dirsy, m, usedtotal, &perimeter_min);

		if (perimeter_min > perimeter_max)
			dirs = dirs_max_copy;
		else {
			dirs = dirs_min;
		}

		// free the unused dirs array
		if (dirs != dirsx) g_free(dirsx);
		if (dirs != dirsy) g_free(dirsy);
		if (dirs != dirs_max_copy) g_free(dirs_max_copy);
	}

	//choose split index
	int splitindex = rtree_split_star_choose_index(dirs, m, usedtotal);
	
	//distribute
	n->used = 0;
	for(i=0; i<usedtotal; i++) {
		if (i < splitindex)
			n->dirs[n->used++] = dirs[i];
		else
			new_dir->dirs[new_dir->used++] = dirs[i];
	}
	
	assert(n->used + new_dir->used == m+1);
	
	n->mbr = rtree_compute_mbr(n);
	new_dir->mbr = rtree_compute_mbr(new_dir);
	
#ifdef DEBUG
	printf("Dirs on existing: ");
	for(i=0; i<n->used; i++) printf("%p ", n->dirs[i]);
	printf("\n");
	
	printf("Dirs on new: ");
	for(i=0; i<new_dir->used; i++) printf("%p ", new_dir->dirs[i]);
	printf("\n");

	double finaloverlap = ENVELOPE_AREA(EnvelopeIntersection(n->mbr, new_dir->mbr));
	printf("Final overlap %f\n\n", finaloverlap);
#endif
	
	g_free(dirs);
	return new_dir;
}

rtree_root *rtree_new_rstar(const int m, const int servers) {
	return rtree_new(
		m,
		servers,
		rtree_choose_subtree_star,
		rtree_split_star,
		rtree_split_star_leaf,
		rstar_reinsert_on_parent,
		NULL, NULL);
}



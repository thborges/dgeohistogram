#include <ogr_core.h>
#include <assert.h>
#include "glibwrap.h"
#include "rtree.h"
#include "ogrext.h"

void rtree_window_search_recursive(rtree_node *node, const int m, const rtree_window *window, GList **results, rtree_window_stat *stats, char onlydirs) {
	char accessed = FALSE;
	if (node->type == DIRECTORY) {
		for(int i = 0; i < node->used; i++) {
			rtree_node *dir = node->dirs[i];
			/*GEOSGeometryH geo = dir->mbrgeo ? dir->mbrgeo : (dir->mbrgeo = EnvelopeToGeometry(dir->mbr));
			if (OGR_G_Intersects(geo, window->geo)) {
				rtree_window_search_recursive(dir, window, results, stats, onlydirs);
				accessed = TRUE;
				(node->server == dir->server) ? stats->lmessages++ : stats->rmessages++;
			}*/
			if (ENVELOPE_INTERSECTS(dir->mbr, window->mbr)) {
				rtree_window_search_recursive(dir, m, window, results, stats, onlydirs);
				accessed = TRUE;
				(node->server == dir->server) ? stats->lmessages++ : stats->rmessages++;
			}
		}
		if (accessed) 
			stats->truedirs++;
		else
			stats->falsedirs++;
	}
	else if (onlydirs == TRUE) {
		for(int i=0; i<node->used; i++)
			*results = g_list_prepend(*results, &node->leaves[i]);
	}
	else {
		for(int i=0; i<node->used; i++) {
			stats->geomchecked++;
			//int chull_intersects = node->leaves[i].chull == NULL ||
			//	GEOSIntersects(node->leaves[i].chull, window->chull);
			if (/*chull_intersects && */GEOSPreparedIntersects(window->pgeo, node->leaves[i].geo)) {
				accessed = TRUE;
				*results = g_list_prepend(*results, &node->leaves[i]);
			}
		}
		if (accessed)
			stats->trueleaves += node->used < m ? 0.5 : 1.0;
		else
			stats->falseleaves += node->used < m ? 0.5 : 1.0;
	}
}

GList *rtree_window_search(const rtree_root *root, const GEOSGeometryH window, const EnvelopeC windowenv, rtree_window_stat *stats) {
	GList *results = NULL;
	
	rtree_window swin;
	swin.geo = window;
	swin.pgeo = GEOSPrepare(swin.geo);
	swin.mbr = windowenv;
	swin.chull = NULL;//OGR_G_ConvexHull(window);

	rtree_window_search_recursive(root->root, root->m, &swin, &results, stats, FALSE);
	
	//OGR_G_DestroyGeometry(swin.chull);

	return results;
}

GList *rtree_window_csearch(const rtree_root *root, const GEOSGeometryH window, const EnvelopeC windowenv, rtree_window_stat *stats) {
	GList *results = NULL;
	
	rtree_window swin;
	swin.geo = window;
	swin.pgeo = GEOSPrepare(swin.geo);
	swin.mbr = windowenv;
	swin.chull = NULL;//OGR_G_ConvexHull(window);

	rtree_window_search_recursive(root->root, root->m, &swin, &results, stats, TRUE);
	
	//OGR_G_DestroyGeometry(swin.chull);

	return results;
}



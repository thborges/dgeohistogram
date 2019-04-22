#include <ogr_core.h>
#include <assert.h>
#include "glibwrap.h"
#include "rtree.h"
#include "ogrext.h"
#include "dataset.h"
#include "join.h"

void rtree_join_inl_print_stats(rtree_join_stat *stats, size_t results) {
	if (stats->lines == 0)
		printf("\nresults   |io_pages     |cache_miss   |mbr compared |geos compared|  mbrs false |  geos false | queue1 msize| queue2 msize|\n");

	fprintf(stderr, "%'10ld|%'13ld|%'12.4f%%|%'13d|%'13d|%'13d|%'13d|%'13d|%'13d|\r",
		results, 
		stats->io_pages, 
		stats->cache_miss*100,
		stats->mbrs_compared,
		stats->geos_compared, 
		stats->mbrs_compared_false,
		stats->geos_compared_false, 
		stats->queue1_size, 
		stats->queue2_size);

	stats->lines++;
}

dataset *rtree_join_inl (dataset *r, const rtree_root *s, rtree_join_stat *stats, enum JoinPredicateCheck pcheck) {

	dataset *results = dataset_create("", 2);
    
	int queue2_size = 0;

	dataset_iter pair;
    dataset_foreach(pair, r) {
		queue2_size++;
		rtree_leaf *lgeo = get_join_pair_leaf(pair.item, pcheck);
		assert(lgeo->geo != NULL);
		
		// window search - breadth first traversal
		GList *wresults = NULL;
		GList *parents = g_list_prepend(NULL, s->root);
		while (parents != NULL) {
			rtree_node *n = parents->data;
			parents = g_list_delete_link(parents, parents);
			if (n->type == DIRECTORY) {
				for(int j=0; j<n->used; j++) {
					stats->mbrs_compared++;
					if (ENVELOPE_INTERSECTS(lgeo->mbr, n->dirs[j]->mbr)) {
						parents = g_list_prepend(parents, n->dirs[j]);
					}
					else
						stats->mbrs_compared_false++;
				}
			}
			else {
				for(int j=0; j<n->used; j++)
					wresults = g_list_prepend(wresults, &n->leaves[j]);
			}
		}

		//printf("(%d: %d) ", queue2_size, g_list_length(wresults));

		// Refinement
		GList *geo;
		for(geo = wresults; geo != NULL; geo = g_list_next(geo)) {
			rtree_leaf *nret = geo->data;
			//lrubuffer_get(lru, nret);

			stats->mbrs_compared++;
			int intersects = ENVELOPE_INTERSECTS(lgeo->mbr, nret->mbr); 
			if (!intersects) {
				stats->mbrs_compared_false++;
				continue;
			}
		
			// get the prepared geometry, or prepare if none prepared
			const GEOSPreparedGeometry *pgeo = lgeo->pgeo ? lgeo->pgeo : nret->pgeo;
			GEOSGeometryH geo = !lgeo->pgeo ? lgeo->geo : nret->geo;
			if (!pgeo) {
				pgeo = lgeo->pgeo = GEOSPrepare(lgeo->geo);
				geo = nret->geo;
			}

			assert(pgeo != NULL);

			stats->geos_compared++;
			if (GEOSPreparedIntersects(pgeo, geo)) {
				rtree_leaf *leaves = dataset_add(results);
				leaves[0] = *lgeo;
				leaves[1] = *nret;

				if (results->metadata.count % 1000 == 0)
					rtree_join_inl_print_stats(stats, results->metadata.count);
			}
			else {
				stats->geos_compared_false++;
			}
		}
	
		g_list_free(wresults);

		stats->queue2_size = queue2_size;
		//stats->io_pages = lru->get_extmem;
		//stats->cache_miss = (double)lru->get_extmem / (lru->get_extmem + lru->get_buffer);
	};

	rtree_join_inl_print_stats(stats, results->metadata.count);
	return results;
}


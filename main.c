
#include <ogr_api.h>
#include <geos_c.h>
#include <cpl_conv.h>
#include <time.h>
#include <arpa/inet.h>
#include <locale.h>
#include "glibwrap.h"
#include "rtree.h"
#include "utils.h"
#include "rtree-test.h"
#include "rtree-star.h"
#include "rtree-gut.h"
#include "rtree-lazy.h"
#include "geosext.h"
#include "dataset.h"

enum IndexType {
	RSTAR,
	RLAZY,
	R0,
	RGUT,
};

rtree_root *build_rtree(char *shpfile, int m_size, int cluster_size, int shapes, enum IndexType rtreetype);
dataset *read_geos(char *shpfile);

void window_search(rtree_root *rtree, char* windowfile);
void gpu_test(rtree_root *rtree);


void cuda_host_buffer_init();
void cuda_host_free_buffer();

void GEOSAuxDestroy(gpointer g) {
	GEOSGeom_destroy((GEOSGeometry*)g);
}

void rtree_join_pair_free(gpointer g) {
	g_free(g);
}

GList *joinplan = NULL;

int main(int argc, char* argv[]) {
	//srand(time(NULL));

	if (argc < 4) {
		printf(
"Use: %s file.shp m_size cluster_size [-rgut] [-rlazy] [-r0] [-s] [-g] [-gpu]\n", argv[0]);
		return 1;
	}

	// fix to print thousand separator
	setlocale(LC_NUMERIC, "");

	int m_size = atoi(argv[2]);
	if (m_size <= 1) {
		fprintf(stderr, "Invalid m_size value.\n");
		return 1;
	}

	int cluster_size = atoi(argv[3]);
	if (cluster_size <= 0) {
		fprintf(stderr, "Invalid cluster size value.\n");
		return 1;
	}

	bool build = false;
	bool search = false;
	bool shapes = false;
	bool gpu = false;
	enum IndexType rtreetype = RSTAR;

	char *joinshp = NULL;
	char *windowfile = NULL;
	bool skip = false;
	bool added_first = false;

	if (argc >= 4) {

		OGRRegisterAll();
		initGEOS(geos_messages, geos_messages);
		//cuda_host_buffer_init();
		lru = lrubuffer_new(512);
	
		printf("%4s|%4s|%10s|%10s|%10s|%10s|%10s|%10s|%10s|%10s| %s\n", "M", "Clus", "Items", "Dirs", "Leaves", "Indexed", "Fill", "OverlapD", "OverlapL", "Time (ms)", "File");

		for(int i = 4; i < argc; i++) {
			if (skip) {
				skip  = false;
				continue;
			}

			if (strcmp(argv[i], "-b") == 0) build = true;
			if (strcmp(argv[i], "-s") == 0) {
				search = true;
				if (argc > i+1) {
					windowfile = argv[i+1];
				}
				else {
					fprintf(stderr, "Please specify the windows file.");
				}
			}

			if (strcmp(argv[i], "-g") == 0) shapes = true;
			if (strcmp(argv[i], "-gpu") == 0) gpu = true;
			
			if (strcmp(argv[i], "-rgut") == 0) rtreetype = RGUT;
			if (strcmp(argv[i], "-rlazy") == 0) rtreetype = RLAZY;
			if (strcmp(argv[i], "-r0") == 0) rtreetype = R0;
		}
	}

	rtree_root *rtree = NULL;

	if ((build || search) && !added_first)
		rtree = build_rtree(argv[1], m_size, cluster_size, shapes, rtreetype);

	//search
	if (search)
		window_search(rtree, windowfile);

	printf("Finishing...\n");
	//getchar();

	//cuda_host_free_buffer();
	lrubuffer_destroy(lru);

	finishGEOS();

	return 0;
}

rtree_root *build_rtree(char *shpfile, int m_size, int cluster_size, int shapes, enum IndexType rtreetype) {
	OGRDataSourceH ds;
	ds = OGROpen(shpfile, false, NULL);
	if (ds == NULL) {
		fprintf(stderr, "Invalid shapefile: %s\n", shpfile);
		return NULL;
	}

	OGRLayerH layer = OGR_DS_GetLayer(ds, 0);
	OGR_L_ResetReading(layer);
	
	// start the r-tree building
	clock_t cs = clock();

	rtree_root *rtree = NULL;
	switch (rtreetype) {
		case RSTAR:
			rtree = rtree_new_rstar(m_size, cluster_size);
			break;

		case RLAZY:
			rtree = rtree_new_rlazy(m_size, cluster_size);
			break;

		case R0:
			rtree = rtree_new_r0(m_size, cluster_size);
			break;
		case RGUT:
			rtree = rtree_new_rtree_gut_quad(m_size, cluster_size);
			break;

	}

	size_t wkb_current_size = 4096;
	unsigned char *wkb = g_new(unsigned char, wkb_current_size);

	
	OGRFeatureH feature;
	OGRGeometryH geometry;
	int i = 0;
	while((feature = OGR_L_GetNextFeature(layer)) != NULL) {

		geometry = OGR_F_GetGeometryRef(feature);
		if (geometry != NULL) {
			size_t wkb_size = OGR_G_WkbSize(geometry);
			if (wkb_size > wkb_current_size)
				wkb = g_renew(unsigned char, wkb, wkb_size);

			OGR_G_ExportToWkb(geometry, (OGRwkbByteOrder)(htonl(1) == 1 ? 0 : 1), wkb);

			GEOSGeometry *ggeo = GEOSGeomFromWKB_buf(wkb, wkb_size);
			//GEOSGeometry *cggeo = GEOSGetCentroid(ggeo);	
			
			if (ggeo) {
				rtree_append(rtree, ggeo);
			}
			i++;


		}

		OGR_F_Destroy(feature);
	}
	g_free(wkb);

	clock_t cf = clock();
	double runtime_diff_ms = (cf-cs) * 1000. / CLOCKS_PER_SEC;

	count_node_struct c = count_nodes(rtree);
	printf("%4d|%4d|%10d|%10d|%10d|%10d| %-2.5f%%|", m_size, cluster_size, i, c.dirs,
		(c.leaves + (int)ceil(c.subleaves)), c.items,
		((double)c.items / rtree->m) / (c.leaves + ceil(c.subleaves)) * 100);

	check_mbrs(rtree);

	print_overlap(rtree, OGR_G_GetSpatialReference(geometry));
	printf("%10.0f| %s\n", runtime_diff_ms, shpfile);

	if (shapes) {
		printf("Generating shape files with MBRs...\n");
		write_dirs_shapefile(rtree);
	}

	OGR_DS_Destroy(ds);
	return rtree;
}

void window_search(rtree_root *rtree, char *windowfile) {
	int wi = 0;
	char wkt[1024];
	FILE *file = fopen(windowfile, "r");

	if (!file) {
		fprintf(stderr, "Cannot open file %s.\n", windowfile);
		return;
	}
	else
		printf("Searching for windows...\n");

	rtree_window_stat stats;

	printf("window|results   |true_dirs |trueleaf  |falsedirs |falseleaf |lmessages |rmessages |geomchecked |        us|\n");

	double runtime_total_time = 0.0;
	while (fgets(wkt, sizeof(wkt), file)) {
	 	char *wktp = &wkt[0];
		GEOSGeometryH win = GEOSGeomFromWKT(wktp);
		
		memset(&stats, 0, sizeof(rtree_window_stat));

		printf("W%05d|", wi);
	
		clock_t cs = clock();
		GList *results = rtree_window_search(rtree, win, &stats);
		clock_t cf = clock();
		double runtime_diff_us = (cf-cs) * 1000000. / CLOCKS_PER_SEC;


		printf("%10d|%10d|%10d|%10d|%10d|%10d|%10d|%12d|%10.1f|\n",
				g_list_length(results), stats.truedirs, (int)ceil(stats.trueleaves),
				stats.falsedirs, (int)ceil(stats.falseleaves), stats.lmessages,
				stats.rmessages, stats.geomchecked, runtime_diff_us);

		runtime_total_time += runtime_diff_us;

		g_list_free(results);
		GEOSGeom_destroy(win);
		wi++;
	}
	fclose(file);

	printf("%96s:%10.1f|\n", "Total time", runtime_total_time);
}

dataset *read_geos(char *shpfile) {
	dataset *results = dataset_create_mem("", 1);

	OGRDataSourceH ds;
	ds = OGROpen(shpfile, false, NULL);
	if (ds == NULL) {
		fprintf(stderr, "Invalid shapefile.\n");
		return results;
	}

	OGRLayerH layer = OGR_DS_GetLayer(ds, 0);
	OGR_L_ResetReading(layer);

	clock_t cs = clock();
	
	OGRFeatureH feature;
	OGRGeometryH geometry;
	while((feature = OGR_L_GetNextFeature(layer)) != NULL) {
		geometry = OGR_F_GetGeometryRef(feature);
		if (geometry != NULL) {
			size_t wkb_size = OGR_G_WkbSize(geometry);
			unsigned char* wkb = g_new(unsigned char, wkb_size);
			OGR_G_ExportToWkb(geometry, (OGRwkbByteOrder)(htonl(1) == 1 ? 0 : 1), wkb);

			GEOSGeometryH ggeo = GEOSGeomFromWKB_buf(wkb, wkb_size);

			if (ggeo) {
				rtree_leaf *leaves = dataset_add(results);
				leaves[0].geo = ggeo;
				GEOSGetEnvelope(ggeo, &leaves[0].mbr);

				// Test OGRGet
				/*OGREnvelope e;
				OGR_G_GetEnvelope(geometry, &e);

				Envelope f;
				GEOSGetEnvelope(ggeo, &f);

				if (!(e.MinX == f.MinX &&
					  e.MinY == f.MinY &&
					  e.MaxX == f.MaxX &&
					  e.MaxY == f.MaxY))
					printf("Diferente!\n(%f %f, %f %f) e \n(%f %f, %f %f)\n", e.MinX, e.MinY, e.MaxX, e.MaxY, f.MinX, f.MinY, f.MaxX, f.MaxY);
*/
			}
		}
		OGR_F_Destroy(feature);
	}
	OGR_DS_Destroy(ds);

	clock_t cf = clock();
	double runtime_diff_us = (cf-cs) * 1000. / CLOCKS_PER_SEC;

	printf("%4s|%4s|%'10d|%10s|%10s|%10s|%10s|%10s|%10s|%10.1f| %s\n", "N", "N", results->metadata.count, "N", "N", "N", "N", "N", "N", runtime_diff_us, shpfile);

	return results;
}



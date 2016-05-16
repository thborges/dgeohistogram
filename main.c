
#include <libgen.h>
#include <ogr_api.h>
#include <geos_c.h>
#include <cpl_conv.h>
#include <time.h>
#include <arpa/inet.h>
#include <locale.h>
#include "dataset.h"
#include "glibwrap.h"
#include "ogrext.h"
#include "rtree.h"
#include "rtree-star.h"
#include "histogram.h"
#include "minskew.h"

dataset *read_geos(char *shpfile);

OGRDataSourceH ogr_ds;

int main(int argc, char* argv[]) {
	//srand(time(NULL));
	OGRRegisterAll();
	initGEOS(geos_messages, geos_messages);

	if (argc < 4) {
		printf("Use: %s [mbrc, centr, areaf, areafs] file.shp queries.shp\n", argv[0]);
		return 1;
	}

	enum HistogramHashMethod hm = HHASH_AREAFRAC;
	if (strcmp(argv[1], "mbrc") == 0)
		hm = HHASH_MBRCENTER;
	else
	if (strcmp(argv[1], "centr") == 0)
		hm = HHASH_CENTROID;
	else
	if (strcmp(argv[1], "areaf") == 0)
		hm = HHASH_AREAFRAC;
	else
	if (strcmp(argv[1], "areafs") == 0)
		hm = HHASH_AREAFRACSPLIT;
	else {
		printf("Method %s does not exists.\n", argv[1]);
		exit(1);
	}
	
	dataset *ds = read_geos(argv[2]);

	// chamar a função que cria o histograma
	histogram_generate(ds, hm, 0);
	histogram_print_geojson(ds);
	histogram_print(ds, CARDIN);

	// create min skew histogram
	GList *minskewh = minskew_generate_hist(ds, 5000);
	minskew_print_hist(ds, minskewh);

	// cria uma r*
	rtree_root *rtree = NULL;
	rtree = rtree_new_rstar(30, 10);

	unsigned read = 0;
	dataset_iter_seg iter;
	dataset_foreach(iter, ds) {
		read++;
		GEOSGeometryH geo = dataset_get_leaf_geo(ds, iter.item);
		rtree_append(rtree, geo);
		print_progress_gauge(read, ds->metadata.count);
	}

	double accuracy = 0.0;
	double sum_ei = 0.0;
	double sum_ri = 0.0;
	double mean = 0.0;
	double M2 = 0.0;
	double sum_error = 0.0;
	int n = 0;

	dataset_histogram *hist = &ds->metadata.hist;
	int cells = hist->xqtd*hist->yqtd;

	rtree_window_stat stats;
	dataset *queries = read_geos(argv[3]);
	dataset_iter it;
	dataset_foreach(it, queries) {
		n++;

		GEOSGeometryH geoquery = dataset_get_leaf_geo(queries, it.item);
		Envelope query = it.item->mbr;

		memset(&stats, 0, sizeof(rtree_window_stat));
		GList *results = rtree_window_search(rtree, geoquery, &stats);

		// real cardinality from rtree
		int riq = g_list_length(results);

		// histogram estimate cardinality
		//int rhq = histogram_search_hist(&ds->metadata.hist, query);
		int rhq = minskew_search_hist(minskewh, query);

		int error = abs(rhq-riq);

		// average relative error
		sum_ei += error;
		sum_ri += riq;

		// stdev of error
		double delta = error - mean;
		mean += delta/(double)n;
		M2 += delta*(error - mean);

		// sum error
		sum_error += error;

		// precision and recall
		int pv, pf, nf;
		if (rhq >= riq) {
			pv = riq;
			pf = rhq-riq;
			nf = 0;
		}
		else {
			pv = rhq;
			pf = 0;
			nf = riq-rhq;
		}
		double p = pv==0 && pf==0 ? 1.0 : pv / (double)(pv+pf);
		double r = pv==0 && nf==0 ? 1.0 : pv / (double)(pv+nf);
		accuracy += (p+r)/2.0;
		//printf(" %.2f", (p+r)/2.0);
		//printf(" %d:%d", rhq, riq);
		//printf(" %d", error);

		if (isnan(p+r))
			printf(": pv%d pf%d nf%d riq%d rhq%d\n", pv, pf, nf, riq, rhq);

		g_list_free(results);

		print_progress_gauge(n, cells);

		//printf("\n");
	}
	
	printf("Average Relative Error: %f\nStdevp error: %f\nError sum: %f\n", 
		sum_ei / (double)sum_ri,
		sqrt(M2/(double)n),
		sum_error);
	
	OGR_DS_Destroy(ogr_ds);
	finishGEOS();

	return 0;
}

dataset *read_geos(char *shpfile) {
	dataset *results = dataset_create_mem(basename(shpfile), 1);

	ogr_ds = OGROpen(shpfile, false, NULL);
	if (ogr_ds == NULL) {
		fprintf(stderr, "Invalid shapefile.\n");
		return results;
	}

	OGRLayerH layer = OGR_DS_GetLayer(ogr_ds, 0);
	OGR_L_ResetReading(layer);
	int layer_count = OGR_L_GetFeatureCount(layer, FALSE);

	results->temp_ogr_layer = layer;

	clock_t cs = clock();
	
	unsigned read = 0;
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
				dataset_leaf *leaf = dataset_add(results);
				long long gid = OGR_F_GetFID(feature);
				Envelope mbr = OGRGetEnvelope(geometry);
				dataset_fill_leaf_id(leaf, 0, gid, &mbr);
				leaf[0].points = OGRGetNumPoints(geometry);
			}
		}
		OGR_F_Destroy(feature);

		read++;
		print_progress_gauge(read, layer_count);
	}

	clock_t cf = clock();
	double runtime_diff_us = (cf-cs) * 1000. / CLOCKS_PER_SEC;

	printf("%4s|%4s|%'10d|%10s|%10s|%10s|%10s|%10s|%10s|%10.1f| %s\n", "N", "N", results->metadata.count, "N", "N", "N", "N", "N", "N", runtime_diff_us, shpfile);

	return results;
}


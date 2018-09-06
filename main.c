
#include <cpl_conv.h>
#include <time.h>
#include <arpa/inet.h>
#include <locale.h>
#include <libgen.h>

#include "dataset.h"
#include "glibwrap.h"
#include "ogrext.h"
#include "rtree.h"
#include "rtree-star.h"
#include "histogram.h"
#include "minskew.h"
#include "ehist.h"
#include "dataset_specs.h"

char *dataset_name;

dataset *read_geos(char *shpfile);

OGRDataSourceH ogr_ds;

enum HistogramType {
	HTGRID,
	HTMINSKEW,
	HTEULER
};

int main(int argc, char* argv[]) {

	//srand(time(NULL));

	OGRRegisterAll();
	initGEOS(geos_messages, geos_messages);

	if (argc < 3) {
		printf("Use: %s [grid minskew euler] [mbrc, centr, areaf, areafs] [fix x y, avg, avgstd] file.shp size%%query\n", argv[0]);
		return 1;
	}

	enum HistogramType ht = HTGRID;
	if (strcmp(argv[1], "grid") == 0)
		ht = HTGRID;
	else if (strcmp(argv[1], "minskew") == 0)
		ht = HTMINSKEW;
	else if(strcmp(argv[1], "euler") == 0)
		ht = HTEULER;
	else {
		printf("Histogram type %s does not exists.\n", argv[1]);
		exit(1);
	}

	enum HistogramHashMethod hm = HHASH_AREAFRAC;
	if (ht == HTGRID || ht == HTMINSKEW) {
		if (strcmp(argv[2], "mbrc") == 0)
			hm = HHASH_MBRCENTER;
		else
		if (strcmp(argv[2], "centr") == 0)
			hm = HHASH_CENTROID;
		else
		if (strcmp(argv[2], "areaf") == 0)
			hm = HHASH_AREAFRAC;
		else
		if (strcmp(argv[2], "areafs") == 0)
			hm = HHASH_AREAFRACSPLIT;
		else {
			printf("Hash method %s does not exists.\n", argv[2]);
			exit(1);
		}
	}

	int argatu = 3;
	enum HistogramSplitMethod sm = HSPLIT_AVG_STD;
	int xqtd = 0, yqtd = 0;
	if (strcmp(argv[argatu], "fix") == 0) {
		sm = HSPLIT_FIX;
		argatu++;
		xqtd = atoi(argv[argatu++]);
		yqtd = atoi(argv[argatu++]);
	}
	else if (strcmp(argv[argatu], "avgstd") == 0) {
		argatu++;
		sm = HSPLIT_AVG_STD;
	}
	else if (strcmp(argv[argatu], "avg") == 0) {
		argatu++;
		sm = HSPLIT_AVG;
	}
	else {
		printf("Method %s does not exists.\n", argv[argatu]);
		exit(1);
	}

	// split disabled
	int split_qtd = 1;
	
	dataset_name = argv[argatu++];
	dataset *ds = read_geos(dataset_name);

	HistogramGenerateSpec spec;
	spec.hm = hm;
	spec.sm = sm;
	spec.xqtd = xqtd;
	spec.yqtd = yqtd;

	// chamar a função que cria o histograma
	histogram_generate(ds, spec, CHECKR, split_qtd);
	histogram_print_geojson(ds);
	histogram_print(ds, CARDIN);

	//print_dataset_specs(&ds->metadata.hist);

	// create min skew histogram
	GList *minskewh = NULL;
	if (ht == HTMINSKEW) {
		minskewh = minskew_generate_hist(ds, 300);
		minskew_print_hist(ds, minskewh);
	}

	// create euler histogram
	euler_histogram *eh = NULL;
	if (ht == HTEULER) {
		eh = eh_generate_hist(ds, spec, CHECKR);
		euler_print_hist(ds, eh);
	}

	// the user specified a query?
	if (argc <= argatu)
		return 0;

	double query_size = atof(argv[argatu++]);

	
	

	// cria uma r*
	rtree_root *rtree = NULL;
	rtree = rtree_new_rstar(30, 10);

	unsigned read = 0;
	dataset_iter_seg iter;
	dataset_foreach(iter, ds) {
		read++;
		GEOSGeometryH geo = dataset_get_leaf_geo(ds, iter.item);
		rtree_append(rtree, geo);
		//print_progress_gauge(read, ds->metadata.count);
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
	double width = ds->metadata.hist.mbr.MaxX - ds->metadata.hist.mbr.MinX;
	double height = ds->metadata.hist.mbr.MaxY - ds->metadata.hist.mbr.MinY;
	double wsize = width * query_size;
	double hsize = height * query_size;
	//int qtd = (width / wsize) * (height/hsize) / 2.0;
	int qtd = 500;

	//printf("Query count: %d, w %f, h %f, w_size %f, h_size %f\n", qtd, width, height, wsize, hsize);
	//print_geojson_header();

	while (n < qtd) {
		n++;

		Envelope query;
		query.MinX = ds->metadata.hist.mbr.MinX;
		query.MinY = ds->metadata.hist.mbr.MinY;
		query.MinX += width * (rand()/(double)RAND_MAX);
		query.MinY += height * (rand()/(double)RAND_MAX);
		query.MaxX = query.MinX + wsize;
		query.MaxY = query.MinY + hsize;
		//print_geojson_mbr(query, "0");

	    char wkt[512];
    	sprintf(wkt, "POLYGON((%e %e, %e %e, %e %e, %e %e, %e %e))",
        	query.MinX, query.MinY,
        	query.MaxX, query.MinY,
        	query.MaxX, query.MaxY,
        	query.MinX, query.MaxY,
        	query.MinX, query.MinY);
		
		GEOSGeometryH geoquery = GEOSGeomFromWKT(wkt);

		memset(&stats, 0, sizeof(rtree_window_stat));
		GList *results = rtree_window_search(rtree, geoquery, &stats);

		// real cardinality from rtree
		int riq = g_list_length(results);

		// histogram estimate cardinality
		int rhq = 0;
		if (ht == HTGRID)
			rhq = histogram_search_hist(&ds->metadata.hist, query);
		else if (ht == HTMINSKEW)
			rhq = minskew_search_hist(minskewh, query);
		else if (ht == HTEULER)
			rhq = euler_search_hist(eh, query);

		//printf("Query %d: r: %5d, e: %5d, %5d\n", n, riq, rhq, rhq - riq);

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

		//print_progress_gauge(n, cells);

		//printf("\n");
	}
	//print_geojson_footer();
	
	printf("\nSize\tARE\t\tSTD\t\tSUM\t\tMethod\tName\n");
	printf("%3.2f\t%f\t%f\t%f\t%s\t%s\n",
		query_size, 
		sum_ei / (double)sum_ri,
		sqrt(M2/(double)n),
		sum_error,
		argv[1],
		ds->metadata.name);

	//print result to csv file data.csv in dgeohistogram
	char filename[100] = "data.csv";
	FILE *file;
	file = fopen(filename, "a");
  	fprintf(file, "%3.2f,%f,%f,%f,%s,%s\n", 
		query_size, 
		sum_ei / (double)sum_ri,
		sqrt(M2/(double)n),
		sum_error,
		argv[1],
		ds->metadata.name);
	fclose(file);
	
finish:	
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
				//printf("%lf %lf %lf %lf\n", mbr.MinX, mbr.MinY, mbr.MaxX, mbr.MaxY);
				dataset_fill_leaf_id(leaf, 0, gid, &mbr);
				leaf[0].points = OGRGetNumPoints(geometry);
			}
		}
		OGR_F_Destroy(feature);

		read++;
		//print_progress_gauge(read, layer_count);
	}

	clock_t cf = clock();
	double runtime_diff_us = (cf-cs) * 1000. / CLOCKS_PER_SEC;

	//printf("%4s|%4s|%'10d|%10s|%10s|%10s|%10s|%10s|%10s|%10.1f| %s\n", "N", "N", results->metadata.count, "N", "N", "N", "N", "N", "N", runtime_diff_us, shpfile);

	return results;
}


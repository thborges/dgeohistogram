
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

dataset *read_geos(char *shpfile, OGRDataSourceH ogr_ds);

OGRDataSourceH ogr_ds1;
OGRDataSourceH ogr_ds2;

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
		printf("Use: %s [grid minskew euler] [mbrc, centr, areaf, areafs] [fix x y, avg, avgstd] file1.shp  file2.shp size%%query\n", argv[0]);
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
	if (strcmp(argv[2], "mbrc") == 0)
		hm = HHASH_MBRCENTER;
	else if (strcmp(argv[2], "centr") == 0)
		hm = HHASH_CENTROID;
	else if(strcmp(argv[2], "areaf") == 0)
		hm = HHASH_AREAFRAC;
	else if(strcmp(argv[2], "areafs") == 0)
		hm = HHASH_AREAFRACSPLIT;
	else {
		printf("Histogram hash %s does not exists.\n", argv[2]);
		exit(1);
	}

	int argatu = 3;
	int xqtd = 0, yqtd = 0;
	enum HistogramSplitMethod sm = HSPLIT_AVG_STD;
    if (strcmp(argv[argatu], "fix") == 0) {
		sm = HSPLIT_FIX;
		argatu++;
		xqtd = atoi(argv[argatu++]);
		yqtd = atoi(argv[argatu++]);
	} else if (strcmp(argv[argatu], "avg") == 0) {
		sm = HSPLIT_AVG;
	} else if (strcmp(argv[argatu], "avgstd") == 0) {
		sm = HSPLIT_AVG_STD;
	} else {
		printf("Histogram split %s does not exists.\n", argv[argatu]);
		exit(1);
	}

    //histograma a
    dataset_name = argv[argatu++];
    printf("Carregando arquivo: %s\n", dataset_name);
    dataset *dsA = read_geos(dataset_name, ogr_ds1);

    //histograma b
    dataset_name = argv[argatu++];
    printf("Carregando arquivo: %s\n", dataset_name);
    dataset *dsB = read_geos(dataset_name, ogr_ds2);

	HistogramGenerateSpec spec;
    spec.hm = hm;
    spec.sm = sm;
    spec.xqtd = xqtd;
    spec.yqtd = yqtd;
	spec.split_quantity = 0; // only for HHASH_AREAFRACSPLIT (Isabella)

	// chamar a função que cria o histograma
	printf("Criando histograma para: %s\n", dsA->metadata.name);
	histogram_generate(dsA, spec, 0);
	histogram_print_geojson(dsA);

	printf("Criando histograma para: %s\n", dsB->metadata.name);
	histogram_generate(dsB, spec, 0);
	histogram_print_geojson(dsB);


	printf("\n\n");
	euler_histogram *ehA = NULL;
	euler_histogram *ehB = NULL;

	if (ht == HTEULER) {
		ehA = eh_generate_hist(dsA, spec, CHECKR);
		euler_print_hist(dsA->metadata.name, ehA);
		ehB = eh_generate_hist(dsB, spec, CHECKR);
		euler_print_hist(dsB->metadata.name, ehB);
	}

    euler_histogram *ehIntermed = NULL;
    ehIntermed = eh_generate_intermed_real(dsA,dsB,ehA,ehB);

    //escolhe dataset A ou para Intermediario
    euler_print_hist("intermediario",ehIntermed);

	OGR_DS_Destroy(ogr_ds1);
	OGR_DS_Destroy(ogr_ds2);
	finishGEOS();

	return 0;
}


dataset *read_geos(char *shpfile, OGRDataSourceH ogr_ds) {
	dataset *results = dataset_create_mem(basename(shpfile), 1);

	ogr_ds = OGROpen(shpfile, false, NULL);
	if (ogr_ds == NULL) {
		fprintf(stderr, "Invalid shapefile.\n");
		return results;
	}

	OGRLayerH layer = OGR_DS_GetLayer(ogr_ds, 0);
	OGR_L_ResetReading(layer);
	//int layer_count = OGR_L_GetFeatureCount(layer, FALSE);

	results->temp_ogr_layer = layer;

	//clock_t cs = clock();

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

	//clock_t cf = clock();
	//double runtime_diff_us = (cf-cs) * 1000. / CLOCKS_PER_SEC;

	//printf("%4s|%4s|%'10d|%10s|%10s|%10s|%10s|%10s|%10s|%10.1f| %s\n", "N", "N", results->metadata.count, "N", "N", "N", "N", "N", "N", runtime_diff_us, shpfile);

	return results;
}



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

char *dataset_a;
char *dataset_b;

dataset *read_geos(char *shpfile);

OGRDataSourceH ogr_ds;


int main (int argc, char* argv[]){
	OGRRegisterAll();
	initGEOS(geos_messages, geos_messages);

    dataset_a = argv[1];
    dataset_b = argv[2];
    dataset* ds_a = read_geos(dataset_a);
    dataset* ds_b = read_geos(dataset_b);

	enum HistogramHashMethod hm = HHASH_AREAFRAC;
	enum HistogramSplitMethod sm = HSPLIT_AVG_STD;
    int xqtd = argv[3];
    int yqtd = argv[4];

	HistogramGenerateSpec spec;
	spec.hm = hm;
	spec.sm = sm;
	spec.xqtd = xqtd;
	spec.yqtd = yqtd;

	int split_qtd = 1;

	// chamar a função que cria o histograma
	histogram_generate(ds_a, spec, CHECKR, split_qtd);
	histogram_print_geojson(ds_a);
	histogram_print(ds_a, CARDIN);
	histogram_generate(ds_b, spec, CHECKR, split_qtd);
	histogram_print_geojson(ds_b);
	histogram_print(ds_b, CARDIN);

	euler_histogram *ehr = NULL;
	euler_histogram *ehs = NULL;
    ehs = eh_generate_hist(ds_a, spec, CHECKR);
	euler_print_hist(ds_a, ehs);
    ehr = eh_generate_hist(ds_b, spec, CHECKR);
	euler_print_hist(ds_b, ehr);


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


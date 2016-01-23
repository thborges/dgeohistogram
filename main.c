
#include <ogr_api.h>
#include <geos_c.h>
#include <cpl_conv.h>
#include <time.h>
#include <arpa/inet.h>
#include <locale.h>
#include "dataset.h"
#include "glibwrap.h"
#include "ogrext.h"

dataset *read_geos(char *shpfile);

int main(int argc, char* argv[]) {
	//srand(time(NULL));
	OGRRegisterAll();
	initGEOS(geos_messages, geos_messages);

	if (argc < 2) {
		printf("Use: %s file.shp\n", argv[0]);
		return 1;
	}

	dataset *ds = read_geos(argv[1]);

	// chamar a função que cria o histograma
	histogram_generate(ds, 0);
	histogram_print_geojson(ds);

	finishGEOS();

	return 0;
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
	int layer_count = OGR_L_GetFeatureCount(layer, FALSE);

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
				dataset_leaf *leaves = dataset_add(results);
				leaves[0].geo = ggeo;
				GEOSGetEnvelope(ggeo, &leaves[0].mbr);

			}
		}
		OGR_F_Destroy(feature);

		read++;
		print_progress_gauge(read, layer_count);
	}
	OGR_DS_Destroy(ds);

	clock_t cf = clock();
	double runtime_diff_us = (cf-cs) * 1000. / CLOCKS_PER_SEC;

	printf("%4s|%4s|%'10d|%10s|%10s|%10s|%10s|%10s|%10s|%10.1f| %s\n", "N", "N", results->metadata.count, "N", "N", "N", "N", "N", "N", runtime_diff_us, shpfile);

	return results;
}


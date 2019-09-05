
#include <cpl_conv.h>
#include <time.h>
#include <arpa/inet.h>
#include <locale.h>
#include <libgen.h>
#include <gdal.h>

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

/*Join a dataset r and an rtree s, using INL*/
dataset*		rtree_join_inl	(dataset *r, const rtree_root *s, rtree_join_stat *stats, enum JoinPredicateCheck pcheck);

void write_to_shp(dataset *dh, char *name);

OGRDataSourceH ogr_ds1;
OGRDataSourceH ogr_ds2;

int main(int argc, char* argv[]) {

	//srand(time(NULL));

	OGRRegisterAll();
	initGEOS(geos_messages, geos_messages);

	if (argc < 4) {
		printf("Use: %s file1.shp  file2.shp result.shp\n", argv[0]);
		return 1;
	}

	int argatu = 1;

	//dataset a
	dataset_name = argv[argatu++];
	printf("Carregando arquivo: %s\n", dataset_name);
	dataset *dsA = read_geos(dataset_name, ogr_ds1);

	//dataset b
	dataset_name = argv[argatu++];
	printf("Carregando arquivo: %s\n", dataset_name);
	dataset *dsB = read_geos(dataset_name, ogr_ds2);

	// cria uma r* para dataset A (para o INL)
	rtree_root *rtree_dsB = NULL;
	rtree_dsB = rtree_new_rstar(30, 10);

	unsigned read = 0;
	dataset_iter_seg iter;
	dataset_foreach(iter, dsB) {
		read++;
		GEOSGeometryH geo = dataset_get_leaf_geo(dsB, iter.item);
		rtree_append(rtree_dsB, geo);
		print_progress_gauge(read, dsB->metadata.count);
	}

	rtree_join_stat stats;
	memset(&stats, 0, sizeof stats);

	dataset *result = rtree_join_inl(dsA, rtree_dsB, &stats, CHECKR);

	/*FILE *f = fopen(argv[3], "w");
	  dataset_write_csv(result, f);
	  fclose(f);*/
	fprintf(stderr, "\nWriting shapefile to %s\n", argv[3]);
	write_to_shp(result, argv[3]);

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
				Envelope mbr = OGRGetEnvelope(geometry);
				dataset_fill_leaf_id(leaf, 0, -1, &mbr);
				//leaf[0].points = OGRGetNumPoints(geometry);
				leaf[0].geo = ggeo;
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

void write_to_shp(dataset *dh, char *name) {

	const char *pszDriverName = "ESRI Shapefile";
	GDALDriverH hDriver;
	GDALDatasetH hDS;
	OGRLayerH hLayer;
	OGRFieldDefnH hFieldDefn;
	GDALAllRegister();
	hDriver = GDALGetDriverByName( pszDriverName );
	if( hDriver == NULL )
	{
		printf( "%s driver not available.\n", pszDriverName );
		exit( 1 );
	}
	hDS = GDALCreate( hDriver, name, 0, 0, 0, GDT_Unknown, NULL );
	if( hDS == NULL )
	{
		printf( "Creation of output file failed.\n" );
		exit( 1 );
	}
	hLayer = GDALDatasetCreateLayer( hDS, name, NULL, wkbLineString, NULL );
	if( hLayer == NULL )
	{
		printf( "Layer creation failed.\n" );
		exit( 1 );
	}
	hFieldDefn = OGR_Fld_Create( "ID", OFTInteger);
	if( OGR_L_CreateField( hLayer, hFieldDefn, TRUE ) != OGRERR_NONE )
	{
		printf( "Creating ID field failed.\n" );
		exit( 1 );
	}
	OGR_Fld_Destroy(hFieldDefn);

	int id = 0;
	dataset_iter feature;
	dataset_foreach(feature, dh) {
		OGRFeatureH hFeature;
		OGRGeometryH hPt;
		hFeature = OGR_F_Create( OGR_L_GetLayerDefn( hLayer ) );
		OGR_F_SetFieldInteger( hFeature, OGR_F_GetFieldIndex(hFeature, "ID"), id);

                dataset_leaf *lgeo = get_join_pair_leaf(feature.item, 0);
		char *wkt = GEOSGeomToWKT(lgeo->geo);
		char *wkt2 = wkt;
		OGR_G_CreateFromWkt(&wkt2, NULL, &hPt);
		free(wkt);
		OGR_F_SetGeometry( hFeature, hPt );
		OGR_G_DestroyGeometry(hPt);
		if( OGR_L_CreateFeature( hLayer, hFeature ) != OGRERR_NONE )
		{
			printf( "Failed to create feature in shapefile.\n" );
			exit( 1 );
		}
		OGR_F_Destroy( hFeature );
		id++;
		print_progress_gauge(id, dh->metadata.count);
	}

	GDALClose( hDS );

}


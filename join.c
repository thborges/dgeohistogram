

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

#define GET_VERT_EDGE_EHR(x, y) ((x == ehr->xqtd) ? (x * (2*ehr->yqtd+1) + y) : (x * (2*ehr->yqtd+1) + 2*y + 1))
#define GET_HORZ_EDGE_EHR(x, y) (x * (2*ehs->yqtd+1) + 2*y)
#define GET_VERT_EDGE_EHS(x, y) ((x == ehs->xqtd) ? (x * (2*ehs->yqtd+1) + y) : (x * (2*ehs->yqtd+1) + 2*y + 1))
#define GET_HORZ_EDGE_EHS(x, y) (x * (2*ehs->yqtd+1) + 2*y)

char *dataset_a;
char *dataset_b;

rtree_window_stat stats1;
rtree_window_stat stats2;
dataset *read_geos(char *shpfile);
int real_spatial_join_cardin(rtree_root* r_a, rtree_root* r_b, euler_histogram* ehs, euler_histogram* ehr);

//int histogram_join_cardinality(dataset *dr, dataset *ds);
OGRDataSourceH ogr_ds;


int main (int argc, char* argv[]){
    OGRRegisterAll();
    initGEOS(geos_messages, geos_messages);

    enum HistogramHashMethod hm = HHASH_AREAFRAC;
    enum HistogramSplitMethod sm = HSPLIT_FIX;

    dataset_a = argv[1];
    dataset_b = argv[2];
    printf("%s, %s\n", dataset_a, dataset_b);
    dataset* ds_a = read_geos(dataset_a);
    dataset* ds_b = read_geos(dataset_b);

    int xqtd = atoi(argv[3]);
    int yqtd = atoi(argv[4]);

    printf("%d, %d\n", xqtd, yqtd);
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

    ehr = eh_generate_hist(ds_a, spec, CHECKR);
    euler_print_hist(ds_a, ehr);

    ehs = eh_generate_hist(ds_b, spec, CHECKR);
    euler_print_hist(ds_b, ehs);

    // cria e preenche as r-trees dos dos datasets ds_a e ds_b
    rtree_root *rtree_a = NULL;
    rtree_a = rtree_new_rstar(30, 10);

    rtree_root *rtree_b = NULL;
    rtree_b = rtree_new_rstar(30, 10);


    dataset_iter_seg iter_a;
    unsigned read = 0;
    dataset_foreach(iter_a, ds_a) {
        read++;
        GEOSGeometryH geo = dataset_get_leaf_geo(ds_a, iter_a.item);
        rtree_append(rtree_a, geo);
        print_progress_gauge(read, ds_a->metadata.count);

    }

    dataset_iter_seg iter_b;
    dataset_foreach(iter_b, ds_b) {
        GEOSGeometryH geo = dataset_get_leaf_geo(ds_b, iter_b.item);
        rtree_append(rtree_b, geo);
    }

    double accuracy = 0.0;
    double sum_ei = 0.0;
    double sum_ri = 0.0;
    double mean = 0.0;
    double M2 = 0.0;
    double sum_error = 0.0;
    int n = 0;


    //double a = real_spatial_join_cardin(rtree_a, rtree_b, ehs, ehr);
    int b = euler_join_cardinality(ds_a, ds_b, ehr, ehs);
        
    printf("%d\n", b);

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

// r_a é a rtree do ehs, r_b é a rtree do ehr
int real_spatial_join_cardin(rtree_root* r_a, rtree_root* r_b, euler_histogram* ehs, euler_histogram* ehr){
    printf("size r_a = %d\n", rtree_height(r_a));
    printf("size r_b = %d\n", rtree_height(r_b));

    if(!ENVELOPE_INTERSECTS(ehr->mbr, ehs->mbr))
        return 0;

    double result = 0;
	int xds_atu = 0;

    double xini = MAX(ehr->xtics[0], ehs->xtics[0]);
    double yini = MAX(ehr->ytics[0], ehs->ytics[0]);
    double xfim = MIN(ehr->mbr.MaxX, ehs->mbr.MaxX);
    double yfim = MIN(ehr->mbr.MaxY, ehs->mbr.MaxY);
    //printf("xini = %e\n\
            yini = %e\n\
            xfim = %e\n\
            yfim = %e\n", xini, yini, xfim, yfim);


    int xdr_start = 0;
    while (xdr_start < ehr->xqtd && ehr->xtics[xdr_start+1] < xini)
        xdr_start++;
    int xdr_end = xdr_start+1;
    while (xdr_end < ehr->xqtd && ehr->xtics[xdr_end] <= xfim)
        xdr_end++;
    if (xdr_start == ehr->xqtd)
        return 0;

    // skip non-intersect area on y
    int ydr_start = 0;
    while (ydr_start < ehr->yqtd && ehr->ytics[ydr_start+1] < yini)
        ydr_start++;
    int ydr_end = ydr_start+1;
    while (ydr_end < ehr->yqtd && ehr->ytics[ydr_end] <= yfim)
        ydr_end++;
    if (ydr_start == ehr->yqtd)
        return 0;
    for(int xr = xdr_start; xr < xdr_end; xr++){

        while(xds_atu < ehs->xqtd && ehs->xtics[xds_atu+1] < ehr->xtics[xr]) // skip when end of s < start of r
            xds_atu++;
        int xds_end = xds_atu+1;
        while(xds_end < ehs->xqtd && ehs->xtics[xds_end] <= ehr->xtics[xr+1]) // increment when end of s < start of r
            xds_end++;

        int yds_atu = 0;

        Envelope er, es;

        er.MinX = ehr->xtics[xr];
        er.MaxX = ehr->xtics[xr+1];
        //printf("xr = %d\n", xr);

        for(int yr = ydr_start; yr < ydr_end; yr++){

            while(yds_atu < ehs->yqtd && ehs->ytics[yds_atu+1] < ehr->ytics[yr]) // skip when end of s < start of r
                yds_atu++;
            int yds_end = yds_atu+1;
            while(yds_end < ehs->yqtd && ehs->ytics[yds_end] <= ehr->ytics[yr+1]) // increment when end of s < start of r
                yds_end++;

            er.MinY = ehr->xtics[yr];
            er.MaxY = ehr->xtics[yr+1];

            for(int xs = xds_atu; xs < xds_end; xs++){
                es.MinX = ehs->xtics[xs];
                es.MaxX = ehs->xtics[xs+1];

                for(int ys = yds_atu; ys < yds_end; ys++){
                    es.MinY = ehs->xtics[ys];
                    es.MaxY = ehs->xtics[ys+1];
                    if(ENVELOPE_INTERSECTS(er, es)){

                        char wkt_es[512];
                        sprintf(wkt_es, "POLYGON((%f %f, %f %f, %f %f, %f %f, %f %f))",
                                es.MinX, es.MinY,
                                es.MaxX, es.MinY,
                                es.MaxX, es.MaxY,
                                es.MinX, es.MaxY,
                                es.MinX, es.MinY);

                        char wkt_er[512];
                        sprintf(wkt_er, "POLYGON((%f %f, %f %f, %f %f, %f %f, %f %f))",
                                er.MinX, er.MinY,
                                er.MaxX, er.MinY,
                                er.MaxX, er.MaxY,
                                er.MinX, er.MaxY,
                                er.MinX, er.MinY);
                     //   printf("%s\n", wkt_es);

                        GEOSGeometryH geo_es = GEOSGeomFromWKT(wkt_es);
                        GEOSGeometryH geo_er = GEOSGeomFromWKT(wkt_er);

                        memset(&stats1, 0, sizeof(rtree_window_stat));
                        memset(&stats2, 0, sizeof(rtree_window_stat));

                        GList *results_a = rtree_window_search(r_a, geo_er, &stats1);
                        GList *results_b = rtree_window_search(r_b, geo_es, &stats2);

                        if(g_list_length(results_a) > 0 || g_list_length(results_b) > 0  ){
                            printf("len(results_a) = %d\n", g_list_length(results_a));
                            printf("len(results_b) = %d\n", g_list_length(results_b));
                        }

                        GList *a;

                        g_list_foreach(a, results_a){
                            rtree_leaf* la = (rtree_leaf*)a->data;

                            GList* b;
                            //printf("TO AQUI2 \n");

                            g_list_foreach(b, results_b){
                                rtree_leaf* lb = (rtree_leaf*)b->data;

                                if(GEOSIntersects(la->geo, lb->geo))
                                    result++;
                                else
                                    printf("%f\n", result);
                            }
                        }
                    }
                }
            }
        }
    }


    return (result);


}

int real_join_cardinality(dataset *dr, dataset *ds) {
    rtree_root *r_a = NULL;
    r_a = rtree_new_rstar(30, 10);

    rtree_root *r_b = NULL;
    r_b = rtree_new_rstar(30, 10);


    dataset_iter_seg iter_a;
    unsigned read = 0;
    dataset_foreach(iter_a, ds) {
        read++;
        GEOSGeometryH geo = dataset_get_leaf_geo(dr, iter_a.item);
        rtree_append(r_a, geo);
        print_progress_gauge(read, dr->metadata.count);

    }

    dataset_iter_seg iter_b;
    dataset_foreach(iter_b, dr) {
        GEOSGeometryH geo = dataset_get_leaf_geo(dr, iter_b.item);
        rtree_append(r_b, geo);
    }
	dataset_histogram *hr = &dr->metadata.hist;
	dataset_histogram *hs = &ds->metadata.hist;
	double xini = MAX(hr->xtics[0], hs->xtics[0]);
	double yini = MAX(hr->ytics[0], hs->ytics[0]);
	double xend = MIN(dr->metadata.hist.mbr.MaxX, ds->metadata.hist.mbr.MaxX);
	double yend = MIN(dr->metadata.hist.mbr.MaxY, ds->metadata.hist.mbr.MaxY);

	// skip non-intersect area on x
	int xdr_start = 0;
	while (xdr_start < hr->xqtd && hr->xtics[xdr_start+1] < xini)
		xdr_start++;
	int xdr_end = xdr_start+1;
	while (xdr_end < hr->xqtd && hr->xtics[xdr_end] <= xend)
		xdr_end++;
	if (xdr_start == hr->xqtd)
		return 0;

	// skip non-intersect area on y
	int ydr_start = 0;
	while (ydr_start < hr->yqtd && hr->ytics[ydr_start+1] < yini)
		ydr_start++;
	int ydr_end = ydr_start+1;
	while (ydr_end < hr->yqtd && hr->ytics[ydr_end] <= yend)
		ydr_end++;
	if (ydr_start == hr->yqtd)
		return 0;

	int xds_atu = 0;
	float result = 0;
	for(int xr = xdr_start; xr < xdr_end; xr++) {

		while(xds_atu < hs->xqtd && hs->xtics[xds_atu+1] < hr->xtics[xr]) // skip when end of s < start of r
			xds_atu++;
		int xds_end = xds_atu+1;
		while(xds_end < hs->xqtd && hs->xtics[xds_end] <= hr->xtics[xr+1]) // increment when end of s < start of r
			xds_end++;

		int yds_atu = 0;

		Envelope er;
		Envelope es;
		er.MinX = hr->xtics[xr];
		er.MaxX = hr->xtics[xr+1];

		for(int yr = ydr_start; yr < ydr_end; yr++) {

			while(yds_atu < hs->yqtd && hs->ytics[yds_atu+1] < hr->ytics[yr]) // skip when end of s < start of r
				yds_atu++;
			int yds_end = yds_atu+1;
			while(yds_end < hs->yqtd && hs->ytics[yds_end] <= hr->ytics[yr+1]) // increment when end of s < start of r
				yds_end++;

			er.MinY = hr->ytics[yr];
			er.MaxY = hr->ytics[yr+1];
			double erarea = ENVELOPE_AREA(er);

			for(int xs = xds_atu; xs < xds_end; xs++) {
				es.MinX = hs->xtics[xs];
				es.MaxX = hs->xtics[xs+1];

				for(int ys = yds_atu; ys < yds_end; ys++) {

					es.MinY = hs->ytics[ys];
					es.MaxY = hs->ytics[ys+1];


                    if(ENVELOPE_INTERSECTS(er, es)){
                        char wkt_es[512];
                        sprintf(wkt_es, "POLYGON((%f %f, %f %f, %f %f, %f %f, %f %f))",
                                es.MinX, es.MinY,
                                es.MaxX, es.MinY,
                                es.MaxX, es.MaxY,
                                es.MinX, es.MaxY,
                                es.MinX, es.MinY);

                        char wkt_er[512];
                        sprintf(wkt_er, "POLYGON((%f %f, %f %f, %f %f, %f %f, %f %f))",
                                er.MinX, er.MinY,
                                er.MaxX, er.MinY,
                                er.MaxX, er.MaxY,
                                er.MinX, er.MaxY,
                                er.MinX, er.MinY);
                     //   printf("%s\n", wkt_es);

                        GEOSGeometryH geo_es = GEOSGeomFromWKT(wkt_es);
                        GEOSGeometryH geo_er = GEOSGeomFromWKT(wkt_er);

                        memset(&stats1, 0, sizeof(rtree_window_stat));
                        memset(&stats2, 0, sizeof(rtree_window_stat));

                        GList *results_a = rtree_window_search(r_a, geo_er, &stats1);
                        GList *results_b = rtree_window_search(r_b, geo_es, &stats2);

                        if(g_list_length(results_a) > 0 || g_list_length(results_b) > 0  ){
                            printf("len(results_a) = %d\n", g_list_length(results_a));
                            printf("len(results_b) = %d\n", g_list_length(results_b));
                        }

                        GList *a;

                        g_list_foreach(a, results_a){
                            rtree_leaf* la = (rtree_leaf*)a->data;

                            GList* b;
                            //printf("TO AQUI2 \n");

                            g_list_foreach(b, results_b){
                                rtree_leaf* lb = (rtree_leaf*)b->data;

                                if(GEOSIntersects(la->geo, lb->geo))
                                    result++;
                                else
                                    printf("%f\n", result);
                            }
                        }
                    }
				}
			}
		}
	}

	return (int)result;
}




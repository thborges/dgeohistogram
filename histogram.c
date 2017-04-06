/*
 * histogram.c
 *
 *  Created on: 09/08/2014
 *      Author: thborges
 */

#include <float.h>
#include "histogram.h"

const char *HistogramHashMethodName[4]  = {
	"mbrc",
	"cent",
	"areaf",
	"areafs",
};

void histogram_alloc(dataset_histogram *dh, int xqtd, int yqtd) {
	assert(xqtd > 0 && yqtd > 0 && "X and Y must be greater than zero.");
	dh->xqtd = xqtd;
	dh->yqtd = yqtd;
	dh->xtics = g_new(double, dh->xqtd+1);
	dh->ytics = g_new(double, dh->yqtd+1);
	dh->hcells = g_new0(histogram_cell, dh->xqtd*dh->yqtd);
}

dataset_histogram histogram_join_io(dataset *dr, dataset *ds, enum JoinPredicateCheck next_pcheck, multiway_histogram_estimate *estimate);

inline __attribute__((always_inline))
void fill_hist_cell_centroid(dataset_leaf *l, dataset *ds, dataset_histogram *dh) {

	GEOSGeometryH geo = dataset_get_leaf_geo(ds, l);
	GEOSGeometryH centroid = GEOSGetCentroid(geo);

	double x, y;
	GEOSGeomGetX(centroid, &x);
	GEOSGeomGetY(centroid, &y);

	int xp = (x - ds->metadata.hist.mbr.MinX) / dh->xsize;
	int yp = (y - ds->metadata.hist.mbr.MinY) / dh->ysize;
	histogram_cell *cell = &ds->metadata.hist.hcells[xp*ds->metadata.hist.yqtd +yp];
	cell->cardin++;
	cell->points += l->points;

	GEOSGeom_destroy(centroid);
	if (l->gid != -1) // free because of the call to dataset_get_leaf_geo
		GEOSGeom_destroy(geo);
}

inline __attribute__((always_inline))
void fill_hist_cell_mbr_center(dataset_leaf *l, dataset *ds, dataset_histogram *dh) {
	// using mbr center
	double x = l->mbr.MinX + (l->mbr.MaxX - l->mbr.MinX) / 2.0;
	double y = l->mbr.MinY + (l->mbr.MaxY - l->mbr.MinY) / 2.0;

	int xp = (x - ds->metadata.hist.mbr.MinX) / dh->xsize;
	int yp = (y - ds->metadata.hist.mbr.MinY) / dh->ysize;
	histogram_cell *cell = &ds->metadata.hist.hcells[xp*ds->metadata.hist.yqtd +yp];
	cell->cardin++;
	cell->points += l->points;

	//TODO: calculate avg online
	double delta_x = MIN(l->mbr.MaxX - l->mbr.MinX, dh->xsize);
	double delta_y = MIN(l->mbr.MaxY - l->mbr.MinY, dh->ysize);
	cell->avgwidth += delta_x;
	cell->avgheight += delta_y;
}

void get_xini_xfim(dataset_histogram *dh, Envelope query, 
	int *xini, int *xfim, int *yini, int *yfim) {

	// prevent values of x and y out of histogram bounds
	query = EnvelopeIntersection2(query, dh->mbr);

	*xini = (query.MinX - dh->mbr.MinX) / dh->xsize;
	*xfim = (query.MaxX - dh->mbr.MinX) / dh->xsize;
	*yini = (query.MinY - dh->mbr.MinY) / dh->ysize;
	*yfim = (query.MaxY - dh->mbr.MinY) / dh->ysize;
   
	if (*xfim == dh->xqtd) (*xfim)--;
	if (*yfim == dh->yqtd) (*yfim)--;

	assert(*xini >= 0 && *xfim < dh->xqtd && "x is out of histogram bounds.");
	assert(*yini >= 0 && *yfim < dh->yqtd && "y is out of histogram bounds.");
}

//inline __attribute__((always_inline))
void hash_envelope_area_fraction(dataset_histogram *dh, Envelope ev, double objarea, double points) {

	int xini, xfim, yini, yfim;
	get_xini_xfim(dh, ev, &xini, &xfim, &yini, &yfim);	

	//printf("%d %d %d %d\n", xini, xfim, yini, yfim);

	for(int x = xini; x <= xfim; x++) {
		Envelope rs;
		rs.MinX = dh->xtics[x];
		rs.MaxX = dh->xtics[x+1];

		for(int y = yini; y <= yfim; y++) {
			rs.MinY = dh->ytics[y];
			rs.MaxY = dh->ytics[y+1];

			Envelope inters = EnvelopeIntersection(ev, rs);
			double intarea = ENVELOPE_AREA(inters);
			double rsarea = ENVELOPE_AREA(rs);

			// for point objects, objarea == 0.0
			double fraction = (objarea == 0.0) ? 1.0 : intarea / objarea;

			histogram_cell *cell = &dh->hcells[x*dh->yqtd +y];
			//cell->cardin += 1.0;//fraction;
			cell->cardin += fraction;
			cell->points += points; //object is replicated

			//TODO: calculate avg online
			cell->avgwidth += (inters.MaxX - inters.MinX);
			cell->avgheight += (inters.MaxY - inters.MinY);
		}
	}
}

inline __attribute__((always_inline))
void fill_hist_cell_area_fraction(dataset_leaf *l, dataset *ds, dataset_histogram *dh) {
	// proportional to cover area

	//printf("GID %lld: ", l->gid);

	double objarea = ENVELOPE_AREA(l->mbr);
	hash_envelope_area_fraction(dh, l->mbr, objarea, l->points);

	//printf("%10.5f <= %10.5f <= %10.5f\n", ds->metadata.hist.xtics[xp], x, ds->metadata.hist.xtics[xp+1]);
	//assert(x >= ds->metadata.hist.xtics[xp] && x <= ds->metadata.hist.xtics[xp+1]);
	//assert(y >= ds->metadata.hist.ytics[yp] && y <= ds->metadata.hist.ytics[yp+1]);
}

void envelope_update(Envelope *e, double X, double Y) {
	e->MinX = MIN(e->MinX, X);
	e->MinY = MIN(e->MinY, Y);
	e->MaxX = MAX(e->MaxX, X);
	e->MaxY = MAX(e->MaxY, Y);
}

int fill_hist_cell_area_fraction_with_split(dataset_leaf *l, dataset *ds, dataset_histogram *dh) {
	// proportional to cover area

	int xini, xfim, yini, yfim;
	get_xini_xfim(dh, l->mbr, &xini, &xfim, &yini, &yfim);	

	double objarea = ENVELOPE_AREA(l->mbr);

	int splitted = 0;

	// is a candidate for split?
	int xspan = xfim - xini;
	int yspan = yfim - yini;
	if (xspan >= 2 || yspan >= 2) { // more than two cells?
		splitted = 1;

		GEOSGeometryH geo = dataset_get_leaf_geo(ds, l);

		Envelope split1 = l->mbr;
		Envelope split2 = l->mbr;
		Envelope split3 = l->mbr;
		double filled_space = DBL_MAX;
		
		const GEOSGeometry *linearRing;
		const GEOSCoordSequence *coordSeq;
		int numGeom = GEOSGetNumGeometries(geo);
		for(int n = 0; n < numGeom; n++) {
			const GEOSGeometry *ngeo = GEOSGetGeometryN(geo, n);
			if (GEOSGeomTypeId(ngeo) == GEOS_POLYGON)
				linearRing = GEOSGetExteriorRing(ngeo);
			else
				linearRing = ngeo;

			int numPoints = GEOSGeomGetNumPoints(linearRing);
	        coordSeq = GEOSGeom_getCoordSeq(linearRing);

			double xCoord, yCoord;
       		GEOSCoordSeq_getX(coordSeq, 0, &xCoord);
	        GEOSCoordSeq_getY(coordSeq, 0, &yCoord);

			Envelope e1 = {DBL_MAX, DBL_MAX, -DBL_MAX, -DBL_MAX};
			envelope_update(&e1, xCoord, yCoord);

			for (int p=1; p < numPoints-1; p++) {
				double xCoord, yCoord;
       			GEOSCoordSeq_getX(coordSeq, p, &xCoord);
	        	GEOSCoordSeq_getY(coordSeq, p, &yCoord);
				envelope_update(&e1, xCoord, yCoord);
	
				Envelope e2 = {DBL_MAX, DBL_MAX, -DBL_MAX, -DBL_MAX};
				envelope_update(&e2, xCoord, yCoord);
				for(int p2=p+1; p2 < numPoints; p2++) {
	       			GEOSCoordSeq_getX(coordSeq, p2, &xCoord);
		        	GEOSCoordSeq_getY(coordSeq, p2, &yCoord);
					envelope_update(&e2, xCoord, yCoord);

					// terceiro ponto
					Envelope e3 = {DBL_MAX, DBL_MAX, -DBL_MAX, -DBL_MAX};
					envelope_update(&e3, xCoord, yCoord);

					for(int p3=p2+1; p3 < numPoints; p3++) {
	       				GEOSCoordSeq_getX(coordSeq, p3, &xCoord);
		        		GEOSCoordSeq_getY(coordSeq, p3, &yCoord);
						envelope_update(&e3, xCoord, yCoord);	
					}

					double fs = (ENVELOPE_AREA(e1) + ENVELOPE_AREA(e2) + ENVELOPE_AREA(e3)) / objarea;
					if (fs < filled_space) {
						split1 = e1;
						split2 = e2;
						split3 = e3;
						filled_space = fs;
					}

				}

			}
		}
		
		objarea = ENVELOPE_AREA(split1) + ENVELOPE_AREA(split2) + ENVELOPE_AREA(split3);
		hash_envelope_area_fraction(dh, split1, objarea, l->points);
		hash_envelope_area_fraction(dh, split2, objarea, l->points);
		hash_envelope_area_fraction(dh, split3, objarea, l->points);

		if (l->gid != -1) // free due to the call to dataset_get_leaf_geo
			GEOSGeom_destroy(geo);

/*		print_geojson_header();
		print_geojson_mbr(l->mbr, "orig");
		print_geojson_mbr(split1, "e1");
		print_geojson_mbr(split2, "e2");
		print_geojson_mbr(split3, "e3");
		print_geojson_footer();*/
	}
	else {
		hash_envelope_area_fraction(dh, l->mbr, objarea, l->points);
	}

	return splitted;
}

/* void fill_hist_cell_area_fraction_with_split(dataset_leaf *l, dataset *ds, dataset_histogram *dh) {
	// proportional to cover area

	int xini = (l->mbr.MinX - dh->mbr.MinX) / dh->xsize;
	int xfim = (l->mbr.MaxX - dh->mbr.MinX) / dh->xsize;
	int yini = (l->mbr.MinY - dh->mbr.MinY) / dh->ysize;
	int yfim = (l->mbr.MaxY - dh->mbr.MinY) / dh->ysize;
	double objarea = ENVELOPE_AREA(l->mbr);

	// is a candidate for split?
	int xspan = xfim - xini;
	int yspan = yfim - yini;
	if (xspan >= 2 || yspan >= 2) { // more than two cells?
		GEOSGeometryH geo = dataset_get_leaf_geo(ds, l);

		Envelope split1 = l->mbr;
		Envelope split2 = l->mbr;
		double filled_space = DBL_MAX;
		
		const GEOSGeometry *linearRing;
		const GEOSCoordSequence *coordSeq;
		int numGeom = GEOSGetNumGeometries(geo);
		for(int n = 0; n < numGeom; n++) {
			const GEOSGeometry *ngeo = GEOSGetGeometryN(geo, n);
			if (GEOSGeomTypeId(ngeo) == GEOS_POLYGON)
				linearRing = GEOSGetExteriorRing(ngeo);
			else
				linearRing = ngeo;

			int numPoints = GEOSGeomGetNumPoints(linearRing);
	        coordSeq = GEOSGeom_getCoordSeq(linearRing);

			double xCoord, yCoord;
       		GEOSCoordSeq_getX(coordSeq, 0, &xCoord);
	        GEOSCoordSeq_getY(coordSeq, 0, &yCoord);

			Envelope e1 = {DBL_MAX, DBL_MAX, -DBL_MAX, -DBL_MAX};
			envelope_update(&e1, xCoord, yCoord);

			for (int p=1; p < numPoints-1; p++) {
				double xCoord, yCoord;
       			GEOSCoordSeq_getX(coordSeq, p, &xCoord);
	        	GEOSCoordSeq_getY(coordSeq, p, &yCoord);
				envelope_update(&e1, xCoord, yCoord);
	
				Envelope e2 = {DBL_MAX, DBL_MAX, -DBL_MAX, -DBL_MAX};
				envelope_update(&e2, xCoord, yCoord);
				for(int p2=p+1; p2 < numPoints; p2++) {
	       			GEOSCoordSeq_getX(coordSeq, p2, &xCoord);
		        	GEOSCoordSeq_getY(coordSeq, p2, &yCoord);
					envelope_update(&e2, xCoord, yCoord);
				}

				double fs = (ENVELOPE_AREA(e1) + ENVELOPE_AREA(e2)) / objarea;
				if (fs < filled_space) {
					split1 = e1;
					split2 = e2;
					filled_space = fs;
				}
			}
		}
		
		objarea = ENVELOPE_AREA(split1) + ENVELOPE_AREA(split2);
		hash_envelope_area_fraction(dh, split1, objarea, l->points);
		hash_envelope_area_fraction(dh, split2, objarea, l->points);

		if (l->gid != -1) // free due to the call to dataset_get_leaf_geo
			GEOSGeom_destroy(geo);

		print_geojson_header();
		print_geojson_mbr(l->mbr, "orig");
		print_geojson_mbr(split1, "e1");
		print_geojson_mbr(split2, "e2");
		print_geojson_footer();
	}
	else {
		hash_envelope_area_fraction(dh, l->mbr, objarea, l->points);
	}
	return splitted;

} */


//inline __attribute__((always_inline))
void fill_hist_cell_area_fraction_with_split_old(dataset_leaf *l, dataset *ds, dataset_histogram *dh) {
	// proportional to cover area

	int xini = (l->mbr.MinX - dh->mbr.MinX) / dh->xsize;
	int xfim = (l->mbr.MaxX - dh->mbr.MinX) / dh->xsize;
	int yini = (l->mbr.MinY - dh->mbr.MinY) / dh->ysize;
	int yfim = (l->mbr.MaxY - dh->mbr.MinY) / dh->ysize;
	double objarea = ENVELOPE_AREA(l->mbr);

	// is a candidate for split?
	int xspan = xfim - xini;
	int yspan = yfim - yini;
	if (xspan >= 2 || yspan >= 2) { // more than two cells?
		int xsplit = xspan > yspan;
		GEOSGeometryH geo = dataset_get_leaf_geo(ds, l);
		
		double split_at;
		if (xsplit)
			split_at = (l->mbr.MaxX + l->mbr.MinX) / 2.0;
		else
			split_at = (l->mbr.MaxY + l->mbr.MinY) / 2.0;

		Envelope e1 = {DBL_MAX, DBL_MAX, -DBL_MAX, -DBL_MAX};
		Envelope e2 = {DBL_MAX, DBL_MAX, -DBL_MAX, -DBL_MAX};

		const GEOSGeometry *linearRing;
		const GEOSCoordSequence *coordSeq;
		int numGeom = GEOSGetNumGeometries(geo);
		for(int n = 0; n < numGeom; n++) {
			const GEOSGeometry *ngeo = GEOSGetGeometryN(geo, n);
			if (GEOSGeomTypeId(ngeo) == GEOS_POLYGON)
				linearRing = GEOSGetExteriorRing(ngeo);
			else
				linearRing = ngeo;

			char splitted = '0';
			double last_x, last_y;
			int numPoints = GEOSGeomGetNumPoints(linearRing);
	        coordSeq = GEOSGeom_getCoordSeq(linearRing);
			for (int p=0; p < numPoints; p++) {
				double xCoord, yCoord;
        		GEOSCoordSeq_getX(coordSeq, p, &xCoord);
		        GEOSCoordSeq_getY(coordSeq, p, &yCoord);
				if ((xsplit && xCoord <= split_at) || (!xsplit && yCoord <= split_at)) {
					if (splitted == '0') 
						splitted = 'x';
					else if (splitted == 'y') {
						splitted = '1';
						envelope_update(&e1, last_x, last_y);
						split_at = xsplit ? last_x: last_y;
					}		
					envelope_update(&e1, xCoord, yCoord);
				}
				else {
					if (splitted == '0') 
						splitted = 'y';
					else if (splitted == 'x') {
						splitted = '1';
						envelope_update(&e2, last_x, last_y);
						split_at = xsplit ? last_x: last_y;
					}		
					envelope_update(&e2, xCoord, yCoord);
				}
				last_x = xCoord;
				last_y = yCoord;
			}
			
		}
		if (xsplit)
			e1.MaxX = e2.MinX = split_at;
		else
			e1.MaxY = e2.MinY = split_at;
		
		objarea = ENVELOPE_AREA(e1) + ENVELOPE_AREA(e2);
		hash_envelope_area_fraction(dh, e1, objarea, l->points);
		hash_envelope_area_fraction(dh, e2, objarea, l->points);

		if (l->gid != -1) // free due to the call to dataset_get_leaf_geo
			GEOSGeom_destroy(geo);

		/*print_geojson_header();
		print_geojson_mbr(l->mbr, "orig");
		print_geojson_mbr(e1, "e1");
		print_geojson_mbr(e2, "e2");
		print_geojson_footer();*/
	}
	else {
		hash_envelope_area_fraction(dh, l->mbr, objarea, l->points);
	}
}


void histogram_generate_cells_fix(dataset *ds, double psizex, double psizey, enum HistogramHashMethod hm, enum JoinPredicateCheck pcheck) {

	dataset_histogram *dh = &ds->metadata.hist;
	dh->xsize = psizex;
	dh->ysize = psizey;
	//printf("Generating histogram of size: %d x %d\n", dh->xqtd, dh->yqtd);

	// X
	double xini = ds->metadata.hist.mbr.MinX;
	for(int i = 0; i < dh->xqtd; i++)
		dh->xtics[i] = xini + (psizex * i);
	dh->xtics[dh->xqtd] = ds->metadata.hist.mbr.MaxX;

	// Y
	double yini = ds->metadata.hist.mbr.MinY;
	for(int i = 0; i < dh->yqtd; i++)
		dh->ytics[i] = yini + (psizey * i);
	dh->ytics[dh->yqtd] = ds->metadata.hist.mbr.MaxY;

	int splitted = 0;
	dataset_iter di;
	dataset_foreach(di, ds) {
		dataset_leaf *l = get_join_pair_leaf(di.item, pcheck);

		switch (hm) {
			case HHASH_CENTROID:
				fill_hist_cell_centroid(l, ds, dh);
				break;
	
			case HHASH_MBRCENTER: 
				fill_hist_cell_mbr_center(l, ds, dh);
				break;

			case HHASH_AREAFRAC:
				fill_hist_cell_area_fraction(l, ds, dh);

				break;

			case HHASH_AREAFRACSPLIT:
				splitted += fill_hist_cell_area_fraction_with_split(l, ds, dh);
				break;

			default:
				printf("Histogram method not defined.\n");
		}
	}

	//TODO: Remove when implement avg online for avgwidth and avgheigt 
	for(int x = 0; x < dh->xqtd; x++) {
		for(int y = 0; y < dh->yqtd; y++) {
			histogram_cell *c = GET_HISTOGRAM_CELL(dh, x, y);
			if (c->cardin > 0.0) {
				c->avgwidth = c->avgwidth / c->cardin;
				c->avgheight = c->avgheight / c->cardin;
			}
		}
	}

	if (hm == HHASH_AREAFRACSPLIT)
		printf("Areafs splitted objects: %d\n", splitted);
}

void histogram_generate_avg(dataset *ds, enum HistogramHashMethod hm, enum JoinPredicateCheck pcheck) {

	double rangex = ds->metadata.hist.mbr.MaxX - ds->metadata.hist.mbr.MinX;
	double rangey = ds->metadata.hist.mbr.MaxY - ds->metadata.hist.mbr.MinY;

	double psizex = ds->metadata.x_average;
	double psizey = ds->metadata.y_average;

	dataset_histogram *dh = &ds->metadata.hist;
	histogram_alloc(dh, ceil(rangex / psizex), ceil(rangey / psizey));

	histogram_generate_cells_fix(ds, psizex, psizey, hm, pcheck);

};

void histogram_generate_hw(dataset *ds, double x, double y, enum HistogramHashMethod hm, enum JoinPredicateCheck pcheck) {

	double rangex = ds->metadata.hist.mbr.MaxX - ds->metadata.hist.mbr.MinX;
	double rangey = ds->metadata.hist.mbr.MaxY - ds->metadata.hist.mbr.MinY;

	double psizex = x;
	double psizey = y;

	const int MAX = 1000;
	
	if ((rangex / psizex) > MAX)
		psizex = rangex / MAX;

	if ((rangey / psizey) > MAX)
		psizey = rangey / MAX;

	dataset_histogram *dh = &ds->metadata.hist;
	histogram_alloc(dh, ceil(rangex / psizex), ceil(rangey / psizey));

	histogram_generate_cells_fix(ds, psizex, psizey, hm, pcheck);

};

void histogram_generate_fix(dataset *ds, int fsizex, int fsizey, enum HistogramHashMethod hm, enum JoinPredicateCheck pcheck) {

	double rangex = ds->metadata.hist.mbr.MaxX - ds->metadata.hist.mbr.MinX;
	double rangey = ds->metadata.hist.mbr.MaxY - ds->metadata.hist.mbr.MinY;

	double psizex = rangex / fsizex;
	double psizey = rangey / fsizey;

	dataset_histogram *dh = &ds->metadata.hist;
	histogram_alloc(dh, fsizex, fsizey);

	histogram_generate_cells_fix(ds, psizex, psizey, hm, pcheck);
	
};

void histogram_build_metadata(dataset *ds, enum JoinPredicateCheck pcheck) {

	char first = 1;
	dataset_iter di;
    dataset_foreach(di, ds) {
		double x, y;
		dataset_leaf *leaf = get_join_pair_leaf(di.item, pcheck);
		
		// metadata: dataset extent
		double new_x = leaf->mbr.MaxX - leaf->mbr.MinX;
		double new_y = leaf->mbr.MaxY - leaf->mbr.MinY;
		//printf("%20.10f\n", new_x);
		double oldmx, oldmy;
		if (first) {
			first = 0;
			ds->metadata.hist.mbr = leaf->mbr;
			ds->metadata.x_average = new_x;
			ds->metadata.y_average = new_y;
			ds->metadata.x_psa = ds->metadata.y_psa = 0.0;
		}
		else {
			ENVELOPE_MERGE(ds->metadata.hist.mbr, leaf->mbr);

			double oldmx = ds->metadata.x_average;
			double oldmy = ds->metadata.y_average;
			ds->metadata.x_average += (new_x - ds->metadata.x_average) / ds->metadata.count;
			ds->metadata.y_average += (new_y - ds->metadata.y_average) / ds->metadata.count;
			ds->metadata.x_psa += (new_x - oldmx) * (new_x - ds->metadata.x_average);
			ds->metadata.y_psa += (new_y - oldmy) * (new_y - ds->metadata.y_average);
		}
	}
}

void histogram_generate(dataset *ds, HistogramGenerateSpec spec, enum JoinPredicateCheck pcheck) {

	histogram_build_metadata(ds, pcheck);

	if (spec.sm == HSPLIT_FIX)
		histogram_generate_fix(ds, spec.xqtd, spec.yqtd, spec.hm, pcheck);
	else if (spec.sm == HSPLIT_AVG)
		histogram_generate_hw(ds, ds->metadata.x_average, ds->metadata.y_average, spec.hm, pcheck);
	else if (spec.sm == HSPLIT_AVG_STD)
		histogram_generate_hw(ds, 
			ds->metadata.x_average + dataset_meta_stddev(ds->metadata, x), 
			ds->metadata.y_average + dataset_meta_stddev(ds->metadata, y),
			spec.hm, pcheck);
	else {
		fprintf(stderr, "Histogram Split Method not found.\n");
	}

	printf("Generated histogram %d x %d, %s.\n", ds->metadata.hist.xqtd,
		ds->metadata.hist.yqtd, HistogramHashMethodName[spec.hm]);
}



int histogram_join_cardinality(dataset *dr, dataset *ds) {
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
					char i = ENVELOPE_INTERSECTS(er, es);
					assert(i != 0);

					/*printf("r[%d][%d] %d x s[%d][%d]%d = %d\n",
						xr, yr, hr->hcells[xr*hr->yqtd + yr],
						xs, ys, hs->hcells[xs*hs->yqtd + ys],
						hr->hcells[xr*hr->yqtd + yr] * hs->hcells[xs*hs->yqtd + ys]);*/

					Envelope inters = EnvelopeIntersection2(er, es);
					double intarea = ENVELOPE_AREA(inters);
					double hrfraction = intarea / erarea;
					double hsfraction = intarea / ENVELOPE_AREA(es);

					histogram_cell *rcell = &hr->hcells[xr*hr->yqtd + yr];
					histogram_cell *scell = &hs->hcells[xs*hs->yqtd + ys];

					double qtdobjr = hrfraction * rcell->cardin;
					double qtdobjs = hsfraction * scell->cardin;
					double intersections = 0;
					if (rcell->cardin >= 1 || scell->cardin >= 1) {
						intersections = qtdobjr * qtdobjs;
						if (intersections < 1.0)
							intersections = 0;
					}
					result += intersections;
				}
			}
		}
	}

	return (int)result;
}

void histogram_distribute_roundrobin(dataset *ds) {
	if (ds->metadata.servers <= 0)
		return;

	printf("Distributing histogram for %d servers using round-robin.\n", ds->metadata.servers);

	dataset_histogram *hist = &ds->metadata.hist;
	int current_server = 0;

	for(int x = 0; x < hist->xqtd; x++) {
		for(int y = 0; y < hist->yqtd; y++) {
			histogram_cell *cell = &hist->hcells[x*hist->yqtd + y];
			if (cell->points > 0.0) {
				cell->copies = 0;
				cell->place = current_server+1;
				SET_IN_PLACE(cell->copies, cell->place);
				current_server = (current_server+1) % ds->metadata.servers;
			}
		}
	}
}

void histogram_distribute_adjacent_grid(dataset *ds) {
	if (ds->metadata.servers <= 0)
		return;

	printf("Distributing histogram for %d servers using adjacent grid.\n", ds->metadata.servers);

	dataset_histogram *hist = &ds->metadata.hist;
	int current_server = 0;

	int qtd = ceil(sqrt(ds->metadata.servers));
	int qtdx = ceil(hist->xqtd / (double)qtd);
	int qtdy = ceil(hist->yqtd / (double)qtd);
	for(int x = 0; x < hist->xqtd; x++) {
		for(int y = 0; y < hist->yqtd; y++) {
			histogram_cell *cell = &hist->hcells[x*hist->yqtd + y];
			if (cell->points > 0) {
				current_server = x/qtdx * qtd + y/qtdy;
				cell->copies = 0;
				cell->place = (current_server % ds->metadata.servers)+1;
				SET_IN_PLACE(cell->copies, cell->place);
			}
		}
	}
}

void histogram_distribute(dataset *ds) {
	histogram_distribute_roundrobin(ds);
	//histogram_distribute_adjacent_grid(ds);
}

void histogram_print(dataset *ds, histogram_type type) {
	char filename[100];
	const char *typenames[] = {"card", "points", "places"};

	char *prefix = getenv("HISTOPREFIX");
	prefix = prefix != NULL ? prefix : "";

	sprintf(filename, "histogram/%shist-%s-%s.dat", prefix, typenames[type], ds->metadata.name);
	FILE *f = fopen(filename, "wb");
	if (f == NULL) {
		perror("Error printing histogram");
		return;
	}

	dataset_histogram *hist = &ds->metadata.hist;
	for(int x = 0; x < hist->xqtd; x++) {
		for(int y = 0; y < hist->yqtd; y++) {
			double value;
			switch (type) {
				case CARDIN: value = hist->hcells[x*hist->yqtd + y].cardin; break;
				case POINTS: value = hist->hcells[x*hist->yqtd + y].points; break;
				case PLACES: value = hist->hcells[x*hist->yqtd + y].place; break;
			}
 
			fprintf(f, "%10.5f %10.5f %lf\n", hist->xtics[x],
				hist->ytics[y], value);
		}
	}
	fclose(f);
}

void histogram_print_geojson(dataset *ds) {
	char filename[100];

	char *prefix = getenv("HISTOPREFIX");
	prefix = prefix != NULL ? prefix : "";

	sprintf(filename, "histogram/%shist-%s.geojson", prefix, ds->metadata.name);
	FILE *f = fopen(filename, "wb");
	if (f == NULL) {
		perror("Error printing histogram");
		return;
	}
	
	fprintf(f, "{'type': 'FeatureCollection', 'features': [\n");

	Envelope e;
	dataset_histogram *hist = &ds->metadata.hist;
	for(int x = 0; x < hist->xqtd; x++) {
		e.MinX = hist->xtics[x];
		e.MaxX = hist->xtics[x+1]; 
		for(int y = 0; y < hist->yqtd; y++) {
			e.MinY = hist->ytics[y];
			e.MaxY = hist->ytics[y+1];

			fprintf(f, "{\"type\": \"Feature\", \"geometry\": {\"type\": \"Polygon\", \"coordinates\": [[");
			fprintf(f, "[%f, %f],", e.MinX, e.MinY);
			fprintf(f, "[%f, %f],", e.MaxX, e.MinY);
			fprintf(f, "[%f, %f],", e.MaxX, e.MaxY);
			fprintf(f, "[%f, %f],", e.MinX, e.MaxY);
			fprintf(f, "[%f, %f]",  e.MinX, e.MinY);
			fprintf(f, "]]}, 'properties': {");
			fprintf(f, "\"name\": \"%d.%d\",", x, y);
			fprintf(f, "\"card\": %f,", hist->hcells[x*hist->yqtd + y].cardin);
			fprintf(f, "\"points\": %f,", hist->hcells[x*hist->yqtd + y].points);
			fprintf(f, "\"place\": %d,", hist->hcells[x*hist->yqtd + y].place);
			fprintf(f, "\"avgwidth\": %f,", hist->hcells[x*hist->yqtd + y].avgwidth);
			fprintf(f, "\"avgheight\": %f,", hist->hcells[x*hist->yqtd + y].avgheight);
			fprintf(f, "}},\n");
		}
	}

	fprintf(f, "]}\n");
	fclose(f);
}

double get_to_pnts(const void *data, const int n) {
	return ((multiway_histogram_estimate*)data)[n].to_pnts;
}

double get_io_pnts(const void *data, const int n) {
	return ((multiway_histogram_estimate*)data)[n].io_pnts;
}

void histogram_print_estimative(char *name, multiway_histogram_estimate *estimate, int servers) {
	int to_pnts = 0, io_pnts = 0;
	printf("%6s     Points  IO Points\n", name);
	printf("------ ---------- ----------\n");
    int max_to_pnts = (int)estimate[1].to_pnts;
    int max_io_pnts = (int)estimate[1].io_pnts;
    int min_to_pnts = (int)estimate[1].to_pnts;
    int min_io_pnts = (int)estimate[1].io_pnts;
	for(int s = 1; s <= servers; s++) {
		printf("%6d %10d %10d\n", s, (int)estimate[s].to_pnts, (int)estimate[s].io_pnts);		
		to_pnts += (int)estimate[s].to_pnts;
		io_pnts += (int)estimate[s].io_pnts;
        if (max_to_pnts < (int)estimate[s].to_pnts)
            max_to_pnts = (int)estimate[s].to_pnts;
        if (max_io_pnts < (int)estimate[s].io_pnts)
            max_io_pnts = (int)estimate[s].io_pnts;
        if (min_to_pnts > (int)estimate[s].to_pnts)
            min_to_pnts = (int)estimate[s].to_pnts;
        if (min_io_pnts > (int)estimate[s].io_pnts)
            min_io_pnts = (int)estimate[s].io_pnts;
	}
	printf("------ ---------- ----------\n");
	printf("total  %10d %10d\n", to_pnts, io_pnts);
	printf("max    %10d %10d\n", max_to_pnts, max_io_pnts);
	printf("min    %10d %10d\n", min_to_pnts, min_io_pnts);

	double to_stdev = stdevd_ex(estimate, 1, servers+1, get_to_pnts);
	double io_stdev = stdevd_ex(estimate, 1, servers+1, get_io_pnts);
	printf("stdev  %10.1f %10.1f\n\n", to_stdev, io_stdev);
}

double histogram_search_hist(dataset_histogram *dh, Envelope query) {

	#define ZU(x,y) (x<y?0:x)
	//#define ZU(x,y) (x)

	// prevent values of x and y out of histogram bounds
	int xini, xfim, yini, yfim;
	get_xini_xfim(dh, query, &xini, &xfim, &yini, &yfim);	

	double result = 0.0;
	int x;
	for(x = xini; x <= xfim; x++) {
		Envelope rs;
		rs.MinX = dh->xtics[x];
		rs.MaxX = dh->xtics[x+1];

		for(int y = yini; y <= yfim; y++) {
			rs.MinY = dh->ytics[y];
			rs.MaxY = dh->ytics[y+1];

			histogram_cell *c = GET_HISTOGRAM_CELL(dh, x, y);
			if (ENVELOPE_INTERSECTS(query, rs)) {

				Envelope inters = EnvelopeIntersection(query, rs);
				double query_x = inters.MaxX - inters.MinX;
				double query_y = inters.MaxY - inters.MinY;
				double univ_x = rs.MaxX - rs.MinX;
				double univ_y = rs.MaxY - rs.MinY;
				double avg_x = c->avgwidth;
				double avg_y = c->avgheight;

				// observing that objects generally doesn't overlap in both axis,
				// 	reduce the probability of intersection in one of them
				if (c->avgwidth > c->avgheight)
					avg_x = MIN(univ_x/c->cardin, avg_x);
				else
					avg_y = MIN(univ_y/c->cardin, avg_y);


				double fraction = MIN(1.0, (ZU(avg_x, query_x) + query_x)/univ_x) 
								* MIN(1.0, (ZU(avg_y, query_y) + query_y)/univ_y);
				/*double fraction = MIN(1.0, (query_x)/univ_x) 
								* MIN(1.0, (query_y)/univ_y);*/
				assert(fraction <= 1.0 && "fraction should not be higher than 1.0");

				result += c->cardin * fraction;

				/*double int_area = ENVELOPE_AREA(inters);
				double bucket_area = ENVELOPE_AREA(rs);
				double fraction = int_area / bucket_area;
				result += fraction * c->cardin;*/
			}
		}
	}
	return result;

}

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
	cell->objcount++;
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
	if (!ENVELOPE_INTERSECTS(query, dh->mbr)) {
		*xini = *yini = 0;
		*xfim = *yfim = -1;
	} else {
		query = EnvelopeIntersection2(query, dh->mbr);

		*xini = (query.MinX - dh->mbr.MinX) / dh->xsize;
		*xfim = (query.MaxX - dh->mbr.MinX) / dh->xsize;
		*yini = (query.MinY - dh->mbr.MinY) / dh->ysize;
		*yfim = (query.MaxY - dh->mbr.MinY) / dh->ysize;

		if (*xfim == dh->xqtd) (*xfim)--;
		if (*yfim == dh->yqtd) (*yfim)--;
   
		while (dh->xtics[*xini] > query.MinX) (*xini)--;
		while (dh->xtics[(*xfim)+1] < query.MaxX) (*xfim)++;
		while (dh->ytics[*yini] > query.MinY) (*yini)--;
		while (dh->ytics[(*yfim)+1] < query.MaxY) (*yfim)++;

		assert(*xini >= 0 && *xfim < dh->xqtd && "x is out of histogram bounds.");
		assert(*yini >= 0 && *yfim < dh->yqtd && "y is out of histogram bounds.");
	}
}

//inline __attribute__((always_inline))
double hash_envelope_area_fraction(dataset_histogram *dh, Envelope ev, double objarea, double points) {

	int xini, xfim, yini, yfim;
	get_xini_xfim(dh, ev, &xini, &xfim, &yini, &yfim);	


	double sum_fraction = 0;
	for(int x = xini; x <= xfim; x++) {
		Envelope rs;
		rs.MinX = dh->xtics[x];
		rs.MaxX = dh->xtics[x+1];

		for(int y = yini; y <= yfim; y++) {
			rs.MinY = dh->ytics[y];
			rs.MaxY = dh->ytics[y+1];

			Envelope inters = EnvelopeIntersection(ev, rs);
			double intarea = ENVELOPE_AREA(inters);
			double fraction;
			if (intarea <= 0.0) { // parallel to one axis
				bool parallel_y = (ev.MaxX - ev.MinX < 1e-30);
				bool parallel_x = (ev.MaxY - ev.MinY < 1e-30);
				if (parallel_x && parallel_y) // point obj
					fraction = 1;
				else if (parallel_x)
					// the part of ev inside rs / the length of rs in X
					fraction = (MIN(rs.MaxX, ev.MaxX) - MAX(rs.MinX, ev.MinX)) / (ev.MaxX - ev.MinX);
				else
					fraction = (MIN(rs.MaxY, ev.MaxY) - MAX(rs.MinY, ev.MinY)) / (ev.MaxY - ev.MinY);
			}
			else {
				fraction = intarea / objarea;
			}
			sum_fraction += fraction;

			histogram_cell *cell = &dh->hcells[x*dh->yqtd +y];
			//cell->cardin += 1.0;//fraction;
			cell->cardin += fraction;
			cell->points += points; //object is replicated

			//TODO: calculate avg online
			double delta_x = (inters.MaxX - inters.MinX);
			double delta_y = (inters.MaxY - inters.MinY);
			cell->avgwidth += delta_x;
			cell->avgheight += delta_y;

			cell->objcount += 1;
		}
	}
	return sum_fraction;
}

inline __attribute__((always_inline))
double fill_hist_cell_area_fraction(dataset_leaf *l, dataset *ds, dataset_histogram *dh) {
	// proportional to cover area

	//printf("GID %lld: ", l->gid);

	double objarea = ENVELOPE_AREA(l->mbr);
	return hash_envelope_area_fraction(dh, l->mbr, objarea, l->points);

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

//inline __attribute__((always_inline))
double fill_hist_cell_area_fraction_with_split(dataset_leaf *l, dataset *ds, dataset_histogram *dh) {
	// proportional to cover area with split

	Envelope eaux = {DBL_MAX, DBL_MAX, -DBL_MAX, -DBL_MAX};
	double sum_fraction = 0;
	bool splitted;

	int xini, xfim, yini, yfim;
	get_xini_xfim(dh, l->mbr, &xini, &xfim, &yini, &yfim);	

	double objarea = ENVELOPE_AREA(l->mbr);

	// is a candidate for split?
	int xspan = xfim - xini + 1;
	int yspan = yfim - yini + 1;
	if (xspan >= 2 || yspan >= 2) { // more than two cells?
		GEOSGeometryH geo = dataset_get_leaf_geo(ds, l);
		
		int envelope_hit[xspan][yspan];
		Envelope envelopes[xspan][yspan];
		Envelope cells_envs[xspan][yspan];
		for(int xi = 0; xi < xspan; xi++) {
			for(int yi = 0; yi < yspan; yi++) {
				envelopes[xi][yi] = eaux;
				envelope_hit[xi][yi] = 0;
				cells_envs[xi][yi].MinX = dh->xtics[xini+xi];
				cells_envs[xi][yi].MaxX = dh->xtics[xini+xi+1]; 
				cells_envs[xi][yi].MinY = dh->ytics[yini+yi];
				cells_envs[xi][yi].MaxY = dh->ytics[yini+yi+1];
			}
		}

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
			assert(numPoints >= 2);

			double xc1, yc1;
       		GEOSCoordSeq_getX(coordSeq, 0, &xc1);
	        GEOSCoordSeq_getY(coordSeq, 0, &yc1);
			for (int p=1; p < numPoints; p++) {
				double xc2, yc2;
        		GEOSCoordSeq_getX(coordSeq, p, &xc2);
		        GEOSCoordSeq_getY(coordSeq, p, &yc2);
				Envelope eseg = {DBL_MAX, DBL_MAX, -DBL_MAX, -DBL_MAX};
				envelope_update(&eseg, xc1, yc1);
				envelope_update(&eseg, xc2, yc2);
		
				bool point_used = false;
				for(int xi = 0; xi < xspan; xi++) {
					for(int yi = 0; yi < yspan; yi++) {
						if (ENVELOPE_INTERSECTS(eseg, cells_envs[xi][yi])) {
							Envelope inters = EnvelopeIntersection2(eseg, cells_envs[xi][yi]);
							envelope_update(&envelopes[xi][yi], inters.MinX, inters.MinY);
							envelope_update(&envelopes[xi][yi], inters.MaxX, inters.MaxY);
							envelope_hit[xi][yi]++;
							point_used = true;
						}
					}
				}
				if (!point_used) {
					assert(point_used && "Point not considered while updating partial envelopes.");
				}

				xc1 = xc2;
				yc1 = yc2;
			}
			
		}
		
		double newobjarea = 0;
		for(int xi = 0; xi < xspan; xi++) {
			for(int yi = 0; yi < yspan; yi++) {
				if (envelope_hit[xi][yi] > 0) {
					double aux = ENVELOPE_AREA(envelopes[xi][yi]);
					if (aux == 0.0) {
						// parallel to x or y axis
						// add a small X or Y difference
						// this may induces error on objects that has large extents parallel to axis
						Envelope ea = envelopes[xi][yi];
						if (ea.MaxX == ea.MinX)
							ea.MinX -= 1e-10;
						else
							ea.MinY -= 1e-10;
						aux = ENVELOPE_AREA(ea);
						envelopes[xi][yi] = ea;
					}
					newobjarea += aux;
				}
			}
		}

		for(int xi = 0; xi < xspan; xi++) {
			for(int yi = 0; yi < yspan; yi++) {
				if (envelope_hit[xi][yi] > 0) {
					histogram_cell *cell = GET_HISTOGRAM_CELL(dh, xini+xi, yini+yi);

					Envelope inters = envelopes[xi][yi];
					double intarea = ENVELOPE_AREA(inters); 
					double fraction = intarea / newobjarea;

					cell->cardin += fraction;
					cell->points += l->points; //object is replicated

					//TODO: calculate avg online
					double delta_x = (inters.MaxX - inters.MinX);
					double delta_y = (inters.MaxY - inters.MinY);
					cell->avgwidth += delta_x;
					cell->avgheight += delta_y;
					cell->objcount += 1;

					sum_fraction += fraction;
				}
			}
		}

		if (l->gid != -1) // free due to the call to dataset_get_leaf_geo
			GEOSGeom_destroy(geo);

		splitted = true;
		double aux = round(sum_fraction*100000.0)/100000.0;
		if (aux < 1.0 || aux > 1.0) {
			printf("Splitted %d, Fraction: %f\n", splitted, sum_fraction);
			printf("%d %d %d %d\n", xini, xfim, yini, yfim);
			printf("%f %f %f %f\n", dh->xtics[xini], dh->xtics[xfim], dh->ytics[yini], dh->ytics[yfim]);

			print_geojson_header();
			print_geojson_mbr(l->mbr, "orig");
			char name[20];
			for(int xi = 0; xi < xspan; xi++) {
				for(int yi = 0; yi < yspan; yi++) {
					if (envelope_hit[xi][yi] > 0) {
						sprintf(name, "e_%d_%d", xi, yi);
						print_geojson_mbr(envelopes[xi][yi], name);
					}
				}
			}
			print_geojson_footer();
		}

	}
	else {
		splitted = false;
		sum_fraction += hash_envelope_area_fraction(dh, l->mbr, objarea, l->points);
	}

	return sum_fraction;
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

	double sum = 0;
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
				sum += fill_hist_cell_area_fraction(l, ds, dh);

				break;

			case HHASH_AREAFRACSPLIT:
				sum += fill_hist_cell_area_fraction_with_split(l, ds, dh);
				break;

			default:
				printf("Histogram method not defined.\n");
		}
	}

	//TODO: Remove when implement avg online for avgwidth and avgheigt 
	for(int x = 0; x < dh->xqtd; x++) {
		for(int y = 0; y < dh->yqtd; y++) {
			histogram_cell *c = GET_HISTOGRAM_CELL(dh, x, y);
			if (c->objcount > 0.0) {
				c->avgwidth = c->avgwidth / c->objcount;
				c->avgheight = c->avgheight / c->objcount;
			}
		}
	}

//	if (hm == HHASH_AREAFRACSPLIT)
//		printf("Areafs splitted objects: %d\n", splitted);
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

double histogram_search_hist_wao(dataset_histogram *dh, Envelope query) {

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
					avg_x = MIN(univ_x/c->objcount, avg_x);
				else
					avg_y = MIN(univ_y/c->objcount, avg_y);

				/*double fraction = MIN(1.0, (ZU(avg_x, query_x) + query_x)/univ_x) 
								* MIN(1.0, (ZU(avg_y, query_y) + query_y)/univ_y);*/
				double fraction = MIN(1.0, (avg_x + query_x)/univ_x) 
								* MIN(1.0, (avg_y + query_y)/univ_y);
				assert(fraction <= 1.0 && "fraction should not be higher than 1.0");

				result += c->cardin * fraction;
			}
		}
	}
	return ceil(result);

}

double histogram_search_hist_mp(dataset_histogram *dh, Envelope query) {

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

				double fraction = MIN(1.0, (avg_x + query_x)/univ_x) 
								* MIN(1.0, (avg_y + query_y)/univ_y);
				assert(fraction <= 1.0 && "fraction should not be higher than 1.0");

				result += c->cardin * fraction;
			}
		}
	}
	return result;

}


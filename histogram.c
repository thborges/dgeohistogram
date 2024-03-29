/*
 * histogram.c
 *
 *  Created on: 09/08/2014
 *      Author: thborges
 */

#include "histogram.h"

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
}

inline __attribute__((always_inline))
void fill_hist_cell_area_fraction(dataset_leaf *l, dataset *ds, dataset_histogram *dh) {
	// proportional to cover area

	int xini = (l->mbr.MinX - ds->metadata.hist.mbr.MinX) / dh->xsize;
	int xfim = (l->mbr.MaxX - ds->metadata.hist.mbr.MinX) / dh->xsize;
	int yini = (l->mbr.MinY - ds->metadata.hist.mbr.MinY) / dh->ysize;
	int yfim = (l->mbr.MaxY - ds->metadata.hist.mbr.MinY) / dh->ysize;
	double objarea = ENVELOPE_AREA(l->mbr);

	for(int x = xini; x <= xfim; x++) {
		Envelope rs;
		rs.MinX = dh->xtics[x];
		rs.MaxX = dh->xtics[x+1];

		for(int y = yini; y <= yfim; y++) {
			rs.MinY = dh->ytics[y];
			rs.MaxY = dh->ytics[y+1];

			Envelope inters = EnvelopeIntersection(l->mbr, rs);
			double intarea = ENVELOPE_AREA(inters);
			double rsarea = ENVELOPE_AREA(rs);

			// for point objects, objarea == 0.0
			double fraction = (objarea == 0.0) ? 1.0 : intarea / objarea;

			histogram_cell *cell = &ds->metadata.hist.hcells[x*ds->metadata.hist.yqtd +y];
			cell->cardin += fraction;
			cell->points += l->points; //object is replicated
		}
	}

	//printf("%10.5f <= %10.5f <= %10.5f\n", ds->metadata.hist.xtics[xp], x, ds->metadata.hist.xtics[xp+1]);
	//assert(x >= ds->metadata.hist.xtics[xp] && x <= ds->metadata.hist.xtics[xp+1]);
	//assert(y >= ds->metadata.hist.ytics[yp] && y <= ds->metadata.hist.ytics[yp+1]);
}

void histogram_generate_cells_fix(dataset *ds, double psizex, double psizey, enum JoinPredicateCheck pcheck) {

	dataset_histogram *dh = &ds->metadata.hist;
	dh->xsize = psizex;
	dh->ysize = psizey;
	printf("Generating histogram of size: %d x %d\n", dh->xqtd, dh->yqtd);

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

	dataset_iter di;
	dataset_foreach(di, ds) {
		dataset_leaf *l = get_join_pair_leaf(di.item, pcheck);

		//fill_hist_cell_centroid(l, ds, dh);
		//fill_hist_cell_mbr_center(l, ds, dh);
		fill_hist_cell_area_fraction(l, ds, dh);
	}
}

void histogram_generate_avg(dataset *ds, enum JoinPredicateCheck pcheck) {

	double rangex = ds->metadata.hist.mbr.MaxX - ds->metadata.hist.mbr.MinX;
	double rangey = ds->metadata.hist.mbr.MaxY - ds->metadata.hist.mbr.MinY;

	double psizex = ds->metadata.x_average;
	double psizey = ds->metadata.y_average;

	dataset_histogram *dh = &ds->metadata.hist;
	histogram_alloc(dh, rangex / psizex + 1, rangey / psizey + 1);

	histogram_generate_cells_fix(ds, psizex, psizey, pcheck);

};

void histogram_generate_hw(dataset *ds, double x, double y, enum JoinPredicateCheck pcheck) {

	double rangex = ds->metadata.hist.mbr.MaxX - ds->metadata.hist.mbr.MinX;
	double rangey = ds->metadata.hist.mbr.MaxY - ds->metadata.hist.mbr.MinY;

	double psizex = x;
	double psizey = y;

	const int MAX = 100;
	
	if ((rangex / psizex) > MAX)
		psizex = rangex / MAX;

	if ((rangey / psizey) > MAX)
		psizey = rangey / MAX;

	dataset_histogram *dh = &ds->metadata.hist;
	histogram_alloc(dh, rangex / psizex + 1, rangey / psizey + 1);

	histogram_generate_cells_fix(ds, psizex, psizey, pcheck);

};

void histogram_generate_fix(dataset *ds, int fsizex, int fsizey, enum JoinPredicateCheck pcheck) {

	double rangex = ds->metadata.hist.mbr.MaxX - ds->metadata.hist.mbr.MinX;
	double rangey = ds->metadata.hist.mbr.MaxY - ds->metadata.hist.mbr.MinY;

	double psizex = rangex / fsizex;
	double psizey = rangey / fsizey;

	dataset_histogram *dh = &ds->metadata.hist;
	histogram_alloc(dh, fsizex, fsizey);

	histogram_generate_cells_fix(ds, psizex, psizey, pcheck);
	
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

void histogram_generate(dataset *ds, enum JoinPredicateCheck pcheck) {

	histogram_build_metadata(ds, pcheck);

	//histogram_generate_fix(ds, 50, 50, pcheck);
	//histogram_generate_hw(ds, ds->metadata.x_average, ds->metadata.y_average, pcheck);
	
	histogram_generate_hw(ds, 
		ds->metadata.x_average + dataset_meta_stddev(ds->metadata, x), 
		ds->metadata.y_average + dataset_meta_stddev(ds->metadata, y),
		pcheck);
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



#include <float.h>
#include "minskew.h"

typedef struct {
	double mean;
	double variance;
	int n;
} variance_result;

void get_ini_fim(dataset_histogram *dh, Envelope ev, int *xini, int *xfim, int *yini, int *yfim) {
	*xini = (ev.MinX - dh->mbr.MinX) / dh->xsize;
	*xfim = (ev.MaxX - dh->mbr.MinX) / dh->xsize;
	*yini = (ev.MinY - dh->mbr.MinY) / dh->ysize;
	*yfim = (ev.MaxY - dh->mbr.MinY) / dh->ysize;

	const double epsilon = 1e-100;
	if (ev.MaxX - dh->xtics[*xfim] < epsilon && *xfim > 0) {
		(*xfim)--;
	}
	if (ev.MaxY - dh->ytics[*yfim] < epsilon && *yfim > 0) {
		(*yfim)--;
	}
	if (*xfim < *xini)
		*xini = *xfim;
	if (*yfim < *yini)
		*yini = *yfim;
}

variance_result calculate_skew_row_col(dataset_histogram *dh, int xini, int xfim, int yini, int yfim) {
	variance_result r;
	r.n = 0;
	r.mean=0.0;
	double M2=0.0;
	for(int i = xini; i <= xfim; i++) {
		for(int j = yini; j <= yfim; j++) {
			r.n++;
			double v = GET_HISTOGRAM_CELL(dh, i, j)->cardin;
			double delta = v - r.mean;
			r.mean += delta/(double)r.n;
			M2 += delta * (v - r.mean);
		}
	}

	r.variance = r.n < 2 ? 0.0 : M2/(double)r.n;
	return r;
}

void calculate_bucket_with_mbr(dataset_histogram *dh, minskew_bucket *b) {
	int xini, xfim, yini, yfim;
	get_ini_fim(dh, b->mbr, &xini, &xfim, &yini, &yfim);
	variance_result r = calculate_skew_row_col(dh, xini, xfim, yini, yfim);
	b->skew = r.n * r.variance;
	b->cardin = r.n * r.mean;
	b->skew_reduction = NAN;
}

void minskew_calculate_skew_reduction(dataset_histogram *dh, minskew_bucket *bucket) {
	int xini, xfim, yini, yfim;
	get_ini_fim(dh, bucket->mbr, &xini, &xfim, &yini, &yfim);

	int xqtd = xfim-xini+1;
	int yqtd = yfim-yini+1;

	bucket->skew_reduction = 0;

	minskew_bucket aux_bucket1, aux_bucket2;
	for(int x=xini; x<xfim; x++) {
		// divide on x
		aux_bucket1.mbr = bucket->mbr;
		aux_bucket2.mbr = bucket->mbr;
		aux_bucket1.mbr.MaxX = aux_bucket2.mbr.MinX = dh->xtics[x+1];
		calculate_bucket_with_mbr(dh, &aux_bucket1);
		calculate_bucket_with_mbr(dh, &aux_bucket2);

		double new_skew = aux_bucket1.skew + aux_bucket2.skew;
		double reduction = bucket->skew - new_skew;
		if (bucket->skew_reduction < reduction) {
			bucket->skew_reduction = reduction;
			bucket->split_axis = 'x';
			bucket->split_point = x+1;
		}
	}

	for(int y=yini; y<yfim; y++) {
		// divide on y
		aux_bucket1.mbr = bucket->mbr;
		aux_bucket2.mbr = bucket->mbr;
		aux_bucket1.mbr.MaxY = aux_bucket2.mbr.MinY = dh->ytics[y+1];
		calculate_bucket_with_mbr(dh, &aux_bucket1);
		calculate_bucket_with_mbr(dh, &aux_bucket2);

		double new_skew = aux_bucket1.skew + aux_bucket2.skew;
		double reduction = bucket->skew - new_skew;
		if (bucket->skew_reduction < reduction) {
			bucket->skew_reduction = reduction;
			bucket->split_axis = 'y';
			bucket->split_point = y+1;
		}
	}
}

GList *minskew_generate_hist(dataset *ds, int buckets_num) {

	GList *minskewhist = NULL;
	dataset_histogram *dh = &ds->metadata.hist;
	
	minskew_bucket *first_bucket = g_new0(minskew_bucket, 1);
	first_bucket->mbr = ds->metadata.hist.mbr;
	calculate_bucket_with_mbr(&ds->metadata.hist, first_bucket);
	
	int buckets = 1;
	minskewhist = g_list_append(minskewhist, first_bucket);
	double global_skew = first_bucket->skew;

	while (buckets < buckets_num) {
		GList *chosen = NULL;
		double skew_reduction = 0;

		GList *item;
		g_list_foreach(item, minskewhist) {
			minskew_bucket *bucket = (minskew_bucket *)item->data;

			if (isnan(bucket->skew_reduction)) {
				minskew_calculate_skew_reduction(dh, bucket);
			}

			if (skew_reduction < bucket->skew_reduction) {
				chosen = item;
				skew_reduction = bucket->skew_reduction;
			}
		}

		// split the chosen bucket

		minskew_bucket *chosen_bucket = (minskew_bucket *)chosen->data;

		minskew_bucket newb1, newb2;
		newb1.mbr = chosen_bucket->mbr;
		newb2.mbr = chosen_bucket->mbr;
		if (chosen_bucket->split_axis == 'x')
			newb1.mbr.MaxX = newb2.mbr.MinX = dh->xtics[chosen_bucket->split_point];
		else
			newb1.mbr.MaxY = newb2.mbr.MinY = dh->ytics[chosen_bucket->split_point];
		calculate_bucket_with_mbr(dh, &newb1);
		calculate_bucket_with_mbr(dh, &newb2);

		global_skew -= chosen_bucket->skew;
		global_skew += newb1.skew + newb2.skew;

		// substitute the old bucket data with b2
		*chosen_bucket = newb1;

		// add the new bucket b2
		minskew_bucket *newb = g_new(minskew_bucket, 1);
		*newb = newb2;
		minskewhist = g_list_append(minskewhist, newb);


		buckets++;
	}

	return minskewhist;
}

double minskew_search_hist(GList *hist, Envelope query) {
	double result = 0.0;
	GList *item;
	g_list_foreach(item, hist) {
		minskew_bucket *b = (minskew_bucket*)item->data;
		if (ENVELOPE_INTERSECTS(query, b->mbr)) {
			Envelope inters = EnvelopeIntersection(query, b->mbr);
			double int_area = ENVELOPE_AREA(inters);
			double bucket_area = ENVELOPE_AREA(b->mbr);
			double fraction = int_area / bucket_area;
			result += fraction * b->cardin;
		}
	}
	return result;
}

void minskew_print_hist(dataset *ds, GList *hist) {
	char filename[100];

	char *prefix = getenv("HISTOPREFIX");
	prefix = prefix != NULL ? prefix : "";

	sprintf(filename, "histogram/%sminskew-%s.geojson", prefix, ds->metadata.name);
	FILE *f = fopen(filename, "wb");
	if (f == NULL) {
		perror("Error printing histogram");
		return;
	}

	fprintf(f, "{'type': 'FeatureCollection', 'features': [\n");

	int i = 0;
	GList *item;
	g_list_foreach(item, hist) {
		minskew_bucket *bucket = (minskew_bucket *)item->data;
		Envelope e = bucket->mbr;

		fprintf(f, "{\"type\": \"Feature\", \"geometry\": {\"type\": \"Polygon\", \"coordinates\": [[");
		fprintf(f, "[%f, %f],", e.MinX, e.MinY);
		fprintf(f, "[%f, %f],", e.MaxX, e.MinY);
		fprintf(f, "[%f, %f],", e.MaxX, e.MaxY);
		fprintf(f, "[%f, %f],", e.MinX, e.MaxY);
		fprintf(f, "[%f, %f]",  e.MinX, e.MinY);
		fprintf(f, "]]}, 'properties': {");
		fprintf(f, "\"name\": \"%d\",", i);
		fprintf(f, "\"skew\": %f,", bucket->skew);
		fprintf(f, "\"card\": %f,", bucket->cardin);
		fprintf(f, "}},\n");

		i++;
	}

	fprintf(f, "]}\n");
	fclose(f);
}

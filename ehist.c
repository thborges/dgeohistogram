
#include "histogram.h"
#include "ehist.h"

void eh_alloc(dataset *ds, euler_histogram *eh, int xqtd, int yqtd, double psizex, double psizey) {
	assert(xqtd > 0 && yqtd > 0 && "X and Y must be greater than zero.");
	eh->xqtd = xqtd;
	eh->yqtd = yqtd;
	eh->xtics = g_new(double, eh->xqtd+1);
	eh->ytics = g_new(double, eh->yqtd+1);
	eh->faces = g_new0(euler_face, eh->xqtd * eh->yqtd);
	eh->edges = g_new0(euler_edge, (xqtd * yqtd) + yqtd + xqtd);
	eh->vertexes = g_new0(euler_vertex, xqtd * yqtd + xqtd + yqtd + 1);

	eh->xsize = psizex;
	eh->ysize = psizey;
	//printf("Generating histogram of size: %d x %d\n", eh->xqtd, eh->yqtd);

	// X tics
	double xini = ds->metadata.hist.mbr.MinX;
	for(int i = 0; i < eh->xqtd; i++)
		eh->xtics[i] = xini + (psizex * i);
	eh->xtics[eh->xqtd] = ds->metadata.hist.mbr.MaxX;

	// Y tics
	double yini = ds->metadata.hist.mbr.MinY;
	for(int i = 0; i < eh->yqtd; i++)
		eh->ytics[i] = yini + (psizey * i);
	eh->ytics[eh->yqtd] = ds->metadata.hist.mbr.MaxY;

	// edges and vertexes
	int v = 0;
	int e = 0;
	for(int i = 0; i <= eh->xqtd; i++) { // x
		for(int j = 0; j <= eh->yqtd; j++) { // y

			// set vetex at (i,j)
			eh->vertexes[v].x = eh->xtics[i];
			eh->vertexes[v].y = eh->ytics[j];
			v++;

			// horizontal edge at vertex v
			if (i < eh->xqtd) {
				eh->edges[e].mbr.MinX = eh->xtics[i];
				eh->edges[e].mbr.MaxX = eh->xtics[i+1];
				eh->edges[e].mbr.MinY = eh->ytics[j];
				eh->edges[e].mbr.MaxY = eh->ytics[j]+1e-10;
				e++;
			}


			// vertical edge at vertex v
			if (j < eh->yqtd) {
				eh->edges[e].mbr.MinY = eh->ytics[j];
				eh->edges[e].mbr.MaxY = eh->ytics[j+1];
				eh->edges[e].mbr.MinX = eh->xtics[i];
				eh->edges[e].mbr.MaxX = eh->xtics[i]+1e-10;
				e++;
			}
		}
	}
}

void eh_hash_ds_objects(dataset *ds, euler_histogram *eh, enum JoinPredicateCheck pcheck) {

	dataset_iter di;
	dataset_foreach(di, ds) {
		dataset_leaf *l = get_join_pair_leaf(di.item, pcheck);
		Envelope ev = l->mbr;
		GEOSGeometryH geo = dataset_get_leaf_geo(ds, l);	

		int xini = (ev.MinX - eh->mbr.MinX) / eh->xsize;
		int xfim = (ev.MaxX - eh->mbr.MinX) / eh->xsize;
		int yini = (ev.MinY - eh->mbr.MinY) / eh->ysize;
		int yfim = (ev.MaxY - eh->mbr.MinY) / eh->ysize;
		if (xfim < eh->xqtd) xfim++;
		if (yfim < eh->yqtd) yfim++; 

		for(int x = xini; x <= xfim; x++) {
			Envelope rs;
			rs.MinX = eh->xtics[x];
			if (x < eh->xqtd)
				rs.MaxX = eh->xtics[x+1];
			else
				rs.MaxX = rs.MinX + 1e-10; // sum a litle fraction to prevent clip error due to empty mbr

			for(int y = yini; y <= yfim; y++) {
				rs.MinY = eh->ytics[y];
				if (y < eh->yqtd)
					rs.MaxY = eh->ytics[y+1];
				else
					rs.MaxY = rs.MinY + 1e-10; // sum a litle fraction to prevent clip error due to empty mbr

				GEOSGeometryH clipped = GEOSClipByRect(geo, rs.MinX, rs.MinY, rs.MaxX, rs.MaxY);
				if (clipped == NULL)
					continue;
				Envelope ev2;
				GEOSEnvelopeGetXY(clipped, &ev2.MinX, &ev2.MaxX, &ev2.MinY, &ev2.MaxY);
				GEOSGeom_destroy(clipped);

				// face
				if (x < eh->xqtd && y < eh->yqtd) {
					if (ENVELOPE_INTERSECTS(ev2, rs)) {
						euler_face *face = &eh->faces[x*eh->yqtd +y];
						face->cardin += 1;
						//double delta_x = ev2.MaxX - ev2.MinX;
						//face->avg_x += (delta_x - face->avg_x) / face->cardin;
					}
				}

				// vertex
				int v = x * (eh->yqtd+1) + y;
				if (ENVELOPE_CONTAINSP(ev2, eh->vertexes[v].x, eh->vertexes[v].y))
					eh->vertexes[v].cardin += 1;

				// horizontal edge
				if (x < eh->xqtd) {
					int e = x * (2*eh->yqtd+1) + 2*y;
					if (ENVELOPE_INTERSECTS(eh->edges[e].mbr, ev2))
						eh->edges[e].cardin += 1;
				}

				// vertical edge
				if (y < eh->yqtd) {
					int e;
					if (x == eh->xqtd) // last column right border
						e = x * (2*eh->yqtd+1) + y;
					else
						e = x * (2*eh->yqtd+1) + 2*y + 1;
					if (ENVELOPE_INTERSECTS(eh->edges[e].mbr, ev2))
						eh->edges[e].cardin += 1;
				}
			}
		}

		if (l->gid != -1) // free due to the call to dataset_get_leaf_geo
			GEOSGeom_destroy(geo);
	}
}

void eh_generate_hw(dataset *ds, euler_histogram *eh, double x, double y, enum JoinPredicateCheck pcheck) {

	double rangex = ds->metadata.hist.mbr.MaxX - ds->metadata.hist.mbr.MinX;
	double rangey = ds->metadata.hist.mbr.MaxY - ds->metadata.hist.mbr.MinY;

	double psizex = x;
	double psizey = y;

	const int MAX = 1000;
	
	if ((rangex / psizex) > MAX)
		psizex = rangex / MAX;

	if ((rangey / psizey) > MAX)
		psizey = rangey / MAX;

	eh_alloc(ds, eh, ceil(rangex / psizex), ceil(rangey / psizey), psizex, psizey);
	eh_hash_ds_objects(ds, eh, pcheck);
};

void eh_generate_fix(dataset *ds, euler_histogram *eh, int fsizex, int fsizey, enum JoinPredicateCheck pcheck) {

	double rangex = ds->metadata.hist.mbr.MaxX - ds->metadata.hist.mbr.MinX;
	double rangey = ds->metadata.hist.mbr.MaxY - ds->metadata.hist.mbr.MinY;

	double psizex = rangex / fsizex;
	double psizey = rangey / fsizey;

	eh_alloc(ds, eh, fsizex, fsizey, psizex, psizey);
	eh_hash_ds_objects(ds, eh, pcheck);
};

euler_histogram *eh_generate_hist(dataset *ds, HistogramGenerateSpec spec, enum JoinPredicateCheck pcheck) {

	euler_histogram *eh = g_new(euler_histogram, 1);
	eh->mbr = ds->metadata.hist.mbr;

	if (spec.sm == HSPLIT_FIX)
		eh_generate_fix(ds, eh, spec.xqtd, spec.yqtd, pcheck);
	else if (spec.sm == HSPLIT_AVG)
		eh_generate_hw(ds, eh, ds->metadata.x_average, ds->metadata.y_average, pcheck);
	else if (spec.sm == HSPLIT_AVG_STD)
		eh_generate_hw(ds, eh,
			ds->metadata.x_average + dataset_meta_stddev(ds->metadata, x), 
			ds->metadata.y_average + dataset_meta_stddev(ds->metadata, y),
			pcheck);
	else {
		fprintf(stderr, "Histogram Split Method not implemented.\n");
	}

	return eh;
//	printf("Generated histogram %d x %d, %s.\n", ds->metadata.hist.xqtd,
//		ds->metadata.hist.yqtd, HistogramHashMethodName[spec.hm]);
}


int euler_search_hist(euler_histogram *eh, Envelope query2) {

	int result = 0;

	Envelope query = EnvelopeIntersection(query2, eh->mbr);

	int xini = (query.MinX - eh->mbr.MinX) / eh->xsize;
	int xfim = (query.MaxX - eh->mbr.MinX) / eh->xsize;
	int yini = (query.MinY - eh->mbr.MinY) / eh->ysize;
	int yfim = (query.MaxY - eh->mbr.MinY) / eh->ysize;

	for(int x = xini; x <= xfim; x++) {
		Envelope rs;
		rs.MinX = eh->xtics[x];
		if (x < eh->xqtd)
			rs.MaxX = eh->xtics[x+1];

		for(int y = yini; y <= yfim; y++) {
			rs.MinY = eh->ytics[y];
			if (y < eh->yqtd)
				rs.MaxY = eh->ytics[y+1];

			// face
			if (x < eh->xqtd && y < eh->yqtd) {
				if (ENVELOPE_INTERSECTS(query, rs)) {
					euler_face *face = &eh->faces[x*eh->yqtd +y];

					Envelope inters = EnvelopeIntersection(query, rs);
					double int_area = ENVELOPE_AREA(inters);
					double face_area = ENVELOPE_AREA(rs);
					double fraction = int_area / face_area;
					result += fraction * face->cardin;
				}
			}

			// vertex
			int v = x * (eh->yqtd+1) + y;
			if (ENVELOPE_CONTAINSP(query, eh->vertexes[v].x, eh->vertexes[v].y))
				result += eh->vertexes[v].cardin;

			// horizontal edge
			if (x < eh->xqtd) {
				int e = x * (2*eh->yqtd+1) + 2*y;
				if (ENVELOPE_INTERSECTS(eh->edges[e].mbr, query)) {
					Envelope inters = EnvelopeIntersection(query, eh->edges[e].mbr);
					double int_length = inters.MaxX - inters.MinX;
					double fraction = int_length / (eh->edges[e].mbr.MaxX - eh->edges[e].mbr.MinX);
					result -= fraction * eh->edges[e].cardin;
				}
			}

			// vertical edge
			if (y < eh->yqtd) {
				int e;
				if (x == eh->xqtd) // last column right border
					e = x * (2*eh->yqtd+1) + y;
				else
					e = x * (2*eh->yqtd+1) + 2*y + 1;
				if (ENVELOPE_INTERSECTS(eh->edges[e].mbr, query)) {
					Envelope inters = EnvelopeIntersection(query, eh->edges[e].mbr);
					double int_length = inters.MaxY - inters.MinY;
					double fraction = int_length / (eh->edges[e].mbr.MaxY - eh->edges[e].mbr.MinY);
					result -= fraction * eh->edges[e].cardin;
				}
			}
		}
	}

	return result;
}

void euler_print_hist(dataset *ds, euler_histogram *eh) {
	char filename[100];

	char *prefix = getenv("HISTOPREFIX");
	prefix = prefix != NULL ? prefix : "";

	sprintf(filename, "histogram/%seuler-%s.geojson", prefix, ds->metadata.name);
	FILE *f = fopen(filename, "wb");
	if (f == NULL) {
		perror("Error printing histogram");
		return;
	}
	
	fprintf(f, "{'type': 'FeatureCollection', 'features': [\n");

	int e = 0;
	int v = 0;
	Envelope env;
	for(int x = 0; x <= eh->xqtd; x++) {
		env.MinX = eh->xtics[x];
		if (x < eh->xqtd)
			env.MaxX = eh->xtics[x+1]; 
	
		for(int y = 0; y <= eh->yqtd; y++) {
			env.MinY = eh->ytics[y];
			if (y < eh->yqtd)
				env.MaxY = eh->ytics[y+1];

			// face
			if (x < eh->xqtd && y < eh->yqtd) {
				fprintf(f, "{\"type\": \"Feature\", \"geometry\": {\"type\": \"Polygon\", \"coordinates\": [[");
				fprintf(f, "[%f, %f],", env.MinX, env.MinY);
				fprintf(f, "[%f, %f],", env.MaxX, env.MinY);
				fprintf(f, "[%f, %f],", env.MaxX, env.MaxY);
				fprintf(f, "[%f, %f],", env.MinX, env.MaxY);
				fprintf(f, "[%f, %f]",  env.MinX, env.MinY);
				fprintf(f, "]]}, 'properties': {");
				fprintf(f, "\"name\": \"f:%d.%d\",", x, y);
				fprintf(f, "\"card\": %f,", eh->faces[x*eh->yqtd + y].cardin);
				fprintf(f, "\"type\": \"face\",");
				fprintf(f, "}},\n");
			}

			// vertex
			fprintf(f, "{\"type\": \"Feature\", \"geometry\": {\"type\": \"Point\", \"coordinates\": [%f, %f]},", 
				eh->vertexes[v].x, eh->vertexes[v].y);
			fprintf(f, "'properties': {");
			fprintf(f, "\"name\": \"v:%d.%d\",", x, y);
			fprintf(f, "\"card\": %f,", eh->vertexes[v].cardin);
			fprintf(f, "\"type\": \"vertex\",");
			fprintf(f, "}},\n");
			v++;

			
			// horizontal edge
			if (x < eh->xqtd) {
				fprintf(f, "{\"type\": \"Feature\", \"geometry\": {\"type\": \"LineString\", \"coordinates\": [[%f, %f], [%f, %f]]},", 
					eh->edges[e].mbr.MinX, eh->edges[e].mbr.MinY, eh->edges[e].mbr.MaxX, eh->edges[e].mbr.MaxY);
				fprintf(f, "'properties': {");
				fprintf(f, "\"name\": \"eh:%d.%d\",", x, y);
				fprintf(f, "\"card\": %f,", eh->edges[e].cardin);
				fprintf(f, "\"type\": \"edgeh\",");
				fprintf(f, "}},\n");
				e++;
			}

			// vertical edge
			if (y < eh->yqtd) {
				fprintf(f, "{\"type\": \"Feature\", \"geometry\": {\"type\": \"LineString\", \"coordinates\": [[%f, %f], [%f, %f]]},", 
					eh->edges[e].mbr.MinX, eh->edges[e].mbr.MinY, eh->edges[e].mbr.MaxX, eh->edges[e].mbr.MaxY);
				fprintf(f, "'properties': {");
				fprintf(f, "\"name\": \"ev:%d.%d\",", x, y);
				fprintf(f, "\"card\": %f,", eh->edges[e].cardin);
				fprintf(f, "\"type\": \"edgev\",");
				fprintf(f, "}},\n");
				e++;
			}

		}
	}

	fprintf(f, "]}\n");
	fclose(f);
}


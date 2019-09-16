
#include "histogram.h"
#include "ehist.h"

#define GET_VERT_EDGE(x, y) ((x == eh->xqtd) ? (x * (2*eh->yqtd+1) + y) : (x * (2*eh->yqtd+1) + 2*y + 1))
#define GET_HORZ_EDGE(x, y) (x * (2*eh->yqtd+1) + 2*y)

#define GET_VERT_EDGE_EHR(x, y) ((x == ehr->xqtd) ? (x * (2*ehr->yqtd+1) + y) : (x * (2*ehr->yqtd+1) + 2*y + 1))
#define GET_HORZ_EDGE_EHR(x, y) (x * (2*ehr->yqtd+1) + 2*y)
#define GET_VERT_EDGE_EHS(x, y) ((x == ehs->xqtd) ? (x * (2*ehs->yqtd+1) + y) : (x * (2*ehs->yqtd+1) + 2*y + 1))
#define GET_HORZ_EDGE_EHS(x, y) (x * (2*ehs->yqtd+1) + 2*y)

void eh_alloc(dataset *ds, euler_histogram *eh, int xqtd, int yqtd, double psizex, double psizey) {
    assert(xqtd > 0 && yqtd > 0 && "X and Y must be greater than zero.");
    eh->xqtd = xqtd;
    eh->yqtd = yqtd;
    eh->xtics = g_new(double, eh->xqtd+1);
    eh->ytics = g_new(double, eh->yqtd+1);
    eh->faces = g_new0(euler_face, eh->xqtd * eh->yqtd);
    eh->edges = g_new0(euler_edge, ((xqtd+1) * yqtd) + ((yqtd+1) * xqtd));
    eh->vertexes = g_new0(euler_vertex, (xqtd+1) * (yqtd+1));

    eh->xsize = psizex;
    eh->ysize = psizey;
    printf("Generating histogram of size: %d x %d\n", eh->xqtd, eh->yqtd);

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
    for(int i = 0; i <= eh->xqtd; i++) { // x
        for(int j = 0; j <= eh->yqtd; j++) { // y

            // set vetex at (i,j)
            eh->vertexes[v].x = eh->xtics[i];
            eh->vertexes[v].y = eh->ytics[j];
            v++;

            // horizontal edge at vertex v
            if (i < eh->xqtd) {
                int e = GET_HORZ_EDGE(i, j);
                eh->edges[e].mbr.MinX = eh->xtics[i];
                eh->edges[e].mbr.MaxX = eh->xtics[i+1];
                eh->edges[e].mbr.MinY = eh->ytics[j];
                eh->edges[e].mbr.MaxY = eh->ytics[j]+1e-10;
            }


            // vertical edge at vertex v
            if (j < eh->yqtd) {
                int e = GET_VERT_EDGE(i, j);
                eh->edges[e].mbr.MinY = eh->ytics[j];
                eh->edges[e].mbr.MaxY = eh->ytics[j+1];
                eh->edges[e].mbr.MinX = eh->xtics[i];
                eh->edges[e].mbr.MaxX = eh->xtics[i]+1e-10;
            }
        }
    }
}

void eh_hash_ds_objects(dataset *ds, euler_histogram *eh, enum JoinPredicateCheck pcheck) {

    dataset_iter di;
    dataset_foreach(di, ds) {
        dataset_leaf *l = get_join_pair_leaf(di.item, pcheck);
        Envelope ev = l->mbr; // pega mbr do objeto
        GEOSGeometryH geo = dataset_get_leaf_geo(ds, l);

        // descobrir a celula de que o objeto intercepta
        int xini = (ev.MinX - eh->mbr.MinX) / eh->xsize; // descobrir a celula do objeto
        int xfim = (ev.MaxX - eh->mbr.MinX) / eh->xsize;
        int yini = (ev.MinY - eh->mbr.MinY) / eh->ysize;
        int yfim = (ev.MaxY - eh->mbr.MinY) / eh->ysize;
        if (xfim < eh->xqtd) xfim++;
        if (yfim < eh->yqtd) yfim++;

    	int xspan = xfim - xini + 1;
    	int yspan = yfim - yini + 1;
        float areasum[xspan][yspan];

        for(int xi = 0; xi < xspan; xi++) {
        	for(int yi = 0; yi < yspan; yi++) {
        		areasum[xi][yi] = 0;
        	}
        }

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

                // face
                if (x < eh->xqtd && y < eh->yqtd) {
                    if (ENVELOPE_INTERSECTS(ev2, rs)) {
                        euler_face *face = &eh->faces[x*eh->yqtd +y];
                        face->cardin += 1;
                        double delta_x = (ev2.MaxX - ev2.MinX);
                        double delta_y = (ev2.MaxY - ev2.MinY);
                        double area = delta_x * delta_y;
                        face->avg_width += (delta_x - face->avg_width) / face->cardin;
                        face->avg_height += (delta_y - face->avg_height) / face->cardin;
                        face->avg_area += (area - face->avg_area) / face->cardin;

                        double aux_area = 0;
                        if (GEOSArea(clipped, &aux_area)) {
                        	areasum[x][y] += aux_area;
                        }

                        face->areasum += areasum[x][y];

                    }

                }

                GEOSGeom_destroy(clipped);

                // vertex
                int v = x * (eh->yqtd+1) + y;
                if (ENVELOPE_CONTAINSP(ev2, eh->vertexes[v].x, eh->vertexes[v].y))
                    eh->vertexes[v].cardin += 1;

                // horizontal edge
                if (x < eh->xqtd) {
                    int e = GET_HORZ_EDGE(x, y);
                    if (ENVELOPE_INTERSECTS(eh->edges[e].mbr, ev2)){
                        double delta_x = ev2.MaxX - ev2.MinX;
                        eh->edges[e].cardin += 1;
                        double edge_size = (eh->edges[e].mbr.MaxX - eh->edges[e].mbr.MinX);
                        eh->edges[e].avg_projection += (delta_x- eh->edges[e].avg_projection ) / eh->edges[e].cardin;
                    }
                    //eh->avg_projection +=
                }

                // vertical edge
                if (y < eh->yqtd) {
                    int e = GET_VERT_EDGE(x, y);
                    if (ENVELOPE_INTERSECTS(eh->edges[e].mbr, ev2)){
                        double delta_y = ev2.MaxY - ev2.MinY;
                        eh->edges[e].cardin += 1;
                        double edge_size = (eh->edges[e].mbr.MaxY - eh->edges[e].mbr.MinY);
                        eh->edges[e].avg_projection += (delta_y - eh->edges[e].avg_projection ) / eh->edges[e].cardin;
                    }
                }
            }
        }

        if (l->gid != -1) // free due to the call to dataset_get_leaf_geo
            GEOSGeom_destroy(geo);
    }
}

euler_histogram *eh_generate_hreal(dataset *ds,dataset *dsb,euler_histogram *eh,enum JoinPredicateCheck pcheck) {

    //gerando um novo histograma
	euler_histogram *eh_real = g_new(euler_histogram, 1);

	eh_real->mbr = ds->metadata.hist.mbr;

    //escolher hist A ou B para copia
	eh_real->xqtd = eh->xqtd;
	eh_real->yqtd = eh->yqtd;
	eh_real->xtics = g_new(double, eh->xqtd+1);
	eh_real->ytics = g_new(double, eh->yqtd+1);
	eh_real->faces = g_new0(euler_face, eh->xqtd * eh->yqtd);
	eh_real->edges = g_new0(euler_edge, ((eh->xqtd+1) * eh->yqtd) + ((eh->yqtd+1) * eh->xqtd));
	eh_real->vertexes = g_new0(euler_vertex, (eh->xqtd+1) * (eh->yqtd+1));

	eh_real->xsize = eh->xsize;
	eh_real->ysize = eh->ysize;
	printf("Generating histogram real of size: %d x %d\n", eh_real->xqtd, eh_real->yqtd);

	// X tics
	for(int i = 0; i < eh->xqtd; i++)
		eh_real->xtics[i] = eh->xtics[i];
	eh_real->xtics[eh->xqtd] = eh->xtics[eh->xqtd];

	// Y tics
	for(int i = 0; i < eh->yqtd; i++)
		eh_real->ytics[i] = eh->ytics[i];
	eh_real->ytics[eh->yqtd] = eh->ytics[eh->yqtd];

	//set edges and vertexes
	int v = 0;
	for(int i = 0; i <= eh_real->xqtd; i++) { // x
		for(int j = 0; j <= eh_real->yqtd; j++) { // y

			// set vetex at (i,j)
			eh_real->vertexes[v].x = eh_real->xtics[i];
			eh_real->vertexes[v].y = eh_real->ytics[j];
			v++;

			// horizontal edge at vertex v
			if (i < eh_real->xqtd) {
				int e = GET_HORZ_EDGE(i, j);
				eh_real->edges[e].mbr.MinX = eh_real->xtics[i];
				eh_real->edges[e].mbr.MaxX = eh_real->xtics[i+1];
				eh_real->edges[e].mbr.MinY = eh_real->ytics[j];
				eh_real->edges[e].mbr.MaxY = eh_real->ytics[j]+1e-10;
			}


			// vertical edge at vertex v
			if (j < eh_real->yqtd) {
				int e = GET_VERT_EDGE(i, j);
				eh_real->edges[e].mbr.MinY = eh_real->ytics[j];
				eh_real->edges[e].mbr.MaxY = eh_real->ytics[j+1];
				eh_real->edges[e].mbr.MinX = eh_real->xtics[i];
				eh_real->edges[e].mbr.MaxX = eh_real->xtics[i]+1e-10;
			}
		}
	}

    eh_hash_ds_objects(dsb, eh_real, pcheck);


    return eh_real;
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

    if (!ENVELOPE_INTERSECTS(query2, eh->mbr))
        return 0;

    double result = 0;
    Envelope query = EnvelopeIntersection2(query2, eh->mbr);

    int xini = (query.MinX - eh->mbr.MinX) / eh->xsize;
    int xfim = (query.MaxX - eh->mbr.MinX) / eh->xsize;
    int yini = (query.MinY - eh->mbr.MinY) / eh->ysize;
    int yfim = (query.MaxY - eh->mbr.MinY) / eh->ysize;
    if (xfim < eh->xqtd) xfim++;
    if (yfim < eh->yqtd) yfim++;

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
            if (ENVELOPE_CONTAINSP(query, eh->vertexes[v].x, eh->vertexes[v].y)){
                result += eh->vertexes[v].cardin;
            }

            // horizontal edge
            if (x < eh->xqtd) {
                int e = GET_HORZ_EDGE(x, y);
                if (ENVELOPE_INTERSECTS(eh->edges[e].mbr, query)){
                    if(eh->edges[e].mbr.MinY != query.MinY && eh->edges[e+1].mbr.MinY != query.MaxY) {
                        Envelope inters = EnvelopeIntersection(query, eh->edges[e].mbr);
                        double int_length = inters.MaxX - inters.MinX;
                        double fraction = int_length / (eh->edges[e].mbr.MaxX - eh->edges[e].mbr.MinX);
                        result -= fraction * eh->edges[e].cardin;
                    }
                }

            }

            // vertical edge
            if (y < eh->yqtd) {
                int e = GET_VERT_EDGE(x, y);
                if (ENVELOPE_INTERSECTS(eh->edges[e].mbr, query)){
                    if (eh->edges[e].mbr.MinX != query.MinX && eh->edges[e+1].mbr.MinX != query.MaxX) {
                        Envelope inters = EnvelopeIntersection(query, eh->edges[e].mbr);
                        double int_length = inters.MaxY - inters.MinY;
                        double fraction = int_length / (eh->edges[e].mbr.MaxY - eh->edges[e].mbr.MinY);
                        result -= fraction * eh->edges[e].cardin;
                    }
                }
            }
        }
    }

    return round(result);
}

void euler_print_hist(char *name, euler_histogram *eh) {
    char filename[100];

    char *prefix = getenv("HISTOPREFIX");
    prefix = prefix != NULL ? prefix : "";

    sprintf(filename, "histogram/%seuler-%s.geojson", prefix, name);
    FILE *f = fopen(filename, "wb");
    if (f == NULL) {
        perror("Error printing histogram");
        return;
    }

    fprintf(f, "{'type': 'FeatureCollection', 'features': [\n");

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
                fprintf(f, "[%.15lf, %.15lf],", env.MinX, env.MinY);
                fprintf(f, "[%.15lf, %.15lf],", env.MaxX, env.MinY);
                fprintf(f, "[%.15lf, %.15lf],", env.MaxX, env.MaxY);
                fprintf(f, "[%.15lf, %.15lf],", env.MinX, env.MaxY);
                fprintf(f, "[%.15lf, %.15lf]",  env.MinX, env.MinY);
                fprintf(f, "]]}, 'properties': {");
                fprintf(f, "\"name\": \"f:%d.%d\",", x, y);
                fprintf(f, "\"card\": %lf,", eh->faces[x*eh->yqtd + y].cardin);
                fprintf(f, "\"avg_heigth\": %lf,", eh->faces[x*eh->yqtd + y].avg_height);
                fprintf(f, "\"avg_width\": %lf,", eh->faces[x*eh->yqtd + y].avg_width);
                fprintf(f, "\"avg_area\": %lf,", eh->faces[x*eh->yqtd + y].avg_area);
                fprintf(f, "\"face_area\": %lf,", eh->xtics[0]*eh->ytics[0] );
                fprintf(f, "\"type\": \"face\",");
                fprintf(f, "}},\n");
            }

            // vertex
            fprintf(f, "{\"type\": \"Feature\", \"geometry\": {\"type\": \"Point\", \"coordinates\": [%.15lf, %.15lf]},",
                    eh->vertexes[v].x, eh->vertexes[v].y);
            fprintf(f, "'properties': {");
            fprintf(f, "\"name\": \"v:%d.%d\",", x, y);
            fprintf(f, "\"card\": %lf,", eh->vertexes[v].cardin);
            fprintf(f, "\"type\": \"vertex\",");
            fprintf(f, "}},\n");
            v++;


            // horizontal edge
            if (x < eh->xqtd) {
                int e = GET_HORZ_EDGE(x, y);
                fprintf(f, "{\"type\": \"Feature\", \"geometry\": {\"type\": \"LineString\", \"coordinates\": [[%.15lf, %.15lf], [%.15lf, %.15lf]]},",
                        eh->edges[e].mbr.MinX, eh->edges[e].mbr.MinY, eh->edges[e].mbr.MaxX, eh->edges[e].mbr.MaxY);
                fprintf(f, "'properties': {");
                fprintf(f, "\"name\": \"eh:%d.%d\",", x, y);
                fprintf(f, "\"card\": %lf,", eh->edges[e].cardin);
                fprintf(f, "\"type\": \"edgeh\",");
                fprintf(f, "}},\n");
            }

            // vertical edge
            if (y < eh->yqtd) {
                int e = GET_VERT_EDGE(x, y);
                fprintf(f, "{\"type\": \"Feature\", \"geometry\": {\"type\": \"LineString\", \"coordinates\": [[%.15lf, %.15lf], [%.15lf, %.15lf]]},",
                        eh->edges[e].mbr.MinX, eh->edges[e].mbr.MinY, eh->edges[e].mbr.MaxX, eh->edges[e].mbr.MaxY);
                fprintf(f, "'properties': {");
                fprintf(f, "\"name\": \"ev:%d.%d\",", x, y);
                fprintf(f, "\"card\": %lf,", eh->edges[e].cardin);
                fprintf(f, "\"type\": \"edgev\",");
                fprintf(f, "}},\n");
            }

        }
    }

    fprintf(f, "]}\n");
    fclose(f);
}


double estimate_intersections_mp_edges_vert(Envelope el, Envelope er, Envelope inters,
        euler_edge *ehr_face, euler_edge *ehs_face) {
    // the code below follows equations (1) and (2) in Mamoulis, Papadias 2001

    // estimate the quantity of objects in LeftDs in the inters window: eqn (1)
    double uy = el.MaxY - el.MinY;

    double avgl_y = ehr_face->avg_projection;
    double wy = inters.MaxY - inters.MinY;
    double qtdobjl = ehr_face->cardin *
        MIN(1,(avgl_y+wy)/uy);

    // estimate the quantity of objects in RightDs in the inters window eqn (1)
    uy = er.MaxY - er.MinY;
    double avgr_y = ehs_face->avg_projection;
    double qtdobjr = ehs_face->cardin *
        MIN(1,(avgr_y+wy)/uy);

    // estimate join result cardinality, eqn (2)
    return qtdobjl * qtdobjr * MIN(1, (avgl_y + avgr_y)/wy);
}

double estimate_intersections_mp_edges_horz(Envelope el, Envelope er, Envelope inters,
        euler_edge *ehr_face, euler_edge *ehs_face) {
    // the code below follows equations (1) and (2) in Mamoulis, Papadias 2001

    // estimate the quantity of objects in LeftDs in the inters window: eqn (1)
    double ux = el.MaxX - el.MinX;

    double avgl_x = ehr_face->avg_projection;
    double wx = inters.MaxX - inters.MinX;
    double qtdobjl = ehr_face->cardin *
        MIN(1,(avgl_x+wx)/ux);

    // estimate the quantity of objects in RightDs in the inters window eqn (1)
    ux = er.MaxX - er.MinX;
    double avgr_x = ehs_face->avg_projection;
    double qtdobjr = ehs_face->cardin *
        MIN(1,(avgr_x+wx)/ux);

    // estimate join result cardinality, eqn (2)
    return qtdobjl * qtdobjr * MIN(1, (avgl_x + avgr_x)/wx);
}

double estimate_intersections_mamoulis_papadias(Envelope el, Envelope er, Envelope inters,
        euler_face *ehr_face, euler_face *ehs_face) {
    // the code below follows equations (1) and (2) in Mamoulis, Papadias 2001

    // estimate the quantity of objects in LeftDs in the inters window: eqn (1)
    double ux = el.MaxX - el.MinX;
    double uy = el.MaxY - el.MinY;
    double avgl_x = ehr_face->avg_width;
    double avgl_y = ehr_face->avg_height;
    double wx = inters.MaxX - inters.MinX;
    double wy = inters.MaxY - inters.MinY;

    if(avgl_x == 0.0)
    	printf("media x zerada \n");
    if(avgl_y == 0.0)
    	printf("media y zerada \n");
    if(ehr_face->cardin == 0.0)
    	printf("face zerada \n");

    double qtdobjl = ehr_face->cardin *
        MIN(1,(avgl_x+wx)/ux) *
        MIN(1,(avgl_y+wy)/uy);
    	assert(MIN(1,(avgl_x+wx)/ux) * MIN(1,(avgl_y+wy)/uy) <= 1.0 && "fraction should not be higher than 1.0");

    // estimate the quantity of objects in RightDs in the inters window eqn (1)
    ux = er.MaxX - er.MinX;
    uy = er.MaxY - er.MinY;
    double avgr_x = ehs_face->avg_width;
    double avgr_y = ehs_face->avg_height;
    //printf("ux = %lf \t uy = %lf\t avgl_x = %lf\t avgl_y = %lf\t wx = %lf\t wy = %lf\t", ux,uy,avgl_x,avgl_y,wx,wy);

    double qtdobjr = ehs_face->cardin *
        MIN(1,(avgr_x+wx)/ux) *
        MIN(1,(avgr_y+wy)/uy);
    	assert(MIN(1,(avgr_x+wx)/ux) *  MIN(1,(avgr_y+wy)/uy) <= 1.0 && "fraction should not be higher than 1.0");

    double minX = MIN(1, (avgl_x + avgr_x)/wx);
    double minY = MIN(1, (avgl_y + avgr_y)/wy);


    // estimate join result cardinality, eqn (2)
    double res = qtdobjl * qtdobjr * minX * minY;

    return res;
}


int euler_join_cardinality(dataset *dr,
        dataset *ds,
        euler_histogram* ehr,
        euler_histogram* ehs,
        rtree_root* rtree_r,
        rtree_root* rtree_s,
        double* stddev) {


    double xini = MAX(ehr->xtics[0], ehs->xtics[0]);
    double yini = MAX(ehr->ytics[0], ehs->ytics[0]);
    double xend = MIN(dr->metadata.hist.mbr.MaxX, ds->metadata.hist.mbr.MaxX);
    double yend = MIN(dr->metadata.hist.mbr.MaxY, ds->metadata.hist.mbr.MaxY);

    unsigned int N = ceil((xend - xini) * (yend - yini));
	double mean = 0.0;
	double M2 = 0.0;
	double sum_error = 0.0;

    int xdr_start = 0;
    while (xdr_start < ehr->xqtd && ehr->xtics[xdr_start+1] < xini)
        xdr_start++;
    int xdr_end = xdr_start+1;
    while (xdr_end < ehr->xqtd && ehr->xtics[xdr_end] <= xend)
        xdr_end++;
    if (xdr_start == ehr->xqtd)
        return 0;

    // skip non-intersect area on y
    int ydr_start = 0;
    while (ydr_start < ehr->yqtd && ehr->ytics[ydr_start+1] < yini)
        ydr_start++;
    int ydr_end = ydr_start+1;
    while (ydr_end < ehr->yqtd && ehr->ytics[ydr_end] <= yend)
        ydr_end++;
    if (ydr_start == ehr->yqtd)
        return 0;

    int xds_atu = 0;
    float result = 0;
    for(int xr = xdr_start; xr < xdr_end; xr++) {

        while(xds_atu < ehs->xqtd && ehs->xtics[xds_atu+1] < ehr->xtics[xr]) // skip when end of s < start of r
            xds_atu++;
        int xds_end = xds_atu+1;
        while(xds_end < ehs->xqtd && ehs->xtics[xds_end] <= ehr->xtics[xr+1]) // increment when end of s < start of r
            xds_end++;

        int yds_atu = 0;

        Envelope er;
        Envelope es;
        er.MinX = ehr->xtics[xr];
        er.MaxX = ehr->xtics[xr+1];

        for(int yr = ydr_start; yr < ydr_end; yr++) {

            while(yds_atu < ehs->yqtd && ehs->ytics[yds_atu+1] < ehr->ytics[yr]) // skip when end of s < start of r
                yds_atu++;
            int yds_end = yds_atu+1;
            while(yds_end < ehs->yqtd && ehs->ytics[yds_end] <= ehr->ytics[yr+1]) // increment when end of s < start of r
                yds_end++;

            er.MinY = ehr->ytics[yr];
            er.MaxY = ehr->ytics[yr+1];
            double erarea = ENVELOPE_AREA(er);

            for(int xs = xds_atu; xs < xds_end; xs++) {
                es.MinX = ehs->xtics[xs];
                es.MaxX = ehs->xtics[xs+1];

                for(int ys = yds_atu; ys < yds_end; ys++) {
                    double estimated_result = 0;

                    es.MinY = ehs->ytics[ys];
                    es.MaxY = ehs->ytics[ys+1];

                    char i = ENVELOPE_INTERSECTS(er, es);
                    assert(i != 0);

                    euler_face *ehr_face = &ehr->faces[xr*ehr->yqtd +yr];
                    euler_face *ehs_face = &ehs->faces[xs*ehs->yqtd +ys];
                    Envelope inters = EnvelopeIntersection2(er, es);
                    double int_area = ENVELOPE_AREA(inters);

                    printf("ehr_face[%d][%d].card = %f, avg_h = %f, avg_w = %f\n"
                            "ehs_face[%d][%d].card = %f, avg_h = %f, avg_w = %f\n", xr, yr, ehr_face->cardin,ehr_face->avg_height, ehr_face->avg_width, xs, ys, ehs_face->cardin, ehs_face->avg_height, ehs_face->avg_width);

                    double intersections = 0;
                    double p = 1;
                    if(ehr_face->cardin >= 1 || ehs_face->cardin >= 1){
                        intersections = estimate_intersections_mamoulis_papadias(er, es, inters, ehr_face, ehs_face);

                        if(intersections< 1.0)
                            intersections = 0;
                    }
                    result +=  intersections;
                    estimated_result += intersections;


                    //vertice
                    int vr = xr * (ehr->yqtd+1) + yr;
                    int vs = xs * (ehs->yqtd+1) + ys;

                    if(ehs->vertexes[vs].x == ehr->vertexes[vr].x && ehs->vertexes[vs].y == ehr->vertexes[vr].y ){
                        result += ehr->vertexes[vr].cardin * ehs->vertexes[vs].cardin;
                        estimated_result += ehr->vertexes[vr].cardin * ehs->vertexes[vs].cardin;
                    }

                    //aresta horizontal
                    int ar = GET_HORZ_EDGE_EHR(xr, yr);
                    int as = GET_HORZ_EDGE_EHS(xs, ys);
                    if (ENVELOPE_INTERSECTS(ehr->edges[ar].mbr, ehs->edges[as].mbr)){
                        printf("ar = %d, as = %d\n", ar, as);
                        Envelope inters = EnvelopeIntersection2(ehr->edges[ar].mbr, ehs->edges[as].mbr);
                        double int_length = inters.MaxX - inters.MinX;
                        double fraction_ar = int_length / (ehr->edges[ar].mbr.MaxX - ehr->edges[ar].mbr.MinX);
                        double fraction_as = int_length / (ehs->edges[as].mbr.MaxX - ehs->edges[as].mbr.MinX);

                        double p = MIN(1, ehr->edges[ar].avg_projection + ehs->edges[as].avg_projection);

                        double cardin_ar = ehr->edges[ar].cardin;
                        double cardin_as = ehr->edges[as].cardin;
                        result -=  cardin_ar * cardin_as * p;
                        estimated_result -=  cardin_ar * cardin_as * p;
                        //result -= estimate_intersections_mp_edges_horz(ehr->edges[ar].mbr, ehs->edges[as].mbr, inters, &ehr->edges[ar], &ehs->edges[as]);
                    }

                    ar = GET_VERT_EDGE_EHR(xr, yr);
                    as = GET_VERT_EDGE_EHS(xs, ys);
                    if (ENVELOPE_INTERSECTS(ehr->edges[ar].mbr, ehs->edges[as].mbr)){
                        Envelope inters = EnvelopeIntersection2(ehr->edges[ar].mbr, ehs->edges[as].mbr);
                        double int_length = inters.MaxY - inters.MinY;
                        double fraction_ar = int_length / (ehr->edges[ar].mbr.MaxY - ehr->edges[ar].mbr.MinY);
                        double fraction_as = int_length / (ehs->edges[as].mbr.MaxY - ehs->edges[as].mbr.MinY);

                        double p = MIN(1, ehr->edges[ar].avg_projection + ehs->edges[as].avg_projection);
                        printf(" edge p = %f\n", ehr->edges[ar].avg_projection + ehs->edges[as].avg_projection);

                        double cardin_ar = ehr->edges[ar].cardin;
                        double cardin_as = ehr->edges[as].cardin;
                        result -=  cardin_ar * cardin_as * p;
                        estimated_result -=  cardin_ar * cardin_as * p;
                        //result -= estimate_intersections_mp_edges_vert(ehr->edges[ar].mbr, ehs->edges[as].mbr, inters, &ehr->edges[ar], &ehs->edges[as]);
                    }
                        //double real_cardin = real_cardin_euler_histogram_cell(rtree_r, rtree_s, inters);
                        double real_cardin = 1;
                        int error = abs(estimated_result - real_cardin);
                        double delta = error - mean;
                        assert(!isnan(delta));

		                mean += delta/(double)N;
                        assert(!isnan(mean));
                        M2 += delta*(error - mean);
                        sum_error += error;


                }
            }
        }
    }


    *stddev = sqrt(M2/(double)N);
    assert(!isnan(*stddev));
    return round(result);
}


int euler_spatial_join(euler_histogram* ehr, euler_histogram* ehs){
    if(!ENVELOPE_INTERSECTS(ehr->mbr, ehs->mbr))
        return 0;

    double result = 0;

    double xini = MAX(ehr->xtics[0], ehs->xtics[0]);
    double yini = MAX(ehr->ytics[0], ehs->ytics[0]);
    double xfim = MIN(ehr->mbr.MaxX, ehs->mbr.MaxX);
    double yfim = MIN(ehr->mbr.MaxY, ehs->mbr.MaxY);

    for(int xr = xini; xr <= xfim; xr++){
        Envelope er, es;
        er.MinX = ehr->xtics[xr];
        er.MaxX = ehr->xtics[xr+1];

        for(int yr = yini; yr <= yfim; yr++){
            er.MinY = ehr->xtics[yr];
            er.MaxY = ehr->xtics[yr+1];

            for(int xs = xini; xs <= xfim; xs++){
                es.MinX = ehs->xtics[xs];
                es.MaxX = ehs->xtics[xs+1];

                for(int ys = yini; ys <= yfim; ys++){
                    es.MinY = ehs->xtics[ys];
                    es.MaxY = ehs->xtics[ys+1];
                    if(ENVELOPE_INTERSECTS(er, es)){
                        //face
                        euler_face *ehr_face = &ehr->faces[xr*ehs->yqtd +yr];
                        euler_face *ehs_face = &ehr->faces[xs*ehs->yqtd +ys];

                        Envelope inters = EnvelopeIntersection(er, es);
                        double int_area = ENVELOPE_AREA(inters);
                        double face_area = ENVELOPE_AREA(es);
                        double fraction = int_area / face_area;
                        //result += fraction * face->cardin;
                        if(ehr_face->avg_height + ehs_face->avg_height >= 1 && ehr_face->avg_width + ehs_face->avg_width >= 1)
                            result += fraction * ehs_face->cardin;
                        else{
                            double p =  ehr_face->avg_area + ehs_face->avg_area + ehr_face->avg_height * ehs_face->avg_width+ehs_face->avg_height * ehr_face->avg_width;

                            result += fraction * ehs_face->cardin * p;
                        }

                        //vertice
                        int vr = xr * (ehr->yqtd+1) + yr;
                        int vs = xs * (ehs->yqtd+1) + ys;
                        if (ENVELOPE_CONTAINSP(es, ehr->vertexes[vr].x, ehr->vertexes[vr].y)){
                            result += ehr->vertexes[vr].cardin * ehs->vertexes[vs].cardin;
                        }

                        //aresta horizontal
                        int ar = GET_HORZ_EDGE_EHR(xr, yr);
                        int as = GET_HORZ_EDGE_EHS(xs, ys);
                        if (ENVELOPE_INTERSECTS(ehs->edges[as].mbr, er)){
                            if(ehs->edges[as].mbr.MinY != es.MinY && ehs->edges[as+1].mbr.MinY != es.MaxY) {
                                Envelope inters = EnvelopeIntersection(es, ehs->edges[as].mbr);
                                double int_length = inters.MaxX - inters.MinX;
                                double fraction_ar = int_length / (ehr->edges[ar].mbr.MaxX - ehr->edges[ar].mbr.MinX);
                                double fraction_as = int_length / (ehs->edges[as].mbr.MaxX - ehs->edges[as].mbr.MinX);
                                double p = MIN(1, fraction_ar * ehs->edges[ar].cardin + fraction_as * ehs->edges[as].cardin);
                                result -=  ehr->edges[ar].cardin * ehs->edges[as].cardin * p;
                            }
                        }

                        ar = GET_VERT_EDGE_EHR(xr, yr);
                        as = GET_VERT_EDGE_EHS(xs, ys);
                        if (ENVELOPE_INTERSECTS(ehs->edges[as].mbr, er)){
                            if (ehs->edges[as].mbr.MinX != es.MinX && ehs->edges[as+1].mbr.MinX != es.MaxX) {
                                Envelope inters = EnvelopeIntersection(es, ehs->edges[as].mbr);
                                double int_length = inters.MaxY - inters.MinY;
                                double fraction_ar = int_length / (ehr->edges[ar].mbr.MaxX - ehr->edges[ar].mbr.MinX);
                                double fraction_as = int_length / (ehs->edges[as].mbr.MaxX - ehs->edges[as].mbr.MinX);
                                double p = MIN(1, fraction_ar * ehs->edges[ar].cardin + fraction_as * ehs->edges[as].cardin);
                                double fraction = int_length / (ehs->edges[as].mbr.MaxY - ehr->edges[as].mbr.MinY);
                                result -=  ehr->edges[ar].cardin * ehs->edges[as].cardin * p;
                            }
                        }
                    }
                }
            }
        }
    }


    return round(result);

}


double estimate_intersections_zu(Envelope el, Envelope er,
Envelope inters, euler_face *ehl_face, euler_face *ehr_face,
dataset *dh_l, dataset *dh_r) {


		//if (!ENVELOPE_INTERSECTS(lcell->usedarea, rcell->usedarea))
		if(!ENVELOPE_INTERSECTS(el, er))
			return 0;

		//Envelope w = EnvelopeIntersection2(el, er);
		double wx = inters.MaxX - inters.MinX;
		double wy = inters.MaxY - inters.MinY;


		 //printf("%c",dh_r->metadata.memory_dataset);


		bool line_to_line = false;
		bool line_to_polygon = true;
		bool left_is_line = false;

		/*if ((dh_l->geom_type == wkbLineString || dh_l->geom_type == wkbMultiLineString) &&
				(dh_r->geom_type == wkbPolygon    || dh_r->geom_type == wkbMultiPolygon)) {
				line_to_polygon = true;
				left_is_line = true;
		}
		else if ((dh_r->geom_type == wkbLineString || dh_r->geom_type == wkbMultiLineString) &&
				(dh_l->geom_type == wkbPolygon    || dh_l->geom_type == wkbMultiPolygon)) {
				line_to_polygon = true;
				left_is_line = false;
		}
		else if ((dh_l->geom_type == wkbLineString || dh_l->geom_type == wkbMultiLineString) &&
		(dh_r->geom_type == wkbLineString || dh_r->geom_type == wkbMultiLineString)) {
			line_to_line = true;
			left_is_line = true;
		}*/

		// estimate the quantity of objects in LeftDs inside inters window
		double ux_l = el.MaxX - el.MinX;
		double uy_l = el.MaxY - el.MinY;
		double avgl_x = ehl_face->avg_width;
		double avgl_y = ehl_face->avg_height;

		// observing that objects generally doesn't overlap in both axis,
		// fix the probability of intersection in one of them
		double usx_l = el.MaxX - el.MinX;
		double usy_l = el.MaxY - el.MinY;

/*
		if(avgl_x == 0)
			printf("avgl_x zerada antes\n");
		if(avgl_y == 0)
			printf("avgl_y zerada antes\n");*/

		if (avgl_x > avgl_y)
			avgl_x = MIN(usx_l/(avgl_y*ehl_face->cardin/usy_l), avgl_x);
		else
			avgl_y = MIN(usy_l/(avgl_x*ehl_face->cardin/usx_l), avgl_y);
		/*
		if(avgl_x == 0)
			printf("avgl_x zerada depois\n");
		if(avgl_y == 0)
			printf("avgl_y zerada depois\n");
		*/
		double qtdobjl = MAX(0.0, ehl_face->cardin) *
				MIN(1.0,(avgl_x + wx)/usx_l) *
				MIN(1.0,(avgl_y + wy)/usy_l);

		// estimate the quantity of objects in RightDs inside inters window
		double ux_r = er.MaxX - er.MinX;
		double uy_r = er.MaxY - er.MinY;
		double avgr_x = ehr_face->avg_width;
		double avgr_y = ehr_face->avg_height;

		// observing that objects generally doesn't overlap in both axis,
		// fix the probability of intersection in one of them
		double usx_r = er.MaxX - er.MinX;
		double usy_r = er.MaxY - er.MinY;

		//printf("avgr_y: %lf\n",avgr_y);
		//printf("avgr_x: %lf\n",avgr_x);

		if(avgr_x == 0)
			printf("avgr_x zerada antes\n");
		if(avgr_y == 0)
			printf("avgr_y zerada antes\n");

		if (avgr_x > avgr_y)
			avgr_x = MIN(usx_r/(avgr_y*ehr_face->cardin/usy_r), avgr_x);
		else
			avgr_y = MIN(usy_r/(avgr_x*ehr_face->cardin/usx_r), avgr_y);

		/*if(avgr_x == 0)
			printf("avgr_x zerada depois\n");
		if(avgr_y == 0)
			printf("avgr_y zerada depois\n\n");*/



		double qtdobjr = MAX(0.0, ehr_face->cardin) *
				MIN(1.0,(avgr_x + wx)/usx_r) *
				MIN(1.0,(avgr_y + wy)/usy_r);

		// estimate join result cardinality
		double result = qtdobjl * qtdobjr *
				MIN(1.0, (avgl_x + avgr_x)/wx) *
				MIN(1.0, (avgl_y + avgr_y)/wy);

		double coef_area = 0;
		if (line_to_polygon) {
			/* disabling this increases the estimation a litle bit for J1,J2
			avgr_x = MIN(wx,avgr_x);
			avgr_y = MIN(wy,avgr_y);
			avgl_x = MIN(wx,avgl_x);
			avgl_y = MIN(wy,avgl_y);*/

			if (left_is_line) {
					double d = ehr_face->areasum/(ux_r*uy_r);
					double f = sqrt(((ux_r*uy_r)/ehr_face->cardin)/(avgr_x*avgr_y));
					//if(ehr_face->cardin == 0)
						//printf("cardinalidade zerada \n");
					//if(avgr_x*avgr_y == 0)
						//printf("avgr_x*avgr_y zerada \n");
					double navgx = avgr_x * f;
					double navgy = avgr_y * f;
					result = qtdobjl * d * MAX(1.0,avgl_x / navgx) * MAX(1.0,avgl_y / navgy);
					//printf("%lf\t",result);
			} else {
				double d = ehl_face->areasum/(ux_l*uy_l);
				double f = sqrt(((ux_l*uy_l)/ehl_face->cardin)/(avgl_x*avgl_y));
				//if(ehl_face->cardin == 0)
					//printf("cardinalidade zerada \n");
				//if(avgl_x*avgl_y == 0)
					//printf("avgl_x*avgl_y zerada \n");
				double navgx = avgl_x * f;
				double navgy = avgl_y * f;
				result = qtdobjr * d * MAX(1.0,avgr_x / navgx) * MAX(1.0,avgr_y / navgy);
			}
		}
		else if (line_to_line) {
			/*
			The probability of two random line segments intersect in unit square = 1/3
			The probability of four random points forms a convex quadrilateral 133/144
			The overall probability then, 1/3 * 133/144 = 133/432 ~= 0.3078
			Source: math.stackexchange.com/questions/134525
			*/

			double line_coef;
			double margina = sqrt(pow(avgl_x,2) + pow(avgl_y,2));
			double marginb = sqrt(pow(avgr_x,2) + pow(avgr_y,2));
			if (margina < marginb)
				line_coef = MIN(133.0/432.0, margina/marginb);
			else
				line_coef = MIN(133.0/432.0, marginb/margina);
			result = qtdobjl * qtdobjr * line_coef *
					MIN(1.0, (avgl_x + avgr_x)/wx) *
					MIN(1.0, (avgl_y + avgr_y)/wy);
		}

		return result;
}

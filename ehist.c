
#include "histogram.h"
#include "ehist.h"


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
                        double delta_x = ev2.MaxX - ev2.MinX;
                        double delta_y = ev2.MaxY - ev2.MinY;
                        double area = delta_x * delta_y;
                        face->avg_width += (delta_x - face->avg_width) / face->cardin;
                        face->avg_heigth += (delta_y - face->avg_heigth) / face->cardin;
                        face->avg_area += (area - face->avg_area) / face->cardin;
                    }

                }

                // vertex
                int v = x * (eh->yqtd+1) + y;
                if (ENVELOPE_CONTAINSP(ev2, eh->vertexes[v].x, eh->vertexes[v].y))
                    eh->vertexes[v].cardin += 1;

                // horizontal edge
                if (x < eh->xqtd) {
                    int e = GET_HORZ_EDGE(x, y);
                    if (ENVELOPE_INTERSECTS(eh->edges[e].mbr, ev2))
                        eh->edges[e].cardin += 1;
                    //eh->avg_projection += 
                }

                // vertical edge
                if (y < eh->yqtd) {
                    int e = GET_VERT_EDGE(x, y);
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

int euler_join_cardinality(dataset *dr, dataset *ds, euler_histogram* ehr, euler_histogram* ehs) {
    double xini = MAX(ehr->xtics[0], ehs->xtics[0]);
    double yini = MAX(ehr->ytics[0], ehs->ytics[0]);
    double xend = MIN(dr->metadata.hist.mbr.MaxX, ds->metadata.hist.mbr.MaxX);
    double yend = MIN(dr->metadata.hist.mbr.MaxY, ds->metadata.hist.mbr.MaxY);

    // skip non-intersect area on x
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

                    es.MinY = ehs->ytics[ys];
                    es.MaxY = ehs->ytics[ys+1];


                    if(ENVELOPE_INTERSECTS(er, es)){
                        //face
                        euler_face *ehr_face = &ehr->faces[xr*ehs->yqtd +yr];
                        euler_face *ehs_face = &ehr->faces[xs*ehs->yqtd +ys];

                        Envelope inters = EnvelopeIntersection2(er, es);
                        double int_area = ENVELOPE_AREA(inters);
                        //double face_area = ENVELOPE_AREA(es);
                        double ehrfraction = int_area / erarea;
                        double ehsfraction = int_area / ENVELOPE_AREA(es);


					    double qtdobjr = ehrfraction * ehr_face->cardin;
					    double qtdobjs = ehsfraction * ehs_face->cardin;

                       printf("ehr_face->avg_heigth + ehs_face->avg_heigth = %f \n", ehr_face->avg_heigth + ehs_face->avg_heigth) ;
                        if(ehr_face->avg_heigth + ehs_face->avg_heigth >= 1 && ehr_face->avg_width + ehs_face->avg_width >= 1){
                            printf("asdfsadf\n");
                            result += qtdobjr * qtdobjs;
                        }
                        
                        else{
                            double p =  ehr_face->avg_area + ehs_face->avg_area + ehr_face->avg_heigth * ehs_face->avg_width+ehs_face->avg_heigth * ehr_face->avg_width;
                            printf("face p = %f\n", p);

                            result +=  (qtdobjr * qtdobjs); 
                        }
                    }

                        /*//vertice
                        int vr = xr * (ehr->yqtd+1) + yr;	
                        int vs = xs * (ehs->yqtd+1) + ys;	
                        if (ENVELOPE_CONTAINSP(er, ehs->vertexes[vs].x, ehs->vertexes[vs].y)){
                            result += ehr->vertexes[vr].cardin * ehs->vertexes[vs].cardin;
                        }

                        //aresta horizontal
                        int ar = GET_HORZ_EDGE_EHR(xr, yr);
                        int as = GET_HORZ_EDGE_EHS(xs, ys);
                        if (ENVELOPE_INTERSECTS(er, ehs->edges[as].mbr)){
                                Envelope inters = EnvelopeIntersection(er, ehs->edges[as].mbr);
                                double int_length = inters.MaxX - inters.MinX;
                                double fraction_ar = int_length / (ehr->edges[ar].mbr.MaxX - ehr->edges[ar].mbr.MinX);
                                double fraction_as = int_length / (ehs->edges[as].mbr.MaxX - ehs->edges[as].mbr.MinX);
                                double p = MIN(1, fraction_ar * ehr->edges[ar].cardin + fraction_as * ehs->edges[as].cardin);
                                //printf("aresta p = %f\n", p);
                                result -=  ehs->edges[ar].cardin * ehr->edges[ar].cardin * p;
                        }

                        ar = GET_VERT_EDGE_EHR(xr, yr);
                        as = GET_VERT_EDGE_EHS(xs, ys);
                        if (ENVELOPE_INTERSECTS(ehs->edges[as].mbr, er)){
                                Envelope inters = EnvelopeIntersection(er, ehs->edges[as].mbr);
                                double int_length = inters.MaxY - inters.MinY;
                                double fraction_ar = int_length / (ehr->edges[ar].mbr.MaxX - ehr->edges[ar].mbr.MinX);
                                double fraction_as = int_length / (ehs->edges[as].mbr.MaxX - ehs->edges[as].mbr.MinX);
                                double p = MIN(1, fraction_ar * ehr->edges[ar].cardin + fraction_as * ehs->edges[as].cardin);
                                double fraction = int_length / (ehs->edges[as].mbr.MaxY - ehr->edges[ar].mbr.MinY);
                                result -= ehr->edges[ar].cardin * ehs->edges[as].cardin * p;
                        }*/
                    
                }
            }
        }
    }




    return (int)result;
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
                        if(ehr_face->avg_heigth + ehs_face->avg_heigth >= 1 && ehr_face->avg_width + ehs_face->avg_width >= 1)
                            result += fraction * ehs_face->cardin;
                        else{
                            double p =  ehr_face->avg_area + ehs_face->avg_area + ehr_face->avg_heigth * ehs_face->avg_width+ehs_face->avg_heigth * ehr_face->avg_width;

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

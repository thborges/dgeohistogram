#include "histogram.h"
#include "ehist.h"

//#define MIN(a,b) (((a)<(b))?(a):(b))
#define GET_VERT_EDGE(x, y) ((x == ehA->xqtd) ? (x * (2*ehA->yqtd+1) + y) : (x * (2*ehA->yqtd+1) + 2*y + 1))
#define GET_HORZ_EDGE(x, y) (x * (2*ehA->yqtd+1) + 2*y)

#define GET_VERT_EDGE_EHR(x, y) ((x == ehA->xqtd) ? (x * (2*ehA->yqtd+1) + y) : (x * (2*ehA->yqtd+1) + 2*y + 1))
#define GET_HORZ_EDGE_EHR(x, y) (x * (2*ehA->yqtd+1) + 2*y)
#define GET_VERT_EDGE_EHS(x, y) ((x == ehB->xqtd) ? (x * (2*ehB->yqtd+1) + y) : (x * (2*ehB->yqtd+1) + 2*y + 1))
#define GET_HORZ_EDGE_EHS(x, y) (x * (2*ehB->yqtd+1) + 2*y)

euler_histogram *eh_generate_intermed(dataset *ds,dataset *dsb,euler_histogram *ehA,euler_histogram *ehB) {

    double intersection = 1;

    //arquivo txt
    FILE *valores;
    valores = fopen("../valores.txt","w");

    if(valores == NULL){
      printf("Erro na abertura do arquivo!");
      return 1;
    }

    fprintf(valores, "%s", "linha\tcoluna\tfaceA\tfaceIntermed");
    fprintf(valores, "%s", "\tarestaVertA\tarestaVertInter\tarestaHorzA\tarestaHorzIntermed\tvertA\tvertIntermed\n");


    //gerando um novo histograma
	euler_histogram *ehIntermed = g_new(euler_histogram, 1);

	ehIntermed->mbr = ds->metadata.hist.mbr;

    //escolher hist A ou B para copia
    ehIntermed->xqtd = ehA->xqtd;
	ehIntermed->yqtd = ehA->yqtd;
	ehIntermed->xtics = g_new(double, ehA->xqtd+1);
	ehIntermed->ytics = g_new(double, ehA->yqtd+1);
	ehIntermed->faces = g_new0(euler_face, ehA->xqtd * ehA->yqtd);
	ehIntermed->edges = g_new0(euler_edge, ((ehA->xqtd+1) * ehA->yqtd) + ((ehA->yqtd+1) * ehA->xqtd));
	ehIntermed->vertexes = g_new0(euler_vertex, (ehA->xqtd+1) * (ehA->yqtd+1));

	ehIntermed->xsize = ehA->xsize;
	ehIntermed->ysize = ehA->ysize;
	printf("Generating histogram Intermed of size: %d x %d\n", ehIntermed->xqtd, ehIntermed->yqtd);

	// X tics
	for(int i = 0; i < ehA->xqtd; i++)
		ehIntermed->xtics[i] = ehA->xtics[i];
	ehIntermed->xtics[ehA->xqtd] = ehA->xtics[ehA->xqtd];

	// Y tics
	for(int i = 0; i < ehA->yqtd; i++)
		ehIntermed->ytics[i] = ehA->ytics[i];
	ehIntermed->ytics[ehA->yqtd] = ehA->ytics[ehA->yqtd];

	//set edges and vertexes
	int v = 0;
	for(int i = 0; i <= ehIntermed->xqtd; i++) { // x
		for(int j = 0; j <= ehIntermed->yqtd; j++) { // y

			// set vetex at (i,j)
			ehIntermed->vertexes[v].x = ehIntermed->xtics[i];
			ehIntermed->vertexes[v].y = ehIntermed->ytics[j];
			v++;

			// horizontal edge at vertex v
			if (i < ehIntermed->xqtd) {
				int e = GET_HORZ_EDGE(i, j);
				ehIntermed->edges[e].mbr.MinX = ehIntermed->xtics[i];
				ehIntermed->edges[e].mbr.MaxX = ehIntermed->xtics[i+1];
				ehIntermed->edges[e].mbr.MinY = ehIntermed->ytics[j];
				ehIntermed->edges[e].mbr.MaxY = ehIntermed->ytics[j]+1e-10;
			}


			// vertical edge at vertex v
			if (j < ehIntermed->yqtd) {
				int e = GET_VERT_EDGE(i, j);
				ehIntermed->edges[e].mbr.MinY = ehIntermed->ytics[j];
				ehIntermed->edges[e].mbr.MaxY = ehIntermed->ytics[j+1];
				ehIntermed->edges[e].mbr.MinX = ehIntermed->xtics[i];
				ehIntermed->edges[e].mbr.MaxX = ehIntermed->xtics[i]+1e-10;
			}
		}
	}



    double xini = MAX(ehA->xtics[0], ehB->xtics[0]);
    double yini = MAX(ehA->ytics[0], ehB->ytics[0]);
    double xend = MIN(ds->metadata.hist.mbr.MaxX, dsb->metadata.hist.mbr.MaxX);
    double yend = MIN(ds->metadata.hist.mbr.MaxY, dsb->metadata.hist.mbr.MaxY);

    unsigned int N = ceil((xend - xini) * (yend - yini));
	double mean = 0.0;
	double M2 = 0.0;
	double sum_error = 0.0;

    int xdr_start = 0;
    while (xdr_start < ehA->xqtd && ehA->xtics[xdr_start+1] < xini)
        xdr_start++;
    int xdr_end = xdr_start+1;
    while (xdr_end < ehA->xqtd && ehA->xtics[xdr_end] <= xend)
        xdr_end++;
    if (xdr_start == ehA->xqtd)
        return 0;

    // skip non-intersect area on y
    int ydr_start = 0;
    while (ydr_start < ehA->yqtd && ehA->ytics[ydr_start+1] < yini)
        ydr_start++;
    int ydr_end = ydr_start+1;
    while (ydr_end < ehA->yqtd && ehA->ytics[ydr_end] <= yend)
        ydr_end++;
    if (ydr_start == ehA->yqtd)
        return 0;

    int cont=0;

    int xds_atu = 0;
    float result = 0;
    for(int xr = xdr_start; xr < xdr_end; xr++) {

        while(xds_atu < ehB->xqtd && ehB->xtics[xds_atu+1] < ehA->xtics[xr]) // skip when end of s < start of r
            xds_atu++;
        int xds_end = xds_atu+1;
        while(xds_end < ehB->xqtd && ehB->xtics[xds_end] <= ehA->xtics[xr+1]) // increment when end of s < start of r
            xds_end++;

        int yds_atu = 0;

        Envelope er;
        Envelope es;
        er.MinX = ehA->xtics[xr];
        er.MaxX = ehA->xtics[xr+1];

        //valores de llinha
        fprintf(valores,"%d",xr);

        for(int yr = ydr_start; yr < ydr_end; yr++) {

            while(yds_atu < ehB->yqtd && ehB->ytics[yds_atu+1] < ehA->ytics[yr]) // skip when end of s < start of r
                yds_atu++;
            int yds_end = yds_atu+1;
            while(yds_end < ehB->yqtd && ehB->ytics[yds_end] <= ehA->ytics[yr+1]) // increment when end of s < start of r
                yds_end++;

            er.MinY = ehA->ytics[yr];
            er.MaxY = ehA->ytics[yr+1];
            double erarea = ENVELOPE_AREA(er);

            euler_face *ehA_face = &ehA->faces[xr*ehA->yqtd +yr];
            euler_face *ehI_face = &ehIntermed->faces[xr*ehIntermed->yqtd +yr];

            //valores coluna
            fprintf(valores,"\t%d",yr);

            //laço histograma B
            for(int xs = xds_atu; xs < xds_end; xs++) {
                es.MinX = ehB->xtics[xs];
                es.MaxX = ehB->xtics[xs+1];

                for(int ys = yds_atu; ys < yds_end; ys++) {
                    double estimated_result = 0;

                    es.MinY = ehB->ytics[ys];
                    es.MaxY = ehB->ytics[ys+1];

                    char i = ENVELOPE_INTERSECTS(er, es);
                    assert(i != 0);

                    euler_face *ehB_face = &ehB->faces[xs*ehB->yqtd +ys];
                    Envelope inters = EnvelopeIntersection2(er, es);
                    double int_area = ENVELOPE_AREA(inters);

                    //printf("ehr_face[%d][%d].card = %f, avg_h = %f, avg_w = %f\n"
                      //      "ehs_face[%d][%d].card = %f, avg_h = %f, avg_w = %f\n", xr, yr, ehr_face->cardin,ehr_face->avg_height, ehr_face->avg_width, xs, ys, ehs_face->cardin, ehs_face->avg_height, ehs_face->avg_width);

                    double intersections = 0;
                    double p = 1;
                    if(ehA_face->cardin >= 1 || ehB_face->cardin >= 1){
                        intersections = estimate_intersections_mamoulis_papadias(er, es, inters, ehA_face, ehB_face);

                        if(intersections > 0)
                        cont++;

                        //if(intersections< 1.0)
                          //  intersections = 0;

                    }
                    result +=  intersections;
                    estimated_result += intersections;

                    //preenche o valor da face do histograma intermediario
                    ehI_face->cardin += intersections;
                    printf("valor face %lf \n",ehI_face->cardin);




                    //vertice
                    int vr = xr * (ehA->yqtd+1) + yr;
                    int vs = xs * (ehB->yqtd+1) + ys;

                    if(ehB->vertexes[vs].x == ehA->vertexes[vr].x && ehB->vertexes[vs].y == ehA->vertexes[vr].y ){
                        result += ehA->vertexes[vr].cardin * ehB->vertexes[vs].cardin;
                        estimated_result += ehA->vertexes[vr].cardin * ehB->vertexes[vs].cardin;
                    }

                    //calculo vertices
                    double objVertFace = 0;
                    if(ehA_face->cardin > 0)
                    	objVertFace = ehA->vertexes[vr].cardin / ehA_face->cardin;

                    //multiplicando a qtd intersecao com objetos que ultrapassam
                    ehIntermed->vertexes[vr].cardin  = objVertFace * ehI_face->cardin;
                    printf("valor vertices %lf \n",ehIntermed->vertexes[vr].cardin );

                    /*
                    //aresta horizontal
                    int ar = GET_HORZ_EDGE_EHR(xr, yr);
                    int as = GET_HORZ_EDGE_EHS(xs, ys);
                    if (ENVELOPE_INTERSECTS(ehA->edges[ar].mbr, ehB->edges[as].mbr)){
                        //printf("ar = %d, as = %d\n", ar, as);
                        Envelope inters = EnvelopeIntersection2(ehA->edges[ar].mbr, ehB->edges[as].mbr);
                        double int_length = inters.MaxX - inters.MinX;

                        double p = MIN(1, ehA->edges[ar].avg_projection + ehB->edges[as].avg_projection);

                        double cardin_ar = ehA->edges[ar].cardin;
                        double cardin_as = ehB->edges[as].cardin;
                        result -=  cardin_ar * cardin_as * p;
                        estimated_result -=  cardin_ar * cardin_as * p;
                        //result -= estimate_intersections_mp_edges_horz(ehr->edges[ar].mbr, ehs->edges[as].mbr, inters, &ehr->edges[ar], &ehs->edges[as]);

                        /*
                        //media dos objetos do hist A
                        double avgObjA = ehA->edges[ar].avg_projection;
                        //printf("media obj A = %lf \n", avgObjA );

                        //intersecao aresta horz
                        double intEdgeHorz;
                        intEdgeHorz = result;

                        //intersecao das faces
                        double intFace = ehI_face->cardin;

                        //cardinalidade da juncao espacial A com B
                        double cardinJS;
                        cardinJS = ehA->edges[ar].cardin * MIN(1, ((avgObjA * intEdgeHorz)  / ehA->xsize ) ) ;

                        intersection = cardinJS * intFace;

                        ehIntermed->edges[ar].cardin = intersection;

                        printf("valor aresta horizontal = %f \n", intersection );



                    }



                    /*

                    //aresta vertical com horizontal
                    ar = GET_HORZ_EDGE_EHR(xr, yr);
                    as = GET_VERT_EDGE_EHS(xs,ys);
                    if(ENVELOPE_INTERSECTS(ehA->edges[ar].mbr, ehB->edges[as].mbr)){
                    	Envelope inters = EnvelopeIntersection2(ehA->edges[ar].mbr, ehB->edges[as].mbr);
                    	double int_length = inters.MaxX - inters.MinX;


                    	//intersetion edge horz with vert
                    	double p = MIN(1, ehA->edges[ar].avg_projection + ehB->edges[as].avg_projection);

                    	double cardin_ar = ehA->edges[ar].cardin;
                    	double cardin_as = ehB->edges[as].cardin;

                    	ehIntermed->edges[ar].cardin += p * cardin_ar * cardin_as;
                    	//ehIntermed->edges[ar].cardin += 1;

                       // printf("valor aresta horz/vert = %f \n", p * cardin_ar * cardin_as );


/*
                    	//cardinalidade da juncao espacial A com B
                    	double cardinJS;
                    	cardinJS += cardin_ar * cardin_as * MIN(1, (( ehA->edges[ar].avg_projection + ehB->edges[as].avg_projection )  / ehA->xsize ) )  * MIN(1, (( ehA->edges[ar].avg_projection + ehB->edges[as].avg_projection )  / ehB->xsize ) );
                    	//cardinJS += p * cardin_ar * cardin_as * MIN(1, (( ehA->edges[ar].avg_projection + ehB->edges[as].avg_projection )  / ehA->xsize ) )  * MIN(1, (( ehA->edges[ar].avg_projection + ehB->edges[as].avg_projection )  / ehB->xsize ) );
                    	ehIntermed->edges[ar].cardin += cardinJS;

                    }


                    ar = GET_VERT_EDGE_EHR(xr, yr);
                    as = GET_HORZ_EDGE_EHS(xs,ys);
                    if(ENVELOPE_INTERSECTS(ehA->edges[ar].mbr, ehB->edges[as].mbr)){
                    	Envelope inters = EnvelopeIntersection2(ehA->edges[ar].mbr, ehB->edges[as].mbr);
                    	double int_length = inters.MaxX - inters.MinX;

                    	//intersetion edge vert with horz
                    	double p = MIN(1, ehA->edges[ar].avg_projection + ehB->edges[as].avg_projection);

                    	double cardin_ar = ehA->edges[ar].cardin;
                    	double cardin_as = ehB->edges[as].cardin;

                    	ehIntermed->edges[ar].cardin += p * cardin_ar * cardin_as;
                    	//ehIntermed->edges[ar].cardin += 1;

                    	// printf("valor aresta vert/horz = %f \n", p * cardin_ar * cardin_as );

                    }
*/
/*

                    //aresta vertical
                    ar = GET_VERT_EDGE_EHR(xr, yr);
                    as = GET_VERT_EDGE_EHS(xs,ys);
                    if (ENVELOPE_INTERSECTS(ehA->edges[ar].mbr, ehB->edges[as].mbr)){
                        Envelope inters = EnvelopeIntersection2(ehA->edges[ar].mbr, ehB->edges[as].mbr);
                        double int_length = inters.MaxY - inters.MinY;
                        double fraction_ar = int_length / (ehA->edges[ar].mbr.MaxY - ehA->edges[ar].mbr.MinY);
                        double fraction_as = int_length / (ehB->edges[as].mbr.MaxY - ehB->edges[as].mbr.MinY);

                        double p = MIN(1, ehA->edges[ar].avg_projection + ehB->edges[as].avg_projection);
                        //printf(" edge p = %f\n", ehA->edges[ar].avg_projection + ehB->edges[as].avg_projection);

                        double cardin_ar = ehA->edges[ar].cardin;
                        double cardin_as = ehB->edges[as].cardin;
                        result -=  cardin_ar * cardin_as * p;
                        estimated_result -=  cardin_ar * cardin_as * p;
                        //result -= estimate_intersections_mp_edges_vert(ehr->edges[ar].mbr, ehs->edges[as].mbr, inters, &ehr->edges[ar], &ehs->edges[as]);
                    }*/
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

            //calculo arestas
            int ahA = GET_HORZ_EDGE_EHR(xr, yr);
            int avA = GET_VERT_EDGE_EHR(xr, yr);

            //objetos que ultrapassam a face e na aresta horizontal
            double objEdgeHorzFace = 0;
            if(ehA_face->cardin > 0)
            	objEdgeHorzFace = ehA->edges[ahA].cardin / ehA_face->cardin;

            //multiplicando a qtd intersecao com objetos que ultrapassam
            ehIntermed->edges[ahA].cardin = objEdgeHorzFace * ehI_face->cardin;
            printf("valor aresta hori %lf \n",ehIntermed->edges[ahA].cardin);

            //objetos que ultrapassam a face e na aresta vertical
            double objEdgeVertFace = 0;
            if(ehA_face->cardin > 0)
            	objEdgeVertFace = ehA->edges[avA].cardin / ehA_face->cardin;

            //multiplicando a qtd intersecao com objetos que ultrapassam
            ehIntermed->edges[avA].cardin = objEdgeVertFace * ehI_face->cardin;
            printf("valor aresta vert %lf \n\n",ehIntermed->edges[avA].cardin);


            fprintf(valores,"\t%lf",ehA_face->cardin);
            fprintf(valores,"\t%lf",ehI_face->cardin);
            fprintf(valores,"\t%lf",ehA->edges[avA].cardin);
            fprintf(valores,"\t%lf",ehIntermed->edges[avA].cardin);
            fprintf(valores,"\t%lf",ehA->edges[ahA].cardin);
            fprintf(valores,"\t%lf",ehIntermed->edges[ahA].cardin);

            int vr = xr * (ehA->yqtd+1) + yr;
            fprintf(valores,"\t%lf",ehA->vertexes[vr].cardin);
            fprintf(valores,"\t%lf",ehIntermed->vertexes[vr].cardin);
            fprintf(valores,"%s","\n");


        }


    }


    printf("cont = %d" , cont);
    printf("result = %f ", result);


    //

	/*
	//verificar valores intersecao
    for(int i = 0; i <= ehIntermed->xqtd; i++){
    	Envelope rsA,rsB;
    	rsA.MinX = ehA->xtics[i];
    	rsB.MinX = ehB->xtics[i];

        if (i < ehA->xqtd)
        	rsA.MaxX = ehA->xtics[i+1];
    	if (i < ehB->xqtd)
    		rsB.MaxX = ehB->xtics[i+1];

        for(int j = 0; j <= ehIntermed->yqtd; j++) {
        	rsA.MinY = ehA->ytics[j];
        	rsB.MinY = ehB->ytics[j];

        	if (j < ehA->yqtd)
   				rsA.MaxY = ehA->ytics[j+1];
           	if (j < ehB->yqtd)
   				rsB.MaxY = ehB->ytics[j+1];

			// face
			if (i < ehIntermed->xqtd && j < ehIntermed->yqtd) {
				//if (ENVELOPE_INTERSECTS(rsA, rsB)) {
				Envelope inters = EnvelopeIntersection2(rsA,rsB);

				euler_face *faceA = &ehA->faces[i*ehA->yqtd +j];
				euler_face *faceB = &ehB->faces[i*ehB->yqtd +j];
				euler_face *faceI = &ehIntermed->faces[i*ehIntermed->yqtd +j];

				if(faceA->cardin >= 1 || faceB->cardin >=1)
					faceI->cardin += estimate_intersections_mamoulis_papadias(rsA,rsB,inters,faceA,faceB);
				printf("cardinar ehInter %lf\n",faceI->cardin );

				//}
			}

			// vertex
			int v = i * (ehA->yqtd+1) + j;
			ehIntermed->vertexes[v].cardin = ehA->vertexes[v].cardin;
			//printf("vertices ehInter %lf vertices ehA %lf\n", ehIntermed->vertexes[v].cardin, ehA->vertexes[v].cardin);
			

    		// aresta horizontal
			if (i < ehIntermed->xqtd) {

				if (i < ehIntermed->xqtd && j < ehIntermed->yqtd) {
					//int e = GET_HORZ_EDGE(i, j);

		           	Envelope inters = EnvelopeIntersection2(rsA,rsB);
                    euler_face *faceA = &ehA->faces[i*ehA->yqtd +j];
                    euler_face *faceB = &ehB->faces[i*ehB->yqtd +j];
                    int e = GET_HORZ_EDGE_EHR(i, j);
                    int d = GET_HORZ_EDGE_EHS(i, j);

					//media dos objetos do hist A
                   	double avgObjA = ehA->edges[e].avg_projection;
					//printf("media obj A = %lf \n", avgObjA );

					//intersecao aresta horz
					double intEdgeHorz;
					intEdgeHorz = estimate_intersections_mp_edges_horz(rsA, rsB, inters, &ehA->edges[e], &ehB->edges[d]);
					//printf("intersecao arest horz = %lf \n", intEdgeHorz);

					//intersecao das faces
					double intFace;
					intFace = estimate_intersections_mamoulis_papadias(rsA, rsB, inters, faceA, faceB);
					//printf("intersecao das faces = %lf \n", intFace);

					//cardinalidade da juncao espacial A com B
					double cardinJS;
					cardinJS = ehA->edges[e].cardin * MIN(1, ((avgObjA * intEdgeHorz)  / ehA->xsize ) ) ;

					//intersection = cardinJS * intFace;

					//ehIntermed->edges[e].cardin = intersection;

					//printf("Intersection = %lf \n", intersection );


				}
			}


			//aresta vertical
			if (j < ehIntermed->yqtd) {

				if (i < ehIntermed->xqtd && j < ehIntermed->yqtd) {
					int e = GET_VERT_EDGE(i, j);

					euler_face *face = &ehA->faces[i*ehA->yqtd +j];
					//printf("aresta cardin ehA %lf\n",ehA->edges[e].cardin);
					if(face->cardin > 0 && ehA->edges[e].cardin > 0)
						ehIntermed->edges[e].cardin = (ehA->edges[e].cardin /face->cardin);;
					else
						ehIntermed->edges[e].cardin = ehA->edges[e].cardin;
					//printf("aresta cardin ehIntermed %lf\n",ehIntermed->edges[e].cardin);
				}

            }

        }
    } */

    fclose(valores);


    return ehIntermed;
}
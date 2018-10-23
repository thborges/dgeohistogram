#include "histogram.h"
#include "ehist.h"

#define GET_VERT_EDGE(x, y) ((x == ehA->xqtd) ? (x * (2*ehA->yqtd+1) + y) : (x * (2*ehA->yqtd+1) + 2*y + 1))
#define GET_HORZ_EDGE(x, y) (x * (2*ehA->yqtd+1) + 2*y)

euler_histogram *eh_generate_intermed(dataset *ds,euler_histogram *ehA,euler_histogram *ehB) {

    int intersecao = 1;

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
	printf("Generating histogram of size: %d x %d\n", ehIntermed->xqtd, ehIntermed->yqtd);


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

	//verificar valores intersecao
    for(int i = 0; i <= ehIntermed->xqtd; i++){
        for(int j = 0; j <= ehIntermed->yqtd; j++) {

			// face
			if (i < ehA->xqtd && j < ehA->yqtd) {
				euler_face *face = &ehA->faces[i*ehA->yqtd +j];
				euler_face *face2 = &ehIntermed->faces[i*ehA->yqtd +j];
				face2->cardin = face->cardin * intersecao;
			}

			// vertex
			int v = i * (ehA->yqtd+1) + j;
			ehIntermed->vertexes[v].cardin = ehA->vertexes[v].cardin * intersecao;
			//printf("vertices ehInter %lf vertices ehA %lf\n", ehIntermed->vertexes[v].cardin, ehA->vertexes[v].cardin);
			

    		// aresta horizontal
			if (i < ehIntermed->xqtd) {

				if (i < ehIntermed->xqtd && j < ehIntermed->yqtd) {
					int e = GET_HORZ_EDGE(i, j);
					euler_face *face = &ehA->faces[i*ehA->yqtd +j];
					//printf("aresta cardin ehA %lf\n",ehA->edges[e].cardin);
		            if(face->cardin > 0 && ehA->edges[e].cardin > 0)
						ehIntermed->edges[e].cardin = (ehA->edges[e].cardin /face->cardin) * intersecao;
					else
						ehIntermed->edges[e].cardin = ehA->edges[e].cardin;
		            //printf("aresta cardin ehIntermed %lf\n",ehIntermed->edges[e].cardin);
				}else{

				}

			}

			//aresta vertical
			if (j < ehIntermed->yqtd) {

				if (i < ehIntermed->xqtd && j < ehIntermed->yqtd) {
					int e = GET_VERT_EDGE(i, j);

					euler_face *face = &ehA->faces[i*ehA->yqtd +j];
					//printf("aresta cardin ehA %lf\n",ehA->edges[e].cardin);
					if(face->cardin > 0 && ehA->edges[e].cardin > 0)
						ehIntermed->edges[e].cardin = (ehA->edges[e].cardin /face->cardin) * intersecao;
					else
						ehIntermed->edges[e].cardin = ehA->edges[e].cardin;
					//printf("aresta cardin ehIntermed %lf\n",ehIntermed->edges[e].cardin);
				}

            }

        }
    }

    return ehIntermed;
}
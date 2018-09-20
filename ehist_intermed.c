#include "histogram.h"
#include "ehist.h"

euler_histogram* eh_generate_intermed(euler_histogram *ehA,euler_histogram *ehB) {

    int intersecao = 1;

    //gerando um novo histograma
	euler_histogram *ehIntermed = g_new(euler_histogram, 1);

    //escolher hist A ou B para copia
    ehIntermed->xqtd = ehA->xqtd;
	ehIntermed->yqtd = ehA->yqtd;
	ehIntermed->xtics = g_new(double, ehA->xqtd+1);
	ehIntermed->ytics = g_new(double, ehA->yqtd+1);
	ehIntermed->faces = g_new0(euler_face, ehA->xqtd * ehA->yqtd);
	ehIntermed->edges = g_new0(euler_edge, ehA->xqtd * (ehA->yqtd+1) + ehA->yqtd * (ehA->xqtd + 1));
	ehIntermed->vertexes = g_new0(euler_vertex, ehA->xqtd * ehA->yqtd+ + ehA->xqtd + ehA->yqtd + 1);

	ehIntermed->xsize = ehA->xsize;
	ehIntermed->ysize = ehA->ysize;

	//verificar valores intersecao

    for(int i; i <= ehIntermed->xqtd; i++){
        for(int j = 0; j <= ehIntermed->yqtd; j++) {

    		// aresta horizontal
			if (i < ehIntermed->xqtd) {
				int e = i * (2*ehA->yqtd+1) + 2*j;
				euler_face *face = &ehA->faces[i*ehA->yqtd +j];
				ehIntermed->edges[e].cardin = (ehA->edges[e].cardin /face->cardin) * intersecao;
			}

			//aresta vertical
			if (j < ehIntermed->yqtd) {
				int e;
				if (i == ehA->xqtd)
					e = i * (2*ehA->yqtd+1) + j;
				else
					e = i * (2*ehA->yqtd+1) + 2*j + 1;

				euler_face *face = &ehA->faces[i*ehA->yqtd +j];
				ehIntermed->edges[e].cardin = (ehA->edges[e].cardin /face->cardin) * intersecao;
            }

        }
    }

    return ehIntermed;
}

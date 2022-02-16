#include "histogram.h"
#include "ehist.h"
#include <math.h>
#include <float.h>

#define GET_VERT_EDGE(eh,x, y) ((x == eh->xqtd) ? (x * (2*eh->yqtd+1) + y) : (x * (2*eh->yqtd+1) + 2*y + 1))
#define GET_HORZ_EDGE(eh,x, y) (x * (2*eh->yqtd+1) + 2*y)

void verifyHistogram(euler_histogram *eh){
	int v = 0;
	for(int i = 0; i <= eh->xqtd; i++) { // x
		for(int j = 0; j <= eh->yqtd; j++) { // y

			if(eh->vertexes[v].cardin > 0){
				printf("\nValor vertice maior que 0");
			}
			v++;


			// horizontal edge at vertex v
			if (i < eh->xqtd) {
				int e = GET_HORZ_EDGE(eh,i, j);
				if(eh->edges[e].cardin > 0)
					printf("\nValor arest horz maior que 0");
			}


			// vertical edge at vertex v
			if (j < eh->yqtd) {
				int e = GET_VERT_EDGE(eh,i, j);
				if(eh->edges[e].cardin > 0)
					printf("\nValor arest vert maior que 0");
			}

		}
	}
}

void salvaArquivo(euler_histogram *ehA, euler_histogram *ehIntermed, euler_histogram *ehC, dataset *dsc){

    //arquivo txt
    FILE *valores,*tabela,*erroMedio;
    FILE *erroCell;
	char nome[200] = "../dados_tabela/values_";

	strcat(nome,dsc->metadata.name);
	strcat(nome,".tese.txt");

    valores = fopen(nome,"w");
    if(valores == NULL){
      printf("Erro na abertura do arquivo!");
    }

    tabela = fopen("../dados_tabela/tabela","a");
    if(tabela == NULL){
      printf("Erro na abertura do arquivo!");
    }

    erroMedio = fopen("../dados_tabela/erroMedio","a");
    if(erroMedio == NULL){
    	printf("Erro na abertura do arquivo!");
    }


    erroCell = fopen("../dados_tabela/erroCell","a");
        if(erroCell == NULL){
        	printf("Erro na abertura do arquivo!");
        }

    fprintf(valores, "%s", "linha \tcoluna \tfaceReal \tfaceIntermed \tarestaVertReal \tarestaVertInter");
    fprintf(valores, "%s", "\tarestaHorzReal \tarestaHorzIntermed  \tvertReal \tvertIntermed ");
    fprintf(valores, "%s", "\t\terroFace \terroArestaVert \terroArestaHorz \terroVert \t\tcardinReal \tcardinEstimado\n");

    if(ftell(tabela) == 0)
    	fprintf(tabela, "%s", "juncao \t cardTotal \t cardFace \t cardAresta \t cardVertice\n");

    if(ftell(erroCell) == 0)
    	fprintf(erroCell, "%s", "Min\tMax\tDesvio\tMedia\t<=5\t<=10\t<=20\t<=50\n");

    double sumFace = 0, sumEdgeVert = 0, sumEdgeHorz = 0,sumVert = 0, sumCardinReal = 0, sumCardinEstimada = 0;
    double sumFaceEstimadaErro = 0, sumEdgeVertEstimadaErro = 0, sumEdgeHorzEstimadaErro = 0,sumVertEstimadaErro= 0;
    double sumFaceEstimada = 0, sumEdgeVertEstimada = 0, sumEdgeHorzEstimada = 0,sumVertEstimada= 0;

    double qtd5 = 0,qtd10 = 0,qtd20 = 0,qtd50 = 0;
    double min = 0, max = 0, media=0, desvio=0;

    for(int xr = 0; xr < ehA->xqtd; xr++){
        //valores de llinha
        fprintf(valores,"%d",xr);

    	for(int yr = 0; yr < ehA->yqtd; yr++){
            //valores coluna
            fprintf(valores,"\t%d",yr);

            euler_face *ehI_face = &ehIntermed->faces[xr*ehIntermed->yqtd +yr];
            euler_face *ehReal_face = &ehC->faces[xr*ehC->yqtd +yr];

            //calculo arestas
            int ahA = GET_HORZ_EDGE(ehIntermed,xr, yr);
            int avA = GET_VERT_EDGE(ehIntermed,xr, yr);


    		fprintf(valores,"\t%lf",ehReal_face->cardin);
			fprintf(valores,"\t%lf",ehI_face->cardin);
			fprintf(valores,"\t%lf",ehC->edges[avA].cardin);
			fprintf(valores,"\t%lf",ehIntermed->edges[avA].cardin);
			fprintf(valores,"\t%lf",ehC->edges[ahA].cardin);
			fprintf(valores,"\t%lf",ehIntermed->edges[ahA].cardin);

			int vr = xr * (ehA->yqtd+1) + yr;
			fprintf(valores,"\t%lf",ehC->vertexes[vr].cardin);
			fprintf(valores,"\t%lf",ehIntermed->vertexes[vr].cardin);

			double erroFace = fabs(ehReal_face->cardin - ehI_face->cardin);
			fprintf(valores,"\t\t%lf",erroFace);
			double erroEdgeVert = fabs(ehC->edges[avA].cardin- ehIntermed->edges[avA].cardin);
			fprintf(valores,"\t%lf",erroEdgeVert);
			double erroEdgeHorz = fabs(ehC->edges[ahA].cardin - ehIntermed->edges[ahA].cardin);
			fprintf(valores,"\t%lf",erroEdgeHorz);
			double erroVert = fabs(ehC->vertexes[vr].cardin - ehIntermed->vertexes[vr].cardin);
			fprintf(valores,"\t%lf",erroVert);

			double cardinReal = ehReal_face->cardin - (ehC->edges[avA].cardin + ehC->edges[ahA].cardin) + ehC->vertexes[vr].cardin;
			fprintf(valores,"\t\t%lf",cardinReal);
			double cardinEstimado = ehI_face->cardin - (ehIntermed->edges[avA].cardin + ehIntermed->edges[ahA].cardin) + ehIntermed->vertexes[vr].cardin;
			fprintf(valores,"\t%lf",cardinEstimado);

			sumFace += ehReal_face->cardin;
			sumEdgeVert += ehC->edges[avA].cardin,
			sumEdgeHorz += ehC->edges[ahA].cardin;
			sumVert += ehC->vertexes[vr].cardin;

			sumFaceEstimada += ehI_face->cardin;
			sumEdgeVertEstimada += ehIntermed->edges[avA].cardin,
			sumEdgeHorzEstimada += ehIntermed->edges[ahA].cardin;
			sumVertEstimada += ehIntermed->vertexes[vr].cardin;

			sumFaceEstimadaErro += erroFace;
			sumEdgeVertEstimadaErro += erroEdgeVert,
			sumEdgeHorzEstimadaErro += erroEdgeHorz;
			sumVertEstimadaErro += erroVert;

			sumCardinReal += cardinReal;
			sumCardinEstimada += cardinEstimado;

			if(min > cardinEstimado)
				min = cardinEstimado;

			if(max < cardinEstimado)
				max = cardinEstimado;

			double cellerror =0;
			if(cardinReal > 0){
				cellerror = fabs(cardinReal - cardinEstimado)/cardinReal;
				//printf("cellerror %f",cellerror);
			}

			if(cellerror <= 0.05)
				qtd5++;
			if(cellerror <= 0.1)
				qtd10++;
			if(cellerror <= 0.2)
				qtd20++;
			if(cellerror <= 0.5)
				qtd50++;


			fprintf(valores,"%s","\n");


    	}
    }

    double totalCell = (ehA->xqtd*ehA->yqtd);
    media = sumCardinEstimada / totalCell;

    double soma = 0;


    for(int i = 0; i <ehA->xqtd;i++){
        for(int j = 0; j <ehA->yqtd;j++){

			int vr = i * (ehA->yqtd+1) + j;
            euler_face *ehI_face = &ehIntermed->faces[i*ehIntermed->yqtd +j];
            //euler_face *ehReal_face = &ehC->faces[i*ehC->yqtd +j];

            //calculo arestas
            int ahA = GET_HORZ_EDGE(ehIntermed,i, j);
            int avA = GET_VERT_EDGE(ehIntermed,i, j);
        	double cardinEstimado = ehI_face->cardin - (ehIntermed->edges[avA].cardin + ehIntermed->edges[ahA].cardin) + ehIntermed->vertexes[vr].cardin;

        	soma += pow(cardinEstimado - media,2);

        }
    }
    desvio = sqrt( soma / totalCell );

    double p5 = 0,p10=0,p20=0,p50=0;


    p5 = ((qtd5*100) / totalCell);
	p10 = ((qtd10*100) / totalCell);
	p20 = ((qtd20*100) / totalCell);
	p50 = ((qtd50*100) / totalCell);

    fprintf(erroCell,"%lf",min);
    fprintf(erroCell,"\t%lf",max);
    fprintf(erroCell,"\t%lf",desvio);
    fprintf(erroCell,"\t%lf",media);
    fprintf(erroCell,"\t%lf",p5);
    fprintf(erroCell,"\t%lf",p10);
    fprintf(erroCell,"\t%lf",p20);
    fprintf(erroCell,"\t%lf\n",p50);



    double erroTotalFace = 0,erroTotalEdgeVert = 0,erroTotalEdgeHorz = 0,erroTotalVertice = 0;
    erroTotalFace = sumFaceEstimadaErro/sumFace;
    erroTotalEdgeVert = sumEdgeVertEstimadaErro/sumEdgeVert;
    erroTotalEdgeHorz = sumEdgeHorzEstimadaErro/sumEdgeHorz;
    erroTotalVertice = sumVertEstimadaErro/sumVert;

    fprintf(valores,"\n\n\t\t%lf",sumFace);
    fprintf(valores,"\t%lf",sumFaceEstimada);
    fprintf(valores,"\t%lf",sumEdgeVert);
    fprintf(valores,"\t%lf",sumEdgeVertEstimada);
    fprintf(valores,"\t%lf",sumEdgeHorz);
    fprintf(valores,"\t%lf",sumEdgeHorzEstimada);
    fprintf(valores,"\t%lf",sumVert);
    fprintf(valores,"\t%lf",sumVertEstimada);
    fprintf(valores,"\t\t%lf",erroTotalFace);
    fprintf(valores,"\t%lf",erroTotalEdgeVert);
    fprintf(valores,"\t%lf",erroTotalEdgeHorz);
    fprintf(valores,"\t%lf",erroTotalVertice);
    fprintf(valores,"\t\t%lf",sumCardinReal);
    fprintf(valores,"\t%lf",sumCardinEstimada);
	fprintf(valores,"%s","\n");

    fprintf(valores,"\n\n\t\t%lf",sumFace);
	fprintf(valores,"\t%lf",sumFaceEstimada);
	fprintf(valores,"\t%lf",sumEdgeVert + sumEdgeHorz);
	fprintf(valores,"\t%lf", (sumEdgeVertEstimada +sumEdgeHorzEstimada));
	fprintf(valores,"\t%lf",sumVert);
	fprintf(valores,"\t%lf",sumVertEstimada);
	fprintf(valores,"\t\t%lf",erroTotalFace);
	fprintf(valores,"\t%lf", ( (sumEdgeVertEstimadaErro + sumEdgeHorzEstimadaErro)/(sumEdgeVert +sumEdgeHorz ) ));
	fprintf(valores,"\t%lf",erroTotalVertice);

   fprintf(tabela,"\t%lf\t",sumCardinEstimada);
   fprintf(tabela,"%lf\t",sumFaceEstimada);
   fprintf(tabela,"%lf\t",sumEdgeVertEstimada +sumEdgeHorzEstimada);
   fprintf(tabela,"%lf\n",sumVertEstimada);

   //erro medio
   fprintf(erroMedio,"%lf",sumFace);
   fprintf(erroMedio,"\t%lf",sumFaceEstimada);
   fprintf(erroMedio,"\t%lf",sumEdgeVert + sumEdgeHorz);
   fprintf(erroMedio,"\t%lf", (sumEdgeVertEstimada +sumEdgeHorzEstimada));
   fprintf(erroMedio,"\t%lf",sumVert);
   fprintf(erroMedio,"\t%lf",sumVertEstimada);
   fprintf(erroMedio,"\t%lf",fabs(sumFace - sumFaceEstimada)/sumFace);
   fprintf(erroMedio,"\t%lf",fabs((sumEdgeVert + sumEdgeHorz) - (sumEdgeVertEstimada +sumEdgeHorzEstimada))/(sumEdgeVert + sumEdgeHorz));
   fprintf(erroMedio,"\t%lf\n",fabs(sumVert - sumVertEstimada)/sumVert);

    fclose(valores);
    fclose(tabela);
    fclose(erroMedio);
    fclose(erroCell);

}


double returnPercent(euler_histogram *ehA, euler_histogram *ehIntermed, int xr, int yr){
	double res = 0.0;

	euler_face *ehI_face = &ehIntermed->faces[xr*ehIntermed->yqtd +yr];
	euler_face *ehA_face = &ehA->faces[xr*ehA->yqtd +yr];

	if(ehA_face->cardin > 0.0)
		res = ehI_face->cardin / ehA_face->cardin;

	return res;

}

double insert_data_edge_vertice(euler_histogram *ehA, euler_histogram *ehIntermed,double result){

	 for(int xr = 0; xr <= ehIntermed->xqtd; xr++){

	    for(int yr = 0; yr <= ehIntermed->yqtd; yr++){

	    	double perface = 0.0;
	    	//calculo vertice
	    	int vr = xr * (ehIntermed->yqtd+1) + yr;
	    	int avA = 0	,ahA = 0;

	    	//calculo arestas
            if (xr < ehIntermed->xqtd){
	    		ahA = GET_HORZ_EDGE(ehIntermed,xr, yr);
	    	}

            if (yr < ehIntermed->yqtd){
            	avA = GET_VERT_EDGE(ehIntermed,xr, yr);
            }

	    	//vertices lado esquerdo
	    	if(xr == 0){

    	    	//calculo porcentagem face
        		perface = returnPercent(ehA,ehIntermed,xr,yr);

    			//1 face
    			if( yr == 0 ){
        			//baixo
    				ehIntermed->vertexes[vr].cardin  = ehA->vertexes[vr].cardin * perface;
    			}else if( yr == ehIntermed->yqtd ){
    				//cima
    				perface = returnPercent(ehA,ehIntermed,xr,yr-1);
    				ehIntermed->vertexes[vr].cardin  = ehA->vertexes[vr].cardin * perface;
    			}else if( yr < ehIntermed->yqtd ){
        		//2 faces
					perface += returnPercent(ehA,ehIntermed,xr,yr-1);
					ehIntermed->vertexes[vr].cardin  = ehA->vertexes[vr].cardin * (perface / 2);

				}

    			result += ehIntermed->vertexes[vr].cardin;

	    	}

    		//vertices lado direito
			if(	xr == ehIntermed->xqtd){

	    		perface = returnPercent(ehA,ehIntermed,xr-1,yr);
				//1 face
				if( yr == 0 ){
					//baixo
					ehIntermed->vertexes[vr].cardin  = ehA->vertexes[vr].cardin * perface;
				}else if( yr == ehIntermed->yqtd ){
					//cima
					perface = returnPercent(ehA,ehIntermed,xr-1,yr-1);
					ehIntermed->vertexes[vr].cardin  = ehA->vertexes[vr].cardin * perface;
				}else if( yr < ehIntermed->yqtd ){
				//2 faces
					perface += returnPercent(ehA,ehIntermed,xr-1,yr-1);
					ehIntermed->vertexes[vr].cardin  = ehA->vertexes[vr].cardin * (perface / 2);

				}

    			result += ehIntermed->vertexes[vr].cardin;
			}

			//vertices inferior
			if(	yr == 0 && xr > 0 && xr < ehIntermed->xqtd ){
	    		perface = returnPercent(ehA,ehIntermed,xr,yr);
				if( xr < ehIntermed->xqtd ){
					perface += returnPercent(ehA,ehIntermed,xr-1,yr);
					ehIntermed->vertexes[vr].cardin  = ehA->vertexes[vr].cardin * (perface / 2);

				}

				result += ehIntermed->vertexes[vr].cardin;
			}
			//vertices superior
			if(	yr == ehIntermed->yqtd && xr > 0 && xr < ehIntermed->xqtd){
	    		perface = returnPercent(ehA,ehIntermed,xr,yr-1);
				perface += returnPercent(ehA,ehIntermed,xr-1,yr-1);
				ehIntermed->vertexes[vr].cardin  = ehA->vertexes[vr].cardin * (perface / 2);


				result += ehIntermed->vertexes[vr].cardin;
			}

	    	//vértices com 4 face
	    	if(	( (xr > 0) && (xr < ehIntermed->xqtd) ) && ( (yr > 0) && (yr < ehIntermed->yqtd) )   ){

		    	//calculo porcentagem face
	    		perface = returnPercent(ehA,ehIntermed,xr,yr);
	    		perface += returnPercent(ehA,ehIntermed,xr-1,yr);
	    		perface += returnPercent(ehA,ehIntermed,xr,yr-1);
	    		perface += returnPercent(ehA,ehIntermed,xr-1,yr-1);

	    		ehIntermed->vertexes[vr].cardin  = ehA->vertexes[vr].cardin * (perface / 4);

    			result += ehIntermed->vertexes[vr].cardin;
	    	}


	    	//arestas infeirores
	    	if(	yr == 0 && xr < ehIntermed->xqtd){
		    	//calculo porcentagem face
	    		perface = returnPercent(ehA,ehIntermed,xr,yr);
		    	ehIntermed->edges[ahA].cardin = ehA->edges[ahA].cardin * perface;

    			result -= ehIntermed->edges[ahA].cardin;

	    	}
	    	//arestas superiores
			if(	yr == ehIntermed->yqtd && xr < ehIntermed->xqtd){

				perface = returnPercent(ehA,ehIntermed,xr,yr-1);
				ehIntermed->edges[ahA].cardin = ehA->edges[ahA].cardin * perface;

    			result -= ehIntermed->edges[ahA].cardin;

			}

			//arestas esquerda
			if(	xr == 0 && yr < ehIntermed->yqtd) {
		    	//calculo porcentagem face
	    		perface = returnPercent(ehA,ehIntermed,xr,yr);
				ehIntermed->edges[avA].cardin = ehA->edges[avA].cardin * perface;

    			result -= ehIntermed->edges[avA].cardin;
			}
			//aresta direita
			if( xr == ehIntermed->xqtd && yr < ehIntermed->yqtd){

				perface = returnPercent(ehA,ehIntermed,xr-1,yr);
				ehIntermed->edges[avA].cardin = ehA->edges[avA].cardin * perface;

    			result -= ehIntermed->edges[avA].cardin;
			}

			//arestas com 2 faces na horizontal
			if(	( (xr > 0) && (xr < ehIntermed->xqtd) && (yr < ehIntermed->yqtd) ) ){
		    	//calculo porcentagem face
	    		perface = returnPercent(ehA,ehIntermed,xr,yr);
				perface += returnPercent(ehA,ehIntermed,xr-1,yr);

		    	ehIntermed->edges[avA].cardin = (ehA->edges[avA].cardin * (perface / 2));


    			result -= ehIntermed->edges[avA].cardin;
			}

			//arestas com 2 faces na vertical
			if(	( (yr > 0) && (yr < ehIntermed->yqtd) && (xr < ehIntermed->xqtd) ) ){
		    	//calculo porcentagem face
	    		perface = returnPercent(ehA,ehIntermed,xr,yr);
				perface += returnPercent(ehA,ehIntermed,xr,yr-1);


		    	ehIntermed->edges[ahA].cardin =  (ehA->edges[ahA].cardin * (perface / 2));

    			result -= ehIntermed->edges[ahA].cardin;
			}

	    }
	 }


	 return result;
}



euler_histogram *eh_generate_intermed_real(dataset *ds,dataset *dsb, euler_histogram *ehA,euler_histogram *ehB) {

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
	//ehIntermed->xtics[ehA->xqtd] = ehA->xtics[ehA->xqtd];
	ehIntermed->xtics[ehA->xqtd] = ds->metadata.hist.mbr.MaxX;

	// Y tics
	for(int i = 0; i < ehA->yqtd; i++)
		ehIntermed->ytics[i] = ehA->ytics[i];
	//ehIntermed->ytics[ehA->yqtd] = ehA->ytics[ehA->yqtd];
	ehIntermed->ytics[ehIntermed->yqtd] = ds->metadata.hist.mbr.MaxY;


	//set edges and vertexes
	int v = 0;
	for(int i = 0; i < ehIntermed->xqtd; i++) { // x
		for(int j = 0; j < ehIntermed->yqtd; j++) { // y

			// set vetex at (i,j)
			v = i * (ehIntermed->yqtd +1 ) + j;
			ehIntermed->vertexes[v].x = ehIntermed->xtics[i];
			ehIntermed->vertexes[v].y = ehIntermed->ytics[j];
			ehIntermed->vertexes[v].cardin = 0.0;
			//v++;

			// horizontal edge at vertex v
			if (i < ehIntermed->xqtd) {
				int e = GET_HORZ_EDGE(ehIntermed, i, j);
				ehIntermed->edges[e].cardin = 0.0;
				ehIntermed->edges[e].mbr.MinX = ehIntermed->xtics[i];
				ehIntermed->edges[e].mbr.MaxX = ehIntermed->xtics[i+1];
				ehIntermed->edges[e].mbr.MinY = ehIntermed->ytics[j];
				ehIntermed->edges[e].mbr.MaxY = ehIntermed->ytics[j]+1e-10;
			}


			// vertical edge at vertex v
			if (j < ehIntermed->yqtd) {
				int e = GET_VERT_EDGE(ehIntermed, i, j);
				ehIntermed->edges[e].cardin = 0.0;
				ehIntermed->edges[e].mbr.MinY = ehIntermed->ytics[j];
				ehIntermed->edges[e].mbr.MaxY = ehIntermed->ytics[j+1];
				ehIntermed->edges[e].mbr.MinX = ehIntermed->xtics[i];
				ehIntermed->edges[e].mbr.MaxX = ehIntermed->xtics[i]+1e-10;
			}

			//usedarea
            if(i < ehIntermed->xqtd && j < ehIntermed->yqtd){
            	euler_face *face = &ehIntermed->faces[i*ehIntermed->yqtd +j];

				face->usedarea.MinX = face->usedarea.MinY = DBL_MAX;
				face->usedarea.MaxX = face->usedarea.MaxY = -DBL_MAX;

            }

			//face
			euler_face *faceA = &ehA->faces[i*ehA->yqtd +j];
			euler_face *faceIntermed = &ehIntermed->faces[i*ehIntermed->yqtd +j];
			faceIntermed->avg_area = faceA->avg_area;
			faceIntermed->avg_height = faceA->avg_height;
			faceIntermed->avg_width = faceA->avg_width;

		}
	}

	verifyHistogram(ehIntermed);


    double xini = MAX(ehA->xtics[0], ehB->xtics[0]);
    double yini = MAX(ehA->ytics[0], ehB->ytics[0]);
    double xend = MIN(ds->metadata.hist.mbr.MaxX, dsb->metadata.hist.mbr.MaxX);
    double yend = MIN(ds->metadata.hist.mbr.MaxY, dsb->metadata.hist.mbr.MaxY);

    //unsigned int N = ceil((xend - xini) * (yend - yini));
	//double mean = 0.0;
	//double M2 = 0.0;
	//double sum_error = 0.0;

    int xdr_start = 0;
    while (xdr_start < ehA->xqtd && ehA->xtics[xdr_start+1] < xini)
        xdr_start++;
    int xdr_end = xdr_start+1;
    while (xdr_end < ehA->xqtd && ehA->xtics[xdr_end] <= xend)
        xdr_end++;
    if (xdr_start == ehA->xqtd)
        printf("xdr_start == ehA->xqtd");

    // skip non-intersect area on y
    int ydr_start = 0;
    while (ydr_start < ehA->yqtd && ehA->ytics[ydr_start+1] < yini)
        ydr_start++;
    int ydr_end = ydr_start+1;
    while (ydr_end < ehA->yqtd && ehA->ytics[ydr_end] <= yend)
        ydr_end++;
    if (ydr_start == ehA->yqtd)
        printf("ydr_start == ehA->yqtd");

    int cont=0;

    int xds_atu = 0;
    double result = 0;

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


        for(int yr = ydr_start; yr < ydr_end; yr++) {

            while(yds_atu < ehB->yqtd && ehB->ytics[yds_atu+1] < ehA->ytics[yr]) // skip when end of s < start of r
                yds_atu++;
            int yds_end = yds_atu+1;
            while(yds_end < ehB->yqtd && ehB->ytics[yds_end] <= ehA->ytics[yr+1]) // increment when end of s < start of r
                yds_end++;

            er.MinY = ehA->ytics[yr];
            er.MaxY = ehA->ytics[yr+1];
            //double erarea = ENVELOPE_AREA(er);

            euler_face *ehA_face = &ehA->faces[xr*ehA->yqtd +yr];
            euler_face *ehI_face = &ehIntermed->faces[xr*ehIntermed->yqtd +yr];


            //laço histograma B
            for(int xs = xds_atu; xs < xds_end; xs++) {
                es.MinX = ehB->xtics[xs];
                es.MaxX = ehB->xtics[xs+1];

                for(int ys = yds_atu; ys < yds_end; ys++) {
                    //double estimated_result = 0;

                    es.MinY = ehB->ytics[ys];
                    es.MaxY = ehB->ytics[ys+1];

                    char i = ENVELOPE_INTERSECTS(er, es);
                    assert(i != 0);

                    euler_face *ehB_face = &ehB->faces[xs*ehB->yqtd +ys];
                    Envelope inters = EnvelopeIntersection2(er, es);
                    //double int_area = ENVELOPE_AREA(inters);


                    double intersections = 0.0;
                    //double p = 1;
                    if(ehA_face->cardin > 0.0 && ehB_face->cardin > 0.0){
                    	//intersections = estimate_intersections_mamoulis_papadias(er, es, inters, ehA_face, ehB_face);

                        intersections = estimate_intersections_zu(er, es, inters, ehA_face, ehB_face,ds,dsb);
                        //intersections *= int_area / ENVELOPE_AREA(es);


                        if(intersections > 0.0)
                        cont++;

                        if (isnan(intersections))
                        	intersections = 0;

                        if (intersections < 0.00001)
                        	intersections = 0.00001;


                    }
                    //preenche o valor da face do histograma intermediario
                    ehI_face->cardin += intersections;

                    result +=  intersections;


                }
            }

        }
    }

    result = insert_data_edge_vertice(ehA,ehIntermed,result);

    printf("\n\ncont = %d" , cont);
    printf("result = %f ", result);

    return ehIntermed;
}



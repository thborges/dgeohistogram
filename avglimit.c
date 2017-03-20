#include <math.h>
#include <float.h>
#include "histogram.h"
#include "Dataset.h"




const double epsilon = 1e-7;
double qtdObjetos = 0;
#define MAX(a,b) ((a) > (b) ? a : b)
#define MIN(a,b) ((a) < (b) ? a : b)



bool doubleEqual(double a, double b) {
	return fabs(a - b) < epsilon;
}

// void ValidaQtdObjetos(histogram *original) {
void ValidaQtdObjetos(dataset_histogram *dh) {


	//GET_HISTOGRAM_CELL(dh, i, j)->cardin;

	if (qtdObjetos == 0) {
		for (int x = 0; x < dh->xqtd; x++) {
					for (int y = 0; y < dh->yqtd; y++) {
				celula *c = GET_CELL(original, x, y);

				// qtdObjetos += c->card;
				qtdObjetos += GET_HISTOGRAM_CELL(dh, i, j)->cardin;


			}

		}

	} else {
		double qtdObjAftMethodo = 0;
		for (int x = 0; x < dh->xqtd; x++) {
					for (int y = 0; y < dh->yqtd; y++) {
				celula *c = GET_CELL(original, x, y);

				// qtdObjAftMethodo += c->card;
				qtdObjAftMethodo += GET_HISTOGRAM_CELL(dh, i, j)->cardin;



			}

		}

		if (doubleEqual(qtdObjetos, qtdObjAftMethodo)) {
			printf("Quantidade de objetos eh a MESMA, %f\n", qtdObjetos);
		} else {
			printf("Quantidade de objetos eh a DIFERENTE, %.15f != %.15f \n",
					qtdObjetos, qtdObjAftMethodo);
		}

	}

}

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

/**
 * Funcao para printar na tela
 * modo == 0 - linha
 * modo == 1 - coluna
 *
 */
histogram* metodoSimples(histogram *original, int modo, int xinicial,
		int yinicial) {
	//declaracao das estruturas auxiliares

	double somaMedias = 0;

	//Modo linha
	if (modo == 0) {
		int contCelzero = 0;
		for (int x = 0; x < original->qtd_colunas; ++x) {
			celula *origi = GET_CELL(original, x, xinicial);

			somaMedias += origi->alturaMedia;
// 			if (origi->alturaMedia < epsilon){
// //				printf("col: %d lin: %d tem alturaZero: %f\n", x,xinicial, origi->alturaMedia);
// 				contCelzero++;
// 			}

		}
//		 printf("Antes somaMedias: %f\n", somaMedias);
		// somaMedias = somaMedias / ((original->qtd_colunas - contCelzero) == 0 ? 1 : (original->qtd_colunas - contCelzero));

		somaMedias = somaMedias / original->qtd_colunas;
		 // printf("Antes somaMedias: %f\n", somaMedias);
		 		

		celula *origi = GET_CELL(original, 0, xinicial);
		 // printf("Limite: %f\n", origi->yfim - origi->yini);
		//Verifico se a soma das medias ficou => 80% da comprimento da dimensão da celula
//		if (somaMedias >= (limitanteAltura * limite) || somaMedias < epsilon) {
		if ( (somaMedias - ((origi->yfim - origi->yini) * limite)) > epsilon ) {
 			// printf("entrou ");

			for (int k = 0; k < original->qtd_colunas; ++k) {

				celula *origi = GET_CELL(original, k, xinicial);
				celula *origiProx = GET_CELL(original, k, xinicial + 1);
//				printf("{(%f, %f);", origi->xini, origi->yini);
//				printf(" (%f, %f)} c: %f \n", origi->xfim, origi->yfim,
//						origi->card);
//
//				printf("{(%f, %f);", origiProx->xini, origiProx->yini);
//				printf(" (%f, %f)} c: %f \n", origiProx->xfim, origiProx->yfim,
//						origiProx->card);
//
//				printf("Antes: %f\n", origi->alturaMedia);

				//como faco uma divisao, se nao tiver card em nenhum das col/lin, da uma divisao por 0.
				if ((origi->card) == 0 || (origiProx->card) == 0) {
					if ((origi->card) == 0) {

//						printf("Original == 0\n");
						origi->alturaMedia = origiProx->alturaMedia;

						origi->larguraMedia = origiProx->larguraMedia;

					} else {
					//	printf("OriginalProx == 0\n");
						//a partir daqui se der merda

					}
				} else {
				//	printf("Nada de zeros\n");

					origi->alturaMedia = ((origi->card * origi->alturaMedia)
							+ (origiProx->card * origiProx->alturaMedia))
							/ (origiProx->card + origi->card);

					origi->larguraMedia = ((origi->card * origi->larguraMedia)
							+ (origiProx->card * origiProx->larguraMedia))
							/ (origiProx->card + origi->card);

				}

				origi->card = origi->card + origiProx->card;
			//	printf("Depos: %f\n", origi->alturaMedia);

				origi->yfim = origiProx->yfim;

			}

			histogram* mod = (histogram*) malloc(sizeof(histogram));
			mod->qtd_colunas = original->qtd_colunas;
			mod->qtd_linhas = original->qtd_linhas - 1;
			mod->matriz = (celula*) malloc(
					sizeof(celula) * MaxLin * MaxCol);

			int cont = 0;

			for (int x = 0; x < original->qtd_colunas; x++) {
				//ignorando a linha que foi mergada com a anterior

				for (int y = 0; y < original->qtd_linhas; y++) {

					if ((xinicial + 1) != y) {
						celula *origi = GET_CELL(original, x, y);

						celula *novo = GET_CELL(mod, x, cont);

						novo->xini = origi->xini;
						novo->yini = origi->yini;
						novo->xfim = origi->xfim;
						novo->yfim = origi->yfim;
						novo->card = origi->card;
						novo->larguraMedia = origi->larguraMedia;
						novo->alturaMedia = origi->alturaMedia;
						cont++;
					}

				}
				cont = 0;

			}

			free(original);

			return (mod);
		}

	}
	//------------------------------------------------------------------
//	Modo coluna
//	-----------------------------------
	else {
		int contCelzero = 0;
		for (int y = 0; y < original->qtd_linhas; ++y) {
			celula *origi = GET_CELL(original, yinicial, y);

			somaMedias += origi->larguraMedia;
			// if (origi->larguraMedia < epsilon) contCelzero++;
		}
	//	printf("Antes somaMedias: %f\n", somaMedias);
		// somaMedias = somaMedias
		// 		/ ((original->qtd_linhas - contCelzero) == 0 ?
		// 				1 : (original->qtd_linhas - contCelzero));

		somaMedias = somaMedias / original->qtd_linhas;

		//preciso calcular assim
		//medo o tamanho fixo inicial, faço uma diferença entre xfim - xini e divido este valor pelo limitanteLargura


		celula *origi = GET_CELL(original, yinicial, 0);
//		printf("Dps somaMedias: %f\nlimitanteLargura * limite: %f\n", somaMedias, limitanteLargura * limite);
//		printf("xfim: %f xini: %f\n", origi->xfim, origi->xini);
//		printf("Limitante: %f\n", (fabs(fabs(origi->xfim) - fabs(origi->xini))) );

		//Verifico se a soma das medias ficou >> 80% da comprimento da dimensão da celula
//		if (somaMedias >= (limitanteLargura * limite) || somaMedias < epsilon) {
		if ((somaMedias - ((origi->xfim - origi->xini) * limite )) > epsilon) {

			for (int k = 0; k < original->qtd_linhas; ++k) {

				celula *origi = GET_CELL(original, yinicial, k);
				celula *origiProx = GET_CELL(original, yinicial + 1, k);
//				printf("{(%f, %f);", origi->xini, origi->yini);
//				printf(" (%f, %f)} c: %f \n", origi->xfim, origi->yfim,
//						origi->card);
//
//				printf("{(%f, %f);", origiProx->xini, origiProx->yini);
//				printf(" (%f, %f)} c: %f \n", origiProx->xfim, origiProx->yfim,
//						origiProx->card);
//
//				printf("Antes: %f\n", origi->larguraMedia);

				//como faco uma divisao, se nao tiver card em nenhum das col/lin, da uma divisao por 0.
				if ((origi->card) == 0 || (origiProx->card) == 0) {
					if ((origi->card) == 0) {

//						printf("Original == 0\n");
						origi->alturaMedia = origiProx->alturaMedia;

						origi->larguraMedia = origiProx->larguraMedia;

					} else {
						// printf("OriginalProx == 0\n");

					}
				} else {
		//			printf("Nada de zeros\n");

					origi->alturaMedia = ((origi->card * origi->alturaMedia)
							+ (origiProx->card * origiProx->alturaMedia))
							/ (origiProx->card + origi->card);

					origi->larguraMedia = ((origi->card * origi->larguraMedia)
							+ (origiProx->card * origiProx->larguraMedia))
							/ (origiProx->card + origi->card);

				}

				origi->card = origi->card + origiProx->card;
		//		printf("Depos: %f\n", origi->larguraMedia);

				origi->xfim = origiProx->xfim;

			}

			histogram* mod = (histogram*) malloc(sizeof(histogram));
			mod->qtd_colunas = original->qtd_colunas - 1;
			mod->qtd_linhas = original->qtd_linhas;
			mod->matriz = (celula*) malloc(
					sizeof(celula) * MaxLin * MaxCol);

			int cont = 0;

			for (int x = 0; x < original->qtd_colunas; x++) {
				//ignorando a coluna que foi mergada com a anterior

				if ((yinicial + 1) == x) {
					x++;
				}
				for (int y = 0; y < original->qtd_linhas; y++) {

					celula *origi = GET_CELL(original, x, y);

					celula *novo = GET_CELL(mod, cont, y);

					novo->xini = origi->xini;
					novo->yini = origi->yini;
					novo->xfim = origi->xfim;
					novo->yfim = origi->yfim;
					novo->card = origi->card;
					novo->larguraMedia = origi->larguraMedia;
					novo->alturaMedia = origi->alturaMedia;

				}
				cont++;


			}

			free(original);

			return (mod);
		}

	}

	return (original);

}



void avglimitRUn(dataset_histogram *dh) {


	histogram* h1 = importarHistograma();



	ValidaQtdObjetos(h1);



/**
Parte para passar por todos
*/
//	//Linha
	int merge = MaxCol;
	int xinicial = 0;
	int qtd_lin = h1->qtd_linhas - 1;

	for (int cont = 0; cont < qtd_lin; cont++) {
		h1 = metodoSimples(h1, 0, xinicial, 0);

//		 printf("qtd: %d \n", h1->qtd_linhas);

		if (h1->qtd_linhas != merge) {

			xinicial--;
			merge--;
		}

		xinicial++;
//		printCelula(h1);
//		printf("Mandando a linha: %d para ser comparada com a proxima \n" , xinicial + 1);
//		printCelulaFile(h1, fp);

	}

	//Coluna
	merge = MaxLin;
	int yinicial = 0;
	// -1 para nao comparar a ultima com ??
	int qtd_col = h1->qtd_colunas - 1;

	for (int cont = 0; cont < qtd_col; cont++) {
		h1 = metodoSimples(h1, 1, 0, yinicial);

		// printf("qtd: %d \n", h1->qtd_colunas);

		if (h1->qtd_colunas != merge) {

			yinicial--;
			merge--;
		}

		yinicial++;
//		printCelula(h1);
//		printf("Mandando a coluna: %d para ser comparada com a proxima \n" , yinicial + 1);
//		printCelulaFile(h1, fp);

	}


	ValidaQtdObjetos(h1);

	printf("%dx%d\n",h1->qtd_colunas,h1->qtd_linhas);


	free(h1);

	return 0;
}



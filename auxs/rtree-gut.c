
#include <float.h>
#include "rtree-star.h"
#include "rtree-reinsert.h"


int rtree_choose_subtree_gut(rtree_root *root, rtree_node *n, const Envelope e) {
	int index = -1;

	double minor_enlargement = DBL_MAX;
	double current_area = DBL_MAX;
	for(int i = 0; i < n->used; i++) {
		Envelope eaux = n->dirs[i]->mbr;
		double initarea = ENVELOPE_AREA(eaux);
		ENVELOPE_MERGE(eaux, e);
		double newarea = ENVELOPE_AREA(eaux);
		double increase = newarea - initarea;
		if (increase < minor_enlargement || (increase == minor_enlargement && newarea < current_area)) {
			index = i;
			minor_enlargement = increase;
			current_area = newarea;
		}
	}
	return index;
}


void pick_seeds(rtree_node *n, int seeds[], int m) {
    int no1 = -1;
    int no2 = -2;
    float menorX = 0.0;
    float menorY = 0.0;
    float maiorX = 0.0;
    float maiorY = 0.0;
    float calculo = 0;
    float dist_maior = 0.0;
    int E1 = 0;
    int E2 = 0;
    int indiceI = 0;
    int indiceJ = 0;


    //IMPRESSÃO DOS NÓS 
    ////printf("NÓ CHEIO!\n");
    for (int i = 0; i < m; i++) {
        ////printf("[%d] - (%.0f, %.0f) (%.0f, %.0f) \n", i + 1, n->dirs[i]->mbr.MinX, n->dirs[i]->mbr.MinY, n->dirs[i]->mbr.MaxX, n->dirs[i]->mbr.MaxY);
    }
    
    
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            
         //ENCONTRA A ÁREA DE E1 E E2
          E1 =  (n->dirs[i]->mbr.MaxX - n->dirs[i]->mbr.MinX) * (n->dirs[i]->mbr.MaxY - n->dirs[i]->mbr.MinY);
          E2 =  (n->dirs[j]->mbr.MaxX - n->dirs[j]->mbr.MinX) * (n->dirs[j]->mbr.MaxY - n->dirs[j]->mbr.MinY);
           
       
          
          //ENCONTRA A ÁREA TOTAL DE (J)
          
          //ENCONTRA O MENOR X E MENOR Y
          if (n->dirs[i]->mbr.MinX <= n->dirs[j]->mbr.MinX){
              menorX = n->dirs[i]->mbr.MinX;
          }else {
              menorX = n->dirs[j]->mbr.MinX;
          }
          
          if (n->dirs[i]->mbr.MinY <= n->dirs[j]->mbr.MinY){
              menorY = n->dirs[i]->mbr.MinY;
          }else {
              menorY = n->dirs[j]->mbr.MinY;
          }
          
          //ENCONTRA O MAIOR X E MAIOR Y          
          if (n->dirs[i]->mbr.MaxX >= n->dirs[j]->mbr.MaxX){
              maiorX = n->dirs[i]->mbr.MaxX;
          }else {
              maiorX = n->dirs[j]->mbr.MaxX;
          }
          
          if (n->dirs[i]->mbr.MaxY >= n->dirs[j]->mbr.MaxY){
              maiorY = n->dirs[i]->mbr.MaxY;
          }else {
              maiorY = n->dirs[j]->mbr.MaxY;
          }
              
          
          //CALCULO DA ÁREA TOTAL MENOS A AREA DOS E1 E E2
          calculo = ((maiorX - menorX) * (maiorY - menorY)) - E1 - E2;
          
          if (calculo >= dist_maior){
              dist_maior = calculo;
              indiceI = i;
              indiceJ = j;
          }
                   
        }
    }

    seeds[0] = indiceI;
    seeds[1] = indiceJ;
    no1 = indiceI;
    no2 = indiceJ;

//    printf("\n---- SEMENTES ENCONTRADAS ----\n");
//    printf("[%d] - (%.0f, %.0f) (%.0f, %.0f) \n", no2 + 1, n->dirs[no2]->mbr.MinX, n->dirs[no2]->mbr.MinY, n->dirs[no2]->mbr.MaxX, n->dirs[no2]->mbr.MaxY);
//    printf("[%d] - (%.0f, %.0f) (%.0f, %.0f) \n", no1 + 1, n->dirs[no1]->mbr.MinX, n->dirs[no1]->mbr.MinY, n->dirs[no1]->mbr.MaxX, n->dirs[no1]->mbr.MaxY);

}

void pick_next(rtree_node *n, rtree_node *g1, rtree_node *g2, int seeds[], int m) {
    double minX1 = 0.0;
    double maxX1 = 0.0;
    double minY1 = 0.0;
    double maxY1 = 0.0;

    double minX2 = 0.0;
    double maxX2 = 0.0;
    double minY2 = 0.0;
    double maxY2 = 0.0;

    double areaG1 = 0.0;
    double areaG2 = 0.0;
    double englobarG1 = 0.0;
    double englobarG2 = 0.0;
    double diferencaG1 = 0.0;
    double diferencaG2 = 0.0;

    Envelope mbr1;
    Envelope mbr2;

    mbr1.MinX = n->dirs[seeds[0]]->mbr.MinX;
    mbr1.MinY = n->dirs[seeds[0]]->mbr.MinY;
    mbr1.MaxX = n->dirs[seeds[0]]->mbr.MaxX;
    mbr1.MaxY = n->dirs[seeds[0]]->mbr.MaxY;

    mbr2.MinX = n->dirs[seeds[1]]->mbr.MinX;
    mbr2.MinY = n->dirs[seeds[1]]->mbr.MinY;
    mbr2.MaxX = n->dirs[seeds[1]]->mbr.MaxX;
    mbr2.MaxY = n->dirs[seeds[1]]->mbr.MaxY;

    areaG1 = ENVELOPE_AREA(mbr1);
    areaG2 = ENVELOPE_AREA(mbr2);
//    printf("AREA INICIAL DE G1: %f \n", areaG1);
//    printf("AREA INICIAL DE G2: %f \n\n", areaG2);

    for (int i = 0; i < n->used; ++i) {
        if (!((seeds[0] == i) || (seeds[1] == i))) {
            minX1 = 0.0;
            maxX1 = 0.0;
            minY1 = 0.0;
            maxY1 = 0.0;
            minX2 = 0.0;
            maxX2 = 0.0;
            minY2 = 0.0;
            maxY2 = 0.0;

            //----VERIFICACAO PARA O PRIMEIRO NÓ--------
            //ENCONTRAR X MENOR E MAIOR
            if (n->dirs[i]->mbr.MinX <= mbr1.MinX) {
                minX1 = n->dirs[i]->mbr.MinX;
            } else {
                minX1 = mbr1.MinX;
            }

            if (n->dirs[i]->mbr.MaxX >= mbr1.MaxX) {
                maxX1 = n->dirs[i]->mbr.MaxX;
            } else {
                maxX1 = mbr1.MaxX;
            }

            //ENCONTRA Y MENOR E MAIOR
            if (n->dirs[i]->mbr.MinY <= mbr1.MinY) {
                minY1 = n->dirs[i]->mbr.MinY;
            } else {
                minY1 = mbr1.MinY;
            }

            if (n->dirs[i]->mbr.MaxY >= mbr1.MaxY) {
                maxY1 = n->dirs[i]->mbr.MaxY;
            } else {
                maxY1 = mbr1.MaxY;
            }

            englobarG1 = (maxX1 - minX1) * (maxY1 - minY1);
            diferencaG1 = englobarG1 - areaG1;

//            printf("Englobar G1: %f \n", englobarG1);
//            printf("DIFERENCA DE G1: %f \n\n", diferencaG1);


            //----- VERIFICAÇÃO PARA O SEGUNDO NÓ ------
            //ENCONTRAR X MENOR E MAIOR
            if (n->dirs[i]->mbr.MinX <= mbr2.MinX) {
                minX2 = n->dirs[i]->mbr.MinX;
            } else {
                minX2 = mbr2.MinX;
            }

            if (n->dirs[i]->mbr.MaxX >= mbr2.MaxX) {
                maxX2 = n->dirs[i]->mbr.MaxX;
            } else {
                maxX2 = mbr2.MaxX;
            }


            //ENCONTRA Y MENOR E MAIOR
            if (n->dirs[i]->mbr.MinY <= mbr2.MinY) {
                minY2 = n->dirs[i]->mbr.MinY;
            } else {
                minY2 = mbr2.MinY;
            }

            if (n->dirs[i]->mbr.MaxY >= mbr2.MaxY) {
                maxY2 = n->dirs[i]->mbr.MaxY;
            } else {
                maxY2 = mbr2.MaxY;
            }


            englobarG2 = (maxX2 - minX2) * (maxY2 - minY2);
            diferencaG2 = englobarG2 - areaG2;

//            printf("Englobar G2: %f \n", englobarG2);
//            printf("DIFERENCA DE G2: %f \n\n", diferencaG2);



            
            

            if (diferencaG1 <= diferencaG2) {
//                printf("G1 é menor: \n\n");
                
                mbr1.MinX = minX1;
                mbr1.MinY = minY1;
                mbr1.MaxX = maxX1;
                mbr1.MaxY = maxY1;
                areaG1 = (mbr1.MaxX - mbr1.MinX) * (mbr1.MaxY - mbr1.MinY);
                minX1 = 0.0;
                maxX1 = 0.0;
                minY1 = 0.0;
                maxY1 = 0.0;  
                
                //ADICIONA O NÓ AO G1, VERIFICANDO m
                if (minm(m) - g2->used - (n->used - i) >= 0) {
					// itens restantes são necessários para atingir minm(M) para g2
                    g2->dirs[g2->used] = n->dirs[i];
                    g2->used++;  
                } else {
                    g1->dirs[g1->used] = n->dirs[i];
                    g1->used++;
                }
                

            } else {
//                printf("G2 é menor: \n\n");
                
                mbr2.MinX = minX2;
                mbr2.MinY = minY2;
                mbr2.MaxX = maxX2;
                mbr2.MaxY = maxY2;
                areaG2 = (mbr2.MaxX - mbr2.MinX) * (mbr2.MaxY - mbr2.MinY);
                minX2 = 0.0;
                maxX2 = 0.0;
                minY2 = 0.0;
                maxY2 = 0.0;
                
                //ADICIONA O NÓ AO G2, VERIFICANDO m
                if (minm(m) - g1->used - (n->used - i) >= 0) {
   					// itens restantes são necessários para atingir minm(M) para g2
	                g1->dirs[g1->used] = n->dirs[i];
                    g1->used++;
                } else {
                    g2->dirs[g2->used] = n->dirs[i];
                    g2->used++;
                }
                

            }
                        
            englobarG1 = 0.0;
            englobarG2 = 0.0;

            
            
        }
    }

}

rtree_node *rtree_split_gut_quad(rtree_root *root, rtree_node *n, rtree_node *newn, rtree_node *parent, int m) {

    int seeds[2];
	rtree_node *g1 = rtree_new_node(root, DIRECTORY);
	rtree_node *new_dir = rtree_new_node(root, DIRECTORY);
	rtree_node *g2 = new_dir;

	n->dirs[n->used] = newn;
	n->used++;

    //PICK SEEDS
    pick_seeds(n, seeds, m);

    //ADICIONA AS SEMENTES AOS NOVOS NÓS
    g1->dirs[g1->used] = n->dirs[seeds[0]];
    g2->dirs[g2->used] = n->dirs[seeds[1]];

    //IMPRIMI OS NOVOS NÓS
//    printf("--------- NOVOS NÓS ----------\n");
//    printf("[%d] - (%.0f, %.0f) (%.0f, %.0f) \n", g1->used + 1, g1->dirs[g1->used]->mbr.MinX, g1->dirs[g1->used]->mbr.MinY, g1->dirs[g1->used]->mbr.MaxX, g1->dirs[g1->used]->mbr.MaxY);
//    printf("[%d] - (%.0f, %.0f) (%.0f, %.0f) \n", g2->used + 1, g2->dirs[g2->used]->mbr.MinX, g2->dirs[g2->used]->mbr.MinY, g2->dirs[g2->used]->mbr.MaxX, g2->dirs[g2->used]->mbr.MaxY);
//    printf("\n");
    g1->used++;
    g2->used++;

    //PICK NEXT
    pick_next(n, g1, g2, seeds, m);

	// finalization
	for(int i = 0; i < g1->used; i++) {
		n->dirs[i] = g1->dirs[i];
	}
	n->used = g1->used;

	n->mbr = rtree_compute_mbr(n);
	g2->mbr = rtree_compute_mbr(g2);


	g_free(g1->dirs);
	g_free(g1);

	return new_dir;
}

void pick_seeds_leaf(rtree_node *n, int seeds[], int m) {
    int no1 = -1;
    int no2 = -2;
    float menorX = 0.0;
    float menorY = 0.0;
    float maiorX = 0.0;
    float maiorY = 0.0;
    float calculo = 0;
    float dist_maior = 0.0;
    int E1 = 0;
    int E2 = 0;
    int indiceI = 0;
    int indiceJ = 0;


    //IMPRESSÃO DOS NÓS 
//    printf("NÓ CHEIO!\n");
    for (int i = 0; i < m; i++) {
//        printf("[%d] - (%.0f, %.0f) (%.0f, %.0f) \n", i + 1, n->leaves[i].mbr.MinX, n->leaves[i].mbr.MinY, n->leaves[i].mbr.MaxX, n->leaves[i].mbr.MaxY);
    }
    
    
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            
         //ENCONTRA A ÁREA DE E1 E E2
          E1 =  (n->leaves[i].mbr.MaxX - n->leaves[i].mbr.MinX) * (n->leaves[i].mbr.MaxY - n->leaves[i].mbr.MinY);
          E2 =  (n->leaves[j].mbr.MaxX - n->leaves[j].mbr.MinX) * (n->leaves[j].mbr.MaxY - n->leaves[j].mbr.MinY);
           
       
          
          //ENCONTRA A ÁREA TOTAL DE (J)
          
          //ENCONTRA O MENOR X E MENOR Y
          if (n->leaves[i].mbr.MinX <= n->leaves[j].mbr.MinX){
              menorX = n->leaves[i].mbr.MinX;
          }else {
              menorX = n->leaves[j].mbr.MinX;
          }
          
          if (n->leaves[i].mbr.MinY <= n->leaves[j].mbr.MinY){
              menorY = n->leaves[i].mbr.MinY;
          }else {
              menorY = n->leaves[j].mbr.MinY;
          }
          
          //ENCONTRA O MAIOR X E MAIOR Y          
          if (n->leaves[i].mbr.MaxX >= n->leaves[j].mbr.MaxX){
              maiorX = n->leaves[i].mbr.MaxX;
          }else {
              maiorX = n->leaves[j].mbr.MaxX;
          }
          
          if (n->leaves[i].mbr.MaxY >= n->leaves[j].mbr.MaxY){
              maiorY = n->leaves[i].mbr.MaxY;
          }else {
              maiorY = n->leaves[j].mbr.MaxY;
          }
              
          
          //CALCULO DA ÁREA TOTAL MENOS A AREA DOS E1 E E2
          calculo = ((maiorX - menorX) * (maiorY - menorY)) - E1 - E2;
          
          if (calculo >= dist_maior){
              dist_maior = calculo;
              indiceI = i;
              indiceJ = j;
          }
                   
        }
    }

    seeds[0] = indiceI;
    seeds[1] = indiceJ;
    no1 = indiceI;
    no2 = indiceJ;

//    printf("\n---- SEMENTES ENCONTRADAS ----\n");
//    printf("[%d] - (%.0f, %.0f) (%.0f, %.0f) \n", no2 + 1, n->leaves[no2].mbr.MinX, n->leaves[no2].mbr.MinY, n->leaves[no2].mbr.MaxX, n->leaves[no2].mbr.MaxY);
//    printf("[%d] - (%.0f, %.0f) (%.0f, %.0f) \n", no1 + 1, n->leaves[no1].mbr.MinX, n->leaves[no1].mbr.MinY, n->leaves[no1].mbr.MaxX, n->leaves[no1].mbr.MaxY);

}

void pick_next_leaf(rtree_node *n, rtree_node *g1, rtree_node *g2, int seeds[], int m) {
    double minX1 = 0.0;
    double maxX1 = 0.0;
    double minY1 = 0.0;
    double maxY1 = 0.0;

    double minX2 = 0.0;
    double maxX2 = 0.0;
    double minY2 = 0.0;
    double maxY2 = 0.0;

    double areaG1 = 0.0;
    double areaG2 = 0.0;
    double englobarG1 = 0.0;
    double englobarG2 = 0.0;
    double diferencaG1 = 0.0;
    double diferencaG2 = 0.0;

    Envelope mbr1;
    Envelope mbr2;

    mbr1.MinX = n->leaves[seeds[0]].mbr.MinX;
    mbr1.MinY = n->leaves[seeds[0]].mbr.MinY;
    mbr1.MaxX = n->leaves[seeds[0]].mbr.MaxX;
    mbr1.MaxY = n->leaves[seeds[0]].mbr.MaxY;

    mbr2.MinX = n->leaves[seeds[1]].mbr.MinX;
    mbr2.MinY = n->leaves[seeds[1]].mbr.MinY;
    mbr2.MaxX = n->leaves[seeds[1]].mbr.MaxX;
    mbr2.MaxY = n->leaves[seeds[1]].mbr.MaxY;

    areaG1 = ENVELOPE_AREA(mbr1);
    areaG2 = ENVELOPE_AREA(mbr2);
//    printf("AREA INICIAL DE G1: %f \n", areaG1);
//    printf("AREA INICIAL DE G2: %f \n\n", areaG2);

    for (int i = 0; i < n->used; ++i) {
        if (!((seeds[0] == i) || (seeds[1] == i))) {
            minX1 = 0.0;
            maxX1 = 0.0;
            minY1 = 0.0;
            maxY1 = 0.0;
            minX2 = 0.0;
            maxX2 = 0.0;
            minY2 = 0.0;
            maxY2 = 0.0;

            //----VERIFICACAO PARA O PRIMEIRO NÓ--------
            //ENCONTRAR X MENOR E MAIOR
            if (n->leaves[i].mbr.MinX <= mbr1.MinX) {
                minX1 = n->leaves[i].mbr.MinX;
            } else {
                minX1 = mbr1.MinX;
            }

            if (n->leaves[i].mbr.MaxX >= mbr1.MaxX) {
                maxX1 = n->leaves[i].mbr.MaxX;
            } else {
                maxX1 = mbr1.MaxX;
            }

            //ENCONTRA Y MENOR E MAIOR
            if (n->leaves[i].mbr.MinY <= mbr1.MinY) {
                minY1 = n->leaves[i].mbr.MinY;
            } else {
                minY1 = mbr1.MinY;
            }

            if (n->leaves[i].mbr.MaxY >= mbr1.MaxY) {
                maxY1 = n->leaves[i].mbr.MaxY;
            } else {
                maxY1 = mbr1.MaxY;
            }

            englobarG1 = (maxX1 - minX1) * (maxY1 - minY1);
            diferencaG1 = englobarG1 - areaG1;

//            printf("Englobar G1: %f \n", englobarG1);
//            printf("DIFERENCA DE G1: %f \n\n", diferencaG1);


            //----- VERIFICAÇÃO PARA O SEGUNDO NÓ ------
            //ENCONTRAR X MENOR E MAIOR
            if (n->leaves[i].mbr.MinX <= mbr2.MinX) {
                minX2 = n->leaves[i].mbr.MinX;
            } else {
                minX2 = mbr2.MinX;
            }

            if (n->leaves[i].mbr.MaxX >= mbr2.MaxX) {
                maxX2 = n->leaves[i].mbr.MaxX;
            } else {
                maxX2 = mbr2.MaxX;
            }


            //ENCONTRA Y MENOR E MAIOR
            if (n->leaves[i].mbr.MinY <= mbr2.MinY) {
                minY2 = n->leaves[i].mbr.MinY;
            } else {
                minY2 = mbr2.MinY;
            }

            if (n->leaves[i].mbr.MaxY >= mbr2.MaxY) {
                maxY2 = n->leaves[i].mbr.MaxY;
            } else {
                maxY2 = mbr2.MaxY;
            }


            englobarG2 = (maxX2 - minX2) * (maxY2 - minY2);
            diferencaG2 = englobarG2 - areaG2;

//            printf("Englobar G2: %f \n", englobarG2);
//            printf("DIFERENCA DE G2: %f \n\n", diferencaG2);



            
            

            if (diferencaG1 <= diferencaG2) {
//                printf("G1 é menor: \n\n");
                
                mbr1.MinX = minX1;
                mbr1.MinY = minY1;
                mbr1.MaxX = maxX1;
                mbr1.MaxY = maxY1;
                areaG1 = (mbr1.MaxX - mbr1.MinX) * (mbr1.MaxY - mbr1.MinY);
                minX1 = 0.0;
                maxX1 = 0.0;
                minY1 = 0.0;
                maxY1 = 0.0;  
                
                //ADICIONA O NÓ AO G1, VERIFICANDO m
                if((g1->used == minm(m)) && (g2->used < minm(m))){
                    g2->leaves[g2->used] = n->leaves[i];
                    g2->used++;  
                }else {
                    g1->leaves[g1->used] = n->leaves[i];
                    g1->used++;
                }
                

            } else {
//                printf("G2 é menor: \n\n");
                
                mbr2.MinX = minX2;
                mbr2.MinY = minY2;
                mbr2.MaxX = maxX2;
                mbr2.MaxY = maxY2;
                areaG2 = (mbr2.MaxX - mbr2.MinX) * (mbr2.MaxY - mbr2.MinY);
                minX2 = 0.0;
                maxX2 = 0.0;
                minY2 = 0.0;
                maxY2 = 0.0;
                
                //ADICIONA O NÓ AO G2, VERIFICANDO m
                if((g2->used == minm(m)) && (g1->used < minm(m))){
                    g1->leaves[g1->used] = n->leaves[i];
                    g1->used++;
                }else {
                    g2->leaves[g2->used] = n->leaves[i];
                    g2->used++;
                }
                

            }
                        
            englobarG1 = 0.0;
            englobarG2 = 0.0;

            
            
        }
    }
}

rtree_node *rtree_split_gut_quad_leaf(rtree_root *root, rtree_node *n, rtree_node *parent, const GEOSGeometryH ngeo, const Envelope e) {

    int seeds[2];
	rtree_node *g1 = rtree_new_node(root, LEAF);
	rtree_node *new_leaf = rtree_new_node(root, LEAF);
	rtree_node *g2 = new_leaf;
	
	n->leaves[n->used].geo = ngeo;
	n->leaves[n->used].mbr = e;
	n->used++;

    //PICK SEEDS
    pick_seeds_leaf(n, seeds, root->m);

    //ADICIONA AS SEMENTES AOS NOVOS NÓS
    g1->leaves[g1->used] = n->leaves[seeds[0]];
    g2->leaves[g2->used] = n->leaves[seeds[1]];

    //IMPRIMI OS NOVOS NÓS
//    printf("--------- NOVOS NÓS ----------\n");
//    printf("[%d] - (%.0f, %.0f) (%.0f, %.0f) \n", g1->used + 1, g1->leaves[g1->used].mbr.MinX, g1->leaves[g1->used].mbr.MinY, g1->leaves[g1->used].mbr.MaxX, g1->leaves[g1->used].mbr.MaxY);
//    printf("[%d] - (%.0f, %.0f) (%.0f, %.0f) \n", g2->used + 1, g2->leaves[g2->used].mbr.MinX, g2->leaves[g2->used].mbr.MinY, g2->leaves[g2->used].mbr.MaxX, g2->leaves[g2->used].mbr.MaxY);
//    printf("\n");
    g1->used++;
    g2->used++;

    //PICK NEXT
    pick_next_leaf(n, g1, g2, seeds, root->m);

	// finalization
	for(int i = 0; i < g1->used; i++) {
		n->leaves[i] = g1->leaves[i];
	}
	n->used = g1->used;

	n->mbr = rtree_compute_mbr(n);
	g2->mbr = rtree_compute_mbr(g2);

	g_free(g1->leaves);
	g_free(g1);

	return new_leaf;
}



rtree_root *rtree_new_rtree_gut_quad(const int m, const int servers) {
	return rtree_new(
		m,
		servers,
		rtree_choose_subtree_gut,
		rtree_split_gut_quad,
		rtree_split_gut_quad_leaf,
		NULL, NULL, NULL);
}




#include <float.h>
#include "rtree-star.h"
#include "rtree-reinsert.h"


int rtree_choose_subtree_linear(rtree_root *root, rtree_node *n, const Envelope e) {
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


void pick_seeds_linear(rtree_node *n, int seeds[], int m) {
    double ultX = 0.0;
    double ultY = 0.0;
    int minX = 0, maxX = 0;
    int minY = 0, maxY = 0;
    double distX = 0.0;
    double distY = 0.0;
    int no1 = -1;
    int no2 = -2;
    float normX = 0.0;
    float normY = 0.0;




    //BUSCA OBJETOS MAIS DISTANTES EM X
    for (int i = 1; i < m; ++i) {
        //ENCONTRA ALCANCE DE X E Y 
        if (n->dirs[i]->mbr.MaxX > ultX) {
            ultX = n->dirs[i]->mbr.MaxX;
        }
        if (n->dirs[i]->mbr.MaxY > ultY) {
            ultY = n->dirs[i]->mbr.MaxY;
        }
        //ECNONTRA OBJETOS MAIS DISTANTES EM X E Y
        if (n->dirs[i]->mbr.MaxX < n->dirs[minX]->mbr.MaxX) {
            minX = i;
        }
        if (n->dirs[i]->mbr.MinX > n->dirs[maxX]->mbr.MinX) {
            maxX = i;
        }
        if (n->dirs[i]->mbr.MaxY < n->dirs[minY]->mbr.MaxY) {
            minY = i;
        }
        if (n->dirs[i]->mbr.MinY > n->dirs[maxY]->mbr.MinY) {
            maxY = i;
        }
    }

    
    
    
    

    //DEFINE O VALOR DA SEPARAÇÃO
    distX = n->dirs[maxX]->mbr.MinX - n->dirs[minX]->mbr.MaxX;
    distY = n->dirs[maxY]->mbr.MinY - n->dirs[minY]->mbr.MaxY;

    //NORMALIZA O VALOR DA SEPARAÇÃO
    normX = (float) distX / (float) ultX;
    normY = (float) distY / (float) ultY;

    //IDENTIFICA A DIMENSAO DO OBJETO MAIS DISTANTE
    if (normX > normY) {
        no1 = maxX;
        no2 = minX;
    } else {
        no1 = maxY;
        no2 = minY;
    }

    seeds[0] = no1;
    seeds[1] = no2;
}

void pick_next_linear(rtree_node *n, rtree_node *g1, rtree_node *g2, int seeds[], int m) {
        for(int i = 0; i < m; i++){
            if(!(i == seeds[0]) || !(i == seeds[1])){
                if(g1->used <= minm(m)){
                    g1->dirs[g1->used] = n->dirs[i];
                    g1->used++;
                }else{
                    g2->dirs[g2->used] = n->dirs[i];
                    g2->used++;
                }
            }        
        }

}

rtree_node *rtree_split_gut_linear(rtree_root *root, rtree_node *n, rtree_node *newn, rtree_node *parent, int m) {

    int seeds[2];
	rtree_node *g1 = rtree_new_node(root, DIRECTORY);
	rtree_node *new_dir = rtree_new_node(root, DIRECTORY);
	rtree_node *g2 = new_dir;

	n->dirs[n->used] = newn;
	n->used++;

    //PICK SEEDS
    pick_seeds_linear(n, seeds, m);

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
    pick_next_linear(n, g1, g2, seeds, m);

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

void pick_seeds_leaf_linear(rtree_node *n, int seeds[], int m) {
    double ultX = 0.0;
    double ultY = 0.0;
    int minX = 0, maxX = 0;
    int minY = 0, maxY = 0;
    double distX = 0.0;
    double distY = 0.0;
    int no1 = -1;
    int no2 = -2;
    float normX = 0.0;
    float normY = 0.0;



    //BUSCA OBJETOS MAIS DISTANTES EM X
    for (int i = 1; i < m; ++i) {
        //ENCONTRA ALCANCE DE X E Y 
        if (n->leaves[i].mbr.MaxX > ultX) {
            ultX = n->leaves[i].mbr.MaxX;
        }
        if (n->leaves[i].mbr.MaxY > ultY) {
            ultY = n->leaves[i].mbr.MaxY;
        }
        //ECNONTRA OBJETOS MAIS DISTANTES EM X E Y
        if (n->leaves[i].mbr.MaxX < n->leaves[minX].mbr.MaxX) {
            minX = i;
        }
        if (n->leaves[i].mbr.MinX > n->leaves[maxX].mbr.MinX) {
            maxX = i;
        }
        if (n->leaves[i].mbr.MaxY < n->leaves[minY].mbr.MaxY) {
            minY = i;
        }
        if (n->leaves[i].mbr.MinY > n->leaves[maxY].mbr.MinY) {
            maxY = i;
        }
    }

    
    
    
    

    //DEFINE O VALOR DA SEPARAÇÃO
    distX = n->leaves[maxX].mbr.MinX - n->leaves[minX].mbr.MaxX;
    distY = n->leaves[maxY].mbr.MinY - n->leaves[minY].mbr.MaxY;

    //NORMALIZA O VALOR DA SEPARAÇÃO
    normX = (float) distX / (float) ultX;
    normY = (float) distY / (float) ultY;

    //IDENTIFICA A DIMENSAO DO OBJETO MAIS DISTANTE
    if (normX > normY) {
        no1 = maxX;
        no2 = minX;
    } else {
        no1 = maxY;
        no2 = minY;
    }

    seeds[0] = no1;
    seeds[1] = no2;

}

void pick_next_leaf_linear(rtree_node *n, rtree_node *g1, rtree_node *g2, int seeds[], int m) {
   
        for(int i = 0; i < m; i++){
            if(!(i == seeds[0]) || !(i == seeds[1])){
                if(g1->used <= minm(m)){
                    g1->leaves[g1->used] = n->leaves[i];
                    g1->used++;
                }else{
                    g2->leaves[g2->used] = n->leaves[i];
                    g2->used++;
                }
            }        
        }
}

rtree_node *rtree_split_gut_linear_leaf(rtree_root *root, rtree_node *n, rtree_node *parent, const GEOSGeometryH ngeo, const Envelope e) {

    int seeds[2];
	rtree_node *g1 = rtree_new_node(root, LEAF);
	rtree_node *new_leaf = rtree_new_node(root, LEAF);
	rtree_node *g2 = new_leaf;
	
	n->leaves[n->used].geo = ngeo;
	n->leaves[n->used].mbr = e;
	n->used++;

    //PICK SEEDS
    pick_seeds_leaf_linear(n, seeds, root->m);

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
    pick_next_leaf_linear(n, g1, g2, seeds, root->m);

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



rtree_root *rtree_new_rtree_gut_linear(const int m, const int servers) {
	return rtree_new(
		m,
		servers,
		rtree_choose_subtree_linear,
		rtree_split_gut_linear,
		rtree_split_gut_linear_leaf,
		NULL, NULL, NULL);
}



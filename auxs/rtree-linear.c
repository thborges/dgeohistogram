
#include <float.h>
#include "rtree-gut.h"
#include "rtree-reinsert.h"


void pick_seeds_linear(rtree_node *n, int seeds[]) {
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
    for (int i = 0; i < n->used; ++i) {
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
	assert(no1 != no2);
	assert(no1 >= 0 && no1 < n->used);
	assert(no2 >= 0 && no2 < n->used);

}

rtree_node *rtree_split_gut_linear(rtree_root *root, rtree_node *n, rtree_node *newn, rtree_node *parent, int m) {

    int seeds[2];
	rtree_node *g1 = rtree_new_node(root, DIRECTORY);
	rtree_node *new_dir = rtree_new_node(root, DIRECTORY);
	rtree_node *g2 = new_dir;

	n->dirs[n->used++] = newn;

    //PICK SEEDS
    pick_seeds_linear(n, seeds);

    //ADICIONA AS SEMENTES AOS NOVOS NÓS
    g1->dirs[g1->used++] = n->dirs[seeds[0]];
    g2->dirs[g2->used++] = n->dirs[seeds[1]];

    //PICK NEXT
	Envelope env_s0 = n->dirs[seeds[0]]->mbr;
	double area_s0 = ENVELOPE_AREA(env_s0);
   	Envelope env_s1 = n->dirs[seeds[1]]->mbr;
	double area_s1 = ENVELOPE_AREA(env_s1);

	int minimum = minm(root->m);
	int remaining = n->used - 2;
	for(int i = 0; i < n->used; i++){
		if (i == seeds[0]) continue;
		if (i == seeds[1]) continue;

		if (g1->used < minimum && remaining <= minimum)
			g1->dirs[g1->used++] = n->dirs[i];
		else
		if (g2->used < minimum && remaining <= minimum)
			g2->dirs[g2->used++] = n->dirs[i];
		else {
			Envelope eaux = env_s0;
			ENVELOPE_MERGE(eaux, n->dirs[i]->mbr);
			double enlarg_s0 = ENVELOPE_AREA(eaux) / area_s0;

			eaux = env_s1;
			ENVELOPE_MERGE(eaux, n->dirs[i]->mbr);
			double enlarg_s1 = ENVELOPE_AREA(eaux) / area_s1;

			if (enlarg_s0 <= enlarg_s1)
            	g1->dirs[g1->used++] = n->dirs[i];
        	else
	            g2->dirs[g2->used++] = n->dirs[i];
        }

		remaining--;
    }

	// finalization
	n->used = 0;
	for(int i = 0; i < g1->used; i++) {
		n->dirs[n->used++] = g1->dirs[i];
	}

	n->mbr = rtree_compute_mbr(n);
	g2->mbr = rtree_compute_mbr(g2);

	assert(n->used + g2->used == m+1);

	g_free(g1->dirs);
	g_free(g1);

	return new_dir;
}

void pick_seeds_leaf_linear(rtree_node *n, int seeds[]) {
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
    for (int i = 0; i < n->used; ++i) {
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
	assert(no1 != no2);
	assert(no1 >= 0 && no1 < n->used);
	assert(no2 >= 0 && no2 < n->used);
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
    pick_seeds_leaf_linear(n, seeds);

    //ADICIONA AS SEMENTES AOS NOVOS NÓS
    g1->leaves[g1->used++] = n->leaves[seeds[0]];
    g2->leaves[g2->used++] = n->leaves[seeds[1]];

    //PICK NEXT
	Envelope env_s0 = n->leaves[seeds[0]].mbr;
	double area_s0 = ENVELOPE_AREA(env_s0);
   	Envelope env_s1 = n->leaves[seeds[1]].mbr;
	double area_s1 = ENVELOPE_AREA(env_s1);

	int minimum = minm(root->m);
	int remaining = n->used - 2;
	for(int i = 0; i < n->used; i++){
		if (i == seeds[0]) continue;
		if (i == seeds[1]) continue;

		if (g1->used < minimum && remaining <= minimum)
			g1->leaves[g1->used++] = n->leaves[i];
		else
		if (g2->used < minimum && remaining <= minimum)
			g2->leaves[g2->used++] = n->leaves[i];
		else {
			Envelope eaux = env_s0;
			ENVELOPE_MERGE(eaux, n->leaves[i].mbr);
			double enlarg_s0 = ENVELOPE_AREA(eaux) / area_s0;

			eaux = env_s1;
			ENVELOPE_MERGE(eaux, n->leaves[i].mbr);
			double enlarg_s1 = ENVELOPE_AREA(eaux) / area_s1;

			if (enlarg_s0 <= enlarg_s1)
            	g1->leaves[g1->used++] = n->leaves[i];
        	else
	            g2->leaves[g2->used++] = n->leaves[i];
        }

		remaining--;
    }

	// finalization
	for(int i = 0; i < g1->used; i++) {
		n->leaves[i] = g1->leaves[i];
	}
	n->used = g1->used;

	n->mbr = rtree_compute_mbr(n);
	g2->mbr = rtree_compute_mbr(g2);

	if (!(n->used + g2->used == root->m+1))
		printf("N: %d, G2: %d, m=%d\n", n->used, g2->used, root->m+1);
	assert(n->used + g2->used == root->m+1);

	g_free(g1->leaves);
	g_free(g1);

	return new_leaf;
}



rtree_root *rtree_new_rtree_gut_linear(const int m, const int servers) {
	return rtree_new(
		m,
		servers,
		rtree_choose_subtree_gut,
		rtree_split_gut_linear,
		rtree_split_gut_linear_leaf,
		NULL, NULL, NULL);
}



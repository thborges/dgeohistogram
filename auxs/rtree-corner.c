

#include <float.h>
#include "rtree.h"
#include "rtree-gut.h"

rtree_node *rtree_split_corner(rtree_root *root, rtree_node *n, rtree_node *newn, rtree_node *parent, int m) {

	rtree_node *g1 = rtree_new_node(root, LEAF);
	rtree_node *new_dir= rtree_new_node(root, LEAF);
	rtree_node *g2 = new_dir;
	
	n->dirs[n->used] = newn;
	n->used++;

    int ultX = 0;
    int ultY = 0;
    int primX = 0;
    int primY = 0;
    float CoverRectX = 0.0;
    float CoverRectY = 0.0;
    float ObjX = 0.0;
    float ObjY = 0.0;
    int C0 = 0;
    int C1 = 0;
    int C2 = 0;
    int C3 = 0;

	typedef struct {
		rtree_node **items;
		int used;
	} node_aux;

	node_aux c0, c1, c2, c3;	
	c0.items = g_new(rtree_node*, m);
   	c1.items = g_new(rtree_node*, m);
	c2.items = g_new(rtree_node*, m);
	c3.items = g_new(rtree_node*, m);
	c0.used = 0;
    c1.used = 0;
    c2.used = 0;
    c3.used = 0;


    //IMPRESSÃO DOS NÓS 
    printf("NÓ CHEIO!\n");
    for (int i = 0; i < n->used; i++) {
        //printf("[%d] - (%.0f, %.0f) (%.0f, %.0f) \n", i + 1, n->dirs[i]->mbr.MinX, n->dirs[i]->mbr.MinY, n->dirs[i]->mbr.MaxX, n->dirs[i]->mbr.MaxY);
    }


    primX = n->dirs[0]->mbr.MinX;
    primY = n->dirs[0]->mbr.MinY;

    for (int i = 1; i < n->used; ++i) {
        //ENCONTRA ALCANCE MAX DE X E Y 
        if (n->dirs[i]->mbr.MaxX > ultX) {
            ultX = n->dirs[i]->mbr.MaxX;
        }
        if (n->dirs[i]->mbr.MaxY > ultY) {
            ultY = n->dirs[i]->mbr.MaxY;
        }

        //ENCONTRA ALCANCE MIN DE X E Y 
        if (n->dirs[i]->mbr.MinX <= primX) {
            primX = n->dirs[i]->mbr.MinX;
        }
        if (n->dirs[i]->mbr.MinY <= primY) {
            primY = n->dirs[i]->mbr.MinY;
        }

    }
    // define o centro da MBR
    CoverRectX = ((ultX + primX) / 2.0);
    CoverRectY = ((ultY + primY) / 2.0);


    //ENCONTRA O CENTRO DE CADA OBJETO e INCREMENTA O CONTADOR DOS CANTOS;
    for (int i = 0; i < n->used; ++i) {
        //ENCONTRA O CENTRO DO OBJETO
        ObjX = ((n->dirs[i]->mbr.MaxX + n->dirs[i]->mbr.MinX) / 2.0);
        ObjY = ((n->dirs[i]->mbr.MaxY + n->dirs[i]->mbr.MinY) / 2.0);


        //DEFINE O CANTO MAIS PROXIMO E INCREMENTA CONTADOR

        if (ObjX > CoverRectX) {
            if (ObjY > CoverRectY) {
                C2++;
                c2.items[c2.used] = n->dirs[i];
                c2.used++;
            } else {
                C3++;
                c3.items[c3.used] = n->dirs[i];
                c3.used++;
            }

        } else {
            if (ObjY > CoverRectY) {
                C1++;
                c1.items[c1.used] = n->dirs[i];
                c1.used++;
            } else {
                C0++;
                c0.items[c0.used] = n->dirs[i];
                c0.used++;
            }
        }
    }


    // DISTRIBUIR OS OBJETOS EM SEUS NOVOS GRUPOS

    if (C0 > C2) {
        //move C0 para N1 e C2 para N2
        for (int i = 0; i < c0.used; i++) {
            g1->dirs[g1->used] = c0.items[i];
            g1->used++;
        }
        for (int i = 0; i < c2.used; i++) {
            g2->dirs[g2->used] = c2.items[i];
            g2->used++;
        }

    } else {
        //move C2 para N1 e C0 para N2
        for (int i = 0; i < c2.used; i++) {
            g1->dirs[g1->used] = c2.items[i];
            g1->used++;
        }
        for (int i = 0; i < c0.used; i++) {
            g2->dirs[g2->used] = c0.items[i];
            g2->used++;
        }

    }

    if (C1 > C3) {
        //move C1 para N2 e C3 para N1
        for (int i = 0; i < c1.used; i++) {
            g2->dirs[g2->used] = c1.items[i];
            g2->used++;
        }
        for (int i = 0; i < c3.used; i++) {
            g1->dirs[g1->used] = c3.items[i];
            g1->used++;
        }
    } else {
        //move C3 para N2 e C1 para N1 
        if (C3 > C1) {
            for (int i = 0; i < c3.used; i++) {
                g2->dirs[g2->used] = c3.items[i];
                g2->used++;
            }
            for (int i = 0; i < c1.used; i++) {
                g1->dirs[g1->used] = c1.items[i];
                g1->used++;
            }
        } 
    }

	// finalization
	for(int i = 0; i < g1->used; i++) {
		n->dirs[i] = g1->dirs[i];
	}
	n->used = g1->used;

	n->mbr = rtree_compute_mbr(n);
	g2->mbr = rtree_compute_mbr(g2);

	g_free(c0.items);
	g_free(c1.items);
	g_free(c2.items);
	g_free(c3.items);

	g_free(g1->dirs);
	g_free(g1);

	return new_dir;
}


rtree_node *rtree_split_corner_leaf(rtree_root *root, rtree_node *n, rtree_node *parent, const GEOSGeometryH ngeo, const Envelope e) {

}

rtree_root *rtree_new_rtree_corner(const int m, const int servers) {
	return rtree_new(
		m,
		servers,
		rtree_choose_subtree_gut,
		rtree_split_corner,
		rtree_split_corner_leaf,
		NULL, NULL, NULL);
}



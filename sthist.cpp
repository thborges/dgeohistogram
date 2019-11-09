
#include "sthist.h"
#include<algorithm>

int computeNBi(STHist &hotSpotTree, int bRoot, int iSon);
void constructionHotSpotsTree(STHist &hotSpotTree, int bRoot, int NB);
int detectedHotSpot(STHist& hotSpotTree, int bRoot, int s, int f);
bool yCoordinateAscending(coordenada a, coordenada b);
bool xCoordinateAscending(coordenada a, coordenada b);
int limitIndexY_binarySearch(const vector<coordenada>& vetor, int begin, int end, double value);

void STHist_generate(dataset* ds, STHist& toFillOut, int buckets_num){
	STBucket root;
	Envelope root_mbr;
	unsigned int cardin=0;

	dataset_iter di;
	dataset_foreach(di, ds) {
		root_mbr = di.item->mbr;
		break;
	}

	dataset_foreach(di, ds) {
		Envelope mbr = di.item->mbr;
		coordenada objeto;
		objeto.x = mbr.MinX;
		objeto.y = mbr.MinY;
		root.O.push_front(objeto);
		cardin++;
		if(mbr.MinX < root_mbr.MinX)  root_mbr.MinX = mbr.MinX;
		if(mbr.MaxX > root_mbr.MaxX)  root_mbr.MaxX = mbr.MaxX;
		if(mbr.MinY < root_mbr.MinY)  root_mbr.MinY = mbr.MinY;
		if(mbr.MaxY > root_mbr.MaxY)  root_mbr.MaxY = mbr.MaxY;
	}
	root.mbr = root_mbr;
	root.F = cardin;

	toFillOut.push_back(root);

	constructionHotSpotsTree(toFillOut,0,buckets_num);

}

int computeNBi(STHist &hotSpotTree, int bRoot, int iSon){
	return 1;
}

void constructionHotSpotsTree(STHist &hotSpotTree, int bRoot, int NB){
	// B : hotSpotTree[i].mbr
	// O : (hotSpotTree[i].Ox, hotSpotTree[i].Oy)
	// NB : NB

	int s = max(NB, 2);
	int f = max(NB, 2);


	int qtd_hotSpot = detectedHotSpot(hotSpotTree, bRoot, s, f);

	for(int i =0; i < qtd_hotSpot; i++){
		hotSpotTree[bRoot].childrean.push_back(bRoot+i);
	}

	for(int i =0; i < hotSpotTree[bRoot].childrean.size(); i++){
		int NBi = computeNBi(hotSpotTree, bRoot, i);
		constructionHotSpotsTree(hotSpotTree, i, NBi);
	}
}




int detectedHotSpot(STHist& hotSpotTree, int bRoot, int s, int f){
	// Entrada:
		// D : hotSpotTree[i].mbr
		// O : (hotSpotTree[i].Ox, hotSpotTree[i].Oy)
		// s : s
		// f : f
	// SAIDA:
		//hotSpotsCount : quantidade de hotSpots encontrados ( que ja são adicionados a arvore de hotSpots dentro dessa propria função).


	int hotSpotsCount = 0;

	forward_list<coordenada> SortedListX(hotSpotTree[bRoot].O.begin(), hotSpotTree[bRoot].O.end());
	forward_list<coordenada> SortedListY(hotSpotTree[bRoot].O.begin(), hotSpotTree[bRoot].O.end());

	SortedListX.sort(xCoordinateAscending);
	SortedListY.sort(yCoordinateAscending);

/*1:*/

/*2:*/	Envelope D = hotSpotTree[bRoot].mbr;
	unsigned int F = hotSpotTree[bRoot].F;

/*3:*/	float w =  (1/pow(s,0.5))* (D.MaxX - D.MinY);
/*4:*/	float h =  (1/pow(s,0.5))* (D.MaxY - D.MinY);

/*5, 6:*/
	for(forward_list<coordenada>::iterator oi= SortedListX.begin(); oi != SortedListX.end(); oi++){
		double W = oi->x + w;
		vector<coordenada> SortedCandiListY;

		for(forward_list<coordenada>::iterator it = oi; it != SortedListX.end(); it++){
			if(it->x <= W){  // verificar no artigo se menor ou menor igual
				SortedCandiListY.push_back(*it);
			}else{
				break;
			}
/*8:*/		}

		unsigned int freqW = SortedCandiListY.size();
/*7:*/		if(freqW >= F/f){
/*8:*/			sort(SortedCandiListY.begin(), SortedCandiListY.end(),yCoordinateAscending);
			double H = SortedCandiListY[0].y + h;

/*9, 10:*/		for(unsigned int oj=0; oj < freqW; oj++){
				unsigned int freqH = -oj + 1 + limitIndexY_binarySearch(SortedCandiListY,oj,SortedCandiListY.size(), H);
				if(freqH >= F/f){
/*12*/
					hotSpot R;
					R.mbr.MinX = oi->x;
					R.mbr.MaxX = oi->x + W;
					R.mbr.MinY = SortedCandiListY[oj].y;
					R.mbr.MaxY =SortedCandiListY[oj].y + H;
					R.F = freqH;
					R.O.insert_after(R.O.begin(),SortedCandiListY.begin()+oj, SortedCandiListY.begin()+oj+freqH);

/*13*/					hotSpotTree.push_back(R);
					hotSpotsCount++;


					auto oY = SortedListY.begin();
					for(unsigned int oC = 0; oC < SortedCandiListY.size(); oC++){
						while(oY != SortedListY.end()){
							if((*oY).x == SortedCandiListY[oC].x  && (*oY).y == SortedCandiListY[oC].y){
								oY = SortedListY.erase_after(oY);
							}else{
								oY++;
							}
						}
					}

					sort(SortedCandiListY.begin(), SortedCandiListY.end(),xCoordinateAscending);

					auto oX = SortedListX.begin();
					for(unsigned int oC = 0; oC < SortedCandiListY.size(); oC++){
						while(oX != SortedListX.end()){
							if((*oX).x == SortedCandiListY[oC].x  && (*oX).y == SortedCandiListY[oC].y){
								oX = SortedListX.erase_after(oX);
							}else{
								oX++;
							}
						}
					}
				}
			}
		}
	}
	return hotSpotsCount;
}


bool yCoordinateAscending(coordenada a, coordenada b){
	return a.y < b.y;
}

bool xCoordinateAscending(coordenada a, coordenada b){
	return a.x < b.x;
}

int limitIndexY_binarySearch(const vector<coordenada>& vetor, int begin, int end, double value){
	int meio = (end-begin)/2 +begin;
	while(end-begin >3){
		if(value < vetor[meio].y){
			end = meio;
		}else{
			begin = meio;
		}
		meio = (end - begin)/2;
	}
	int x = begin;
	while(x < end) {
		if(vetor[x+1].y > value )
			return x;
	}
	return x;
}

double STHist_search(STHist &hist, Envelope query);

void STHist_print(STHist &hist);

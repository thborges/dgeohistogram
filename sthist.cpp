
#include "sthist.h"
#include<algorithm>

STHist *STHist_generate(dataset *ds, int buckets_num);
int computeNBi(STHist &hotSpotTree, int bRoot, int iSon);
void constructionHotSpotsTree(STHist &hotSpotTree, int bRoot, int NB);
void detectedHotSpot(STHist &hotSpotTree, int bRoot, int s, int f);
bool yCoordinateAscending(coordenada a, coordenada b);
bool xCoordinateAscending(coordenada a, coordenada b);
int limitIndex_binarySearch(const vector<coordenada> &vetor, int begin, int end, double value);

STHist* STHist_generate(dataset* ds, int buckets_num){
	STBucket root;
	Envelope root_mbr;
	unsigned int cardin=0;
	double Ox, Oy;

	dataset_iter di;
	dataset_foreach(di, ds) {
		coordenada objeto;
		objeto.x = di.item->mbr.MinX;
		objeto.y = di.item->mbr.MinY;
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

	STHist hist;
	hist.push_back(root);

	constructionHotSpotsTree(hist,0,buckets_num);
}



void constructionHotSpotsTree(STHist& hotSpotTree, int bRoot, int NB){
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



int computeNBi(STHist &hotSpotTree, int bRoot, int iSon){
	return 1;
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

/*2:*/	Envelope D = hotSpotTree[i].mbr;
	unsigned int F = hotSpotTree[i].F;

/*3:*/	float w =  (1/pow(s,0.5))* (D.MaxX - D.MinY);
/*4:*/	float h =  (1/pow(s,0.5))* (D.MaxY - D.MinY);

/*5, 6:*/

	for(coordenada& oi= SortedListX.begin(); oi != SortedListX.end(); oi++){
		double W = SortedListX[oi] + w;
		vector<coordenada> SortedCandiListY;

		for(coordenada& it = oi; it != SortedListX.end(); it++){
			if(it.x <= W){  // verificar no artigo se menor ou menor igual
				SortedCandiListY.push_back(it);
			}else{
				break;
			}
/*8:*/		}

		usigned int freqW = SortedCandiListY.size();
/*7:*/		if(freqW >= F/f){
/*8:*/			sort(SortedCandiListY.begin(), SortedCandiListY.end(),yCoordinateAscending);
			double H = SortedCandiListY[0].y + h;

/*9, 10:*/		for(int oj=0; oj < freqW; oj++){
				usigned int freqH = -oj + 1 + limitIndex_binarySearch(SortedCandiListY,oj,SortedCandiListY.size(), H);
				if(freqH >= F/f){
/*12*/
					hotSpot R;
					R.mbr.MinX = SortedListX[oi].x;
					R.mbr.MinY = SortedCandiListY[oj].y;
					R.mbr.MaxX = SortedListX[io].x + w;
					R.mbr.MaxY = SortedCandiListY[oj].y +h;
					R.F = freqH;
					R.O.insert(R.O.begin(),SortedCandiListY.begin()+oj, SortedCandiListY.begin()+oj+freqH);

/*13*/					hotSpotTree.push_back(R);
					hotSpotsCount++;


					auto oY = SortedListY.begin();
					for(usingnet int oC = 0; oC < SortedCandiListY.size(); oC++){
						while(oY != SortedListY.end()){
							if((*oY).x == SortedCandiListY[oC].x  && (*oY).y == SortedCandiListY[oC].y){
								oY = SortedListY.remove_after(oY);
							}else{
								oY++;
							}
						}
					}

					sort(SortedCandiListY.begin(), SortedCandiListY.end(),xCoordinateAscending);

					auto oX = SortedListX.begin();
					for(usingnet int oC = 0; oC < SortedCandiListY.size(); oC++){
						while(oX != SortedListX.end)){
							if((*oX).x == SortedCandiListY[oC].x  && (*oX).y == SortedCandiListY[oC].y){
								oX = SortedListX.remove_after(oX);
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

bool xCoordinateAscendinge(coordenada a, coordenada b){
	return a.x < b.x;
}


bool idCoordinateAscending(coordenada a, coordenada b){
	return a.id < b.id;
}

int limitIndex_binarySearch(const vector<double> &vetor; int begin; int end, double value){
	int meio = (end-begin)/2 +begin;
	while(end-begin >3){
		if(value < vetor[meio]){
			end = meio;
		}else{
			begin = meio;
		}
		meio = (end - begin)/2;
	}
	int x = begin;
	while(x < end) {
		if(vector[x+1] > value )
			return x;
	}
	return x;
}

double STHist_search(STHist *hist, Envelope query);

void STHist_print(STHist *hist);

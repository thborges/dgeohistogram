
#include <float.h>
#include "minskew.h"
#include "utils.h"

int partition(int start, int end, int n, double A[n]) {
	double x = A[end];
	int i = start-1;
	for(int j=start; j<end; j++) {
		if (A[j] <= x) {
			i++;
			double aux = A[i];
			A[i] = A[j];
			A[j] = aux;
		}
	}
	double aux = A[i+1];
	A[i+1] = x;
	A[end] = aux;
	return i+1;
}

double selectnth_rec(int start, int end, int n, double values[n], int nth) {
	if (start == end)
		return values[start];
	else {
		int q = partition(start, end, n, values);
		if (q+1 == nth)
			return values[q];
		else if (q > nth)
			return selectnth_rec(start, q-1, n, values, nth);
		else
			return selectnth_rec(q, end, n, values, nth);
	}
}

double selectnth(int n, double values[n], int nth) {
	return selectnth_rec(0, n-1, n, values, nth);
}

void print_dataset_specs(dataset_histogram *dh) {

	//variance_result r = calculate_skew_row_col(dh, 
	//	0, dh->xqtd-1, 0, dh->yqtd-1);

	// find values
	int qtd = 0;
	double total = 0.0;
	double max = 0.0;
	double min = DBL_MAX;
	double values[dh->xqtd * dh->yqtd];
	for(int x = 0; x < dh->xqtd; x++) {
		for(int y = 0; y < dh->yqtd; y++) {
			double aux = GET_HISTOGRAM_CELL(dh, x, y)->cardin;
			//if (aux > 0) {
				values[qtd] = aux;
				total += aux;
				qtd++;
			//}
			if (max < aux)
				max = aux;
			if (min > aux)
				min = aux;
		}
	}

	double mean = total / (double)qtd;
	double stdev = stdevd(values, 0, qtd);
	double median = selectnth(qtd, values, qtd/2);

	double skewness = (mean - median) / stdev;
	printf("Dataset specs (avg %f, min %f, max %f, median %f, stdev %f)\n\tLocality skewness: %f\n", 
		mean, min, max, median, stdev, skewness);
	printf("\tLocality skewness 2: %f\n", stdev/(max-min));
	printf("\tCV: %f\n", stdev/mean);
	printf("\tRSD: %f\n", pow(stdev,2)/mean);
}


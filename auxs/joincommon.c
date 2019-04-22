
#include <stdio.h>
#include <utils.h>
#include "joincommon.h"

double get_multiway_estimate_to_pnts(const void *data, const int n) {
	return ((multiway_histogram_estimate*)data)[n].to_pnts;
}

double get_multiway_estimate_io_pnts(const void *data, const int n) {
	return ((multiway_histogram_estimate*)data)[n].io_pnts;
}

void multiway_totalize_estimate(multiway_histogram_estimate *estimate, int servers,
	double *totalpnts, double *totalcomm, 
	double *mkspan, double *max_comm,
	double *stdev_mkspan, double *stdev_comm,
	double *mkspan_gap) {

	double to_pnts = 0, io_pnts = 0;
    double max_to_pnts = estimate[1].to_pnts;
    double max_io_pnts = estimate[1].io_pnts;
	for(int s = 1; s <= servers; s++) {
		to_pnts += estimate[s].to_pnts;
		io_pnts += estimate[s].io_pnts;
        if (max_to_pnts < estimate[s].to_pnts)
            max_to_pnts = estimate[s].to_pnts;
        if (max_io_pnts < estimate[s].io_pnts)
            max_io_pnts = estimate[s].io_pnts;
	}

	if (totalpnts)
		*totalpnts = to_pnts;
	if (totalcomm)
		*totalcomm = io_pnts;
	if (mkspan)
		*mkspan = max_to_pnts;
	if (max_comm)
		*max_comm = max_io_pnts;
	if (stdev_mkspan)
		*stdev_mkspan = stdevd_ex(estimate, 1, servers+1, get_multiway_estimate_to_pnts);
	if (stdev_comm)
		*stdev_comm = stdevd_ex(estimate, 1, servers+1, get_multiway_estimate_io_pnts);
	if (mkspan_gap) {
		double min_mksp = (to_pnts / servers);
		*mkspan_gap = (max_to_pnts - min_mksp) / min_mksp * 100;
	}
}
	
void multiway_print_estimate(char *name, multiway_histogram_estimate *estimate, int servers) {

	double to_pnts = 0, io_pnts = 0;
	printf("%6s     Points  IO Points\n", name);
	printf("------ ---------- ----------\n");
    double max_to_pnts = estimate[1].to_pnts;
    double max_io_pnts = estimate[1].io_pnts;
    double min_to_pnts = estimate[1].to_pnts;
    double min_io_pnts = estimate[1].io_pnts;
	for(int s = 1; s <= servers; s++) {
		printf("%6d %10.0f %10.0f\n", s, estimate[s].to_pnts, estimate[s].io_pnts);		
		to_pnts += estimate[s].to_pnts;
		io_pnts += estimate[s].io_pnts;
        if (max_to_pnts < estimate[s].to_pnts)
            max_to_pnts = estimate[s].to_pnts;
        if (max_io_pnts < estimate[s].io_pnts)
            max_io_pnts = estimate[s].io_pnts;
        if (min_to_pnts > estimate[s].to_pnts)
            min_to_pnts = estimate[s].to_pnts;
        if (min_io_pnts > estimate[s].io_pnts)
            min_io_pnts = estimate[s].io_pnts;
	}
	printf("------ ---------- ----------\n");
	printf("total  %10.0f %10.0f\n", to_pnts, io_pnts);
	printf("max    %10.0f %10.0f\n", max_to_pnts, max_io_pnts);
	printf("min    %10.0f %10.0f\n", min_to_pnts, min_io_pnts);

	double to_stdev = stdevd_ex(estimate, 1, servers+1, get_multiway_estimate_to_pnts);
	double io_stdev = stdevd_ex(estimate, 1, servers+1, get_multiway_estimate_io_pnts);
	printf("stdev  %10.1f %10.1f\n\n", to_stdev, io_stdev);
}



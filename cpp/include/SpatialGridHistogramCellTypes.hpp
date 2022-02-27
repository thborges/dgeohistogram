/*
 * Cell types for grid histograms
 *
 *  Created on: 2022-02-16
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#pragma once

struct SpatialGridHistogramCellDefault {
	size_t cardin;
	double avg_x;
	double avg_y;

	SpatialGridHistogramCellDefault() {
		cardin = 0;
		avg_x = 0;
		avg_y = 0;
	}
};

struct SpatialGridHistogramCellImproved {
	double cardin;
	double avg_x;
	double avg_y;
	size_t points;
	double objcount;
	double areasum;
	Envelope usedarea;

	SpatialGridHistogramCellImproved() {
		cardin = 0;
		avg_x = 0;
		avg_y = 0;
        points = 0;
		objcount = 0;
		areasum = 0;
	}
};

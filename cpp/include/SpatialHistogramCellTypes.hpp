/*
 * Cell types for grid histograms
 *
 *  Created on: 2022-02-16
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#pragma once

#include "Envelope.hpp"

class SpatialHistogramCellDefault {
public:
	double cardin;
	double avg_x;
	double avg_y;

	SpatialHistogramCellDefault() {
		cardin = 0;
		avg_x = 0;
		avg_y = 0;
	}

	virtual void _empty() {};
};

class SpatialHistogramCellImproved: public SpatialHistogramCellDefault {
public:
	size_t points;
	double objcount;
	double areasum;
	Envelope usedarea;
	struct {
		double x;
		double y;
		double weightsumx;
		double weightsumy;
	} centerofmass;

	SpatialHistogramCellImproved() {
        points = 0;
		objcount = 0;
		areasum = 0;
		centerofmass.x = 0;
		centerofmass.y = 0;
		centerofmass.weightsumx = 0;
		centerofmass.weightsumy = 0;
 	}
};

/*
 * Improved Grid Spatial Histogram based on the work
 * de Oliveira, T. B. (2017). Efficient Processing of Multiway Spatial Join Queries in 
 * Distributed Systems. PhD thesis, Instituto de Informática, Universidade Federal de Goiás,
 * Goiânia, GO, Brasil.
 *
 *  Created on: 2022-02-16
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#pragma once

#include <string>
#include "Dataset.hpp"
#include "Envelope.hpp"
#include "SpatialGridHistogram.hpp"
#include "SpatialHistogramCellTypes.hpp"

class SpatialGridHistogramIHWAF: public SpatialGridHistogram {
public:
	SpatialGridHistogramIHWAF(Dataset& ds);
	~SpatialGridHistogramIHWAF();

	virtual double estimateWQuery(const Envelope& wquery) override;

	double getSize() const override 
    {
        size_t bytes = (xqtd+1)*(yqtd+1)*sizeof(double);
        bytes += (xqtd*yqtd) * sizeof(SpatialHistogramCellImproved);
        return (bytes/(1000.0));
    }
	
	virtual const std::string name() override {
		return "IHWAF";
	}

protected:
	SpatialHistogramCellImproved* hcells;
	SpatialHistogramCellImproved* getHistogramCell(int x, int y) override;
	virtual void allocCells(int size) override;

private:
	void fillHistogramProportionalOverlap(Dataset& ds);
};



/*
 * Grid Spatial Histogram based on the work
 * Nikos Mamoulis and Dimitris Papadias. “Multiway Spatial Joins”. In: 
 * ACM Transactions on Database Systems 26.4 (2001), pp. 424–475
 *
 *  Created on: 2022-02-15
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#pragma once

#include <string>
#include "Dataset.hpp"
#include "Envelope.hpp"
#include "SpatialGridHistogram.hpp"
#include "SpatialHistogramCellTypes.hpp"

class SpatialGridHistogramMP: public SpatialGridHistogram {
public:
	SpatialGridHistogramMP(Dataset& ds, int xqtd, int yqtd);
	~SpatialGridHistogramMP();

	virtual double estimateWQuery(const Envelope& wquery) override;
	
	virtual const std::string name() override {
		return "MP";
	}
	
protected:
	SpatialHistogramCellDefault *hcells;
	virtual SpatialHistogramCellDefault* getHistogramCell(int x, int y) override;
	virtual void allocCells(int size) override;

private:
	void fillHistogramMbrCenter(Dataset& ds);
};



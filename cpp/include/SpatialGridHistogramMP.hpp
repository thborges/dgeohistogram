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
#include "SpatialGridHistogramCellTypes.hpp"

class SpatialGridHistogramMP: public SpatialGridHistogramTemplate<SpatialGridHistogramCellDefault> {
	public:
		SpatialGridHistogramMP(Dataset& ds, int xqtd, int yqtd);
		virtual double estimateWQuery(const Envelope& wquery);
		
		virtual const std::string name() {
			return "MP";
		}

	private:
		void fillHistogramMbrCenter(Dataset& ds);
};



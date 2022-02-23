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
#include "SpatialGridHistogramCellTypes.hpp"

class SpatialGridHistogramIHWAF: public SpatialGridHistogramTemplate<SpatialGridHistogramCellImproved> {
	public:
		SpatialGridHistogramIHWAF(Dataset& ds);
		virtual double estimateWQuery(const Envelope& wquery);
		
		virtual const std::string name() {
			return "IHWAF";
		}

	private:
		void fillHistogramProportionalOverlap(Dataset& ds);
        void getIntersectionIdxs(const DatasetMetadata& meta, 
            const Envelope& query, int *xini, int *xfim, int *yini, int *yfim);
};



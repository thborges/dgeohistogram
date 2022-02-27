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
		virtual double estimateWQuery(const Envelope& wquery) override;
		
		virtual const std::string name() override {
			return "IHWAF";
		}
		
		virtual double getCellCardin(int x, int y) override {
			return getHistogramCell(x, y)->cardin;
		}

		virtual double getCellObjCount(int x, int y) override {
			return getHistogramCell(x, y)->objcount;
		}

		virtual double getCellAvgX(int x, int y) override {
			return getHistogramCell(x, y)->avg_x;
		}

		virtual double getCellAvgY(int x, int y) override {
			return getHistogramCell(x, y)->avg_y;
		}

		virtual const Envelope getUsedArea(int x, int y) override {
			auto cell = getHistogramCell(x, y);
			return cell->usedarea;
		}

	private:
		void fillHistogramProportionalOverlap(Dataset& ds);
};



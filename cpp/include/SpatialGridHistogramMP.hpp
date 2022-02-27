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
		virtual double estimateWQuery(const Envelope& wquery) override;
		
		virtual const std::string name() override {
			return "MP";
		}

		virtual double getCellCardin(int x, int y) override {
			return getHistogramCell(x, y)->cardin;
		}

		virtual double getCellObjCount(int x, int y) override {
			return getCellCardin(x, y);
		}

		virtual double getCellAvgX(int x, int y) override {
			return getHistogramCell(x, y)->avg_x;
		}

		virtual double getCellAvgY(int x, int y) override {
			return getHistogramCell(x, y)->avg_y;
		}

		virtual const Envelope getUsedArea(int x, int y) override {
			return Envelope(xtics[x], ytics[y], xtics[x+1], ytics[y+1]);
		}
		
	private:
		void fillHistogramMbrCenter(Dataset& ds);
};



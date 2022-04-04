/*
 * AB Spatial Histogram based on the work:
 * Cheng, et al. "The Generic Annular Bucket Histogram for Estimating
 * the Seletivity of Spatial Selection and Spatial Join"
 * Geo-spatial Information Science 14.4 (2011): 262-273.
 *
 *  Created on: 2022-03-01
 *      Author: Gabriel Portela Macedo Souza <gportelams@gmail.com>
 */

#pragma once

#include <string>
#include <vector>
#include "Dataset.hpp"
#include "Envelope.hpp"
#include "SpatialGridHistogram.hpp"

struct ABBucket {
    int id;
    double cardin;

    Envelope outer;
    Envelope inner;

    bool contains(const Envelope& mbr) const
    {
        return mbr.MinX >= outer.MinX && mbr.MinX <= inner.MinX &&
               mbr.MaxX >= inner.MaxX && mbr.MaxX <= outer.MaxX &&
               mbr.MinY >= outer.MinY && mbr.MinY <= inner.MinY &&
               mbr.MaxY >= inner.MaxY && mbr.MaxY <= outer.MaxY;
    }
};

class SpatialHistogramAB: public SpatialHistogram {
public:
	SpatialHistogramAB(Dataset& ds, int columns, int rows);
	~SpatialHistogramAB();

	virtual double estimateWQuery(const Envelope& wquery) override;

	virtual const std::string name() override {
		return "AB";
	}

	virtual void printGeoJson(const std::string& filename) override;

    int getcolumns() const { return columns;}
    int getrows() const { return rows;}

private:
    std::vector<ABBucket> buckets;
    int columns;
    int rows;
    double cellWidth;
    double cellHeight;
    double bottomLeftX;
    double bottomLeftY;

    ABBucket getMBRBucket(const Envelope& mbr);
    double intersectionEstimation(const Envelope& wquery);
    double withinEstimation(const Envelope& wquery);
};



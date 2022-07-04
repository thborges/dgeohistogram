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

struct Cell {
    int X;
    int Y;
};

struct Pos {
    double X;
    double Y;
};

struct ABBucket {
    int ID;
    double Cardinality;

    Cell MinCell;
    Cell MaxCell;

    ABBucket() = default;

    bool contains(Cell minCell, Cell maxCell) const
    {
        return MinCell.X == minCell.X &&
               MinCell.Y == minCell.Y &&
               MaxCell.X == maxCell.X &&
               MaxCell.Y == maxCell.Y;
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
    int getNumBuckets() const { return buckets.size(); }
    double getSize() const override {
        size_t bytes = buckets.size() * sizeof(ABBucket);
        return (bytes/(1000.0));
    }

private:
    std::vector<ABBucket> buckets;
    int columns;
    int rows;
    double cellWidth;
    double cellHeight;
    Pos bottomLeft;

    Cell getCell(double x, double y);
    Pos getCellBottomLeft(Cell cell);
    Pos getCellTopRight(Cell cell);
    Envelope getOuterRect(ABBucket);
    Envelope getInnerRect(ABBucket);
    double intersectionEstimation(const Envelope& wquery);
    double withinEstimation(const Envelope& wquery);
};



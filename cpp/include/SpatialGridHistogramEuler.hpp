/*
 * Euler Spatial Histogram based on the work:
 * Sun, Chengyu, et al. "Exploring spatial datasets with histograms." 
 * Distributed and Parallel Databases 20.1 (2006): 57-88.
 *
 *  Created on: 2022-02-27
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#pragma once

#include <string>
#include "Dataset.hpp"
#include "Envelope.hpp"
#include "SpatialGridHistogram.hpp"
#include "SpatialHistogramCellTypes.hpp"

class EulerHistogramEdge {
public:
    double cardin;
    double avg_length;
    Envelope mbr;

    EulerHistogramEdge() {
        cardin = avg_length = 0.0;
    }

    bool isVertical() {
        return mbr.width() < mbr.length();
    }
};

class EulerHistogramVertex {
public:
    double x;
    double y;
    double cardin;
    EulerHistogramVertex() {
        x = y = cardin = 0.0;
    }
};

class SpatialGridHistogramEuler: public SpatialGridHistogram {
public:
	SpatialGridHistogramEuler(Dataset& ds, int xqtd, int yqtd);
	~SpatialGridHistogramEuler();

	virtual double estimateWQuery(const Envelope& wquery) override;

    double getSize() const override 
    {
        size_t bytes = (xqtd+1)*(yqtd+1)*sizeof(double);
        bytes += (xqtd*yqtd) * sizeof(SpatialHistogramCellImproved);
        bytes += ((xqtd+1) * yqtd) + ((yqtd+1) * xqtd) * sizeof(EulerHistogramEdge);
        bytes += (xqtd+1) * (yqtd+1) * sizeof(EulerHistogramVertex);
        return (bytes/(1000.0));
    }
	
	virtual const std::string name() override {
		return "Euler";
	}

	virtual void printGeoJson(const std::string& filename) override;

protected:
	SpatialHistogramCellImproved *faces;
    EulerHistogramEdge *edges;
    EulerHistogramVertex *vertexes;
	virtual SpatialHistogramCellImproved* getHistogramCell(int x, int y) override;
	virtual void allocCells(int size) override;

private:
	void fillHistogram(Dataset& ds);
};



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
#include <fstream>
#include <iomanip>
#include "Dataset.hpp"
#include "Envelope.hpp"
#include "SpatialHistogram.hpp"
#include "SpatialHistogramCellTypes.hpp"

class SpatialGridHistogram: public SpatialHistogram {
public:
        int columns() { return xqtd; };
        int rows() { return yqtd; };
        virtual double getColumnX(int i) const { return xtics[i]; };
        virtual double getRowY(int j) const { return ytics[j]; };
        virtual double getXSize() const { return xsize; };
        virtual double getYSize() const { return ysize; };
        virtual const std::string name() = 0;
        virtual void getIntersectionIdxs(const Envelope& query, int *xini, int *xfim, int *yini, int *yfim);
        virtual SpatialHistogramCellDefault* getHistogramCell(int x, int y) = 0;
        virtual void printGeoJson(const std::string& filename);
        virtual ~SpatialGridHistogram();

protected:
        int xqtd;
        int yqtd;
        double xsize;
        double ysize;
        double *xtics;
        double *ytics;
        
        virtual void histogramAlloc(Dataset& ds, int xqtd, int yqtd);
        virtual void allocCells(int size) = 0;
};

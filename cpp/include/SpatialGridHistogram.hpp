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

class SpatialGridHistogram: public SpatialHistogram {
	public:
        int columns() { return xqtd; };
        int rows() { return yqtd; };

	protected:
        int xqtd;
        int yqtd;
        double xsize;
        double ysize;
        double *xtics;
        double *ytics;

};

template<typename celltype>
class SpatialGridHistogramTemplate: public SpatialGridHistogram {
protected:
    celltype *hcells;
    virtual void histogramAlloc(Dataset& ds, int xqtd, int yqtd) {
        this->xqtd = xqtd;
        this->yqtd = yqtd;
        this->xtics = new double[xqtd+1];
        this->ytics = new double[yqtd+1];
        this->hcells = new celltype[xqtd*yqtd];

        const DatasetMetadata& meta = ds.metadata();
        double rangex = meta.mbr.MaxX - meta.mbr.MinX;
        double rangey = meta.mbr.MaxY - meta.mbr.MinY;
        this->xsize = rangex / xqtd;
        this->ysize = rangey / yqtd;
        
        // X tics
        double xini = meta.mbr.MinX;
        for(int i = 0; i < xqtd; i++)
            xtics[i] = xini + (xsize * i);
        xtics[xqtd] = meta.mbr.MaxX;

        // Y tics
        double yini = meta.mbr.MinY;
        for(int i = 0; i < yqtd; i++)
            ytics[i] = yini + (ysize * i);
        ytics[yqtd] = meta.mbr.MaxY;

        mbr = meta.mbr;
    };
    
    virtual celltype *getHistogramCell(int x, int y) {
        return &hcells[x * yqtd + y];
    }

    virtual void printGeoJson(const std::string& filename) override {
        std::ofstream output;
        output.open(filename);
        output << std::fixed;
        output.precision(15);
        output << "{\"type\": \"FeatureCollection\", \"features\": [\n";

        Envelope e;
        for(int x = 0; x < xqtd; x++) {
            e.MinX = xtics[x];
            e.MaxX = xtics[x+1];
            for(int y = 0; y < yqtd; y++) {
                e.MinY = ytics[y];
                e.MaxY = ytics[y+1];

                output << "{\"type\": \"Feature\", \"geometry\": {\"type\": \"Polygon\", \"coordinates\": [[";
                output << "[" << e.MinX << "," << e.MinY << "],";
                output << "[" << e.MaxX << "," << e.MinY << "],";
                output << "[" << e.MaxX << "," << e.MaxY << "],";
                output << "[" << e.MinX << "," << e.MaxY << "],";
                output << "[" << e.MinX << "," << e.MinY << "]";
                output << "]]}, \"properties\": {";
                output << "\"name\": \"" << x << "." << y << "\",";
                output << "\"card\": " << hcells[x*yqtd + y].cardin << ",";
                output << "\"avg_x\": " << hcells[x*yqtd + y].avg_x << ",";
                output << "\"avg_y\": " << hcells[x*yqtd + y].avg_y;
                
                if (x == xqtd-1 && y == yqtd-1)
                    output << "}}\n";
                else
                    output << "}},\n";
            }
        }

        output << "]}\n";
        output.close();
    }
};

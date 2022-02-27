/*
 * Base Grid Histogram
 *
 *  Created on: 2022-02-27
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#include "../include/SpatialGridHistogram.hpp"

SpatialGridHistogram::~SpatialGridHistogram() {
    delete[] xtics;
    delete[] ytics;
}

void SpatialGridHistogram::histogramAlloc(Dataset& ds, int xqtd, int yqtd) {
    this->xqtd = xqtd;
    this->yqtd = yqtd;
    this->xtics = new double[xqtd+1];
    this->ytics = new double[yqtd+1];
    allocCells(xqtd*yqtd);

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

    hmbr = meta.mbr;
};

void SpatialGridHistogram::getIntersectionIdxs(const Envelope& query, 
    int *xini, int *xfim, int *yini, int *yfim) {

    // prevent values of x and y out of histogram bounds
    if (!query.intersects(mbr())) {
        *xini = *yini = 0;
        *xfim = *yfim = -1;
    } else {
        // prevent values of x and y out of histogram bounds
        Envelope nquery = query.intersection(mbr());

        *xini = (nquery.MinX - hmbr.MinX) / xsize;
        *xfim = (nquery.MaxX - hmbr.MinX) / xsize;
        *yini = (nquery.MinY - hmbr.MinY) / ysize;
        *yfim = (nquery.MaxY - hmbr.MinY) / ysize;

        // turn back one cell, due to rouding errors
        if (*xini > 0) (*xini)--;
        if (*xfim > 0) (*xfim)--;
        if (*yini > 0) (*yini)--;
        if (*yfim > 0) (*yfim)--;

        while (xtics[*xini+1]   < query.MinX && *xini < xqtd) 
            (*xini)++;
        while (xtics[(*xfim)+1] < query.MaxX && *xfim < xqtd-1)
            (*xfim)++;
        while (ytics[*yini+1]   < query.MinY && *yini < yqtd)
            (*yini)++;
        while (ytics[(*yfim)+1] < query.MaxY && *yfim < yqtd-1)
            (*yfim)++;

        if (*xini < 0 || *xfim >= xqtd)
            throw std::runtime_error("x is out of histogram bounds.");

        if (*yini < 0 || *yfim >= yqtd)
            throw std::runtime_error("y is out of histogram bounds.");
    }
}

void SpatialGridHistogram::printGeoJson(const std::string& filename) {
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
            auto *cell = getHistogramCell(x, y);
            output << "\"card\": " << cell->cardin << ",";
            output << "\"avg_x\": " << cell->avg_x << ",";
            output << "\"avg_y\": " << cell->avg_y;
            
            if (x == xqtd-1 && y == yqtd-1)
                output << "}}\n";
            else
                output << "}},\n";
        }
    }

    output << "]}\n";
    output.close();
}

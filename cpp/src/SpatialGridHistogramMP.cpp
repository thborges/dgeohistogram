/*
 * Grid Spatial Histogram based on the work
 * Nikos Mamoulis and Dimitris Papadias. “Multiway Spatial Joins”. In: 
 * ACM Transactions on Database Systems 26.4 (2001), pp. 424–475
 *
 *  Created on: 2022-02-15
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#include "../include/SpatialGridHistogramMP.hpp"

void SpatialGridHistogramMP::allocCells(int size) {
    hcells = new SpatialHistogramCellDefault[size];
}

SpatialGridHistogramMP::~SpatialGridHistogramMP() {
    delete[] hcells;
}

SpatialGridHistogramMP::SpatialGridHistogramMP(Dataset& ds, int xqtd, int yqtd) {
    histogramAlloc(ds, xqtd, yqtd);
    fillHistogramMbrCenter(ds);
}

SpatialHistogramCellDefault* SpatialGridHistogramMP::getHistogramCell(int x, int y) {
    assert(x < xqtd);
    assert(y < yqtd);
    return &hcells[x * yqtd + y];
}

void SpatialGridHistogramMP::fillHistogramMbrCenter(Dataset& ds) {
    // hash dataset objects using mbr center
    const DatasetMetadata& meta = ds.metadata();
    for(const DatasetEntry& de : ds.geoms()) {
        double x = de.mbr.MinX + (de.mbr.MaxX - de.mbr.MinX) / 2.0;
        double y = de.mbr.MinY + (de.mbr.MaxY - de.mbr.MinY) / 2.0;

        int xp = (x - meta.mbr.MinX) / this->xsize;
        int yp = (y - meta.mbr.MinY) / this->ysize;
        SpatialHistogramCellDefault *cell = getHistogramCell(xp, yp);
        cell->cardin++;

        Envelope rs;
        rs.MinX = xtics[xp];
        rs.MaxX = xtics[xp+1];
        rs.MinY = ytics[yp];
        rs.MaxY = ytics[yp+1];
        Envelope inters = de.mbr.intersection(rs);
        cell->avg_x += (inters.width() - cell->avg_x) / cell->cardin;
        cell->avg_y += (inters.length() - cell->avg_y) / cell->cardin;
        
    }
}

double SpatialGridHistogramMP::estimateWQuery(const Envelope& query) {
    double result = 0.0;
    if (!query.intersects(mbr()))
        return result;
    Envelope nquery = query.intersection(mbr());
    int xini, xfim, yini, yfim;
    getIntersectionIdxs(nquery, &xini, &xfim, &yini, &yfim);

    assert(getColumnX(xini) <= nquery.MinX);
    assert(getColumnX(xfim+1) >= nquery.MaxX);
    assert(getRowY(yini) <= nquery.MinY);
    assert(getRowY(yfim+1) >= nquery.MaxY);

    for(int x = xini; x <= xfim; x++) {
        Envelope rs;
        rs.MinX = xtics[x];
        rs.MaxX = xtics[x+1];

        for(int y = yini; y <= yfim; y++) {
            rs.MinY = ytics[y];
            rs.MaxY = ytics[y+1];

            auto *cell = getHistogramCell(x, y);
            if (query.intersects(rs)) {
                /* this formula came from the following paper of same authors:
                   Nikos Mamoulis and Dimitris Papadias. “Advances in Spatial and Temporal Databases”. 
                   In: ed. by Christian S. Jensen et al. Vol. 2121. Lecture Notes in Computer Science. 
                   Springer, 2001. Chap. Selectivity Estimation of Complex Spatial Queries, pp. 155–174.
                   See also Equation 2.2 in de Oliveira, T.B. thesis */
                
                Envelope inters = query.intersection(rs);

                // uniformity assumption! objects aren't uniform located at the cell
                double xprob = std::min(1.0, (cell->avg_x + inters.width()) / rs.width());
                double yprob = std::min(1.0, (cell->avg_y + inters.length()) / rs.length());
                result += cell->cardin * xprob * yprob;
            }
        }
    }
    return result;
}


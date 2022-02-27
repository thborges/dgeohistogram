/*
 * Improved Grid Spatial Histogram based on the work
 * de Oliveira, T. B. (2017). Efficient Processing of Multiway Spatial Join Queries in 
 * Distributed Systems. PhD thesis, Instituto de Informática, Universidade Federal de Goiás,
 * Goiânia, GO, Brasil.
 *
 *  Created on: 2022-02-16
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#include "../include/SpatialGridHistogramIHWAF.hpp"

void SpatialGridHistogramIHWAF::allocCells(int size) {
    hcells = new SpatialHistogramCellImproved[size];
}

SpatialGridHistogramIHWAF::~SpatialGridHistogramIHWAF() {
    delete[] hcells;
}

SpatialHistogramCellImproved* SpatialGridHistogramIHWAF::getHistogramCell(int x, int y) {
    assert(x < xqtd);
    assert(y < yqtd);
    return &hcells[x * yqtd + y];
}

SpatialGridHistogramIHWAF::SpatialGridHistogramIHWAF(Dataset& ds) {

    // IHWAF uses average length of objets to decide the
    // histogram granularity (quantity of cells or split method)
    const DatasetMetadata& metadata = ds.metadata();
    double psizex = metadata.x_average + 6*metadata.x_stdev();
    double psizey = metadata.y_average + 6*metadata.y_stdev();

    double rangex = metadata.mbr.MaxX - metadata.mbr.MinX;
    double rangey = metadata.mbr.MaxY - metadata.mbr.MinY;

    int qtd_x = ceil(rangex / psizex);
    int qtd_y = ceil(rangey / psizey);

	// at least 5x5
	qtd_x = std::max(5, qtd_x);
	qtd_y = std::max(5, qtd_y);

	// at most 2MB per histogram
	const double upper_bound_cell_number = 2*1024*1024 / sizeof(SpatialHistogramCellImproved);
	double adjust = std::sqrt((qtd_x * qtd_y) / upper_bound_cell_number);
	if (adjust > 1.0) {
		qtd_x = ceil(rangex / (psizex * adjust));
		qtd_y = ceil(rangey / (psizey * adjust));
	}

    histogramAlloc(ds, qtd_x, qtd_y);
    fillHistogramProportionalOverlap(ds);
}

void SpatialGridHistogramIHWAF::fillHistogramProportionalOverlap(Dataset& ds) {
    
    for(DatasetEntry de : ds.geoms()) {
        double objarea = de.mbr.area();
        
        int xini, xfim, yini, yfim;
        getIntersectionIdxs(de.mbr, &xini, &xfim, &yini, &yfim);	

        for(int x = xini; x <= xfim; x++) {
            Envelope rs;
            rs.MinX = xtics[x];
            rs.MaxX = xtics[x+1];

            for(int y = yini; y <= yfim; y++) {
                rs.MinY = ytics[y];
                rs.MaxY = ytics[y+1];

                Envelope inters = de.mbr.intersection(rs);
                double intarea = inters.area();
                double fraction;
                if (intarea <= 0.0) { // parallel to one axis
                    bool parallel_y = (de.mbr.MaxX - de.mbr.MinX < 1e-30);
                    bool parallel_x = (de.mbr.MaxY - de.mbr.MinY < 1e-30);
                    if (parallel_x && parallel_y) // point obj
                        fraction = 1;
                    else if (parallel_x) {
                        // the part of de.mbr inside rs / the length of rs in X
                        double length_inside_cell = std::min(rs.MaxX, de.mbr.MaxX) - std::max(rs.MinX, de.mbr.MinX);
                        fraction = length_inside_cell / (de.mbr.MaxX - de.mbr.MinX);
                    } else {
                        double length_inside_cell = std::min(rs.MaxY, de.mbr.MaxY) - std::max(rs.MinY, de.mbr.MinY);
                        fraction = length_inside_cell / (de.mbr.MaxY - de.mbr.MinY);
                    }
                }
                else {
                    fraction = intarea / objarea;
                }

                SpatialHistogramCellImproved *cell = getHistogramCell(x, y);
                cell->cardin += fraction;
                
                // used area in the cell
                cell->usedarea.merge(inters);

                cell->objcount += 1.0;

                // average length, online average calculation
                cell->avg_x += (inters.width() - cell->avg_x) / cell->objcount;
                cell->avg_y += (inters.length() - cell->avg_y) / cell->objcount;
            }
        }
    }
}

double SpatialGridHistogramIHWAF::estimateWQuery(const Envelope& query) {
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
            if (nquery.intersects(cell->usedarea)) {
                /* See de OLIVEIRA, T. B. Efficient Processing of Multiway Spatial Join Queries 
                 * in Distributed Systems. 152 p. Tese (Doutorado) — Instituto de Informática, 
                 * Universidade Federal de Goiás, Goiânia, GO, Brasil, 2017.
                 */
                
                double avg_x = cell->avg_x;
                double avg_y = cell->avg_y;

                // observing that objects generally doesn't overlap in both axis,
                // (e.g. political land limits)
	            // fix the probability of intersection in one of them
	            double ux = cell->usedarea.MaxX - cell->usedarea.MinX;
	            double uy = cell->usedarea.MaxY - cell->usedarea.MinY;
	            if (avg_x > avg_y)
		            avg_x = std::min(avg_x, avg_x/(avg_y * cell->objcount / uy));
	            else
		            avg_y = std::min(avg_y, avg_y/(avg_x * cell->objcount / ux));

                Envelope wqueryadj = nquery.intersection(cell->usedarea);

                // uniformity assumption! objects aren't uniform located at the cell
                double xprob = std::min(1.0, (avg_x + wqueryadj.width()) / cell->usedarea.width());
                double yprob = std::min(1.0, (avg_y + wqueryadj.length()) / cell->usedarea.length());

                result += cell->cardin * xprob * yprob;
            }
        }
    }
    return result;
}


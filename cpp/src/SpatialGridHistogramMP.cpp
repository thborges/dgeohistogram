
#include "SpatialGridHistogramMP.hpp"

SpatialGridHistogramMP::SpatialGridHistogramMP(Dataset& ds, int xqtd, int yqtd) {
    histogramAlloc(ds, xqtd, yqtd);
    fillHistogramMbrCenter(ds);
}

void SpatialGridHistogramMP::fillHistogramMbrCenter(Dataset& ds) {
    // hash dataset objects using mbr center
    const DatasetMetadata& meta = ds.metadata();
    for(const DatasetEntry& de : ds.geoms()) {
        double x = de.mbr.MinX + (de.mbr.MaxX - de.mbr.MinX) / 2.0;
        double y = de.mbr.MinY + (de.mbr.MaxY - de.mbr.MinY) / 2.0;

        int xp = (x - meta.mbr.MinX) / this->xsize;
        int yp = (y - meta.mbr.MinY) / this->ysize;
        SpatialGridHistogramCellDefault *cell = &hcells[xp * yqtd + yp];
        cell->cardin++;

        Envelope rs;
        rs.MinX = xtics[xp];
        rs.MaxX = xtics[xp+1];
        rs.MinY = ytics[yp];
        rs.MaxY = ytics[yp+1];
        Envelope inters = de.mbr.intersection(rs);
        double delta_x = abs(inters.MaxX - inters.MinX);
        double delta_y = abs(inters.MaxY - inters.MinY);
        cell->avg_x += delta_x;
        cell->avg_y += delta_y;
    }

    for(int x = 0; x < xqtd; x++) {
        for(int y = 0; y < yqtd; y++) {
            auto *cell = getHistogramCell(x, y);
            if (cell->cardin > 0.0) {
                cell->avg_x /= cell->cardin;
                cell->avg_y /= cell->cardin;
            }
        }
    }
}

double SpatialGridHistogramMP::estimateWQuery(const Envelope& query) {
    double result = 0.0;

    int xini = (query.MinX - mbr->MinX) / xsize;
    int xfim = (query.MaxX - mbr->MinX) / xsize;
    int yini = (query.MinY - mbr->MinY) / ysize;
    int yfim = (query.MaxY - mbr->MinY) / ysize;

    xfim = std::min(xfim, xqtd-1);
    yfim = std::min(yfim, yqtd-1);
    xini = std::max(xini, 0);
    yini = std::max(yini, 0);

    const double epsilon = 1e-100;
    if (query.MaxX - xtics[xfim] < epsilon && xfim > 0) {
        xfim--;
    }
    if (query.MaxY - ytics[yfim] < epsilon && yfim > 0) {
        yfim--;
    }
    if (xfim < xini)
        xini = xfim;
    if (yfim < yini)
        yini = yfim;

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


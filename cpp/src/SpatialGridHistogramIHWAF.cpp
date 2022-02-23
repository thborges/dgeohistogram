
#include "SpatialGridHistogramIHWAF.hpp"

SpatialGridHistogramIHWAF::SpatialGridHistogramIHWAF(Dataset& ds) {

    // IHWAF uses average length of objets to decide the
    // histogram granularity (quantity of cells or split method)
    const DatasetMetadata& metadata = ds.metadata();
    double psizex = metadata.x_average + 4*metadata.x_stdev();
    double psizey = metadata.y_average + 4*metadata.y_stdev();

    double rangex = metadata.mbr.MaxX - metadata.mbr.MinX;
    double rangey = metadata.mbr.MaxY - metadata.mbr.MinY;

    int qtd_x = ceil(rangex / psizex);
    int qtd_y = ceil(rangey / psizey);

	// at least 5x5
	qtd_x = std::max(5, qtd_x);
	qtd_y = std::max(5, qtd_y);

	// at most 2MB per histogram
	const double upper_bound_cell_number = 2*1024*1024 / sizeof(SpatialGridHistogramCellImproved);
	double adjust = std::sqrt((qtd_x * qtd_y) / upper_bound_cell_number);
	if (adjust > 1.0) {
		qtd_x = ceil(rangex / (psizex * adjust));
		qtd_y = ceil(rangey / (psizey * adjust));
	}

    histogramAlloc(ds, qtd_x, qtd_y);
    fillHistogramProportionalOverlap(ds);
}

void SpatialGridHistogramIHWAF::getIntersectionIdxs(const DatasetMetadata& meta, const Envelope& query, 
	int *xini, int *xfim, int *yini, int *yfim) {

	// prevent values of x and y out of histogram bounds
	if (!query.intersects(meta.mbr)) {
		*xini = *yini = 0;
		*xfim = *yfim = -1;
	} else {
		// prevent values of x and y out of histogram bounds
		Envelope nquery = query.intersection(meta.mbr);

		*xini = (nquery.MinX - meta.mbr.MinX) / xsize;
		*xfim = (nquery.MaxX - meta.mbr.MinX) / xsize;
		*yini = (nquery.MinY - meta.mbr.MinY) / ysize;
		*yfim = (nquery.MaxY - meta.mbr.MinY) / ysize;

		if (*xfim == xqtd) (*xfim)--;
		if (*yfim == yqtd) (*yfim)--;
   
		while (xtics[*xini]     > query.MinX) (*xini)--;
		while (xtics[(*xfim)+1] < query.MaxX) (*xfim)++;
		while (ytics[*yini]     > query.MinY) (*yini)--;
		while (ytics[(*yfim)+1] < query.MaxY) (*yfim)++;

		if (*xini < 0 && *xfim >= xqtd)
            throw std::runtime_error("x is out of histogram bounds.");

		if (*yini < 0 && *yfim >= yqtd)
            throw std::runtime_error("y is out of histogram bounds.");
	}
}

void SpatialGridHistogramIHWAF::fillHistogramProportionalOverlap(Dataset& ds) {
    
    const DatasetMetadata& meta = ds.metadata();
    
    for(DatasetEntry de : ds.geoms()) {
        double objarea = de.mbr.area();
        
        int xini, xfim, yini, yfim;
        getIntersectionIdxs(meta, de.mbr, &xini, &xfim, &yini, &yfim);	

        double sum_fraction = 0;
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
                sum_fraction += fraction;

                auto *cell = getHistogramCell(x, y);
                cell->cardin += fraction;
                cell->objcount += 1.0;

                // average length, online average calculation
                double delta_x = (inters.MaxX - inters.MinX);
                double delta_y = (inters.MaxY - inters.MinY);
                cell->avg_x += (delta_x - cell->avg_x) / cell->objcount;
                cell->avg_y += (delta_y - cell->avg_y) / cell->objcount;
            }
        }
    }
}

double SpatialGridHistogramIHWAF::estimateWQuery(const Envelope& query) {
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

            SpatialGridHistogramCellImproved *cell = &hcells[x * yqtd + y];
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


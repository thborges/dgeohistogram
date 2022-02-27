/*
 * Euler Spatial Histogram based on the work:
 * Sun, Chengyu, et al. "Exploring spatial datasets with histograms." 
 * Distributed and Parallel Databases 20.1 (2006): 57-88.
 *
 *  Created on: 2022-02-27
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#define USE_UNSTABLE_GEOS_CPP_API
#include <geos/geom/Geometry.h>
#include <geos/geom/Envelope.h>
#include "../include/SpatialGridHistogramEuler.hpp"

#define GET_VERT_EDGE(x, y) ((x == xqtd) ? (x * (2*yqtd+1) + y) : (x * (2*yqtd+1) + 2*y + 1))
#define GET_HORZ_EDGE(x, y) (x * (2*yqtd+1) + 2*y)
const double eulerh_epsylon = 1e-10;

void SpatialGridHistogramEuler::allocCells(int size) {
    faces = new SpatialHistogramCellImproved[size];
    edges = new EulerHistogramEdge[((xqtd+1) * yqtd) + ((yqtd+1) * xqtd)];
    vertexes = new EulerHistogramVertex[(xqtd+1) * (yqtd+1)];
}

SpatialGridHistogramEuler::~SpatialGridHistogramEuler() {
    delete[] faces;
    delete[] edges;
    delete[] vertexes;
}

SpatialGridHistogramEuler::SpatialGridHistogramEuler(Dataset& ds, int xqtd, int yqtd) {
    histogramAlloc(ds, xqtd, yqtd);

    // set edges and vertexes positions
    int v = 0;
    for(int i = 0; i <= xqtd; i++) { // x
        for(int j = 0; j <= yqtd; j++) { // y

            // set vetex at (i,j)
            vertexes[v].x = xtics[i];
            vertexes[v].y = ytics[j];
            v++;

            // horizontal edge at vertex v
            if (i < xqtd) {
                int e = GET_HORZ_EDGE(i, j);
                edges[e].mbr.MinX = xtics[i];
                edges[e].mbr.MaxX = xtics[i+1];
                edges[e].mbr.MinY = ytics[j];
                edges[e].mbr.MaxY = ytics[j] + eulerh_epsylon;
            }

            // vertical edge at vertex v
            if (j < yqtd) {
                int e = GET_VERT_EDGE(i, j);
                edges[e].mbr.MinY = ytics[j];
                edges[e].mbr.MaxY = ytics[j+1];
                edges[e].mbr.MinX = xtics[i];
                edges[e].mbr.MaxX = xtics[i] + eulerh_epsylon;
            }
        }
    }

    fillHistogram(ds);
}

SpatialHistogramCellImproved* SpatialGridHistogramEuler::getHistogramCell(int x, int y) {
    assert(x < xqtd);
    assert(y < yqtd);
    return &faces[x * yqtd + y];
}

void SpatialGridHistogramEuler::fillHistogram(Dataset& ds) {

    const DatasetMetadata& meta = ds.metadata();
    for(const DatasetEntry& de : ds.geoms()) {
        int xini, xfim, yini, yfim;
        getIntersectionIdxs(de.mbr, &xini, &xfim, &yini, &yfim);

        // we use +1 on xfim and yfim to iterate over the edges 
        // and vertexes in the histogram top and right borders.
        for(int x = xini; x <= xfim+1; x++) {
            Envelope facembr;
            facembr.MinX = xtics[x];
            if (x < xqtd)
                facembr.MaxX = xtics[x+1];
            else
                facembr.MaxX = facembr.MinX + eulerh_epsylon; // sum a litle fraction to prevent clip error due to empty mbr

            for(int y = yini; y <= yfim+1; y++) {
                facembr.MinY = ytics[y];
                if (y < yqtd)
                    facembr.MaxY = ytics[y+1];
                else
                    facembr.MaxY = facembr.MinY + eulerh_epsylon; // sum a litle fraction to prevent clip error due to empty mbr

                GEOSGeometry *clipped_geo = GEOSClipByRect(de.geo, facembr.MinX, facembr.MinY, 
                    facembr.MaxX, facembr.MaxY);

                if (clipped_geo != NULL){
					Envelope clipped_mbr;
                    const geos::geom::Envelope *ev = ((const geos::geom::Geometry*)clipped_geo)->getEnvelopeInternal();
                    clipped_mbr.MinX = ev->getMinX();
                    clipped_mbr.MinY = ev->getMinY();
                    clipped_mbr.MaxX = ev->getMaxX();
                    clipped_mbr.MaxY = ev->getMaxY();

                    // face
                    if (x < xqtd && y < yqtd && facembr.intersects(clipped_mbr)) {
                        auto *face = getHistogramCell(x, y);
                        face->cardin += 1.0;
                        face->avg_x += (clipped_mbr.width() - face->avg_x) / face->cardin;
                        face->avg_y += (clipped_mbr.length() - face->avg_y) / face->cardin;
                        face->usedarea.merge(clipped_mbr);
                    }

                    // vertex
                    int v = x * (yqtd+1) + y;
                    if (clipped_mbr.contains(vertexes[v].x, vertexes[v].y))
                        vertexes[v].cardin += 1.0;
                    
                    auto compute_edge = [&](auto& edge) {
						if (edge.mbr.intersects(clipped_mbr)) {
							edge.cardin += 1.0;
                            const Envelope& inters = edge.mbr.intersection(clipped_mbr);
                            if (edge.isVertical())
                                edge.avg_length += (inters.length() - edge.avg_length) / edge.cardin;
                            else
                                edge.avg_length += (inters.width() - edge.avg_length) / edge.cardin;
						}
                    };

                    // horizontal edge
                    if (x < xqtd) {
                        int e = GET_HORZ_EDGE(x, y);
                        compute_edge(edges[e]);
                    }

                    // vertical edge
                    if (y < yqtd) {
                        int e = GET_VERT_EDGE(x, y);
                        compute_edge(edges[e]);
                    }
                }
                GEOSGeom_destroy(clipped_geo);
            }
        }
    }
}

double SpatialGridHistogramEuler::estimateWQuery(const Envelope& query) {
    double result = 0.0;
    if (!query.intersects(hmbr))
        return result;
    Envelope nquery = query.intersection(hmbr);
    int xini, xfim, yini, yfim;
    getIntersectionIdxs(nquery, &xini, &xfim, &yini, &yfim);

    assert(getColumnX(xini) <= nquery.MinX);
    assert(getColumnX(xfim+1) >= nquery.MaxX);
    assert(getRowY(yini) <= nquery.MinY);
    assert(getRowY(yfim+1) >= nquery.MaxY);

    Envelope facembr;
    for(int x = xini; x <= xfim; x++) {
        facembr.MinX = xtics[x];
        if (x < xqtd)
            facembr.MaxX = xtics[x+1];

        for(int y = yini; y <= yfim; y++) {
            facembr.MinY = ytics[y];
            if (y < yqtd)
                facembr.MaxY = ytics[y+1];

            // face
            if (x < xqtd && y < yqtd) {
                auto *cell = getHistogramCell(x, y);
                if (query.intersects(cell->usedarea)) {
                    double avg_x = cell->avg_x;
			        double avg_y = cell->avg_y;
                    double ux = cell->usedarea.width();
                    double uy = cell->usedarea.length();
                    if (avg_x > avg_y)
                        avg_x = std::min(avg_x, avg_x/(avg_y * cell->cardin / uy));
                    else
                        avg_y = std::min(avg_y, avg_y/(avg_x * cell->cardin / ux));

                    Envelope inters = query.intersection(cell->usedarea);
                    double xprob = std::min(1.0, (avg_x + inters.width()) / cell->usedarea.width());
                    double yprob = std::min(1.0, (avg_y + inters.length()) / cell->usedarea.length());
                    result += cell->cardin * xprob * yprob;
                }
            }

            // vertex
            int v = x * (yqtd+1) + y;
            if (query.contains(vertexes[v].x, vertexes[v].y)) {
                result += vertexes[v].cardin;
            }

            auto subtract_edge = [&](auto& edge) {
                if (edge.mbr.intersects(query)) {
                    const Envelope& inters = edge.mbr.intersection(query);
                    double prob;
                    if (edge.isVertical())
                        prob = std::min(1.0, (edge.avg_length + inters.length()) / edge.mbr.length());
                    else
                        prob = std::min(1.0, (edge.avg_length + inters.width()) / edge.mbr.width());
                    // worst!
                    // result -= edge.cardin * prob;
                    result -= edge.cardin;
                }
            };

            // horizontal edge
            if (x < xqtd) {
                int e = GET_HORZ_EDGE(x, y);
                if (edges[e].mbr.MinY != query.MinY && edges[e].mbr.MinY != query.MaxY)
                    subtract_edge(edges[e]);
            }

            // horizontal edge
            if (y < yqtd) {
                int e = GET_VERT_EDGE(x, y);
                if (edges[e].mbr.MinX != query.MinX && edges[e].mbr.MinX != query.MaxX)
                    subtract_edge(edges[e]);
            }
        }
    }
    return result;
}

void SpatialGridHistogramEuler::printGeoJson(const std::string& filename) {
    std::ofstream output;
    output.open(filename);
    output << std::fixed;
    output.precision(15);
    output << "{\"type\": \"FeatureCollection\", \"features\": [\n";

    int v = 0;
    Envelope e;
    for(int x = 0; x <= xqtd; x++) {
        e.MinX = xtics[x];
        if (x < xqtd)
            e.MaxX = xtics[x+1];
        
        for(int y = 0; y <= yqtd; y++) {
            e.MinY = ytics[y];
            if (y < yqtd)
                e.MaxY = ytics[y+1];

            if (x < xqtd && y < yqtd) {
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
                output << "\"avg_y\": " << cell->avg_y << ",";
                output << "\"type\": \"face\"}},\n";
            }

            // vertex
            output << "{\"type\": \"Feature\", \"geometry\": {\"type\": \"Point\", \"coordinates\": [" << 
                vertexes[v].x << ", " << vertexes[v].y << "]},";
            output << "\"properties\": {" <<
                      "\"name\": \"v:" << x << "." << y << "\", " <<
                      "\"card\": " << vertexes[v].cardin << ", " <<
                      "\"type\": \"vertex\"}}";
            if (x == xqtd && y == yqtd)
                output << "\n";
            else
                output << ",\n";
            v++;

            auto printedge = [&](auto& edge, const std::string& type) {
                output << "{\"type\": \"Feature\", \"geometry\": {\"type\": \"LineString\", \"coordinates\": [[" <<
                    edge.mbr.MinX << ", " << edge.mbr.MinY << "], [" << 
                    edge.mbr.MaxX << ", " << edge.mbr.MaxY << "]]},";
                output << "\"properties\": {" <<
                          "\"name\": \"eh:" << x << "." << y << "\"," <<
                          "\"card\": " << edge.cardin << "," <<
                          "\"type\": \"" << type << "\"}},\n";
            };

            // horizontal edge
            if (x < xqtd) {
                int e = GET_HORZ_EDGE(x, y);
                printedge(edges[e], "edgeh");
            }
            
            // vertical edge
            if (y < yqtd) {
                int e = GET_VERT_EDGE(x, y);
                printedge(edges[e], "edgev");
            }
        }
    }

    output << "]}\n";
    output.close();
}
/*
 * Eulerskew Spatial Histogram based on the work of
 * Swarup Acharya, Viswanath Poosala, and Sridhar Ramaswamy. Selectivity Estimation 
 * in Spatial Databases”. In: SIGMOD Record 28.2 (1999), pp. 13–24.
 * and
 * Sun, Chengyu, et al. "Exploring spatial datasets with histograms." 
 * Distributed and Parallel Databases 20.1 (2006): 57-88.
 *
 * Created on: 2022-07-07
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#include <geos/geom/Geometry.h>
#include "../include/SpatialHistogramEulerskew.hpp"
#include "rtree.h"

//#define DEBUG_QUERYNO
#ifdef DEBUG_QUERYNO
extern int global_qryno;
extern double global_query_size;
extern int global_debug_qryno;
extern double global_debug_qrysize;
extern rtree_root *global_search_rtree;
#endif

SpatialHistogramEulerskew::SpatialHistogramEulerskew(SpatialGridHistogram &hist, 
	Dataset &ds, int bucketsNeeded): basehist(hist) {
    this->bucketsNeeded = bucketsNeeded;
	this->hmbr = basehist.mbr();
	generateBuckets();
    generateEdgesAndVertices();
    hashSpatialObjects(ds);
}

double SpatialHistogramEulerskew::estimateWQuery(const Envelope& wquery) {
	double result = 0.0;
	for(const EulerskewBucket *b : buckets.getIntersections(wquery)) {
		if (wquery.intersects(b->usedarea)) {
			#define IMPROVED_IHWAF
			//#define DEFAULT_MP
			//#define AREA_BASED
			
			#ifdef AREA_BASED
			// TODO: This is a very simple uniformity assumption.
			// One can improve this by using the avg_x and avg_y
			// form the original histogram used to create the buckets
			Envelope inters = wquery.intersection(b->mbr);
			double fraction = inters.area() / b->mbr.area();
			result += fraction * b->cardin;
			#endif

			#ifdef DEFAULT_MP
			/* this formula came from the following paper:
				Nikos Mamoulis and Dimitris Papadias. “Advances in Spatial and Temporal Databases”. 
				In: ed. by Christian S. Jensen et al. Vol. 2121. Lecture Notes in Computer Science. 
				Springer, 2001. Chap. Selectivity Estimation of Complex Spatial Queries, pp. 155–174.
				See also Equation 2.2 in de Oliveira, T.B. thesis */
			// uniformity assumption! objects aren't uniform located at the cell
			Envelope inters = wquery.intersection(b->mbr);
			double xprob = std::min(1.0, (b->avg_x + inters.width()) / b->mbr.width());
			double yprob = std::min(1.0, (b->avg_y + inters.length()) / b->mbr.length());
			result += b.objcount * xprob * yprob;
			#endif

			#ifdef IMPROVED_IHWAF
			/* See de OLIVEIRA, T. B. Efficient Processing of Multiway Spatial Join Queries 
 			 * in Distributed Systems. 152 p. Tese (Doutorado) — Instituto de Informática, 
			 * Universidade Federal de Goiás, Goiânia, GO, Brasil, 2017.
			 */
			double avg_x = b->avg_x;
			double avg_y = b->avg_y;

			// observing that objects generally doesn't overlap in both axis,
			// (e.g. political land limits)
			// fix the probability of intersection in one of them
			double ux = b->usedarea.width();
			double uy = b->usedarea.length();
			if (avg_x > avg_y)
				avg_x = std::min(avg_x, avg_x/(avg_y * b->objcount / uy));
			else
				avg_y = std::min(avg_y, avg_y/(avg_x * b->objcount / ux));

			// uniformity assumption! objects aren't uniform located at the cell
			Envelope inters = wquery.intersection(b->usedarea);
			double xprob = std::min(1.0, (avg_x + inters.width()) / b->usedarea.width());
			double yprob = std::min(1.0, (avg_y + inters.length()) / b->usedarea.length());
			result += b->objcount * xprob * yprob;
			#endif

			#ifdef DEBUG_QUERYNO
            if (global_qryno == global_debug_qryno && global_query_size == global_debug_qrysize) {
				if (!b->usedarea.isEmpty()) {
					rtree_window_stat stats;
					memset(&stats, 0, sizeof(rtree_window_stat));
					Envelope inter = b->mbr.intersection(wquery);
					GEOSGeometry *geoquery = inter.toGEOSGeom();
					EnvelopeC queryc;
					queryc.MinX = inter.MinX;
					queryc.MaxX = inter.MaxX;
					queryc.MinY = inter.MinY;
					queryc.MaxY = inter.MaxY;
					GList *results = rtree_window_search(global_search_rtree, geoquery, queryc, &stats);
                	printf("Face %10.2lf %10u %10.2lf: %d\n", b->objcount * xprob * yprob, 
						g_list_length(results),
						(g_list_length(results) - b->objcount * xprob * yprob) * -1.0,
						b->id);
					g_list_free(results);
					GEOSGeom_destroy(geoquery);
				}
            }
			#endif
		}
	}

    // edges
	for(const EulerskewEdge *edge : edges.getIntersections(wquery)) {
		if (wquery.intersects(edge->mbr)) {
			double int_length, fraction;
			if (edge->isVertical()) {
				// ignora edges coincidentes com a borda da consulta
				if ((fabs(edge->mbr.MinX - wquery.MinX) <= 1e-10) ||
					(fabs(edge->mbr.MinX - wquery.MaxX) <= 1e-10) ||
                    (fabs(edge->mbr.MaxY - wquery.MinY) <= 1e-10) ||
                    (fabs(edge->mbr.MinY - wquery.MaxY) <= 1e-10))
					continue;

				Envelope inters = wquery.intersection(edge->mbr);
				int_length = inters.MaxY - inters.MinY;
				fraction = int_length / edge->mbr.length();
			} else {
				// ignora edges coincidentes com a borda da consulta
				if ((fabs(edge->mbr.MinY - wquery.MinY) <= 1e-10) ||
					(fabs(edge->mbr.MinY - wquery.MaxY) <= 1e-10) ||
                    (fabs(edge->mbr.MaxX - wquery.MinX) <= 1e-10) ||
                    (fabs(edge->mbr.MinX - wquery.MaxX) <= 1e-10))
					continue;

				Envelope inters = wquery.intersection(edge->mbr);
				int_length = inters.MaxX - inters.MinX;
				fraction = int_length / edge->mbr.width();
			}
			
			double aux = fraction * edge->cardin;
			result -= aux;

			#ifdef DEBUG_QUERYNO
            if (global_qryno == global_debug_qryno && global_query_size == global_debug_qrysize) {
                printf("Edge %d: %lf\n", edge->id, -aux);
            }
			#endif
		}
    }

    // vertices
    for(EulerskewVertex *vertex: vertices.getIntersectionsSL(wquery)) {
		//if (wquery.contains(vertex.x, vertex.y)) {
			// ignora vertices coincidentes com a borda da consulta
			if ((fabs(vertex->x - wquery.MinX) <= 1e-10) ||
				(fabs(vertex->x - wquery.MaxX) <= 1e-10) ||
				(fabs(vertex->y - wquery.MinY) <= 1e-10) ||
				(fabs(vertex->y - wquery.MaxY) <= 1e-10))
				continue;
			result += vertex->cardin;
			
			#ifdef DEBUG_QUERYNO
			if (global_qryno == global_debug_qryno && global_query_size == global_debug_qrysize) {
                printf("Vertex %d: %lf\n", vertex->id, vertex->cardin);
            }
			#endif
		//}
    }

	return result;
}

VarianceResult SpatialHistogramEulerskew::calculateSkewRowCol(int xini, int xfim, 
	int yini, int yfim) {
	
	VarianceResult r;
	r.n = 0;
	r.mean = 0.0;
	r.cardin = 0.0;
	r.objcount = 0;
	r.avg_x = 0.0;
	r.avg_y = 0.0;
	double M2 = 0.0;
	for(int i = xini; i <= xfim; i++) {
		for(int j = yini; j <= yfim; j++) {
			r.n++;
			
			SpatialHistogramCellDefault *cell = basehist.getHistogramCell(i, j);
			double v = cell->cardin;
			double delta = v - r.mean;
			r.mean += delta/(double)r.n;
			M2 += delta * (v - r.mean);

			// MinSkew doesn't have average lengths. This is an improvement.
			double nqtd = (r.cardin+v);
			if (nqtd > 0.0) {
				r.avg_x = ((r.avg_x * r.cardin) + (v * cell->avg_x)) / nqtd;
				r.avg_y = ((r.avg_y * r.cardin) + (v * cell->avg_y)) / nqtd;
				r.cardin = nqtd;
			}

			SpatialHistogramCellImproved *celli = dynamic_cast<SpatialHistogramCellImproved*>(cell);
			if (celli) {
				r.objcount += celli->objcount;
				r.usedarea.merge(celli->usedarea);
			}
		}
	}

	r.variance = r.n < 2 ? 0.0 : M2/(double)r.n;
	return r;
}

void SpatialHistogramEulerskew::calculateBucketWithMbr(EulerskewBucket &bucket) {
	
	int xini, xfim, yini, yfim;
	basehist.getIntersectionIdxs(bucket.mbr, &xini, &xfim, &yini, &yfim);
	VarianceResult vr = calculateSkewRowCol(xini, xfim, yini, yfim);
	bucket.skew = vr.n * vr.variance;
	bucket.cardin = vr.cardin;
	bucket.objcount = vr.objcount;
	bucket.avg_x = vr.avg_x;
	bucket.avg_y = vr.avg_y;
	bucket.skew_reduction = NAN;
	bucket.usedarea = vr.usedarea;
}

void SpatialHistogramEulerskew::minskewCalculateSkewReduction(EulerskewBucket &bucket) {

	int xini, xfim, yini, yfim;
	basehist.getIntersectionIdxs(bucket.mbr, &xini, &xfim, &yini, &yfim);

	bucket.skew_reduction = 0;

	EulerskewBucket aux_bucket1, aux_bucket2;
	aux_bucket1.mbr = bucket.mbr;
	aux_bucket2.mbr = bucket.mbr;

	for(int x=xini; x<xfim; x++) {
		// divide on x
		aux_bucket1.mbr.MaxX = aux_bucket2.mbr.MinX = basehist.getColumnX(x+1);
		calculateBucketWithMbr(aux_bucket1);
		calculateBucketWithMbr(aux_bucket2);

		double new_skew = aux_bucket1.skew + aux_bucket2.skew;
		double reduction = bucket.skew - new_skew;
		if (bucket.skew_reduction < reduction) {
			bucket.skew_reduction = reduction;
			bucket.split_axis = 'x';
			bucket.split_point = x+1;
		}
	}

	aux_bucket1.mbr = bucket.mbr;
	aux_bucket2.mbr = bucket.mbr;
	for(int y=yini; y<yfim; y++) {
		// divide on y
		aux_bucket1.mbr.MaxY = aux_bucket2.mbr.MinY = basehist.getRowY(y+1);
		calculateBucketWithMbr(aux_bucket1);
		calculateBucketWithMbr(aux_bucket2);

		double new_skew = aux_bucket1.skew + aux_bucket2.skew;
		double reduction = bucket.skew - new_skew;
		if (bucket.skew_reduction < reduction) {
			bucket.skew_reduction = reduction;
			bucket.split_axis = 'y';
			bucket.split_point = y+1;
		}
	}
}

void SpatialHistogramEulerskew::generateBuckets() {
	
	EulerskewBucket firstBucket;
	firstBucket.id = 0;
	firstBucket.mbr = basehist.mbr();
	calculateBucketWithMbr(firstBucket);
	
	buckets.push_back(firstBucket);
	double global_skew = firstBucket.skew;

	while (bucketsNeeded > buckets.size()) {
		EulerskewBucket *chosen = NULL;
		double skew_reduction = 0;

		for(EulerskewBucket& bucket: buckets) {
			if (isnan(bucket.skew_reduction)) {
				minskewCalculateSkewReduction(bucket);
			}

			if (skew_reduction < bucket.skew_reduction) {
				chosen = &bucket;
				skew_reduction = bucket.skew_reduction;
			}
		}

		// If chosen == NULL, no existing bucket can be
		// splitted in a way to improve the skew. In this
		// case, the histogram will not generate the quantity
		// of buckets determined in bucketsNeeded 
		if (!chosen)
			break;

		// split the chosen bucket
		EulerskewBucket newb1, newb2;
		newb1.mbr = chosen->mbr;
		newb2.mbr = chosen->mbr;
		if (chosen->split_axis == 'x')
			newb1.mbr.MaxX = newb2.mbr.MinX = basehist.getColumnX(chosen->split_point);
		else
			newb1.mbr.MaxY = newb2.mbr.MinY = basehist.getRowY(chosen->split_point);
		calculateBucketWithMbr(newb1);
		calculateBucketWithMbr(newb2);

		global_skew -= chosen->skew;
		global_skew += newb1.skew + newb2.skew;

 		assert(abs(chosen->cardin - newb1.cardin + newb2.cardin) < 1e100 && "Cardinality differs when bucket was spplited.");

		// substitute the old bucket data with b2
		newb1.id = chosen->id;
		*chosen = newb1;

		// add the new bucket b2
		newb2.id = buckets.size();
		buckets.push_back(newb2);
	}

	// recompute cardinality

}

void SpatialHistogramEulerskew::addEdgeIfNotExists(const Envelope& newEdgeMbr) {
    
    // It can be improved if edges are sorted (linesweep) or if
    // a R-tree like container is used.
    bool edgeExists = false;
    for(const EulerskewEdge& edge : edges) {
        if (edge.mbr.contains(newEdgeMbr) || newEdgeMbr.contains(edge.mbr)) {
            edgeExists = true;
            break;
        }
    }

	static int edge_count = 0;
    if (!edgeExists)
    {
        EulerskewEdge edge;
		edge.id = ++edge_count;
        edge.mbr = newEdgeMbr;
        edges.push_back(edge);
    }
}

void SpatialHistogramEulerskew::addVertexIfNotExists(double X, double Y) {
    
    // It can be improved if vertices are sorted (linesweep) or if
    // a R-tree like container is used.
    bool vertexExists = false;
    for(const EulerskewVertex& vertex : vertices) {
		if (vertex.x == X && vertex.y == Y) {
			vertexExists = true;
            break;
		}
    }

	static int vertex_count = 0; 
    if (!vertexExists) {
        EulerskewVertex v;
		v.id = ++vertex_count;
        v.x = X;
        v.y = Y;
        vertices.push_back(v);
    }
}

void SpatialHistogramEulerskew::generateEdgesAndVertices() {
    
	for(EulerskewBucket& bucket : buckets) {
		Envelope newEdge;

		// horizontal edge 1
		newEdge.MinX = bucket.mbr.MinX;
		newEdge.MinY = bucket.mbr.MinY;
		newEdge.MaxX = bucket.mbr.MaxX;
		newEdge.MaxY = bucket.mbr.MinY + 1e-10;
        addEdgeIfNotExists(newEdge);

		// horizontal edge 2
		newEdge.MinX = bucket.mbr.MinX;
		newEdge.MinY = bucket.mbr.MaxY;
		newEdge.MaxX = bucket.mbr.MaxX;
		newEdge.MaxY = bucket.mbr.MaxY + 1e-10;
        addEdgeIfNotExists(newEdge);

		// vertical edge 1
		newEdge.MinX = bucket.mbr.MinX;
		newEdge.MinY = bucket.mbr.MinY;
		newEdge.MaxX = bucket.mbr.MinX + 1e-10;
		newEdge.MaxY = bucket.mbr.MaxY;
        addEdgeIfNotExists(newEdge);

		// vertical edge 2
		newEdge.MinX = bucket.mbr.MaxX;
		newEdge.MinY = bucket.mbr.MinY;
		newEdge.MaxX = bucket.mbr.MaxX + 1e-10;
		newEdge.MaxY = bucket.mbr.MaxY;
        addEdgeIfNotExists(newEdge);

		// vertex: below left
        addVertexIfNotExists(bucket.mbr.MinX, bucket.mbr.MinY);

		// vertex: above left
		addVertexIfNotExists(bucket.mbr.MinX, bucket.mbr.MaxY);
        
		// vertex: below right
		addVertexIfNotExists(bucket.mbr.MaxX, bucket.mbr.MinY);
        
		// vertex: above right
		addVertexIfNotExists(bucket.mbr.MaxX, bucket.mbr.MaxY);
	}
}

void SpatialHistogramEulerskew::hashSpatialObjects(Dataset &ds) {
	buckets.buildRTree();
	edges.buildRTree();
	vertices.sortForSL();

	for(EulerskewBucket &b : buckets) {
		b.avg_x = b.avg_y = b.objcount = 0;
		b.usedarea = Envelope();
	}

    for(const DatasetEntry& de : ds.geoms()) {
		// recompute buckets cardinality
		for(EulerskewBucket *b : buckets.getIntersections(de.mbr)) {
			Envelope clipped_mbr;
			if (b->mbr.contains(de.mbr)) {
				clipped_mbr = de.mbr;
			} else {
				GEOSGeometry *clipped_geo = GEOSClipByRect(de.geo, b->mbr.MinX, b->mbr.MinY, 
					b->mbr.MaxX, b->mbr.MaxY);
				if (clipped_geo != NULL) {
					const geos::geom::Envelope *ev = ((const geos::geom::Geometry*)clipped_geo)->getEnvelopeInternal();
					if (!ev->isNull())
						clipped_mbr = ev;
					GEOSGeom_destroy(clipped_geo);
				}
			}
			if (!clipped_mbr.isEmpty()) {
				b->usedarea.merge(clipped_mbr);
				b->objcount += 1.0;

				// average length, online average calculation
				b->avg_x += (clipped_mbr.width() - b->avg_x) / b->objcount;
				b->avg_y += (clipped_mbr.length() - b->avg_y) / b->objcount;
			}
		}

		// edges
		for(EulerskewEdge *edge : edges.getIntersections(de.mbr)) {
            //if (edge.mbr.intersects(de.mbr)) {
                edge->cardin += 1;
                if (edge->isVertical()) {
                    double delta_y = de.mbr.MaxY - de.mbr.MinY;
                    edge->avg_length += (delta_y - edge->avg_length) / edge->cardin;
                } else {
                    double delta_x = de.mbr.MaxX - de.mbr.MinX;
                    edge->avg_length += (delta_x - edge->avg_length) / edge->cardin;
				}
            //}
		}

        // vertices
        for(EulerskewVertex *v : vertices.getIntersectionsSL(de.mbr)) {
            //if (de.mbr.contains(v->x, v->y)) {
				if (!v->aux_geometry)
					v->aux_geometry = GEOSGeom_createPointFromXY(v->x, v->y);
				if (GEOSContains(de.geo, v->aux_geometry)) {
					v->cardin += 1;				
				}
			//}
		}
    }

	for(EulerskewVertex &v : vertices) {
		GEOSGeom_destroy(v.aux_geometry);
		v.aux_geometry = NULL;
	}
}

void SpatialHistogramEulerskew::printGeoJson(const std::string& filename) {
	std::ofstream output;
	output.open(filename);
	output << std::fixed;
	output.precision(15);
	output << "{\"type\": \"FeatureCollection\", \"features\": [\n";

	for(EulerskewBucket& b : buckets) {
		if (!b.usedarea.isEmpty()) {
			output << "{\"type\": \"Feature\", \"geometry\": {\"type\": \"Polygon\", \"coordinates\": [[";
			output << "[" << b.usedarea.MinX << "," << b.usedarea.MinY << "],";
			output << "[" << b.usedarea.MaxX << "," << b.usedarea.MinY << "],";
			output << "[" << b.usedarea.MaxX << "," << b.usedarea.MaxY << "],";
			output << "[" << b.usedarea.MinX << "," << b.usedarea.MaxY << "],";
			output << "[" << b.usedarea.MinX << "," << b.usedarea.MinY << "]";
			output << "]]}, \"properties\": {";
			output << "\"name\": \"" << b.id << "\",";
			output << "\"card\": " << b.cardin << ",";
            output << "\"objects\": " << b.objcount << ",";
			output << "\"avg_x\": " << b.avg_x << ",";
			output << "\"avg_y\": " << b.avg_y << ",";
			output << "\"skew\": " << b.skew << "}},\n";
		}
	}

	for(EulerskewEdge& ee : edges) {
		output << "{\"type\": \"Feature\", \"geometry\": {\"type\": \"LineString\", \"coordinates\": [";
        output << "[" << ee.mbr.MinX << "," << ee.mbr.MinY << "],";
        output << "[" << ee.mbr.MaxX << "," << ee.mbr.MaxY << "]";
        output << "]},\"properties\": {";
		output << "\"name\": " << ee.id << ",";
		output << "\"card\": " << ee.cardin << ",";
        if (ee.isVertical())
            output << "\"avg_y\": " << ee.avg_length << ",";
        else
            output << "\"avg_x\": " << ee.avg_length << ",";
		output << "\"type\": \"edge\"";
		output << "}},\n";
	}

	// vertex
	bool first = true;
    for(EulerskewVertex &v : vertices) {
        if (!first)
            output << ",\n";
		first = false;

		output << "{\"type\": \"Feature\", \"geometry\": {\"type\": \"Point\", \"coordinates\": ";
        output << "[" << v.x << "," << v.y << "]},";
		output << "\"properties\": {";
		output << "\"name\": " << v.id << ",";
		output << "\"card\": " << v.cardin << ",";
		output << "\"type\": \"vertex\"}}";
	}

	output << "]}\n";
	output.close();
}

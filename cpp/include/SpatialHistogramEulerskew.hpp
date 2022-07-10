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

#pragma once

#include <string>
#include <list>
#include "Dataset.hpp"
#include "Envelope.hpp"
#include "SpatialGridHistogram.hpp"
#include "SpatialHistogramMinskew.hpp"
#include "DequeBinarySearchPoint.hpp"

struct EulerskewBucket {
	int id;
	double skew;
	char split_axis;
	short split_point;
	double skew_reduction;

	Envelope mbr;
	double cardin;
	double objcount;
	double avg_x;
	double avg_y;
	Envelope usedarea;
};

class EulerskewEdge {
public:
	int id;
    double cardin;
    double avg_length;
    Envelope mbr;

    EulerskewEdge() {
        id = 0;
		cardin = avg_length = 0.0;
    }

    bool isVertical() const {
        return mbr.width() < mbr.length();
    }
};

class EulerskewVertex {
public:
    int id;
	double x;
    double y;
    double cardin;
	GEOSGeometry *aux_geometry;
    EulerskewVertex() {
        id = 0;
		x = y = cardin = 0.0;
		aux_geometry = NULL;
    }
};

class SpatialHistogramEulerskew: public SpatialHistogram {
	public:
		SpatialHistogramEulerskew(SpatialGridHistogram &hist, Dataset &ds, int bucket_num);
		virtual double estimateWQuery(const Envelope& wquery) override;
		
		virtual const std::string name() override {
			return "Eulerskew";
		};

		virtual void printGeoJson(const std::string& filename) override;

		int bucketCount() {
			return buckets.size();
		}

	private:
        int bucketsNeeded;
		SpatialGridHistogram& basehist;
        SimpleEnvelopeRTree<EulerskewBucket, 30> buckets;
        SimpleEnvelopeRTree<EulerskewEdge, 30> edges;
        DequeBinarySearchPoint<EulerskewVertex> vertices;
		void generateBuckets();
        void generateEdgesAndVertices();
		void hashSpatialObjects(Dataset &ds);
		void addEdgeIfNotExists(const Envelope& newEdgeMbr);
		void addVertexIfNotExists(double X, double Y);
        void minskewCalculateSkewReduction(EulerskewBucket &bucket);
		void calculateBucketWithMbr(EulerskewBucket &bucket);
		VarianceResult calculateSkewRowCol(int xini, int xfim, int yini, int yfim);
};

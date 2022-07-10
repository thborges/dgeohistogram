/*
 * Minskew Spatial Histogram based on the work
 * Swarup Acharya, Viswanath Poosala, and Sridhar Ramaswamy. Selectivity Estimation 
 * in Spatial Databases”. In: SIGMOD Record 28.2 (1999), pp. 13–24.
 *
 *  Created on: 2022-02-23
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#pragma once

#include <string>
#include <list>
#include "Dataset.hpp"
#include "Envelope.hpp"
#include "SpatialGridHistogram.hpp"
#include "SimpleEnvelopeRTree.hpp"

struct MinskewBucket {
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

struct VarianceResult {
	double mean;
	double variance;
	int n;
	double cardin;
	double objcount;
	double avg_x;
	double avg_y;
	Envelope usedarea;
};

class SpatialHistogramMinskew: public SpatialHistogram {
	public:
		SpatialHistogramMinskew(SpatialGridHistogram &hist, int bucket_num);
		virtual double estimateWQuery(const Envelope& wquery) override;
		
		virtual const std::string name() override {
			return "MinSkew";
		};

		virtual void printGeoJson(const std::string& filename) override;

		int bucketCount() {
			return buckets.size();
		}

	private:
        int bucketsNeeded;
		SpatialGridHistogram& basehist;
        SimpleEnvelopeRTree<MinskewBucket, 30> buckets;
		void generateBuckets();
        void minskewCalculateSkewReduction(MinskewBucket &bucket);
		void calculateBucketWithMbr(MinskewBucket &bucket);
		VarianceResult calculateSkewRowCol(int xini, int xfim, int yini, int yfim);
};


/*
 * Minskew Spatial Histogram based on the work
 * Swarup Acharya, Viswanath Poosala, and Sridhar Ramaswamy. Selectivity Estimation 
 * in Spatial Databases”. In: SIGMOD Record 28.2 (1999), pp. 13–24.
 *
 *  Created on: 2022-02-23
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#include "../include/SpatialHistogramMinskew.hpp"

SpatialHistogramMinskew::SpatialHistogramMinskew(SpatialGridHistogram &hist, 
	int bucketsNeeded): basehist(hist) {
    this->bucketsNeeded = bucketsNeeded;
	this->hmbr = basehist.mbr();
	generateBuckets();
	buckets.buildRTree();
}

double SpatialHistogramMinskew::estimateWQuery(const Envelope& wquery) {
	double result = 0.0;
	// A simple loop through buckets. Can be improved if a RTree is used to
	// store them; or a sweep line algorithm
	for(const MinskewBucket* b : buckets.getIntersections(wquery)) {
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
			/* this formula came from the following paper of same authors:
				Nikos Mamoulis and Dimitris Papadias. “Advances in Spatial and Temporal Databases”. 
				In: ed. by Christian S. Jensen et al. Vol. 2121. Lecture Notes in Computer Science. 
				Springer, 2001. Chap. Selectivity Estimation of Complex Spatial Queries, pp. 155–174.
				See also Equation 2.2 in de Oliveira, T.B. thesis */
			// uniformity assumption! objects aren't uniform located at the cell
			Envelope inters = wquery.intersection(b->mbr);
			double xprob = std::min(1.0, (b->avg_x + inters.width()) / b->mbr.width());
			double yprob = std::min(1.0, (b->avg_y + inters.length()) / b->mbr.length());
			result += b->cardin * xprob * yprob;
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
			result += b->cardin * xprob * yprob;
			#endif
		}
	}
	return result;
}

VarianceResult SpatialHistogramMinskew::calculateSkewRowCol(int xini, int xfim, 
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
			} else {
				// if not using improved cell, use the mbr as usedarea.
				// e.g. when building minskew from MP histogram
				Envelope e(basehist.getColumnX(i), basehist.getRowY(j),
					basehist.getColumnX(i+1), basehist.getRowY(j+1));
				r.usedarea.merge(e);
			}
		}
	}

	r.variance = r.n < 2 ? 0.0 : M2/(double)r.n;
	return r;
}

void SpatialHistogramMinskew::calculateBucketWithMbr(MinskewBucket &bucket) {
	
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

void SpatialHistogramMinskew::minskewCalculateSkewReduction(MinskewBucket &bucket) {

	int xini, xfim, yini, yfim;
	basehist.getIntersectionIdxs(bucket.mbr, &xini, &xfim, &yini, &yfim);

	bucket.skew_reduction = 0;

	MinskewBucket aux_bucket1, aux_bucket2;
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

void SpatialHistogramMinskew::generateBuckets() {
	
	MinskewBucket firstBucket;
	firstBucket.id = 0;
	firstBucket.mbr = basehist.mbr();
	calculateBucketWithMbr(firstBucket);
	
	buckets.push_back(firstBucket);
	double global_skew = firstBucket.skew;

	while (bucketsNeeded > buckets.size()) {
		MinskewBucket *chosen = NULL;
		double skew_reduction = 0;

		for(MinskewBucket& bucket: buckets) {
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
		MinskewBucket newb1, newb2;
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
}

void SpatialHistogramMinskew::printGeoJson(const std::string& filename) {
	std::ofstream output;
	output.open(filename);
	output << std::fixed;
	output.precision(15);
	output << "{\"type\": \"FeatureCollection\", \"features\": [\n";

	bool first = true;
	for(MinskewBucket& b : buckets) {
		if (!b.usedarea.isEmpty()) {
			if (!first)
				output << ",\n";
			first = false;
			output << "{\"type\": \"Feature\", \"geometry\": {\"type\": \"Polygon\", \"coordinates\": [[";
			output << "[" << b.usedarea.MinX << "," << b.usedarea.MinY << "],";
			output << "[" << b.usedarea.MaxX << "," << b.usedarea.MinY << "],";
			output << "[" << b.usedarea.MaxX << "," << b.usedarea.MaxY << "],";
			output << "[" << b.usedarea.MinX << "," << b.usedarea.MaxY << "],";
			output << "[" << b.usedarea.MinX << "," << b.usedarea.MinY << "]";
			output << "]]}, \"properties\": {";
			output << "\"name\": \"" << b.id << "\",";
			output << "\"card\": " << b.cardin << ",";
			output << "\"avg_x\": " << b.avg_x << ",";
			output << "\"avg_y\": " << b.avg_y << ",";
			output << "\"skew\": " << b.skew << "}}";
			//output << "\"avg_x\": " << hcells[x*yqtd + y].avg_x << ",";
			//output << "\"avg_y\": " << hcells[x*yqtd + y].avg_y;
		}
	}

	output << "]}\n";
	output.close();
}
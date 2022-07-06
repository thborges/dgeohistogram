/*
 * A simple experiment using aleatory window queries
 *
 *  Created on: 2022-02-15
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#include <iostream>
#include "Dataset.hpp"
#include "SpatialGridHistogramMP.hpp"
#include "SpatialGridHistogramIHWAF.hpp"
#include "SpatialGridHistogramEuler.hpp"
#include "SpatialHistogramMinskew.hpp"
#include "UniformWQueryExperiment.hpp"

void printMessageAndGenGeoJson(SpatialHistogramMinskew& hist, const std::string& filename);
void printMessageAndGenGeoJson(SpatialGridHistogram& hist, const std::string& filename);
extern "C" void geos_messages(const char *fmt, ...);

int main(int argc, char *argv[]) {
	initGEOS(geos_messages, geos_messages);

	if (argc <= 1) {
		std::cout << argv[0] << " filename.[shp,geojson]\n";
		return 1;
	}
	std::string filename(argv[1]);

	try {
		Dataset ds(filename);
		std::cout << "Loaded " << filename << ". Objects: " << ds.size() << std::endl;

		// IHWAF histogram
		SpatialGridHistogramIHWAF histIHWAF(ds);
		printMessageAndGenGeoJson(histIHWAF, filename);

		int cols = histIHWAF.columns();
		int rows = histIHWAF.rows();

		// Mamoulis/Papadias histogram
		SpatialGridHistogramMP histMP(ds, cols, rows);
		printMessageAndGenGeoJson(histMP, filename);

		// MinSkew histogram
		//SpatialHistogramMinskew histMinSkew(histMP, 0.1 * cols * rows);
		//SpatialHistogramMinskew histMinSkew(histMP, 1000);
		SpatialHistogramMinskew histMinSkew(histIHWAF, 0.8 * cols * rows);
		printMessageAndGenGeoJson(histMinSkew, filename);

		// Euler histogram
		SpatialGridHistogramEuler histEuler(ds, cols, rows);
		printMessageAndGenGeoJson(histEuler, filename);

		// Histogram list to experiment with
		std::vector<SpatialHistogram*> hists;
		hists.push_back(&histMP);
		hists.push_back(&histIHWAF);
		hists.push_back(&histMinSkew);
		hists.push_back(&histEuler);

		// Query sizes
		std::vector<double> qsizes;
		//qsizes.push_back(0.01);
		qsizes.push_back(0.05);
		qsizes.push_back(0.10);
		qsizes.push_back(0.15);
		qsizes.push_back(0.20);
		qsizes.push_back(0.25);
		qsizes.push_back(0.30);

		UniformWQueryExperiment exp(ds, hists, qsizes);
		exp.run();
	}
	catch (const std::exception& e) {
		std::cerr << e.what();
	}

	return 0;
}

void printMessageAndGenGeoJson(SpatialHistogramMinskew& hist, const std::string& filename) {
	std::cout << hist.name() << "\t" << hist.bucketCount() << "\t";
	std::string histname(filename);
	histname += "." + hist.name() + ".geojson";
	hist.printGeoJson(histname);
	std::cout << histname << std::endl;
}

void printMessageAndGenGeoJson(SpatialGridHistogram& hist, const std::string& filename) {
	std::cout << hist.name() << "\t" << hist.columns() << "x" << hist.rows() << "\t";
	std::string histname(filename);
	histname += "." + hist.name() + ".geojson";
	hist.printGeoJson(histname);
	std::cout << histname << std::endl;
}

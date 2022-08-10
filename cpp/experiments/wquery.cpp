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
#include "SpatialHistogramEulerskew.hpp"
#include "UniformWQueryExperiment.hpp"

void printMessageAndGenGeoJson(SpatialHistogramMinskew& hist, const std::string& filename);
void printMessageAndGenGeoJson(SpatialHistogramEulerskew& hist, const std::string& filename);
void printMessageAndGenGeoJson(SpatialGridHistogram& hist, const std::string& filename);
extern "C" void geos_messages(const char *fmt, ...);

int main(int argc, char *argv[]) {
	initGEOS(geos_messages, geos_messages);

	if (argc <= 1 || argc == 3) {
		std::cout << argv[0] << " filename.[shp,geojson] [cols rows]\n";
		return 1;
	}
	std::string filename(argv[1]);
	int cols, rows;
	bool input_cols_rows = argc > 2;
	if (input_cols_rows) {
		cols = atoi(argv[2]);
		rows = atoi(argv[3]);
	}

	try {
		Dataset ds(filename);
		std::cout << "Loaded " << filename << ". Objects: " << ds.size() << std::endl;

		// IHWAF histogram
		SpatialGridHistogramIHWAF *histIHWAF;
		if (input_cols_rows)
			histIHWAF = new SpatialGridHistogramIHWAF(ds, cols, rows);
		else
			histIHWAF = new SpatialGridHistogramIHWAF(ds);
		printMessageAndGenGeoJson(*histIHWAF, filename);

		int cols = histIHWAF->columns();
		int rows = histIHWAF->rows();

		// Mamoulis/Papadias histogram
		SpatialGridHistogramMP histMP(ds, cols, rows);
		printMessageAndGenGeoJson(histMP, filename);

		// MinSkew histogram
		SpatialHistogramMinskew histMinSkew(histMP, 0.5 * cols * rows);
		//SpatialHistogramMinskew histMinSkew(*histIHWAF, 0.3 * (cols * rows));
		printMessageAndGenGeoJson(histMinSkew, filename);

		// Eulerskew histogram
		SpatialHistogramEulerskew histEulerskew(*histIHWAF, ds, 0.5 * (cols * rows));
		printMessageAndGenGeoJson(histEulerskew, filename);

		// Euler histogram
		SpatialGridHistogramEuler histEuler(ds, cols, rows);
		printMessageAndGenGeoJson(histEuler, filename);

		// Histogram list to experiment with
		std::vector<SpatialHistogram*> hists;
		hists.push_back(&histMP);
		hists.push_back(histIHWAF);
		hists.push_back(&histMinSkew);
		hists.push_back(&histEulerskew);
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

		delete histIHWAF;
	}
	catch (const std::exception& e) {
		std::cerr << e.what();
	}

	return 0;
}

void printMessageAndGenGeoJson(SpatialHistogramMinskew& hist, const std::string& filename) {
	std::cout << std::setw(10) << hist.name() << "\t" << hist.bucketCount() << "\t";
	std::string histname(filename);
	histname += "." + hist.name() + ".geojson";
	hist.printGeoJson(histname);
	std::cout << histname << std::endl;
}

void printMessageAndGenGeoJson(SpatialHistogramEulerskew& hist, const std::string& filename) {
	std::cout << std::setw(10) << hist.name() << "\t" << hist.bucketCount() << "\t";
	std::string histname(filename);
	histname += "." + hist.name() + ".geojson";
	hist.printGeoJson(histname);
	std::cout << histname << std::endl;
}

void printMessageAndGenGeoJson(SpatialGridHistogram& hist, const std::string& filename) {
	std::cout << std::setw(10) << hist.name() << "\t" << hist.columns() << "x" << hist.rows() << "\t";
	std::string histname(filename);
	histname += "." + hist.name() + ".geojson";
	hist.printGeoJson(histname);
	std::cout << histname << std::endl;
}

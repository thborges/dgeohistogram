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
#include "RandomWQueryExperiment.hpp"

void printMessageAndGenGeoJson(SpatialGridHistogram& hist, const std::string& filename) {
	std::cout << "Built grid histogram " << hist.name() << " with " << 
	hist.columns() << "x" << hist.rows() << " cols/rows" << std::endl;
	std::string histname(filename);
	histname += "." + hist.name() + ".geojson";
	hist.printGeoJson(histname);
	std::cout << "Generated " << histname << " with the histogram data." << std::endl;
}

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

		SpatialGridHistogramIHWAF histIHWAF(ds);
		printMessageAndGenGeoJson(histIHWAF, filename);

		// Mamoulis/Papadias histogram
		SpatialGridHistogramMP histMP(ds, histIHWAF.columns(), histIHWAF.rows());
		printMessageAndGenGeoJson(histMP, filename);

		std::vector<SpatialHistogram*> hists;
		hists.push_back(&histMP);
		hists.push_back(&histIHWAF);

		std::vector<double> qsizes;
		qsizes.push_back(0.01);
		qsizes.push_back(0.05);
		qsizes.push_back(0.10);
		qsizes.push_back(0.15);
		qsizes.push_back(0.20);
		qsizes.push_back(0.25);
		qsizes.push_back(0.30);

		RandomWQueryExperiment exp(ds, hists, qsizes);
		exp.run();
	}
	catch (const std::exception& e) {
		std::cerr << e.what();
	}

	return 0;
}


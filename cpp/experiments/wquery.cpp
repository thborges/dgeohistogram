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
#include "SpatialHistogramAB.hpp"
#include "UniformWQueryExperiment.hpp"

void printMessageAndGenGeoJson(SpatialHistogramMinskew& hist, const std::string& filename);
void printMessageAndGenGeoJson(SpatialHistogramAB& hist, const std::string& filename);
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

		const DatasetMetadata& metadata = ds.metadata();
    	double dsw = metadata.mbr.MaxX - metadata.mbr.MinX;
    	double dsh = metadata.mbr.MaxY - metadata.mbr.MinY;
		double ratio = dsw / dsh;
		int cols = ceil(dsw / ratio);
		int rows = ceil(dsh / ratio);
		printf("%dx%d\n", cols, rows);
		// Query sizes
		std::vector<double> qsizes;
		//qsizes.push_back(0.01);
		qsizes.push_back(0.05);
		qsizes.push_back(0.10);
		qsizes.push_back(0.15);
		qsizes.push_back(0.20);
		qsizes.push_back(0.25);
		qsizes.push_back(0.30);

		std::vector<int> gridSizes;
		gridSizes.push_back(5);
		gridSizes.push_back(15);
		gridSizes.push_back(25);
		gridSizes.push_back(50);
		gridSizes.push_back(100);
		gridSizes.push_back(125);
		gridSizes.push_back(150);
		gridSizes.push_back(175);
		gridSizes.push_back(200);
		gridSizes.push_back(225);

		// MinSkew histogram
		//SpatialHistogramMinskew histMinSkew(histMP, 0.1 * histMP.columns() * histMP.rows());
		//SpatialHistogramMinskew histMinSkew(histMP, 1000);
		//SpatialHistogramMinskew histMinSkew(histIHWAF, 0.8 * histIHWAF.columns() * histIHWAF.rows());
		//printMessageAndGenGeoJson(histMinSkew, filename);

		for (auto size : gridSizes) {
			// IHWAF histogram
			auto start = std::chrono::system_clock::now();
			SpatialGridHistogramIHWAF histIHWAF(ds);
			std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - start;
			std::cout << elapsed.count() * 1000.0 << "\t";
			printMessageAndGenGeoJson(histIHWAF, filename);

			// Mamoulis/Papadias histogram
			start = std::chrono::system_clock::now();
			SpatialGridHistogramMP histMP(ds, size, size);
			elapsed = std::chrono::system_clock::now() - start;
			std::cout << elapsed.count() * 1000.0 << "\t";
			printMessageAndGenGeoJson(histMP, filename);

			// Euler histogram
			start = std::chrono::system_clock::now();
			SpatialGridHistogramEuler histEuler(ds, size, size);
			elapsed = std::chrono::system_clock::now() - start;
			std::cout << elapsed.count() * 1000.0 << "\t";
			printMessageAndGenGeoJson(histEuler, filename);

			// AB histogram
			start = std::chrono::system_clock::now();
			SpatialHistogramAB histAB(ds, size, size);
			elapsed = std::chrono::system_clock::now() - start;
			std::cout << elapsed.count() * 1000.0 << "\t";
			printMessageAndGenGeoJson(histAB, filename);

			// Histogram list to experiment with
			std::vector<SpatialHistogram*> hists;
			hists.push_back(&histMP);
			hists.push_back(&histIHWAF);
			//hists.push_back(&histMinSkew);
			hists.push_back(&histEuler);
			hists.push_back(&histAB);

			UniformWQueryExperiment exp(ds, hists, qsizes);
			exp.run(filename);
		}

	}
	catch (const std::exception& e) {
		std::cerr << e.what();
	}

	return 0;
}

void printMessageAndGenGeoJson(SpatialHistogramMinskew& hist, const std::string& filename) {
    /*
	std::cout << hist.name() << "\t"
			  //<< hist.bucketCount() << "\t\t"
			  << hist.getSize() << " KB\t";
    */
	std::string histname(filename);
	histname += "." + hist.name() + ".geojson";
	hist.printGeoJson(histname);
	//std::cout << histname << std::endl;
}

void printMessageAndGenGeoJson(SpatialHistogramAB& hist, const std::string& filename) {
    /*
	std::cout << hist.name() << "\t"
			  << hist.getcolumns() << "x" << hist.getrows() << " (" << hist.getNumBuckets() << ")\t"
			  << hist.getSize() << " KB\t";
    */
	std::cout << hist.name() << "\t"
			  << hist.getcolumns() << "x" << hist.getrows() << "\t"
			  << hist.getNumBuckets() << "\t"
			  << hist.getSize() << "\n";
	std::string histname(filename);
	histname += "." + std::to_string(hist.getcolumns()) + "." + hist.name() + ".geojson";
	hist.printGeoJson(histname);
	//std::cout << histname << std::endl;
}

void printMessageAndGenGeoJson(SpatialGridHistogram& hist, const std::string& filename) {
    /*
	std::cout << hist.name() << "\t"
			  //<< hist.columns() << "x" << hist.rows() << "\t\t"
			  << hist.getSize() << " KB\t";
    */
	std::cout << hist.name() << "\t"
			  << hist.columns() << "x" << hist.rows() << "\t"
			  << hist.getSize() << "\n";
	std::string histname(filename);
	histname += "." + std::to_string(hist.columns()) + "." + hist.name() + ".geojson";
	hist.printGeoJson(histname);
	//std::cout << histname << std::endl;
}

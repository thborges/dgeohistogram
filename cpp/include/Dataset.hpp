/*
 * Dataset Container
 *
 *  Created on: 2022-02-15
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#pragma once

#include <string>
#include <list>
#include <geos_c.h>
#include <ogrsf_frmts.h>
#include "Envelope.hpp"
#include "Property.hpp"

struct DatasetEntry {
	GEOSGeometry* geo;
	Envelope mbr;
};

struct DatasetMetadata {
	std::string name;
	Envelope mbr;
	double x_average;
	double y_average;
	double x_psa; //x power sum average for stddev
	double y_psa; //y power sum average for stddev
	unsigned count;
	//unsigned geocount;
	OGRwkbGeometryType geomtype;

	double x_stdev() const { return std::sqrt(x_psa / (count-1)); }
	double y_stdev() const { return std::sqrt(y_psa / (count-1)); }
	
	DatasetMetadata() {
		x_average = y_average = x_psa = y_psa = count = 0;
	}
};

class Dataset {
	public:
		Dataset(const std::string filename);
		
		size_t size() {
			return geometries.size();
		}

		const std::list<DatasetEntry>& geoms() {
			return geometries;
		}

		const DatasetMetadata& metadata() {
			return meta;
		}

	private:
		std::list<DatasetEntry> geometries;
		DatasetMetadata meta;
		
		void addGeomAndBuildMeta(const DatasetEntry& e);
};



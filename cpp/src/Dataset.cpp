/*
 * Dataset Container
 *
 *  Created on: 2022-02-15
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#include "../include/Dataset.hpp"
#include <iostream>

Dataset::Dataset(const std::string filename) {
	
	GDALAllRegister();
	GDALDatasetUniquePtr gdalds(GDALDataset::Open(filename.c_str(), GDAL_OF_VECTOR));
	if (gdalds == NULL)
		throw std::runtime_error("Invalid shapefile.");
	
	size_t pos = filename.find_last_of("/");
	meta.name = filename.substr(pos != std::string::npos ? pos+1 : 0);

	if (gdalds->GetLayerCount() > 1)
		throw std::runtime_error("Currently, datasets with more than one layer are not supported.");

	size_t wkb_max = 1024;
	unsigned char *wkb = new unsigned char[wkb_max];
	for(OGRLayer* layer: gdalds->GetLayers()) {
		
		// if the dataset has more than one layer, metadata could be inconsistent
		meta.geomtype = OGR_GT_Flatten(layer->GetGeomType());
        
		for(const auto& feature: *layer) {
			const OGRGeometry *geometry = feature->GetGeometryRef();
			if (geometry != NULL) {
				size_t wkb_size = geometry->WkbSize();

				// increase memory if needed
				if (wkb_size > wkb_max) {
					delete[] wkb;
					wkb = new unsigned char[wkb_size];
					wkb_max = wkb_size;
				}
				
				OGREnvelope env;
				geometry->getEnvelope(&env);
				meta.mbr.merge(env.MinX, env.MinY, env.MaxX, env.MaxY);

				geometry->exportToWkb(wkbXDR, wkb);

				GEOSGeometry *ggeo = GEOSGeomFromWKB_buf(wkb, wkb_size);
				Envelope mbr(env.MinX, env.MinY, env.MaxX, env.MaxY);
				addGeomAndBuildMeta({ggeo, mbr});
			}
		}
	}

	std::cout << "Dataset object average: x: " << meta.x_average << ", y: " << meta.y_average << std::endl;
	delete[] wkb;
}

void Dataset::addGeomAndBuildMeta(const DatasetEntry& e) {

	geometries.push_back(e);
	meta.count++;
	
	// metadata: object extent
	double new_x = e.mbr.MaxX - e.mbr.MinX;
	double new_y = e.mbr.MaxY - e.mbr.MinY;

	if (meta.count == 1) {
		meta.x_average = new_x;
		meta.y_average = new_y;
		meta.x_psa = meta.y_psa = 0.0;
	} else {
		double oldmx = meta.x_average;
		double oldmy = meta.y_average;
		meta.x_average += (new_x - meta.x_average) / meta.count;
		meta.y_average += (new_y - meta.y_average) / meta.count;
		meta.x_psa += (new_x - oldmx) * (new_x - meta.x_average);
		meta.y_psa += (new_y - oldmy) * (new_y - meta.y_average);
	}
}
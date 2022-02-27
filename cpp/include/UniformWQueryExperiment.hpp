/*
 * Random window query experiment
 *
 *  Created on: 2022-02-16
 *      Author: Thiago Borges de Oliveira <thborges@gmail.com>
 */

#pragma once

#include "GenericExperiment.hpp"
#include "SpatialHistogram.hpp"
#include "Dataset.hpp"

#include "rtree.h"
#include "rtree-star.h"

class UniformWQueryExperiment: public GenericExperiment {
public:
    UniformWQueryExperiment(Dataset& d, std::vector<SpatialHistogram*>& hs, 
		std::vector<double> qsizes): ds(d), hists(hs), 
		query_sizes(qsizes) {}

    const virtual void run() {
		// cria uma r*
		rtree_root *rtree = NULL;
		rtree = rtree_new_rstar(30, 10);

		std::printf("Building R*-Tree\n");
		for(const DatasetEntry& de: ds.geoms()) {
			EnvelopeC env{de.mbr.MinX, de.mbr.MinY, de.mbr.MaxX, de.mbr.MaxY};
			rtree_append(rtree, de.geo, env);
		}

		std::printf("\nSize\tARE\t\tSTD\t\tSUM\t\tMethod\tName\n");
		const DatasetMetadata& meta = ds.metadata();
		rtree_window_stat stats;
		double width = meta.mbr.MaxX - meta.mbr.MinX;
		double height = meta.mbr.MaxY - meta.mbr.MinY;

		for(double query_size: query_sizes) {
			// metrics for each histogram
			int hs = hists.size();
			double sum_ei[hs];
			double sum_ri[hs];
			double mean[hs];
			double M2[hs];
			double sum_error[hs];
			for(int i = 0; i < hs; i++) {
				sum_ei[i] = 0.0;
				sum_ri[i] = 0.0;
				mean[i] = 0.0;
				M2[i] = 0.0;
				sum_error[i] = 0.0;
			}

			double wsize = width * query_size;
			double hsize = height * query_size;

			// quantidade de consultas para cobrir o dataset, considerando
			// uma distribuicao uniforme
			int qtd = ceil((width / wsize) * (height/hsize));
			int qtdqx = (width / wsize);

			//std::printf("Query count: %d, w %f, h %f, w_size %f, h_size %f\n", qtd, width, height, wsize, hsize);

			int n = 0;
			int qryno = 0;
			while (n < qtd) {
				n++;

				EnvelopeC queryc;
				queryc.MinX = meta.mbr.MinX + (qryno / qtdqx) * wsize;
				queryc.MinY = meta.mbr.MinY + (qryno % qtdqx) * hsize;
				queryc.MaxX = queryc.MinX + wsize;
				queryc.MaxY = queryc.MinY + hsize;
				Envelope query{queryc.MinX, queryc.MinY, queryc.MaxX, queryc.MaxY};

				char wkt[512];
				std::sprintf(wkt, "POLYGON((%lf %lf, %lf %lf, %lf %lf, %lf %lf, %lf %lf))",
					query.MinX, query.MinY,
					query.MaxX, query.MinY,
					query.MaxX, query.MaxY,
					query.MinX, query.MaxY,
					query.MinX, query.MinY);
			
				GEOSGeometryH geoquery = GEOSGeomFromWKT(wkt);

				memset(&stats, 0, sizeof(rtree_window_stat));
				GList *results = rtree_window_search(rtree, geoquery, queryc, &stats);

				// real cardinality from rtree
				int riq = g_list_length(results);

				for(int i = 0; i < hs; i++) {
					// histogram estimate cardinality
					int rhq = hists[i]->estimateWQuery(query);
					//printf("Query %d: r: %5d, e: %5d, %5d\n", n, riq, rhq, rhq - riq);

					int error = abs(rhq-riq);

					// average relative error
					sum_ei[i] += error;
					sum_ri[i] += riq;

					// stdev of error
					double delta = error - mean[i];
					mean[i] += delta/(double)n;
					M2[i] += delta*(error - mean[i]);

					// sum error
					sum_error[i] += error;
				}

				g_list_free(results);
				qryno++;
			}
		
			for(int i = 0; i < hs; i++) {
				std::printf("%3.2lf\t%lf\t%lf\t%lf\t%s\t%s\n",
					query_size, 
					sum_ei[i] / (double)sum_ri[i],
					sqrt(M2[i]/(double)n),
					sum_error[i],
					hists[i]->name().c_str(),
					ds.metadata().name.c_str());
			}
			std::printf("\n");
		}
	}

private:
	Dataset& ds;
    std::vector<SpatialHistogram*>& hists;
	std::vector<double> query_sizes;
};
